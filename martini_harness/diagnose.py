"""Finding out why a step failed from the program output it left behind.

GROMACS ends a failed run with a block like::

    Program:     gmx mdrun, version 2024.2
    Source file: src/gromacs/mdlib/sim_util.cpp (line 297)

    Fatal error:
    At step 2 detected non-finite forces on 2 atoms

    For more information and tips for troubleshooting, ...

The message is extracted and classified (LINCS, GPU, non-finite forces from
``-pforce``, ...). Warnings printed before the error, such as LINCS warnings
or the large forces listed by ``-pforce``, are counted as notes.
"""

import re
from dataclasses import dataclass, field
from pathlib import Path

from .config import COORDINATES_STEP

_TAIL_BYTES = 256 * 1024
_MAX_MESSAGE = 300
_HELP_LINE = "For more information and tips for troubleshooting"
_BLOCK_KEYS = ("Program:", "Source file:", "Function:", "MPI rank:", "---")

# (reason, pattern) pairs, checked in order against the error message.
_REASONS = [
    ("GPU", re.compile(r"\b(GPUs?|CUDA|OpenCL|SYCL|HIP)\b|cuda\w*Error")),
    ("LINCS", re.compile(r"LINCS")),
    ("non-finite forces (pforce)", re.compile(r"non-finite forces")),
    ("non-finite energy", re.compile(r"not finite|\bnan\b", re.IGNORECASE)),
    ("too many grompp warnings", re.compile(r"Too many warnings")),
    (
        "atoms moved too far",
        re.compile(r"moved more than|domain decomposition|cut-?off length is longer", re.IGNORECASE),
    ),
    ("out of memory", re.compile(r"out of memory|bad_alloc|not enough memory", re.IGNORECASE)),
    ("topology", re.compile(r"moleculetype|atom name|No such|does not match|Incorrect number")),
]
_LARGE_FORCE_LINE = re.compile(r"^step\s+\d+\s+atom\s+\d+\s+x\s", re.MULTILINE)
_SIGNAL = re.compile(r"exited with code (-\d+)")


@dataclass
class Diagnosis:
    reason: str
    message: str = ""
    notes: list[str] = field(default_factory=list)

    def __str__(self) -> str:
        text = self.reason
        if self.notes:
            text += f" ({'; '.join(self.notes)})"
        if self.message:
            text += f": {self.message}"
        return text


def gromacs_error(text: str) -> str | None:
    """The message of the last GROMACS error block in ``text``, on one line."""
    end = text.rfind(_HELP_LINE)
    if end == -1:
        return None
    start = text.rfind("Program:", 0, end)
    lines = [
        line.strip()
        for line in text[max(start, 0) : end].splitlines()
        if line.strip() and not line.strip().startswith(_BLOCK_KEYS)
    ]
    if lines and lines[0].endswith(":"):  # "Fatal error:", "Inconsistency in user input:", ...
        lines = lines[1:]
    message = " ".join(lines)
    return message[: _MAX_MESSAGE - 3] + "..." if len(message) > _MAX_MESSAGE else message


def diagnose_text(text: str) -> Diagnosis | None:
    """Classify the program output of a failed GROMACS run."""
    notes = []
    lincs_warnings = text.count("LINCS WARNING")
    large_forces = len(_LARGE_FORCE_LINE.findall(text))

    message = gromacs_error(text)
    if message is not None:
        reason = next((name for name, pattern in _REASONS if pattern.search(message)), "other")
    elif "Segmentation fault" in text:
        reason, message = "crashed (segmentation fault)", ""
    elif lincs_warnings:
        reason, message = "LINCS", "LINCS warnings but no error message; the run was interrupted"
    else:
        return None

    if lincs_warnings and reason == "non-finite energy":
        # The usual way a system blows up: constraints fail first, then energies diverge.
        reason = "LINCS, then non-finite energy"
        notes.append(f"{lincs_warnings} LINCS warning(s)")
    elif lincs_warnings and reason != "LINCS":
        notes.append(f"after {lincs_warnings} LINCS warning(s)")
    if large_forces and reason != "non-finite forces (pforce)":
        notes.append(f"{large_forces} large force(s) listed by -pforce")
    return Diagnosis(reason, message, notes)


def diagnose_step(workdir: Path, step: str, error: str | None = None) -> Diagnosis | None:
    """Work out why ``step`` failed from the files it wrote in ``workdir``.

    Args:
        workdir: The simulation directory.
        step: Step name (a stage name, or ``coords``).
        error: The exception message recorded for the failure, if any.
    """
    if error:
        signal = _SIGNAL.search(error)
        if signal:
            number = -int(signal.group(1))
            hint = " (out of memory or wall-time limit?)" if number == 9 else ""
            return Diagnosis(f"killed by signal {number}{hint}")
    if step == COORDINATES_STEP:
        return Diagnosis("coordinate generation (polyply)", (error or "")[:_MAX_MESSAGE])

    grompp_out = workdir / f"{step}.grompp.out"
    candidates = [workdir / f"{step}.mdrun.out", workdir / f"{step}.log"]
    if grompp_out.is_file():
        # Output of mdrun older than this attempt's grompp is from an earlier attempt.
        started = grompp_out.stat().st_mtime
        candidates = [path for path in candidates if path.is_file() and path.stat().st_mtime >= started]
        candidates.append(grompp_out)
    for path in candidates:
        if path.is_file():
            diagnosis = diagnose_text(_tail(path))
            if diagnosis is not None:
                return diagnosis
    return None


def _tail(path: Path) -> str:
    with open(path, "rb") as handle:
        handle.seek(0, 2)
        handle.seek(max(0, handle.tell() - _TAIL_BYTES))
        return handle.read().decode("utf-8", errors="replace")
