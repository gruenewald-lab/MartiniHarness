"""Checking whether the production stage of each simulation finished.

A production run counts as complete when the last frame of its dH/dλ output
(``<stage>.xvg``, written by mdrun with ``-deffnm <stage>``) has been
written, i.e. its last time reaches the end time of the run. This looks at
the data itself, so it also works when a state file is missing or a run was
done outside the harness.
"""

import os
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from .mdp import get_mdp_param
from .simulation import STATUS_DONE, STATUS_FAILED, Simulation

COMPLETE = "complete"
INCOMPLETE = "incomplete"
FAILED = "failed"
NOT_STARTED = "not-started"
CATEGORIES = (COMPLETE, INCOMPLETE, FAILED, NOT_STARTED)

# GROMACS defaults for parameters that may be missing from an MDP file
_MDP_DEFAULTS = {"tinit": "0", "dt": "0.001", "nsteps": "0", "nstdhdl": "50"}
_TAIL_BYTES = 64 * 1024
_TIME_TOLERANCE = 1e-6


@dataclass
class ProductionCheck:
    name: str
    category: str
    last_time: float | None = None
    end_time: float | None = None
    failed_step: str | None = None
    error: str | None = None
    next_step: str | None = None

    @property
    def progress(self) -> float | None:
        """Fraction of the run written to the xvg file, if known."""
        if self.last_time is None or not self.end_time:
            return None
        return min(self.last_time / self.end_time, 1.0)


def last_xvg_time(path: Path) -> float | None:
    """Time (first column) of the last complete data line of an xvg file.

    Only the end of the file is read. A partially written last line (from a
    run that is still writing or was killed) is ignored.
    """
    try:
        with open(path, "rb") as handle:
            handle.seek(0, os.SEEK_END)
            size = handle.tell()
            handle.seek(max(0, size - _TAIL_BYTES))
            tail = handle.read().decode("ascii", errors="replace")
    except OSError:
        return None

    lines = [
        line for line in tail.splitlines()
        if line.strip() and not line.lstrip().startswith(("#", "@"))
    ]
    if size > _TAIL_BYTES and lines:
        lines = lines[1:]  # the first line of the tail may be cut off
    if not tail.endswith("\n") and lines:
        lines = lines[:-1]  # still being written
    for line in reversed(lines):
        try:
            return float(line.split()[0])
        except ValueError:
            continue
    return None


def expected_end_time(mdp_text: str, nsteps_override: int | None = None) -> float | None:
    """Time of the last dH/dλ frame a run with these MDP settings writes.

    Returns None if it cannot be known (infinite run or no dH/dλ output).
    """

    def value(key: str) -> str:
        return get_mdp_param(mdp_text, key) or _MDP_DEFAULTS[key]

    try:
        tinit, dt = float(value("tinit")), float(value("dt"))
        nsteps = int(nsteps_override if nsteps_override is not None else value("nsteps"))
        nstdhdl = int(value("nstdhdl"))
    except ValueError:
        return None
    if nsteps < 0 or nstdhdl <= 0:
        return None
    return tinit + (nsteps // nstdhdl) * nstdhdl * dt


def check_production(simulation: Simulation, stage: str, nsteps_override: int | None = None) -> ProductionCheck:
    """Classify one simulation by the state of its production stage ``stage``."""
    workdir = simulation.spec.workdir
    check = ProductionCheck(simulation.name, NOT_STARTED)

    mdout = workdir / f"{stage}_mdout.mdp"
    step = next(step for step in simulation.steps if step.name == stage)
    mdp_text = mdout.read_text() if mdout.is_file() else step.signature["mdp"]
    if nsteps_override is None and "nsteps" in step.signature["mdrun"]:
        nsteps_override = int(step.signature["mdrun"]["nsteps"])
    check.end_time = expected_end_time(mdp_text, nsteps_override)

    xvg = workdir / f"{stage}.xvg"
    if xvg.is_file():
        check.last_time = last_xvg_time(xvg)

    statuses = simulation.status()
    if (
        check.end_time is not None
        and check.last_time is not None
        and check.last_time >= check.end_time - _TIME_TOLERANCE
    ):
        check.category = COMPLETE
        return check

    failed = next((s for s in statuses if s.status == STATUS_FAILED), None)
    if failed is not None:
        check.category = FAILED
        check.failed_step = failed.name
        check.error = failed.error
    elif xvg.is_file():
        check.category = INCOMPLETE
    check.next_step = next(
        (s.name for s in statuses if s.status != STATUS_DONE), f"none ({stage}.xvg missing)"
    )
    return check


def summarize(checks: list[ProductionCheck], stage: str) -> str:
    """Human-readable summary of production checks."""
    total = len(checks)
    by_category = Counter(check.category for check in checks)

    def share(count: int) -> str:
        return f"{count:>6}  ({100 * count / total:5.1f}%)" if total else f"{count:>6}"

    lines = [f"Production stage {stage!r}: {total} simulation(s)"]
    lines.append(f"  complete     {share(by_category[COMPLETE])}")

    incomplete = [c for c in checks if c.category == INCOMPLETE]
    detail = ""
    progress = sorted(p for p in (c.progress for c in incomplete) if p is not None)
    if progress:
        detail = f"  running or interrupted, median progress {100 * progress[len(progress) // 2]:.0f}%"
    lines.append(f"  incomplete   {share(len(incomplete))}{detail}")

    failed = Counter(c.failed_step for c in checks if c.category == FAILED)
    detail = f"  at {_format_counts(failed)}" if failed else ""
    lines.append(f"  failed       {share(by_category[FAILED])}{detail}")

    waiting = Counter(c.next_step or "-" for c in checks if c.category == NOT_STARTED)
    detail = f"  next step: {_format_counts(waiting)}" if waiting else ""
    lines.append(f"  not started  {share(by_category[NOT_STARTED])}{detail}")

    unknown = sum(1 for c in checks if c.end_time is None)
    if unknown:
        lines.append(
            f"  ({unknown} simulation(s) have no known end time: nsteps < 0 or nstdhdl = 0; "
            "they can never count as complete)"
        )
    return "\n".join(lines)


def _format_counts(counts: Counter) -> str:
    return ", ".join(f"{key} {count}" for key, count in counts.most_common())
