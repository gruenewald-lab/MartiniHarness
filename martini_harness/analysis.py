"""Running an analysis program over every simulation of a campaign.

An analysis runs a command (a GROMACS tool such as ``gmx awh`` or any other
program) inside the directory of each simulation whose stage finished, then
takes the last data line of one output file as that simulation's result, or the
signed difference of one column between its first and last data line (e.g. the
free energy difference of an AWH PMF). The results of all simulations are
collected into one table.
"""

import csv
import logging
import re
import shlex
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .check import expected_end_time, first_xvg_row, last_xvg_row, stage_run_settings
from .config import AnalysisConfig
from .mdp import TemplateError, render_value
from .simulation import Simulation

logger = logging.getLogger(__name__)

OK = "ok"
NOT_FINISHED = "not-finished"
FAILED = "failed"

_LEGEND = re.compile(r'^@\s+s(\d+)\s+legend\s+"(.*)"')
_XAXIS_LABEL = re.compile(r'^@\s+xaxis\s+label\s+"(.*)"')
# xmgrace markup such as \f{Symbol}, \S, \N in axis labels and legends
_XVG_MARKUP = re.compile(r"\\f\{[^}]*\}|\\[A-Za-z]")


class AnalysisError(RuntimeError):
    """Raised when an analysis fails for one simulation."""


@dataclass
class AnalysisResult:
    name: str
    status: str
    values: dict[str, str] = field(default_factory=dict)
    error: str | None = None


def run_analysis(
    analysis: AnalysisConfig, simulation: Simulation, nsteps_override: int | None = None
) -> AnalysisResult:
    """Run ``analysis`` on one simulation; errors are reported in the result.

    Args:
        analysis: The analysis to run.
        simulation: The simulation whose directory is analysed.
        nsteps_override: Campaign-wide mdrun ``-nsteps``, if any, for ``{end_time}``.
    """
    workdir = simulation.spec.workdir
    if not (workdir / f"{analysis.stage}.gro").is_file():
        return AnalysisResult(simulation.name, NOT_FINISHED)

    log = workdir / f"{analysis.name}.analysis.out"
    try:
        # GROMACS takes everything after the last dot of a path as the file
        # extension, so the scratch directory is passed relative to workdir
        # (the command's cwd) and its name has no dots.
        with tempfile.TemporaryDirectory(prefix="analysis_tmp_", dir=workdir) as tmpdir:
            variables = _variables(analysis, simulation, Path(tmpdir), nsteps_override)
            command = _render(analysis.command, variables)
            argv = shlex.split(command) if isinstance(command, str) else [str(arg) for arg in command]
            stdin = _render(analysis.stdin, variables)
            output_name = _render(analysis.output, variables)
            output = workdir / output_name

            logger.debug("[%s] %s: %s", simulation.name, analysis.name, shlex.join(argv))
            with open(log, "w") as handle:
                try:
                    process = subprocess.run(
                        argv, cwd=workdir, input=stdin, text=True,
                        stdout=handle, stderr=subprocess.STDOUT,
                    )
                except OSError as exc:
                    raise AnalysisError(f"cannot run {argv[0]}: {exc}") from exc
            if process.returncode != 0:
                raise AnalysisError(f"{argv[0]} exited with code {process.returncode}; see {log}")

            last = last_xvg_row(output) if output.is_file() else None
            if last is None:
                found = "not written" if not output.exists() else "has no data lines"
                raise AnalysisError(f"output {output_name} {found}; see {log}")
            if analysis.result == "last":
                names = analysis.columns or _xvg_column_names(output)
                values = {
                    names[index] if index < len(names) else f"col{index}": text
                    for index, text in enumerate(last)
                }
            else:
                values = {
                    (analysis.columns or [analysis.name])[0]: _difference(analysis, output, last)
                }
    except (AnalysisError, TemplateError) as exc:
        return AnalysisResult(simulation.name, FAILED, error=str(exc))
    return AnalysisResult(simulation.name, OK, values)


def _difference(analysis: AnalysisConfig, output: Path, last: list[str]) -> str:
    """Signed difference of ``analysis.column`` between the first and the last data line."""
    first = first_xvg_row(output)
    try:
        first_value, last_value = float(first[analysis.column]), float(last[analysis.column])
    except (IndexError, ValueError, TypeError) as exc:
        raise AnalysisError(f"{output.name} has no numeric column {analysis.column}") from exc
    value = first_value - last_value if analysis.result == "first-last" else last_value - first_value
    return f"{value:.6g}"


def run_analyses(
    analysis: AnalysisConfig, simulations: list[Simulation], nsteps_override: int | None = None
) -> list[AnalysisResult]:
    """Run ``analysis`` on every simulation, one after another."""
    results = []
    for index, simulation in enumerate(simulations, start=1):
        result = run_analysis(analysis, simulation, nsteps_override)
        if result.status == FAILED:
            logger.error("[%s] %s failed: %s", simulation.name, analysis.name, result.error)
        else:
            logger.debug("(%d/%d) %s: %s", index, len(simulations), simulation.name, result.status)
        results.append(result)
    return results


def write_table(path: Path, simulations: list[Simulation], results: list[AnalysisResult]) -> None:
    """Write one tab-separated row per simulation: its variables, status and values."""
    rows = []
    for simulation, result in zip(simulations, results):
        spec = simulation.spec
        variables = {"target": spec.target, "solvent": spec.solvent, **spec.variables}
        variables.pop("name", None)
        rows.append(
            {"name": result.name, **variables, "status": result.status, **result.values,
             "error": result.error or ""}
        )
    fields = list(dict.fromkeys(key for row in rows for key in row if key != "error"))
    fields.append("error")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", restval="", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _variables(
    analysis: AnalysisConfig, simulation: Simulation, tmpdir: Path, nsteps_override: int | None
) -> dict[str, Any]:
    variables = {
        **simulation.spec.variables,
        "stage": analysis.stage,
        "workdir": str(simulation.spec.workdir),
        "tmpdir": tmpdir.name,
    }
    mdp_text, nsteps = stage_run_settings(simulation, analysis.stage, nsteps_override)
    end_time = expected_end_time(mdp_text, nsteps, interval=None)
    if end_time is not None:
        variables["end_time"] = f"{end_time:g}"  # as in GROMACS file names (awh_t10000.xvg)
    return variables


def _render(value: Any, variables: dict[str, Any]) -> Any:
    if "{end_time}" in str(value) and "end_time" not in variables:
        raise TemplateError(f"end time of the run is unknown (nsteps < 0?) in {value!r}")
    return render_value(value, variables)


def _xvg_column_names(path: Path) -> list[str]:
    """Column names from the xvg header: the x-axis label, then the legends."""
    xlabel, legends = "x", {}
    with open(path, errors="replace") as handle:
        for line in handle:
            if not line.startswith(("#", "@")) and line.strip():
                break
            if match := _XAXIS_LABEL.match(line):
                xlabel = _clean_label(match.group(1)) or xlabel
            elif match := _LEGEND.match(line):
                legends[int(match.group(1))] = _clean_label(match.group(2))
    names = [xlabel]
    for index in range(max(legends, default=-1) + 1):
        names.append(legends.get(index) or f"col{index + 1}")
    return names


def _clean_label(label: str) -> str:
    return " ".join(_XVG_MARKUP.sub("", label).split())
