"""Turning a campaign into simulations and running them step by step.

Every simulation consists of a fixed ``topology`` and ``coords`` step followed
by the stages from the campaign file. Each step has a cache key derived from
its inputs and the key of the step before it; a step is skipped when the state
file says it finished with the same key and its outputs still exist. Once one
step runs, all later steps run too.
"""

import dataclasses
import functools
import hashlib
import json
import logging
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .backend import Backend, BackendError, GromacsWrapperBackend
from .config import (
    COORDINATES_STEP,
    TOPOLOGY_STEP,
    CampaignConfig,
    ConfigError,
    SimulationDefaults,
    StageConfig,
    SAFE_NAME,
)
from .diagnose import diagnose_step
from .library import LibraryError, MoleculeLibrary
from .matrix import expand_matrix
from .mdp import TemplateError, is_minimization, render_mdp, render_value
from .state import DONE, FAILED, StateError, StateFile
from .topology import count_atoms, read_moleculetypes, render_topology

logger = logging.getLogger(__name__)

SPEC_FILE_NAME = "simulation.json"

# Step statuses reported by Simulation.status()
STATUS_DONE = "done"
STATUS_PENDING = "pending"
STATUS_STALE = "stale"
STATUS_FAILED = "failed"

# SimulationOutcome.failed when a simulation fails outside of a step
SETUP_FAILURE = "setup"


@dataclass
class SimulationSpec:
    """One fully resolved matrix combination."""

    name: str
    target: str
    solvent: str
    n_target: int
    n_solvent: int
    density: float
    # (molecule, count) pairs: the solvent itself, or the components of a mixture
    solvent_molecules: list[tuple[str, int]]
    includes: list[Path]
    workdir: Path
    variables: dict[str, Any] = field(default_factory=dict)
    target_atoms: int | None = None


@dataclass
class Step:
    name: str
    signature: dict[str, Any]
    outputs: list[Path]
    action: Callable[[Backend], None]
    key: str = ""


@dataclass
class StepStatus:
    name: str
    status: str
    error: str | None = None
    reason: str | None = None


@dataclass
class SimulationOutcome:
    name: str
    ran: list[str] = field(default_factory=list)
    skipped: list[str] = field(default_factory=list)
    failed: str | None = None
    error: str | None = None
    reason: str | None = None

    @property
    def success(self) -> bool:
        return self.failed is None


class Simulation:
    """A single simulation: its spec plus the ordered steps that produce it."""

    def __init__(self, spec: SimulationSpec, config: CampaignConfig):
        self.spec = spec
        self.steps = plan_steps(spec, config)
        key = ""
        for step in self.steps:
            key = _hash({"step": step.name, "signature": step.signature, "upstream": key})
            step.key = key

    @property
    def name(self) -> str:
        return self.spec.name

    def _is_done(self, state: StateFile, step: Step) -> bool:
        return state.is_done(step.name, step.key) and all(path.exists() for path in step.outputs)

    def status(self) -> list[StepStatus]:
        """Report the status of every step without running anything."""
        try:
            state = StateFile(self.spec.workdir)
        except StateError as exc:
            return [StepStatus(step.name, STATUS_FAILED, str(exc)) for step in self.steps[:1]] + [
                StepStatus(step.name, STATUS_PENDING) for step in self.steps[1:]
            ]
        statuses = []
        blocked = False
        for step in self.steps:
            entry = state.get(step.name)
            if blocked:
                statuses.append(StepStatus(step.name, STATUS_PENDING))
            elif self._is_done(state, step):
                statuses.append(StepStatus(step.name, STATUS_DONE))
            else:
                blocked = True
                if entry is None:
                    statuses.append(StepStatus(step.name, STATUS_PENDING))
                elif entry["status"] == FAILED:
                    statuses.append(
                        StepStatus(step.name, STATUS_FAILED, entry.get("error"), entry.get("reason"))
                    )
                else:
                    statuses.append(StepStatus(step.name, STATUS_STALE))
        return statuses

    def run(self, backend: Backend, force: bool = False) -> SimulationOutcome:
        """Run all steps that are not done yet.

        Args:
            backend: Backend that executes polyply and GROMACS.
            force: Rerun every step, ignoring the state file.

        Returns:
            What ran, what was skipped and what failed. Errors are caught and
            reported in the outcome rather than raised.
        """
        outcome = SimulationOutcome(self.name)
        workdir = self.spec.workdir
        try:
            workdir.mkdir(parents=True, exist_ok=True)
            (workdir / SPEC_FILE_NAME).write_text(
                json.dumps(dataclasses.asdict(self.spec), indent=2, default=str)
            )
            state = StateFile(workdir)
        except (OSError, StateError) as exc:
            logger.error("[%s] cannot prepare %s: %s", self.name, workdir, exc)
            outcome.failed = SETUP_FAILURE
            outcome.error = str(exc)
            return outcome

        names = [step.name for step in self.steps]
        rerun = force
        for index, step in enumerate(self.steps):
            if not rerun and self._is_done(state, step):
                outcome.skipped.append(step.name)
                continue
            rerun = True
            state.discard(names[index:])
            for path in step.outputs:
                path.unlink(missing_ok=True)

            logger.info("[%s] running %s", self.name, step.name)
            try:
                step.action(backend)
                missing = [path.name for path in step.outputs if not path.exists()]
                if missing:
                    raise BackendError(f"expected output(s) missing: {', '.join(missing)}")
            except (Exception, SystemExit) as exc:
                diagnosis = diagnose_step(workdir, step.name, str(exc))
                reason = str(diagnosis) if diagnosis else None
                logger.error("[%s] %s failed: %s", self.name, step.name, reason or exc)
                state.record(step.name, step.key, FAILED, str(exc), reason)
                outcome.failed = step.name
                outcome.error = str(exc)
                outcome.reason = reason
                return outcome
            state.record(step.name, step.key, DONE)
            outcome.ran.append(step.name)
        return outcome


def plan_steps(spec: SimulationSpec, config: CampaignConfig) -> list[Step]:
    """Build the ordered steps of a simulation.

    Raises:
        ConfigError: If an MDP template is missing or a template fails to render.
    """
    workdir = spec.workdir
    top_file = workdir / "system.top"
    coords_file = workdir / "init.gro"

    top_text = render_topology(
        spec.includes,
        spec.name,
        [(spec.target, spec.n_target), *spec.solvent_molecules],
    )

    def write_topology(backend: Backend) -> None:
        top_file.write_text(top_text)

    def generate_coordinates(backend: Backend) -> None:
        backend.generate_coordinates(
            top_file, coords_file, spec.density, spec.name, workdir / "coords.log"
        )

    steps = [
        Step(
            TOPOLOGY_STEP,
            {"topology": top_text, "includes": {str(p): _file_hash(p) for p in spec.includes}},
            [top_file],
            write_topology,
        ),
        Step(COORDINATES_STEP, {"density": spec.density}, [coords_file], generate_coordinates),
    ]

    structures = {COORDINATES_STEP: coords_file}
    previous = COORDINATES_STEP
    for stage in config.stages:
        steps.append(_stage_step(stage, spec, config, top_file, structures, previous))
        structures[stage.name] = workdir / f"{stage.name}.gro"
        previous = stage.name
    return steps


def _stage_step(
    stage: StageConfig,
    spec: SimulationSpec,
    config: CampaignConfig,
    top_file: Path,
    structures: dict[str, Path],
    previous: str,
) -> Step:
    workdir = spec.workdir
    mdp_template = config.mdp_path(stage)
    if not mdp_template.is_file():
        raise ConfigError(f"stage {stage.name!r}: MDP file not found: {mdp_template}")
    try:
        mdp_text = render_mdp(
            mdp_template.read_text(),
            substitute=render_value(stage.substitute, spec.variables),
            params=render_value(stage.mdp_params, spec.variables),
        )
        grompp_flags = render_value(stage.grompp, spec.variables)
        mdrun_flags = render_value(stage.mdrun, spec.variables)
        # Campaign-wide mdrun flags (threads, GPU offload) are about the
        # machine, not the result, so they stay out of the step signature.
        global_mdrun_flags = (
            {} if is_minimization(mdp_text) else render_value(config.mdrun, spec.variables)
        )
    except TemplateError as exc:
        raise ConfigError(f"stage {stage.name!r}: {exc}") from exc

    run_mdrun_flags = {**global_mdrun_flags, **mdrun_flags}
    if spec.target_atoms == 1 and str(run_mdrun_flags.get("bonded", "")).lower() == "gpu":
        # mdrun refuses -bonded gpu when there is nothing bonded to compute.
        del run_mdrun_flags["bonded"]
        logger.debug("[%s] %s: single-bead target, not passing -bonded gpu", spec.name, stage.name)

    input_step = stage.input or previous
    input_structure = structures[input_step]
    mdp_file = workdir / f"{stage.name}.mdp"
    tpr_file = workdir / f"{stage.name}.tpr"

    def run_stage(backend: Backend) -> None:
        mdp_file.write_text(mdp_text)
        backend.grompp(
            {
                "f": mdp_file,
                "c": input_structure,
                "p": top_file,
                "o": tpr_file,
                "po": workdir / f"{stage.name}_mdout.mdp",
                "maxwarn": stage.maxwarn,
                **grompp_flags,
            },
            log=workdir / f"{stage.name}.grompp.out",
        )
        backend.mdrun(
            {"s": tpr_file, "deffnm": workdir / stage.name, **run_mdrun_flags},
            log=workdir / f"{stage.name}.mdrun.out",
        )

    signature = {
        "mdp": mdp_text,
        "input": input_step,
        "maxwarn": stage.maxwarn,
        "grompp": grompp_flags,
        "mdrun": mdrun_flags,
    }
    return Step(stage.name, signature, [tpr_file, workdir / f"{stage.name}.gro"], run_stage)


def build_simulations(
    config: CampaignConfig, output_dir: Path | None = None
) -> list[Simulation]:
    """Expand the build matrix and resolve every combination into a Simulation.

    Args:
        config: The campaign configuration.
        output_dir: Overrides ``config.output_dir``.

    Raises:
        ConfigError: Listing every problem found across all combinations.
    """
    output_dir = Path(output_dir).resolve() if output_dir else config.output_dir
    errors = []

    defined_in: dict[str, Path] = {}
    for include in config.topology.includes:
        if not include.is_file():
            errors.append(f"topology include not found: {include}")
            continue
        for molname in read_moleculetypes(include):
            defined_in.setdefault(molname, include)

    molecule_itps: dict[str, Path] = {}
    for molname, itp in config.molecules.items():
        if not itp.is_file():
            errors.append(f"molecule {molname!r}: ITP file not found: {itp}")
        elif molname not in read_moleculetypes(itp):
            errors.append(f"molecule {molname!r}: no [ moleculetype ] {molname} in {itp}")
        else:
            molecule_itps[molname] = itp

    library = None
    if config.molecule_library is not None:
        try:
            library = MoleculeLibrary(config.molecule_library)
        except LibraryError as exc:
            errors.append(str(exc))
    if errors:
        raise ConfigError("\n".join(errors))

    combinations = expand_matrix(
        config.matrix.axes, config.matrix.include, config.matrix.exclude
    )
    param_names = set(SimulationDefaults.model_fields)
    simulations = []
    seen_names: set[str] = set()
    for index, combo in enumerate(combinations, start=1):
        label = f"matrix combination {index} {combo}"
        missing = [key for key in ("target", "solvent") if key not in combo]
        if missing:
            errors.append(f"{label}: missing {', '.join(missing)}")
            continue

        try:
            params = SimulationDefaults.model_validate(
                {
                    **config.simulation.model_dump(),
                    **{key: value for key, value in combo.items() if key in param_names},
                }
            )
        except ValueError as exc:
            errors.append(f"{label}: {exc}")
            continue

        variables = {**params.model_dump(exclude={"name"}), **combo}
        try:
            name = render_value(params.name, variables)
        except TemplateError as exc:
            errors.append(f"{label}: simulation name: {exc}")
            continue
        if not SAFE_NAME.match(name):
            errors.append(f"{label}: simulation name {name!r} is not a valid directory name")
            continue
        if name in seen_names:
            errors.append(f"{label}: duplicate simulation name {name!r}")
            continue
        seen_names.add(name)
        variables["name"] = name

        solvent = str(combo["solvent"])
        if solvent in config.solvents:
            try:
                solvent_molecules = split_count(config.solvents[solvent], params.n_solvent)
            except ValueError as exc:
                errors.append(f"{name}: solvent {solvent!r}: {exc}")
                continue
        else:
            solvent_molecules = [(solvent, params.n_solvent)]

        includes = list(config.topology.includes)
        target_itp = None
        roles = [("target", str(combo["target"]))]
        roles += [("solvent", molname) for molname, _ in solvent_molecules]
        for role, molname in roles:
            if molname in molecule_itps:
                itp = molecule_itps[molname]
            elif molname in defined_in:
                if role == "target":
                    target_itp = defined_in[molname]
                continue
            elif library is not None and molname in library:
                itp = library.path(molname)
                if not itp.is_file():
                    errors.append(f"{name}: {role} {molname!r}: library file missing: {itp}")
                    continue
            else:
                errors.append(
                    f"{name}: {role} {molname!r} is not defined in the topology "
                    "includes, the molecules section or the molecule library"
                )
                continue
            if role == "target":
                target_itp = itp
            if itp not in includes:
                includes.append(itp)

        spec = SimulationSpec(
            name=name,
            target=str(combo["target"]),
            solvent=solvent,
            n_target=params.n_target,
            n_solvent=params.n_solvent,
            density=params.density,
            solvent_molecules=solvent_molecules,
            includes=includes,
            workdir=output_dir / name,
            variables=variables,
            target_atoms=count_atoms(target_itp, str(combo["target"])) if target_itp else None,
        )
        try:
            simulations.append(Simulation(spec, config))
        except ConfigError as exc:
            errors.append(f"{name}: {exc}")

    if errors:
        raise ConfigError("\n".join(errors))
    return simulations


def split_count(amounts: dict[str, float], total: int) -> list[tuple[str, int]]:
    """Divide ``total`` molecules over a mixture in proportion to ``amounts``.

    Rounds by largest remainder, so the counts always add up to ``total``.

    Raises:
        ValueError: If a component would get no molecules.
    """
    scale = total / sum(amounts.values())
    exact = {molname: amount * scale for molname, amount in amounts.items()}
    counts = {molname: int(value) for molname, value in exact.items()}
    by_remainder = sorted(exact, key=lambda molname: exact[molname] - counts[molname], reverse=True)
    for molname in by_remainder[: total - sum(counts.values())]:
        counts[molname] += 1
    empty = [molname for molname, count in counts.items() if count == 0]
    if empty:
        raise ValueError(f"n_solvent {total} leaves no {', '.join(empty)} molecules")
    return list(counts.items())


def select_simulations(simulations: list[Simulation], names: list[str]) -> list[Simulation]:
    """Return the simulations with the given names, in campaign order."""
    known = {sim.name for sim in simulations}
    unknown = [name for name in names if name not in known]
    if unknown:
        raise ConfigError(
            f"unknown simulation(s): {', '.join(unknown)}; available: {', '.join(sorted(known))}"
        )
    return [sim for sim in simulations if sim.name in names]


def run_campaign(
    simulations: list[Simulation],
    backend: Backend | None = None,
    force: bool = False,
) -> list[SimulationOutcome]:
    """Run simulations one after another.

    A simulation that fails, for whatever reason, is reported in its outcome
    and skipped; the loop continues with the next one. Only an interrupt
    (Ctrl-C) stops the campaign.
    """
    backend = backend or GromacsWrapperBackend()
    outcomes = []
    for index, simulation in enumerate(simulations, start=1):
        logger.info("(%d/%d) simulation %s", index, len(simulations), simulation.name)
        try:
            outcome = simulation.run(backend, force=force)
        except (Exception, SystemExit) as exc:
            logger.exception("[%s] unexpected error, skipping", simulation.name)
            outcome = SimulationOutcome(
                simulation.name, failed=SETUP_FAILURE, error=f"{type(exc).__name__}: {exc}"
            )
        outcomes.append(outcome)
    return outcomes


def _hash(data: Any) -> str:
    return hashlib.sha256(json.dumps(data, sort_keys=True, default=str).encode()).hexdigest()


def _file_hash(path: Path) -> str | None:
    try:
        stat = path.stat()
    except OSError:
        return None
    return _cached_file_hash(path, stat.st_size, stat.st_mtime_ns)


@functools.lru_cache(maxsize=4096)
def _cached_file_hash(path: Path, size: int, mtime_ns: int) -> str:
    # size and mtime are part of the cache key, so an edited file is rehashed.
    return hashlib.sha256(path.read_bytes()).hexdigest()
