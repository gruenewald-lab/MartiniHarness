"""Campaign configuration: YAML schema and loading.

Relative paths in a campaign file are resolved against the directory of that
file, except ``stages[].mdp`` (relative to ``mdp_dir``) and
``topology.includes`` (relative to ``forcefield_dir``).

A campaign file may list other YAML files under ``imports``. Each imported
file is a partial campaign (and may import further files); they are
deep-merged in order, then the importing file's own keys are merged on top.
Mappings merge key by key, everything else (lists included) is replaced.
Paths in an imported file are resolved against that file's directory.
"""

import re
from pathlib import Path
from typing import Any

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    ValidationError,
    ValidationInfo,
    field_validator,
    model_validator,
)
from ruamel.yaml import YAML

PACKAGE_DATA = Path(__file__).parent / "data"
DEFAULT_FORCEFIELD_DIR = PACKAGE_DATA / "gmx_files"
DEFAULT_MDP_DIR = PACKAGE_DATA / "mdp_files"
DEFAULT_INCLUDES = ["martini_v3.0.0.itp", "martini_v3.0.0_solvents_v1.itp"]

# Steps that run before the user-defined stages; stage names may not clash.
TOPOLOGY_STEP = "topology"
COORDINATES_STEP = "coords"
RESERVED_STEP_NAMES = {TOPOLOGY_STEP, COORDINATES_STEP}

IMPORTS_KEY = "imports"
# Top-level keys holding a path relative to the file that sets them.
FILE_RELATIVE_KEYS = ("output_dir", "forcefield_dir", "mdp_dir", "molecule_library")

SAFE_NAME = re.compile(r"^[A-Za-z0-9._+-]+$")

Scalar = str | int | float | bool


class ConfigError(ValueError):
    """Raised when a campaign configuration is invalid."""


class _Model(BaseModel):
    model_config = ConfigDict(extra="forbid")


class StageConfig(_Model):
    """One grompp + mdrun step of the simulation pipeline.

    String values in ``substitute``, ``mdp_params``, ``grompp`` and ``mdrun``
    may reference matrix variables, e.g. ``"{target}"``.
    """

    name: str
    mdp: str
    input: str | None = Field(
        default=None,
        description="Step whose .gro file is the input structure (default: previous step).",
    )
    maxwarn: int = Field(default=0, ge=0)
    substitute: dict[str, str] = Field(
        default_factory=dict, description="Literal placeholder replacements in the MDP."
    )
    mdp_params: dict[str, Scalar | list[Scalar]] = Field(
        default_factory=dict, description="MDP parameters to set or override."
    )
    grompp: dict[str, Any] = Field(default_factory=dict, description="Extra grompp flags.")
    mdrun: dict[str, Any] = Field(default_factory=dict, description="Extra mdrun flags.")

    @field_validator("name")
    @classmethod
    def _check_name(cls, name: str) -> str:
        if not SAFE_NAME.match(name):
            raise ValueError(f"stage name {name!r} may only contain letters, digits and ._+-")
        if name in RESERVED_STEP_NAMES:
            raise ValueError(f"stage name {name!r} is reserved")
        return name


class AnalysisConfig(_Model):
    """A program run in the directory of every simulation whose ``stage`` finished.

    ``command``, ``output`` and ``stdin`` may reference matrix variables plus
    ``{stage}``, ``{workdir}``, ``{tmpdir}`` (a scratch directory inside the
    simulation directory, where the command runs; removed afterwards) and ``{end_time}`` (simulated time of the stage in ps,
    formatted like GROMACS file names, e.g. ``10000``).
    """

    name: str
    stage: str | None = Field(default=None, description="Stage analysed (default: the last stage).")
    command: str | list[str] = Field(
        description="Program and arguments; a string is split like a shell command line."
    )
    output: str = Field(
        description="File (relative to the simulation directory) whose last data line is the result."
    )
    stdin: str | None = Field(
        default=None, description="Text sent to the program's standard input, e.g. a group selection."
    )
    columns: list[str] = Field(
        default_factory=list,
        description="Names of the result columns (default: taken from the xvg legends).",
    )

    @field_validator("name")
    @classmethod
    def _check_name(cls, name: str) -> str:
        if not SAFE_NAME.match(name):
            raise ValueError(f"analysis name {name!r} may only contain letters, digits and ._+-")
        return name


class SimulationDefaults(_Model):
    """Per-simulation parameters; every field can be overridden by the matrix."""

    name: str = "{target}_in_{solvent}"
    n_target: int = Field(default=1, ge=1)
    n_solvent: int = Field(default=8000, ge=1)
    density: float = Field(default=900.0, gt=0, description="Target density in kg/m^3.")


class MatrixConfig(BaseModel):
    """Build matrix: every key other than include/exclude is an axis."""

    model_config = ConfigDict(extra="allow")

    include: list[dict[str, Scalar]] = Field(default_factory=list)
    exclude: list[dict[str, Scalar]] = Field(default_factory=list)

    @model_validator(mode="after")
    def _check_axes(self) -> "MatrixConfig":
        for axis, values in self.axes.items():
            if not isinstance(values, list) or not values:
                raise ValueError(f"matrix axis {axis!r} must be a non-empty list")
            if not all(isinstance(value, Scalar) for value in values):
                raise ValueError(f"matrix axis {axis!r} may only contain scalar values")
        if not self.axes and not self.include:
            raise ValueError("matrix needs at least one axis or include entry")
        return self

    @property
    def axes(self) -> dict[str, list[Any]]:
        return dict(self.model_extra or {})


class TopologyConfig(_Model):
    includes: list[Path] = Field(
        default_factory=lambda: [Path(name) for name in DEFAULT_INCLUDES],
        description="Force-field ITPs included in every topology.",
    )


def default_stages() -> list[StageConfig]:
    """Minimization, equilibration and AWH solvation free energy."""
    return [
        StageConfig(name="min", mdp="min.mdp"),
        StageConfig(name="eq", mdp="eq.mdp", maxwarn=2),
        StageConfig(name="awh", mdp="awh.mdp", substitute={"TOCOUPLE": "{target}"}),
    ]


class CampaignConfig(_Model):
    """A full simulation campaign."""

    title: str = "Unnamed campaign"
    output_dir: Path = Path("campaign")
    forcefield_dir: Path = DEFAULT_FORCEFIELD_DIR
    mdp_dir: Path = DEFAULT_MDP_DIR
    topology: TopologyConfig = Field(default_factory=TopologyConfig)
    molecules: dict[str, Path] = Field(
        default_factory=dict, description="Target molecule name -> ITP file."
    )
    molecule_library: Path | None = Field(
        default=None,
        description="Library from `martini_harness split-itp`, searched after `molecules`.",
    )
    matrix: MatrixConfig
    simulation: SimulationDefaults = Field(default_factory=SimulationDefaults)
    stages: list[StageConfig] = Field(default_factory=default_stages)
    stage_overrides: dict[str, dict[str, Any]] = Field(
        default_factory=dict,
        description="Stage name -> keys deep-merged into that stage, e.g. of an imported pipeline.",
    )
    mdrun: dict[str, Any] = Field(
        default_factory=dict,
        description=(
            "mdrun flags for every stage except energy minimizations, e.g. "
            "{bonded: gpu, ntomp: 8, ntmpi: 1}. A stage's own mdrun flags take precedence."
        ),
    )
    analyses: list[AnalysisConfig] = Field(
        default_factory=list, description="Analyses run by `martini_harness analyze`."
    )

    @model_validator(mode="after")
    def _resolve_paths(self, info: ValidationInfo) -> "CampaignConfig":
        base_dir = Path((info.context or {}).get("base_dir", Path.cwd()))
        self.output_dir = _resolve(self.output_dir, base_dir)
        self.forcefield_dir = _resolve(self.forcefield_dir, base_dir)
        self.mdp_dir = _resolve(self.mdp_dir, base_dir)
        if self.molecule_library is not None:
            self.molecule_library = _resolve(self.molecule_library, base_dir)
        self.topology.includes = [
            _resolve(path, self.forcefield_dir) for path in self.topology.includes
        ]
        self.molecules = {
            name: _resolve(path, base_dir) for name, path in self.molecules.items()
        }
        return self

    @model_validator(mode="after")
    def _apply_stage_overrides(self) -> "CampaignConfig":
        stages = {stage.name: index for index, stage in enumerate(self.stages)}
        for name, override in self.stage_overrides.items():
            if name not in stages:
                raise ValueError(
                    f"stage_overrides: no stage named {name!r} (stages: {', '.join(stages)})"
                )
            if "name" in override:
                raise ValueError(f"stage_overrides.{name}: a stage cannot be renamed")
            index = stages[name]
            merged = _deep_merge(self.stages[index].model_dump(exclude_unset=True), override)
            try:
                self.stages[index] = StageConfig.model_validate(merged)
            except ValidationError as exc:
                raise ValueError(
                    "; ".join(
                        ".".join(str(part) for part in ("stage_overrides", name, *error["loc"]))
                        + f": {error['msg']}"
                        for error in exc.errors()
                    )
                ) from exc
        return self

    @model_validator(mode="after")
    def _check_stages(self) -> "CampaignConfig":
        if not self.stages:
            raise ValueError("at least one stage is required")
        seen = [COORDINATES_STEP]
        for stage in self.stages:
            if stage.name in seen:
                raise ValueError(f"duplicate stage name {stage.name!r}")
            if stage.input is not None and stage.input not in seen:
                raise ValueError(
                    f"stage {stage.name!r}: input {stage.input!r} must be "
                    f"{COORDINATES_STEP!r} or an earlier stage"
                )
            seen.append(stage.name)
        return self

    @model_validator(mode="after")
    def _check_analyses(self) -> "CampaignConfig":
        stage_names = [stage.name for stage in self.stages]
        names = set()
        for analysis in self.analyses:
            if analysis.name in names:
                raise ValueError(f"duplicate analysis name {analysis.name!r}")
            names.add(analysis.name)
            if analysis.stage is None:
                analysis.stage = stage_names[-1]
            elif analysis.stage not in stage_names:
                raise ValueError(
                    f"analysis {analysis.name!r}: no stage named {analysis.stage!r} "
                    f"(stages: {', '.join(stage_names)})"
                )
        return self

    def mdp_path(self, stage: StageConfig) -> Path:
        return _resolve(Path(stage.mdp), self.mdp_dir)


def _resolve(path: Path, base_dir: Path) -> Path:
    path = Path(path).expanduser()
    return path if path.is_absolute() else (base_dir / path).resolve()


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """Merge ``override`` into a copy of ``base``; only mappings merge recursively."""
    merged = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def _read_mapping(path: Path) -> dict[str, Any]:
    try:
        data = YAML(typ="safe").load(path.read_text()) or {}
    except Exception as exc:
        raise ConfigError(f"cannot read {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise ConfigError(f"{path}: top level must be a mapping")
    return data


def _resolve_file_paths(data: dict[str, Any], base_dir: Path) -> None:
    """Make the file-relative paths in one file's data absolute, in place."""
    for key in FILE_RELATIVE_KEYS:
        if isinstance(data.get(key), str):
            data[key] = str(_resolve(Path(data[key]), base_dir))
    molecules = data.get("molecules")
    if isinstance(molecules, dict):
        for name, itp in molecules.items():
            if isinstance(itp, str):
                molecules[name] = str(_resolve(Path(itp), base_dir))


def _load_with_imports(path: Path, chain: tuple[Path, ...] = ()) -> dict[str, Any]:
    """Read ``path`` and merge its imports (recursively) underneath it.

    Args:
        path: Absolute path of the file to load.
        chain: Files currently being imported, outermost first; used to
            detect cycles.

    Raises:
        ConfigError: On unreadable or missing files, malformed ``imports``
            or an import cycle.
    """
    data = _read_mapping(path)
    imports = data.pop(IMPORTS_KEY, None) or []
    if isinstance(imports, str):
        imports = [imports]
    if not isinstance(imports, list) or not all(isinstance(entry, str) for entry in imports):
        raise ConfigError(f"{path}: {IMPORTS_KEY} must be a file name or a list of file names")
    _resolve_file_paths(data, path.parent)

    merged: dict[str, Any] = {}
    for entry in imports:
        imported = _resolve(Path(entry), path.parent)
        if imported in (*chain, path):
            cycle = " -> ".join(str(file) for file in (*chain, path, imported))
            raise ConfigError(f"import cycle: {cycle}")
        if not imported.is_file():
            raise ConfigError(f"{path}: imported file not found: {entry} ({imported})")
        merged = _deep_merge(merged, _load_with_imports(imported, (*chain, path)))
    return _deep_merge(merged, data)


def load_campaign(config_file: str | Path) -> CampaignConfig:
    """Load and validate a campaign YAML file, including its imports.

    Raises:
        ConfigError: If a file cannot be read or parsed, an import is
            missing or cyclic, or the merged campaign fails validation.
    """
    config_file = Path(config_file)
    path = config_file.expanduser().resolve()
    if not path.is_file():
        raise ConfigError(f"cannot read {config_file}: file not found")
    data = _load_with_imports(path)
    try:
        return CampaignConfig.model_validate(data, context={"base_dir": path.parent})
    except ValueError as exc:
        raise ConfigError(f"{config_file}: {exc}") from exc
