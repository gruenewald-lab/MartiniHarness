"""Campaign configuration: YAML schema and loading.

Relative paths in a campaign file are resolved against the directory of that
file, except ``stages[].mdp`` (relative to ``mdp_dir``) and
``topology.includes`` (relative to ``forcefield_dir``).
"""

import re
from pathlib import Path
from typing import Any

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
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
    matrix: MatrixConfig
    simulation: SimulationDefaults = Field(default_factory=SimulationDefaults)
    stages: list[StageConfig] = Field(default_factory=default_stages)

    @model_validator(mode="after")
    def _resolve_paths(self, info: ValidationInfo) -> "CampaignConfig":
        base_dir = Path((info.context or {}).get("base_dir", Path.cwd()))
        self.output_dir = _resolve(self.output_dir, base_dir)
        self.forcefield_dir = _resolve(self.forcefield_dir, base_dir)
        self.mdp_dir = _resolve(self.mdp_dir, base_dir)
        self.topology.includes = [
            _resolve(path, self.forcefield_dir) for path in self.topology.includes
        ]
        self.molecules = {
            name: _resolve(path, base_dir) for name, path in self.molecules.items()
        }
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

    def mdp_path(self, stage: StageConfig) -> Path:
        return _resolve(Path(stage.mdp), self.mdp_dir)


def _resolve(path: Path, base_dir: Path) -> Path:
    path = Path(path).expanduser()
    return path if path.is_absolute() else (base_dir / path).resolve()


def load_campaign(config_file: str | Path) -> CampaignConfig:
    """Load and validate a campaign YAML file.

    Raises:
        ConfigError: If the file cannot be parsed or fails validation.
    """
    config_file = Path(config_file)
    try:
        data = YAML(typ="safe").load(config_file.read_text()) or {}
    except Exception as exc:
        raise ConfigError(f"cannot read {config_file}: {exc}") from exc
    if not isinstance(data, dict):
        raise ConfigError(f"{config_file}: top level must be a mapping")
    try:
        return CampaignConfig.model_validate(
            data, context={"base_dir": config_file.resolve().parent}
        )
    except ValueError as exc:
        raise ConfigError(f"{config_file}: {exc}") from exc
