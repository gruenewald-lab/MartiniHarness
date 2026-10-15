"""Configuration handling for simulation campaigns."""

from pathlib import Path
from typing import Annotated, Any

import pydantic
from pydantic import BaseModel, Field, ConfigDict, model_validator
from ruamel.yaml import YAML


class SimulationStep(BaseModel):
    """A single simulation step configuration.

    Attributes:
        name: Name of the simulation step (e.g., "minimization", "equilibration").
        mdp_file: Path to the MDP configuration file.
        input_structure: Input structure file pattern or specific file.
        output_prefix: Prefix for output files (without extension).
        maxwarn: Maximum number of GROMACS warnings to ignore.
        extra_args: Extra arguments to pass to grompp/mdrun.
    """

    name: str
    mdp_file: str
    input_structure: str | None = None  # Will default based on step order
    output_prefix: str | None = None
    maxwarn: int = Field(default=4, ge=0)
    extra_args: dict[str, Any] = Field(default_factory=dict)

    @model_validator(mode="after")
    def set_defaults(self) -> "SimulationStep":
        """Set default values if not provided."""
        if self.output_prefix is None:
            self.output_prefix = self.name
        return self


class SimulationConfig(BaseModel):
    """Configuration for a single simulation run.

    Attributes:
        name: Name of the simulation (e.g., "target_in_solvent").
        target_molecule: Name or identifier of the target molecule.
        solvent: Type of solvent to use.
        n_target: Number of target molecules.
        n_solvent: Number of solvent molecules.
        density: Target density in kg/m³.
        pforce: Positional restraint force constant (kJ/mol/nm²).
        mdp_files: List of MDP files for each step.
    """

    name: str
    target_molecule: str
    solvent: str
    n_target: int = Field(default=1, ge=1)
    n_solvent: int = Field(default=8000, ge=1)
    density: float = Field(default=900.0, gt=0)
    pforce: float = Field(default=800_000.0, ge=0)
    mdp_files: list[str] = []

    model_config = pydantic.ConfigDict(extra="ignore")


class CampaignConfig(BaseModel):
    """Configuration for a full simulation campaign.

    Attributes:
        title: Title/description of the campaign.
        output_dir: Directory for all simulation outputs.
        simulations: List of individual simulations to run.
        gromacs_path: Path to GROMACS force field directory.
        mdp_directory: Directory containing MDP files.
    """

    model_config = ConfigDict(
        env_prefix="MARTINI_",
        extra="ignore",
        populate_by_name=True,
    )

    title: str = "Unnamed Campaign"
    output_base: Path = Field(default=Path("./campaigns"), description="Base output directory")
    simulations: list[SimulationConfig] = []

    # Paths - will be resolved relative to package or absolute
    gromacs_path: str | None = None
    mdp_directory: str | None = None

    @model_validator(mode="after")
    def resolve_paths(self) -> "CampaignConfig":
        """Resolve relative paths to absolute paths."""
        if self.gromacs_path and not Path(self.gromacs_path).is_absolute():
            # Resolve relative to current working directory
            self.gromacs_path = str(Path(self.gromacs_path).resolve())

        if self.mdp_directory and not Path(self.mdp_directory).is_absolute():
            self.mdp_directory = str(Path(self.mdp_directory).resolve())

        return self

    def get_gmx_files_dir(self) -> Path:
        """Get the GROMACS force field files directory."""
        if self.gromacs_path:
            return Path(self.gromacs_path)
        # Default to package data location
        from martini_harness import get_package_data_path
        return get_package_data_path("gmx_files")

    def get_mdp_dir(self) -> Path:
        """Get the MDP files directory."""
        if self.mdp_directory:
            return Path(self.mdp_directory)
        # Default to package data location
        from martini_harness import get_package_data_path
        return get_package_data_path("mdp_files")


def load_yaml_config(config_file: str | Path) -> CampaignConfig:
    """Load a YAML configuration file and parse it into a CampaignConfig.

    Args:
        config_file: Path to the YAML configuration file.

    Returns:
        Parsed CampaignConfig instance.
    """
    yaml = YAML()
    yaml.preserve_quotes = True

    with open(config_file) as f:
        data = yaml.load(f) or {}

    # Parse simulations separately to handle potential list format
    if "simulations" in data and isinstance(data["simulations"], list):
        sim_configs = []
        for sim_data in data["simulations"]:
            if isinstance(sim_data, dict):
                sim_configs.append(SimulationConfig(**sim_data))
        data["simulations"] = sim_configs

    return CampaignConfig(**data)


def dump_config_to_yaml(config: CampaignConfig, output_file: str | Path) -> None:
    """Dump a CampaignConfig to a YAML file.

    Args:
        config: The configuration to dump.
        output_file: Path to the output YAML file.
    """
    import json

    yaml = YAML()
    yaml.preserve_quotes = True

    with open(output_file, "w") as f:
        # Convert to dict and write
        config_dict = config.model_dump(exclude_none=True)

        # Handle Path objects for JSON serialization
        def serialize_value(val):
            if isinstance(val, Path):
                return str(val)
            elif hasattr(val, "__dict__"):
                return val.__dict__
            else:
                return val

        config_dict = json.loads(json.dumps(config_dict, default=lambda o: str(o)))

        yaml.dump(config_dict, f)
