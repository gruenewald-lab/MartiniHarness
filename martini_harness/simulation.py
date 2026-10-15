"""Simulation workflow management."""

import logging
import os
from contextlib import contextmanager
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any

import gromacs as gmx
import numpy as np
from polyply.src.topology import Topology
from polyply import gen_coords
from vermouth.forcefield import ForceField

from .config import CampaignConfig, SimulationConfig

logger = logging.getLogger(__name__)


# Type alias for stage type classification
StageType = str  # Can be one of: 'SETUP', 'COORDINATE_GENERATION', 'MINIMIZATION', 'EQUILIBRATION', 'AWH', 'PRODUCTION'

STAGE_TYPE_MAPPING = {
    "SETUP": "SETUP",
    "COORDINATE_GENERATION": "COORDINATE_GENERATION", 
    "MINIMIZATION": "MINIMIZATION",
    "EQUILIBRATION": "EQUILIBRATION",
    "AWH": "AWH",
    "PRODUCTION": "PRODUCTION",  # General production stage (can be AWH or MD)
}

# Mapping from MDP filename to default stage type
MDP_TO_STAGE_TYPE = {
    "setup.mdp": "SETUP",
    "coordinate_generation.mdp": "COORDINATE_GENERATION",
    "min.mdp": "MINIMIZATION",
    "eq.mdp": "EQUILIBRATION",
    "awh.mdp": "AWH",
    "prod.mdp": "PRODUCTION",
}

# Expected output files per stage type (relative paths from workdir)
EXPECTED_FILES_PER_STAGE_TYPE = {
    "SETUP": ["system.top"],
    "COORDINATE_GENERATION": ["init.gro", "coords.gro"],
    "MINIMIZATION": ["min.gro", "min.log", "min.cpt"],
    "EQUILIBRATION": ["eq.gro", "eq.log", "eq.cpt"],
    "AWH": ["awh.xvg", "awh.tpr"],
    "PRODUCTION": ["prod.gro", "prod.log", "prod.cpt"],  # XPM is optional, xvg files common too
}

# Stage types that are checkpointable (can be restarted)
CHECKPOINTABLE_STAGE_TYPES = {"MINIMIZATION", "EQUILIBRATION", "AWH", "PRODUCTION"}


def get_stage_type_from_mdp_filename(mdp_filename: str) -> str:
    """Determine stage type from MDP filename.
    
    Args:
        mdp_filename: Name of the MDP file (e.g., 'min.mdp', 'prod.mdp')
        
    Returns:
        Stage type string (SETUP, COORDINATE_GENERATION, MINIMIZATION, etc.)
    """
    # Extract base name without .mdp extension
    base_name = mdp_filename.lower()
    if not base_name.endswith(".mdp"):
        base_name += ".mdp"
    
    return MDP_TO_STAGE_TYPE.get(base_name, "PRODUCTION")


def get_expected_output_files(stage_type: str) -> list[str]:
    """Get expected output files for a given stage type.
    
    Args:
        stage_type: The stage type (SETUP, MINIMIZATION, EQUILIBRATION, etc.)
        
    Returns:
        List of expected output file names (relative paths from workdir)
    """
    return EXPECTED_FILES_PER_STAGE_TYPE.get(stage_type, [])


def check_files_exist(workdir: Path, file_patterns: list[str]) -> list[Path]:
    """Check if files exist in a directory, supporting wildcard patterns.
    
    Args:
        workdir: Working directory to search in
        file_patterns: List of file patterns (e.g., ['prod.gro', '*.xvg'])
        
    Returns:
        List of Path objects for found files
    """
    import glob as glob_module
    
    found_files = []
    for pattern in file_patterns:
        if "*" in pattern or "?" in pattern:
            # Handle wildcard patterns using glob
            full_pattern = workdir / pattern
            matches = list(workdir.glob(pattern))
            found_files.extend(matches)
        else:
            # Exact file name
            file_path = workdir / pattern
            if file_path.exists():
                found_files.append(file_path)
    
    return found_files


def check_stage_complete_by_type(workdir: Path, stage_config: dict) -> bool:
    """Check if a simulation stage is complete by examining expected output files.
    
    Supports wildcard patterns for PRODUCTION stages which may generate multiple xvg files.
    
    Args:
        workdir: Working directory path (absolute).
        stage_config: Dictionary containing stage configuration with keys:
            - stage_type: One of SETUP, COORDINATE_GENERATION, MINIMIZATION, EQUILIBRATION, AWH, PRODUCTION
            - output_prefix: Optional prefix for stage-specific files (e.g., "prod" for prod.gro)
            
    Returns:
        True if the stage appears to be complete, False otherwise.
    """
    import os
    
    stage_type = stage_config.get("stage_type", "PRODUCTION")
    output_prefix = stage_config.get("output_prefix", "")
    
    # Get expected files for this stage type
    expected_files = get_expected_output_files(stage_type)
    
    if not expected_files:
        logger.warning(f"Unknown stage type {stage_type}, cannot verify completion")
        return False
    
    # Adjust file names based on output_prefix if provided
    adjusted_files = []
    for fname in expected_files:
        # Replace common prefixes with custom prefix where applicable
        if output_prefix and ("prod" in fname or "awh" in fname or "eq" in fname or "min" in fname):
            # Replace the stage-specific prefix with provided one
            import re
            adjusted = re.sub(r'^(prod|awh|eq|min)\.', f'{output_prefix}.', fname)
            adjusted_files.append(adjusted)
        else:
            adjusted_files.append(fname)
    
    # Check if all expected files exist (with wildcard support for PRODUCTION stages)
    found_files = check_files_exist(workdir, adjusted_files)
    
    # For checkpointable stages, we need most critical files; for setup/coord_gen, strict check
    is_checkpointable = stage_type in CHECKPOINTABLE_STAGE_TYPES
    
    if len(found_files) == len(adjusted_files):
        return True  # All files found
    elif is_checkpointable and len(found_files) >= int(len(adjusted_files) * 0.75):
        # For checkpointable stages, partial completion counts (need ~75%)
        logger.debug(f"Stage {stage_type}: Found {len(found_files)}/{len(adjusted_files)} expected files")
        return True
    
    logger.debug(f"Stage {stage_type}: Only found {len(found_files)}/{len(adjusted_files)} expected files")
    return False


class SimulationStage(Enum):
    """Enumeration of simulation stages (deprecated, kept for backward compatibility)."""

    SETUP = "setup"
    COORDINATE_GENERATION = "coordinate_generation"
    MINIMIZATION = "minimization"
    EQUILIBRATION = "equilibration"
    AWH = "awh"
    MD_PRODUCTION = "md_production"


@dataclass
class SimulationResult:
    """Result of a simulation stage.

    Attributes:
        stage: The simulation stage that was completed.
        success: Whether the stage completed successfully.
        message: Human-readable status message.
        output_files: List of output files generated.
        errors: Any errors encountered during execution.
    """

    stage: SimulationStage
    success: bool = False
    message: str = ""
    output_files: list[Path] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)


class CoordinateGenerator:
    """Generate simulation coordinates using Polyply.

    Attributes:
        top_path: Path to the topology file.
    """

    def __init__(self, top_path: str | Path):
        self.top_path = Path(top_path)

    def generate(
        self,
        output_name: str = "init.gro",
        density: float = 900.0,
        name: str = "coords",
    ) -> Path:
        """Generate initial coordinates.

        Args:
            output_name: Name of the output .gro file (relative to workdir).
            density: Target density in kg/m³.
            name: Name for the generated structure.

        Returns:
            Path object pointing to the generated coordinate file (relative to workdir).
        """
        logger.info("Generating coordinates with Polyply (density=%.1f kg/m³)", density)

        gen_coords(
            toppath=str(self.top_path), outpath=output_name, density=density, name=name
        )
        return Path(output_name)


def write_topology(
    solvent: str,
    target: str,
    n_target: int,
    n_solvent: int,
    forcefield_itp: Path,
    solvent_itp: Path | None = None,
    molecule_itp_dir: Path | None = None,
    output_file: str | Path | None = None,
) -> Path:
    """Write a simulation topology file.

    Args:
        solvent: Solvent name/identifier.
        target: Target molecule name/identifier.
        n_target: Number of target molecules.
        n_solvent: Number of solvent molecules.
        forcefield_itp: Path to the force field ITP file.
        solvent_itp: Optional path to solvent-specific ITP file.
        molecule_itp_dir: Path to molecule definition ITP directory/file.
        output_file: Output topology file path (default: system.top).

    Returns:
        Path to the written topology file.
    """
    from martini_harness import get_package_data_path

    lines = [f'#include "{forcefield_itp}"']

    if solvent_itp:
        lines.append(f'#include "{solvent_itp}"')

    # Get molecule directory or use default location
    mol_dir = molecule_itp_dir
    if not mol_dir:
        from martini_harness import get_package_data_path
        mol_dir = get_package_data_path("row42.itp")

    lines.append(f'#include "{mol_dir}"')

    # Add system and molecules sections
    lines += [
        "",
        "[ system ]",
        f"{target}_in_{solvent}",
        "",
        "[ molecules ]",
        f"{target}  {n_target}",
        f"{solvent}  {n_solvent}",
    ]

    output_file = Path(output_file or "system.top")
    output_file.write_text("\n".join(lines))
    return output_file


class SimulationRunner:
    """Run GROMACS simulation steps.

    Attributes:
        workdir: Working directory for this simulation.
        forcefield_itp: Force field ITP file path.
        solvent_itp: Solvent ITP file path (optional).
        molecule_itp_dir: Molecule definition directory/file.
    """

    def __init__(
        self,
        workdir: str | Path,
        forcefield_itp: Path,
        solvent_itp: Path | None = None,
        molecule_itp_dir: Path | None = None,
    ):
        self.workdir = Path(workdir)
        self.forcefield_itp = forcefield_itp
        self.solvent_itp = solvent_itp
        self.molecule_itp_dir = molecule_itp_dir

        # Ensure work directory exists
        self.workdir.mkdir(parents=True, exist_ok=True)

    @contextmanager
    def _working_directory(self):
        """Context manager to temporarily change into the simulation workdir."""
        original_cwd = os.getcwd()
        try:
            os.chdir(self.workdir)
            yield
        finally:
            os.chdir(original_cwd)

    def run_grompp(
        self,
        mdp_file: str | Path,
        structure: str | Path,
        topology: str | Path,
        output_prefix: str | Path,
        maxwarn: int = 4,
        extra_args: dict[str, Any] | None = None,
    ) -> Path:
        """Run GROMACS grompp to create a TPR file.

        Args:
            mdp_file: Path to the MDP parameter file.
            structure: Input structure file (e.g., .gro).
            topology: Topology file (.top).
            output_prefix: Output prefix for the TPR file.
            maxwarn: Maximum number of GROMACS warnings to ignore.
            extra_args: Extra arguments for grompp.

        Returns:
            Path to the generated TPR file.
        """
        mdp_file = Path(mdp_file)
        structure = Path(structure)
        topology = Path(topology)
        output_prefix = Path(output_prefix)

        logger.info("Running grompp with MDP: %s", mdp_file.name)

        args = {
            "f": mdp_file,
            "c": structure,
            "p": topology,
            "o": f"{output_prefix}.tpr",
            "maxwarn": maxwarn,
            **(extra_args or {}),
        }

        gmx.grompp(**args)

        return Path(f"{output_prefix}.tpr")

    def run_mdrun(
        self,
        tpr_file: str | Path,
        output_prefix: str | Path,
        extra_args: dict[str, Any] | None = None,
    ) -> list[Path]:
        """Run GROMACS mdrun to execute a simulation.

        Args:
            tpr_file: Input TPR file.
            output_prefix: Output prefix (without extension).
            extra_args: Extra arguments for mdrun.

        Returns:
            List of generated output files.
        """
        tpr_file = Path(tpr_file)
        output_prefix = Path(output_prefix)

        logger.info("Running mdrun with TPR: %s", tpr_file.name)

        args = {"s": tpr_file, "deffnm": output_prefix}

        if extra_args:
            args.update(extra_args)

        # Add verbose mode by default
        args.setdefault("v", True)

        gmx.mdrun(**args)

        # Return common output files (relative paths from workdir)
        output_files = [
            Path(f"{output_prefix}.gro"),
            Path(f"{output_prefix}.log"),
            Path(f"{output_prefix}.ener"),
            Path(f"{output_prefix}.trr"),
        ]
        return [f for f in output_files if f.exists()]

    def run_stage(
        self,
        stage: SimulationStage,
        mdp_file: str | Path,
        structure: str | Path,
        topology: str | Path,
        output_prefix: str | None = None,
        maxwarn: int = 4,
        extra_args: dict[str, Any] | None = None,
    ) -> SimulationResult:
        """Run a complete simulation stage (grompp + mdrun).

        Args:
            stage: The simulation stage enum.
            mdp_file: Path to MDP file.
            structure: Input structure file.
            topology: Topology file.
            output_prefix: Output prefix for files.
            maxwarn: Maximum GROMACS warnings to ignore.
            extra_args: Extra arguments for both grompp and mdrun.

        Returns:
            SimulationResult with status information.
        """
        output_prefix = output_prefix or stage.value

        try:
            # Run grompp
            tpr_file = self.run_grompp(
                mdp_file=mdp_file,
                structure=structure,
                topology=topology,
                output_prefix=output_prefix,
                maxwarn=maxwarn,
                extra_args=extra_args,
            )

            # Run mdrun
            outputs = self.run_mdrun(tpr_file=tpr_file, output_prefix=output_prefix)

            logger.info(
                "Stage %s completed successfully. Output files: %s",
                stage.value,
                [f.name for f in outputs],
            )

            return SimulationResult(
                stage=stage, success=True, message="Completed successfully"
            )

        except Exception as e:
            error_msg = f"{stage.value} failed: {e}"
            logger.warning(error_msg)
            return SimulationResult(stage=stage, success=False, message=str(e))


class SimulationCampaign:
    """Manage a complete simulation campaign with multiple configurations.

    Attributes:
        config: The campaign configuration.
        campaign_dir: Base directory for all outputs.
    """

    def __init__(self, config: CampaignConfig):
        self.config = config
        self.campaign_dir = Path(config.output_base) / config.title.replace(" ", "_")

        # Ensure campaign directory exists
        self.campaign_dir.mkdir(parents=True, exist_ok=True)

        logger.info(
            "Campaign '%s' initialized. Output will be in: %s",
            self.config.title,
            self.campaign_dir,
        )

    def run_simulation(self, sim_config: SimulationConfig) -> list[SimulationResult]:
        """Run a single simulation with all stages.

        Args:
            sim_config: Configuration for the specific simulation.

        Returns:
            List of results from each simulation stage.
        """
        logger.info("Running simulation: %s", sim_config.name)

        # Create work directory for this simulation
        workdir = self.campaign_dir / sim_config.name
        workdir.mkdir(parents=True, exist_ok=True)

        # Initialize components
        gmx_files_dir = self.config.get_gmx_files_dir()
        mdp_dir = self.config.get_mdp_dir()

        runner = SimulationRunner(
            workdir=workdir,
            forcefield_itp=gmx_files_dir / "martini_v3.0.0.itp",
            solvent_itp=gmx_files_dir / "martini_v3.0.0_solvents_v1.itp",
            molecule_itp_dir=None,  # Will be handled separately
        )

        results: list[SimulationResult] = []

        try:
            # Write topology file
            top_file = write_topology(
                solvent=sim_config.solvent,
                target=sim_config.target_molecule,
                n_target=sim_config.n_target,
                n_solvent=sim_config.n_solvent,
                forcefield_itp=runner.forcefield_itp,
                solvent_itp=runner.solvent_itp,
            )

            # Coordinate generation stage
            logger.info("Generating initial coordinates")
            coords_file = CoordinateGenerator(str(top_file)).generate(
                output_path=workdir / "init.gro",
                density=sim_config.density,
                name="coords",
            )

            results.append(SimulationResult(stage=SimulationStage.COORDINATE_GENERATION))

            # Minimization
            min_mdp = mdp_dir / "min.mdp"

            results.append(
                runner.run_stage(
                    stage=SimulationStage.MINIMIZATION,
                    mdp_file=min_mdp,
                    structure=workdir / "init.gro",
                    topology=top_file,
                    output_prefix="min",
                    maxwarn=4,
                )
            )

            # Check if minimization was successful
            min_gro = workdir / "min.gro"
            if not min_gro.exists():
                logger.warning("min.gro not found, falling back to init.gro")
                min_gro = workdir / "init.gro"

            # Equilibration
            eq_mdp = mdp_dir / "eq.mdp"
            results.append(
                runner.run_stage(
                    stage=SimulationStage.EQUILIBRATION,
                    mdp_file=eq_mdp,
                    structure=min_gro,
                    topology=top_file,
                    output_prefix="eq",
                    maxwarn=4,
                    extra_args={"pforce": sim_config.pforce},
                )
            )

            # Check if equilibration was successful
            eq_gro = workdir / "eq.gro"
            if not eq_gro.exists():
                raise IOError("Equilibration failed: eq.gro not found")

            # AWH (Adaptive Biasing Force) simulation
            awh_mdp_src = mdp_dir / "awh.mdp"
            awh_mdp_dst = workdir / "awh.mdp"

            with open(awh_mdp_src, "r") as f:
                awh_content = f.read()

            # Replace TOCOUPLE placeholder with target molecule name
            awh_content = awh_content.replace("TOCOUPLE", sim_config.target_molecule)

            with open(awh_mdp_dst, "w") as f:
                f.write(awh_content)

            results.append(
                runner.run_stage(
                    stage=SimulationStage.AWH,
                    mdp_file=awh_mdp_dst,
                    structure=eq_gro,
                    topology=top_file,
                    output_prefix="awh",
                    maxwarn=4,
                    extra_args={"pforce": sim_config.pforce},
                )
            )

        except Exception as e:
            logger.error("Simulation failed for %s: %s", sim_config.name, e)
            raise

        finally:
            # Clean up temporary AWH MDP file if it exists
            awh_mdp_dst = workdir / "awh.mdp"
            if awh_mdp_dst.exists():
                try:
                    awh_mdp_dst.unlink()
                except OSError as e:
                    logger.warning("Could not remove %s: %s", awh_mdp_dst, e)

        return results

    def run_all(self) -> list[SimulationResult]:
        """Run all simulations in the campaign.

        Returns:
            List of results from all simulation runs.
        """
        logger.info("Starting campaign '%s' with %d simulations", self.config.title, len(self.config.simulations))

        all_results: list[SimulationResult] = []

        for sim_config in self.config.simulations:
            logger.info("\n%s Running simulation: %s", "=" * 50, sim_config.name)
            results = self.run_simulation(sim_config)
            all_results.extend(results)

            # Log summary
            success_count = sum(1 for r in results if r.success)
            total_count = len(results)
            logger.info("Completed %s/%d stages", success_count, total_count)

        return all_results


def run_simulation_workflow(
    target: str,
    solvent: str,
    n_target: int = 1,
    n_solvent: int = 8000,
    density: float = 900.0,
    pforce: float = 800_000.0,
    workdir: Path | None = None,
    gmx_files_dir: Path | None = None,
    mdp_directory: Path | None = None,
    stop_after: SimulationStage | None = None,
) -> list[SimulationResult]:
    """Run a complete simulation workflow for a target-solvent pair.

    This is the main entry point for running simulations from scripts or CLI.

    Args:
        target: Target molecule name/identifier.
        solvent: Solvent type (e.g., "OCO", "W", "CLF").
        n_target: Number of target molecules.
        n_solvent: Number of solvent molecules.
        density: Target density in kg/m³.
        pforce: Positional restraint force constant.
        workdir: Working directory (default: current directory / "{target}_in_{solvent}").
        gmx_files_dir: Directory containing GROMACS files.
        mdp_directory: Directory containing MDP files.
        stop_after: Optional stage to run up to (inclusive). If specified, the workflow
            will complete this stage and then return without running subsequent stages.

    Returns:
        List of SimulationResult objects for each stage.
    
    Checkpoints: Automatically checks for existing workdir and completed stages,
    skipping already finished steps to allow restart from previous state.
    
    Early Exit: If stop_after is specified (e.g., SimulationStage.AWH), the workflow will
        run through that stage and then return without running MD_PRODUCTION or later stages.
    """


def check_simulated_stage(workdir: Path, stage: SimulationStage) -> bool:
    """Check if a simulation stage has been completed by verifying expected output files.

    Args:
        workdir: Working directory path (absolute).
        stage: The simulation stage to check.

    Returns:
        True if the stage appears to be complete, False otherwise.
    """
    import os
    
    # Map stages to their expected output files (relative paths from workdir)
    expected_files = {
        SimulationStage.SETUP: ["system.top"],
        SimulationStage.COORDINATE_GENERATION: ["init.gro", "coords.gro"],
        SimulationStage.MINIMIZATION: ["min.gro", "min.log", "min.cpt"],
        SimulationStage.EQUILIBRATION: ["eq.gro", "eq.log", "eq.cpt"],
        SimulationStage.AWH: ["awh.xvg", "awh.tpr"],
        SimulationStage.MD_PRODUCTION: ["prod.gro", "prod.log", "prod.cpt"],  # XPM is optional
    }
    
    expected = expected_files.get(stage)
    if not expected:
        logger.warning(f"Unknown stage {stage.name}, cannot verify completion")
        return False
    
    # Check all required files exist
    missing = []
    for file_pattern in expected:
        if (workdir / file_pattern).exists():
            continue
        else:
            missing.append(file_pattern)
    
    if missing:
        logger.debug(f"Stage {stage.value} incomplete - missing: {', '.join(missing)}")
        return False
    
    logger.info(f"Stage {stage.value} appears to be complete")
    return True


def check_complete_simulation(workdir: Path) -> list[SimulationResult]:
    """Check which simulation stages are already complete.

    Args:
        workdir: Working directory path (absolute).

    Returns:
        List of SimulationResult objects for all completed stages.
    """
    import os
    
    # Check if workdir exists at all
    if not workdir.exists():
        logger.info(f"Working directory {workdir} does not exist, starting fresh")
        return []
    
    if not workdir.is_dir():
        raise ValueError(f"Expected {workdir} to be a directory, but it is not")
    
    # Check each stage and collect completed ones
    completed = []
    
    for stage in SimulationStage:
        if check_simulated_stage(workdir, stage):
            completed.append(SimulationResult(stage=stage, success=True))
        else:
            logger.info(f"Stage {stage.value} not yet complete")
            break  # Stop checking as soon as we find an incomplete stage
    
    return completed


def run_simulation_workflow(
    target: str,
    solvent: str,
    n_target: int = 1,
    n_solvent: int = 8000,
    density: float = 900.0,
    pforce: float = 800_000.0,
    workdir: Path | None = None,
    gmx_files_dir: Path | None = None,
    mdp_directory: Path | None = None,
) -> list[SimulationResult]:
    # Import here to avoid circular imports
    from martini_harness import get_package_data_path

    workdir = workdir or Path(f"{target}_in_{solvent}")
    
    # Ensure workdir is absolute so we can use it in _working_directory context manager
    if not workdir.is_absolute():
        workdir = Path(".").resolve() / workdir
    
    logger.info(f"Working directory: {workdir}")
    
    # Check for existing simulation state and collect completed stages
    completed_stages = check_complete_simulation(workdir)
    
    if completed_stages:
        logger.info(f"Found {len(completed_stages)} completed stage(s)")
        all_done = any(stage == SimulationStage.MD_PRODUCTION for stage in [r.stage for r in completed_stages])
        
        if all_done:
            # All stages complete - return existing results without re-running
            logger.info("All simulation stages already completed. Returning existing results.")
            
            # Create fresh result objects pointing to the workdir
            final_results = []
            for result in completed_stages:
                final_results.append(SimulationResult(stage=result.stage, success=True))
            
            return final_results
        else:
            # Partial completion - will continue from where it left off
            incomplete_stage = next((stage for stage in SimulationStage if not check_simulated_stage(workdir, stage)), None)
            if incomplete_stage:
                logger.info(f"Resuming from incomplete stage: {incomplete_stage.value}")
    
    gmx_files_dir = gmx_files_dir or get_package_data_path("gmx_files")
    mdp_directory = mdp_directory or get_package_data_path("mdp_files")

    # Create campaign config for this single simulation - use empty string as title 
    # to avoid nested directory structure, workdir itself becomes the campaign dir
    sim_config = SimulationConfig(
        name=f"{target}_in_{solvent}",
        target_molecule=target,
        solvent=solvent,
        n_target=n_target,
        n_solvent=n_solvent,
        density=density,
        pforce=pforce,
    )

    campaign_config = CampaignConfig(
        title="",  # Empty to avoid nested directory structure - we use workdir directly
        output_base=Path("."),
        simulations=[sim_config],
        gromacs_path=str(gmx_files_dir),
        mdp_directory=str(mdp_directory),
    )

    campaign = SimulationCampaign(campaign_config)
    
    # Run simulation within the workdir context manager to change CWD and reset afterwards
    from martini_harness.simulation import SimulationRunner
    
    runner = SimulationRunner(
        workdir=workdir,
        forcefield_itp=gmx_files_dir / "martini_v3.0.0.itp",
        solvent_itp=gmx_files_dir / "martini_v3.0.0_solvents_v1.itp",
        molecule_itp_dir=None,  # Will be handled separately in write_topology
    )
    
    results: list[SimulationResult] = []

    try:
        with runner._working_directory():
            # All paths here should use relative paths since we're now in workdir
            mdp_dir = Path(mdp_directory)
            
            # Check if SETUP stage is already complete before writing topology
            if check_simulated_stage(workdir, SimulationStage.SETUP):
                logger.info("SETUP stage already complete, skipping")
                # Find existing .top file in workdir
                top_file = next(Path(".").glob("*.top"), None)
                if not top_file:
                    raise FileNotFoundError("No topology file found for SETUP stage")
            else:
                # Write topology file (uses absolute paths for library files, outputs to local dir)
                top_file = write_topology(
                    solvent=sim_config.solvent,
                    target=sim_config.target_molecule,
                    n_target=sim_config.n_target,
                    n_solvent=sim_config.n_solvent,
                    forcefield_itp=runner.forcefield_itp,
                    solvent_itp=runner.solvent_itp,
                )

            # Check if coordinate generation is complete before running it
            if not check_simulated_stage(workdir, SimulationStage.COORDINATE_GENERATION):
                logger.info("Generating initial coordinates")
                coords_gen = CoordinateGenerator(top_file)
                coords_file = coords_gen.generate(
                    output_name="init.gro",
                    density=sim_config.density,
                    name="coords",
                )

            results.append(SimulationResult(stage=SimulationStage.COORDINATE_GENERATION))

            # Minimization - skip if already complete
            min_mdp = mdp_dir / "min.mdp"
            
            if check_simulated_stage(workdir, SimulationStage.MINIMIZATION):
                logger.info("MINIMIZATION stage already complete, skipping")
            else:
                results.append(
                    runner.run_stage(
                        stage=SimulationStage.MINIMIZATION,
                        mdp_file=min_mdp,
                        structure=Path("init.gro"),  # Relative path from workdir
                        topology=top_file,
                        output_prefix="min",
                        maxwarn=4,
                    )
                )

            # Check if minimization was successful
            min_gro = Path("min.gro")
            if not min_gro.exists():
                logger.warning("min.gro not found, falling back to init.gro")
                min_gro = Path("init.gro")

            # Ensure init.gro exists for subsequent stages
            init_gro = Path("init.gro")
            if not init_gro.exists():
                raise FileNotFoundError("Coordinate file (init.gro) not found after COORDINATE_GENERATION stage")

            # Equilibration - skip if already complete
            eq_mdp = mdp_dir / "eq.mdp"
            
            if check_simulated_stage(workdir, SimulationStage.EQUILIBRATION):
                logger.info("EQUILIBRATION stage already complete, skipping")
            else:
                results.append(
                    runner.run_stage(
                        stage=SimulationStage.EQUILIBRATION,
                        mdp_file=eq_mdp,
                        structure=min_gro,
                        topology=top_file,
                        output_prefix="eq",
                        maxwarn=4,
                    )
                )

            # Check if equilibration was successful
            eq_gro = Path("eq.gro")
            if not eq_gro.exists():
                logger.warning("eq.gro not found, falling back to min.gro")
                eq_gro = Path("min.gro")

            # AWH - skip if already complete
            awh_mdp = mdp_dir / "awh.mdp"
            
            if check_simulated_stage(workdir, SimulationStage.AWH):
                logger.info("AWH stage already complete, skipping")
            else:
                awh_mdp = mdp_dir / "awh.mdp"
                results.append(
                    runner.run_stage(
                        stage=SimulationStage.AWH,
                        mdp_file=awh_mdp,
                        structure=eq_gro,
                        topology=top_file,
                        output_prefix="awh",
                        maxwarn=4,
                    )
                )

            # MD_PRODUCTION - skip if already complete
            prod_mdp = mdp_dir / "prod.mdp"
            
            if check_simulated_stage(workdir, SimulationStage.MD_PRODUCTION):
                logger.info("MD_PRODUCTION stage already complete, skipping")
            else:
                results.append(
                    runner.run_stage(
                        stage=SimulationStage.MD_PRODUCTION,
                        mdp_file=prod_mdp,
                        structure=eq_gro,
                        topology=top_file,
                        output_prefix="prod",
                        maxwarn=4,
                    )
                )

    finally:
        # Ensure we're back in original directory even if an error occurs
        pass  # _working_directory already handles this via context manager

    return results


if __name__ == "__main__":
    # Example usage from script (similar to original run_simulation.py)
    import sys

    if len(sys.argv) < 3:
        print("Usage: python -m martini_harness.simulation <target> <solvent> [n_target] [n_solvent] [--stop-after STAGE]")
        print("")
        print("Arguments:")
        print("  target       Target molecule name/identifier (e.g., 'row42')")
        print("  solvent      Solvent type (e.g., 'OCO', 'W', 'CLF')")
        print("  n_target     Number of target molecules (default: 1)")
        print("  n_solvent    Number of solvent molecules (default: 8000)")
        print("")
        print("Optional:")
        print("  --stop-after STAGE   Stop after completing this stage (e.g., 'SETUP', 'COORDINATE_GENERATION',")
        print("                       'MINIMIZATION', 'EQUILIBRATION', 'AWH'). Defaults to running all stages.")
        sys.exit(1)

    target = sys.argv[1]
    solvent = sys.argv[2]
    n_target = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    n_solvent = int(sys.argv[4]) if len(sys.argv) > 4 else 8000

    # Parse optional --stop-after argument
    stop_after = None
    i = 5
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == "--stop-after":
            if i + 1 >= len(sys.argv):
                print("Error: --stop-after requires a stage name argument")
                sys.exit(1)
            stop_after_str = sys.argv[i + 1]
            try:
                stop_after = SimulationStage[stop_after_str.upper()]
            except KeyError:
                valid_stages = [stage.name for stage in SimulationStage]
                print(f"Error: Invalid stage '{stop_after_str}'. Valid stages: {', '.join(valid_stages)}")
                sys.exit(1)
            i += 2
        else:
            # Unknown argument, skip it
            i += 1

    results = run_simulation_workflow(
        target=target, 
        solvent=solvent, 
        n_target=n_target, 
        n_solvent=n_solvent,
        stop_after=stop_after
    )

