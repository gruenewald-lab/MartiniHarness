"""Command-line interface for Martini Harness."""

import argparse
import logging
from pathlib import Path

from .config import CampaignConfig, load_yaml_config
from .simulation import SimulationCampaign


def setup_logging(level: int = logging.INFO) -> None:
    """Setup logging configuration.

    Args:
        level: Logging level (e.g., logging.DEBUG, logging.INFO).
    """
    handler = logging.StreamHandler()
    handler.setLevel(level)
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    handler.setFormatter(formatter)

    logger = logging.getLogger("martini_harness")
    logger.addHandler(handler)
    logger.setLevel(level)


def cmd_run(args: argparse.Namespace) -> int:
    """Run a simulation campaign from YAML configuration.

    Args:
        args: Parsed command line arguments.

    Returns:
        Exit code (0 for success, non-zero for failure).
    """
    setup_logging(logging.DEBUG if args.verbose else logging.INFO)

    logger = logging.getLogger("martini_harness")

    # Load configuration
    try:
        config_path = Path(args.config)
        if not config_path.exists():
            logger.error(f"Configuration file not found: {config_path}")
            return 1

        config = load_yaml_config(config_path)
        logger.info("Loaded configuration from %s", args.config)

        # Override output directory if specified
        if args.output_dir:
            config.output_base = Path(args.output_dir).resolve()

        # Create and run campaign
        campaign = SimulationCampaign(config)

        try:
            results = campaign.run_all()

            # Summary
            total_stages = len(results)
            successful_stages = sum(1 for r in results if r.success)

            logger.info(
                "Campaign completed: %d/%d stages successful",
                successful_stages,
                total_stages,
            )

            return 0 if successful_stages == total_stages else 1

        except Exception as e:
            logger.error("Campaign execution failed: %s", e)
            return 1

    except Exception as e:
        logger.error("Failed to load configuration: %s", e)
        import traceback
        traceback.print_exc()
        return 1


def cmd_list_molecules(args: argparse.Namespace) -> int:
    """List available molecules from the force field.

    Args:
        args: Parsed command line arguments.

    Returns:
        Exit code (0 for success).
    """
    setup_logging(logging.WARNING)

    try:
        import sys

        # Read row42.itp or similar molecule definition file
        from martini_harness import get_package_data_path

        mol_file = get_package_data_path("row42.itp")

        if not mol_file.exists():
            print(f"Molecule definition file not found: {mol_file}")
            return 1

        # Parse the ITP file to extract molecule names
        try:
            from polyply.src.topology import Topology
            from vermouth.forcefield import ForceField

            ff = ForceField("test")
            top = Topology.from_gmx_topfile(str(mol_file), name="test")
            mol_names = list(top.force_field.blocks.keys())

            print("Available molecules:")
            for mol in sorted(mol_names):
                print(f"  {mol}")

            return 0

        except ImportError:
            # Fallback: just read the file and try to extract names manually
            print(f"Molecule definitions in {mol_file}:")
            with open(mol_file, "r") as f:
                content = f.read()

            # Try to find [ molecules ] section
            if "[ moleculetype ]" in content or "[ molecules ]" in content:
                print("  (Could not parse ITP file without polyply/vermouth)")
                return 1

        return 0

    except Exception as e:
        logging.getLogger("martini_harness").error(f"Error listing molecules: {e}")
        return 1


def main() -> int:
    """Main entry point for the CLI.

    Returns:
        Exit code.
    """
    parser = argparse.ArgumentParser(
        prog="martini-harness",
        description="A framework for running MD simulation campaigns with GROMACS and Polyply.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Run a campaign from YAML configuration:
    martini-harness run --config campaign.yaml

  Run a single simulation directly:
    martini-harness run-single target W

  List available molecules:
    martini-harness list-molecules

  Create a new campaign configuration (template):
    martini-harness new-config my_campaign
        """,
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Run command
    run_parser = subparsers.add_parser("run", help="Run a simulation campaign from YAML config")
    run_parser.add_argument(
        "config", type=Path, nargs="?", help="YAML configuration file path"
    )
    run_parser.add_argument(
        "--output-dir", "-o", type=Path, help="Override output directory"
    )
    run_parser.add_argument(
        "--verbose", "-v", action="store_true", help="Enable verbose logging"
    )
    run_parser.set_defaults(func=cmd_run)

    # List molecules command
    list_molecules_parser = subparsers.add_parser(
        "list-molecules", help="List available molecules from the force field"
    )
    list_molecules_parser.set_defaults(func=cmd_list_molecules)

    # New config command (placeholder for future)
    new_config_parser = subparsers.add_parser(
        "new-config", help="Create a new campaign configuration template"
    )
    new_config_parser.add_argument("name", help="Name of the new configuration")
    new_config_parser.set_defaults(func=lambda args: 1 if hasattr(args, 'func') else None)

    # Parse arguments
    parsed_args = parser.parse_args()

    if not parsed_args.command:
        parser.print_help()
        return 0

    # Execute command function
    func = getattr(parsed_args, "func", None)
    if func and callable(func):
        try:
            exit_code = func(parsed_args)
            return exit_code if isinstance(exit_code, int) else 0
        except Exception as e:
            logging.getLogger("martini_harness").error(f"Command failed: {e}")
            import traceback

            traceback.print_exc()
            return 1
    else:
        parser.print_help()
        return 0


if __name__ == "__main__":
    exit(main())
