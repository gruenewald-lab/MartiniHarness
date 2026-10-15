"""Command-line interface for Martini Harness."""

import argparse
import logging
import sys
from pathlib import Path

from .config import ConfigError, load_campaign
from .simulation import (
    STATUS_DONE,
    STATUS_FAILED,
    Simulation,
    build_simulations,
    run_campaign,
    select_simulations,
)
from .topology import read_moleculetypes

logger = logging.getLogger("martini_harness")


def setup_logging(verbose: bool = False) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)-7s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )


def _load_simulations(args: argparse.Namespace) -> list[Simulation]:
    config = load_campaign(args.config)
    simulations = build_simulations(config, output_dir=args.output_dir)
    if args.only:
        simulations = select_simulations(simulations, args.only)
    return simulations


def print_status(simulations: list[Simulation]) -> None:
    """Print one line per simulation with the status of each step."""
    if not simulations:
        print("No simulations.")
        return
    width = max(len(sim.name) for sim in simulations)
    for sim in simulations:
        statuses = sim.status()
        steps = "  ".join(f"{s.name}:{s.status}" for s in statuses)
        print(f"{sim.name:<{width}}  {steps}")
        for s in statuses:
            if s.status == STATUS_FAILED and s.error:
                print(f"{'':<{width}}    {s.name} error: {s.error}")
    print(f"\nOutput directory: {simulations[0].spec.workdir.parent}")


def cmd_run(args: argparse.Namespace) -> int:
    simulations = _load_simulations(args)
    if args.dry_run:
        print_status(simulations)
        pending = sum(
            any(s.status != STATUS_DONE for s in sim.status()) for sim in simulations
        )
        print(f"{pending}/{len(simulations)} simulation(s) have steps to run.")
        return 0

    outcomes = run_campaign(simulations, force=args.force)
    failed = [o for o in outcomes if not o.success]
    logger.info("Campaign finished: %d/%d simulations succeeded", len(outcomes) - len(failed), len(outcomes))
    for outcome in failed:
        logger.error("  %s failed at %s: %s", outcome.name, outcome.failed, outcome.error)
    return 1 if failed else 0


def cmd_status(args: argparse.Namespace) -> int:
    print_status(_load_simulations(args))
    return 0


def cmd_list_molecules(args: argparse.Namespace) -> int:
    for itp in args.itp_files:
        for name in read_moleculetypes(itp):
            print(f"{name}\t{itp}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="martini_harness",
        description="Run Martini simulation campaigns defined in YAML files.",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="enable debug logging")
    subparsers = parser.add_subparsers(dest="command", required=True)

    def add_campaign_args(sub: argparse.ArgumentParser) -> None:
        sub.add_argument("config", type=Path, help="campaign YAML file")
        sub.add_argument(
            "-o", "--output-dir", type=Path, help="override the campaign's output_dir"
        )
        sub.add_argument(
            "--only", action="append", metavar="NAME",
            help="restrict to this simulation (repeatable)",
        )

    run_parser = subparsers.add_parser("run", help="run a campaign, resuming where it stopped")
    add_campaign_args(run_parser)
    run_parser.add_argument(
        "--dry-run", action="store_true", help="show what would run without running it"
    )
    run_parser.add_argument(
        "--force", action="store_true", help="rerun all steps, ignoring saved state"
    )
    run_parser.set_defaults(func=cmd_run)

    status_parser = subparsers.add_parser("status", help="show the progress of a campaign")
    add_campaign_args(status_parser)
    status_parser.set_defaults(func=cmd_status)

    list_parser = subparsers.add_parser(
        "list-molecules", help="list the molecule types defined in ITP files"
    )
    list_parser.add_argument("itp_files", type=Path, nargs="+")
    list_parser.set_defaults(func=cmd_list_molecules)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    setup_logging(args.verbose)
    try:
        return args.func(args)
    except ConfigError as exc:
        logger.error("Invalid configuration:\n%s", exc)
        return 2
    except KeyboardInterrupt:
        logger.error("Interrupted; rerun the same command to resume.")
        return 130


if __name__ == "__main__":
    sys.exit(main())
