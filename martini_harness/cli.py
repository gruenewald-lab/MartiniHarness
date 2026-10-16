"""Command-line interface for Martini Harness."""

import argparse
import json
import logging
import os
import random
import sys
import time
from pathlib import Path

from .config import ConfigError, load_campaign
from .library import (
    LibraryError,
    MoleculeLibrary,
    sample_molecules,
    sample_per_group,
    split_itp,
)
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


def cmd_split_itp(args: argparse.Namespace) -> int:
    start = time.perf_counter()
    result = split_itp(args.source, args.output, force=args.force)
    elapsed = time.perf_counter() - start
    for warning in result.warnings[:20]:
        logger.warning(warning)
    if len(result.warnings) > 20:
        logger.warning("... and %d more warnings", len(result.warnings) - 20)
    if result.duplicates:
        logger.warning(
            "%d duplicate molecule name(s) skipped (first definition kept), e.g. %s",
            len(result.duplicates), ", ".join(result.duplicates[:5]),
        )
    logger.info("Wrote %d molecules to %s in %.1f s", result.written, args.output, elapsed)
    return 0


def cmd_sample(args: argparse.Namespace) -> int:
    library = MoleculeLibrary(args.library)
    seed = args.seed if args.seed is not None else random.SystemRandom().randrange(2**31)
    if args.per_group is not None:
        names, sizes = sample_per_group(library, args.per_group, seed)
        how = f"{args.per_group} from each of {len(sizes)} groups"
        logger.info("Group sizes: %s", ", ".join(f"{g}: {n}" for g, n in sizes.items()))
    else:
        names = sample_molecules(library, args.n, seed)
        how = "uniformly"

    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    library_path = os.path.relpath(library.library_dir.resolve(), output.parent)
    lines = [
        f"# {len(names)} of {len(library)} molecules sampled {how} from {library_path} (seed {seed}).",
        "# Import this file from a campaign: `imports: [" + output.name + "]`.",
        f"molecule_library: {json.dumps(library_path)}",
        "matrix:",
        "  target:",
    ]
    lines += [f"    - {json.dumps(name)}" for name in names]
    output.write_text("\n".join(lines) + "\n")
    logger.info("Sampled %d of %d molecules (seed %d) into %s", len(names), len(library), seed, output)
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

    split_parser = subparsers.add_parser(
        "split-itp", help="split a multi-molecule ITP into a library of per-molecule ITPs"
    )
    split_parser.add_argument("source", type=Path, help="ITP file with many [ moleculetype ]s")
    split_parser.add_argument("-o", "--output", type=Path, required=True, help="library directory")
    split_parser.add_argument("--force", action="store_true", help="replace an existing library")
    split_parser.set_defaults(func=cmd_split_itp)

    sample_parser = subparsers.add_parser(
        "sample", help="randomly pick molecules from a library into an importable YAML file"
    )
    sample_parser.add_argument("library", type=Path, help="library directory from split-itp")
    how_many = sample_parser.add_mutually_exclusive_group(required=True)
    how_many.add_argument("-n", type=int, help="number of molecules, drawn uniformly")
    how_many.add_argument(
        "--per-group", type=int, metavar="K",
        help="draw K molecules from every group; the group is the name part before the "
        "first '_' (e.g. the length in <length>_<molid>)",
    )
    sample_parser.add_argument("--seed", type=int, help="random seed (default: random, recorded in the file)")
    sample_parser.add_argument("-o", "--output", type=Path, required=True, help="YAML file to write")
    sample_parser.set_defaults(func=cmd_sample)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    setup_logging(args.verbose)
    try:
        return args.func(args)
    except ConfigError as exc:
        logger.error("Invalid configuration:\n%s", exc)
        return 2
    except LibraryError as exc:
        logger.error("%s", exc)
        return 2
    except KeyboardInterrupt:
        logger.error("Interrupted; rerun the same command to resume.")
        return 130


if __name__ == "__main__":
    sys.exit(main())
