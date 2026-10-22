#!/usr/bin/env python3
"""Add a CGSmiles column (coarse-grained resolution only) to an `analyze` table.

For every row, the molecule named in the ``target`` column is looked up in a
molecule library (the directory made by ``martini_harness split-itp``), its
ITP is read and the bead graph is written as a CGSmiles string such as
``{[#TC3]1[#TC3][#TC3]1}``. Nodes are the bead types of ``[ atoms ]``, edges
are the pairs listed in ``[ bonds ]`` and ``[ constraints ]`` (everything
inside ``#ifdef``/``#ifndef`` is taken, so the FLEXIBLE and constrained
variants give the same graph).

Virtual sites are usually not bonded to anything, but CGSmiles needs a
connected graph. By default such a site is attached to the first bead it is
constructed from, which adds no ring; use ``--vsites all`` to attach it to
every constructing bead. A molecule whose graph is still disconnected is
reported and gets an empty cell.

The table is updated in place unless ``-o`` is given. Rerunning replaces the
column. Note that ``martini_harness analyze`` rewrites its table from scratch,
so run this script again afterwards.

    python scripts/add_cgsmiles.py campaigns/screen/dG.tsv --library library/
"""

import argparse
import csv
import os
import re
import sys
import tempfile
from pathlib import Path

import networkx as nx
from cgsmiles.write_cgsmiles import write_cgsmiles_graph

INDEX_FILE_NAME = "index.tsv"
_SECTION = re.compile(r"^\[\s*(\w+)\s*\]")
_VSITE_PARENTS = {"virtual_sites2": 2, "virtual_sites3": 3, "virtual_sites4": 4}


def read_library_index(library: Path) -> dict[str, Path]:
    """Molecule name -> ITP path, from the index of a library made by split-itp."""
    index_file = library / INDEX_FILE_NAME
    if not index_file.is_file():
        raise SystemExit(f"no molecule library at {library} ({INDEX_FILE_NAME} missing)")
    index = {}
    for line in index_file.read_text(encoding="utf-8", errors="surrogateescape").splitlines():
        if line and not line.startswith("#"):
            name, _, relative = line.partition("\t")
            index[name] = library / relative
    return index


def bead_graph(itp: Path, molecule: str, vsites: str = "first") -> nx.Graph:
    """Graph of the beads of ``molecule`` in ``itp``; nodes carry their bead type as ``fragname``."""
    atoms: dict[int, str] = {}
    edges: set[frozenset[int]] = set()
    parents: dict[int, list[int]] = {}
    in_molecule = False
    section = None

    for raw in itp.read_text(errors="replace").splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line or line.startswith("#"):
            continue
        if header := _SECTION.match(line):
            section = header.group(1).lower()
            continue
        fields = line.split()
        if section == "moleculetype":
            in_molecule = fields[0] == molecule
        elif not in_molecule:
            continue
        elif section == "atoms":
            atoms[int(fields[0])] = fields[1]
        elif section in ("bonds", "constraints"):
            edges.add(frozenset((int(fields[0]), int(fields[1]))))
        elif section in _VSITE_PARENTS:
            parents[int(fields[0])] = [int(f) for f in fields[1 : 1 + _VSITE_PARENTS[section]]]
        elif section == "virtual_sitesn":
            constructing = fields[2::2] if fields[1] == "3" else fields[2:]  # 3: atom weight pairs
            parents[int(fields[0])] = [int(f) for f in constructing]

    if not atoms:
        raise ValueError(f"no atoms found for molecule {molecule!r} in {itp}")

    graph = nx.Graph()
    for index in sorted(atoms):
        graph.add_node(index, fragname=atoms[index])
    for edge in edges:
        if not edge <= atoms.keys() or len(edge) != 2:
            raise ValueError(f"bond {sorted(edge)} refers to an unknown atom in {itp}")
        graph.add_edge(*edge, order=1)

    for site, constructing in parents.items():
        if site in graph and graph.degree(site) == 0:
            for parent in constructing[: 1 if vsites == "first" else None]:
                graph.add_edge(site, parent, order=1)
    if not nx.is_connected(graph):
        raise ValueError(f"the beads of {molecule!r} in {itp} do not form one connected graph")
    return graph


def cgsmiles(itp: Path, molecule: str, vsites: str = "first") -> str:
    return write_cgsmiles_graph(bead_graph(itp, molecule, vsites))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Add a CG-resolution CGSmiles column to a results table (tab-separated).",
    )
    parser.add_argument("table", type=Path, help="TSV file, e.g. campaigns/screen/dG.tsv")
    parser.add_argument(
        "--library", type=Path, required=True, help="molecule library directory (from split-itp)"
    )
    parser.add_argument("-o", "--output", type=Path, help="write here instead of updating the table")
    parser.add_argument("--target-column", default="target", help="column with the molecule name")
    parser.add_argument("--cgsmiles-column", default="cgsmiles", help="column to write")
    parser.add_argument(
        "--vsites", choices=("first", "all"), default="first",
        help="connect a virtual site without bonds to its first or to all constructing beads "
        "(default: first)",
    )
    args = parser.parse_args(argv)

    with open(args.table, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    if args.target_column not in fields:
        raise SystemExit(f"{args.table}: no column {args.target_column!r} (columns: {', '.join(fields)})")
    if args.cgsmiles_column not in fields:
        fields.append(args.cgsmiles_column)

    index = read_library_index(args.library)
    cache: dict[str, str] = {}
    problems: dict[str, str] = {}
    for row in rows:
        molecule = row[args.target_column]
        if molecule not in cache:
            try:
                if molecule not in index:
                    raise ValueError("not in the molecule library")
                cache[molecule] = cgsmiles(index[molecule], molecule, args.vsites)
            except (ValueError, OSError) as exc:
                cache[molecule] = ""
                problems[molecule] = str(exc)
        row[args.cgsmiles_column] = cache[molecule]

    output = args.output or args.table
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(dir=output.parent, prefix=".cgsmiles_", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n", restval="")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(tmp_name, output)
    except BaseException:
        Path(tmp_name).unlink(missing_ok=True)
        raise

    for molecule, problem in problems.items():
        print(f"warning: {molecule}: {problem}", file=sys.stderr)
    print(
        f"{len(rows)} row(s), {len(cache) - len(problems)} of {len(cache)} molecule(s) -> {output}",
        file=sys.stderr,
    )
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
