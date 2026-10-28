"""Molecule libraries: one ITP file per molecule plus an index.

A large multi-molecule ITP is split once with :func:`split_itp`; afterwards
campaigns only read the index and the few small ITP files they use.

Layout of a library directory::

    index.tsv           # "name<TAB>relative/path.itp", one line per molecule
    3f/mol_a.itp        # files are sharded into 256 subdirectories by a hash
    a1/mol_b.itp        # of the name, keeping directories small
    _preamble.itp       # only if the source had content before the first molecule
"""

import hashlib
import os
import random
import re
import shutil
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

INDEX_FILE_NAME = "index.tsv"
PREAMBLE_FILE_NAME = "_preamble.itp"

_ENCODING = {"encoding": "utf-8", "errors": "surrogateescape"}
_MOLECULETYPE = re.compile(r"^[ \t]*\[[ \t]*moleculetype[ \t]*\]", re.MULTILINE | re.IGNORECASE)
_TOPOLOGY_ONLY = re.compile(r"^[ \t]*\[[ \t]*(?:system|molecules)[ \t]*\]", re.MULTILINE | re.IGNORECASE)
_UNSAFE_CHARS = re.compile(r"[^A-Za-z0-9._+-]")


class LibraryError(ValueError):
    """Raised for a missing or malformed molecule library."""


@dataclass
class SplitResult:
    library_dir: Path
    written: int = 0
    duplicates: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    preamble: Path | None = None


def split_itp(source: str | Path, library_dir: str | Path, force: bool = False) -> SplitResult:
    """Split a multi-molecule ITP file into a molecule library.

    The source is read in one go and cut at every ``[ moleculetype ]``
    header, so even files with tens of thousands of molecules split in
    seconds. Every block up to the next header becomes one ITP file.

    Args:
        source: The multi-molecule ITP file.
        library_dir: Directory to create; must not contain a library yet.
        force: Replace an existing library in ``library_dir``.

    Returns:
        Counts and warnings. Duplicate molecule names keep the first
        definition. Blocks with unbalanced ``#if``/``#endif`` are written but
        reported, since a conditional spanning molecules cannot be split.
    """
    source = Path(source)
    library_dir = Path(library_dir)
    if (library_dir / INDEX_FILE_NAME).exists():
        if not force:
            raise LibraryError(f"{library_dir} already contains a library (use force to replace it)")
        shutil.rmtree(library_dir)

    text = source.read_text(**_ENCODING)
    starts = _moleculetype_offsets(text)
    if not starts:
        raise LibraryError(f"no [ moleculetype ] found in {source}")

    library_dir.mkdir(parents=True, exist_ok=True)
    result = SplitResult(library_dir)

    preamble = text[: starts[0]]
    if _has_content(preamble):
        result.preamble = library_dir / PREAMBLE_FILE_NAME
        result.preamble.write_text(preamble, **_ENCODING)
        result.warnings.append(
            f"content before the first molecule was saved to {PREAMBLE_FILE_NAME}; "
            "it is not included automatically"
        )

    index: dict[str, str] = {}
    used_paths: set[str] = set()
    files: list[tuple[str, str]] = []
    ends = starts[1:] + [len(text)]
    for start, end in zip(starts, ends):
        block = text[start:end]
        if end == len(text):
            trailing = _TOPOLOGY_ONLY.search(block)
            if trailing:
                block = block[: trailing.start()]
                result.warnings.append("ignored [ system ]/[ molecules ] sections at the end of the file")

        name = _molecule_name(block)
        if name is None:
            result.warnings.append(f"skipped a [ moleculetype ] without a name at offset {start}")
            continue
        if name in index:
            result.duplicates.append(name)
            continue
        if block.count("#if") != block.count("#endif"):
            result.warnings.append(f"{name}: unbalanced #if/#endif, check the written file")

        relative = _relative_path(name, used_paths)
        index[name] = relative
        files.append((relative, block))

    for shard in {relative.split("/", 1)[0] for relative, _ in files}:
        (library_dir / shard).mkdir(exist_ok=True)
    # Creating tens of thousands of small files is dominated by system calls,
    # which run in parallel across threads.
    with ThreadPoolExecutor(max_workers=8) as pool:
        for _ in pool.map(lambda item: _write_file(library_dir / item[0], item[1]), files):
            pass
    result.written = len(files)

    _write_index(library_dir, index, source)
    return result


def _moleculetype_offsets(text: str) -> list[int]:
    """Offsets of the lines holding a ``[ moleculetype ]`` header."""
    lowered = text.lower()
    if len(lowered) != len(text):  # a few Unicode characters change length when lowered
        return [match.start() for match in _MOLECULETYPE.finditer(text)]
    offsets = []
    position = lowered.find("moleculetype")
    while position != -1:
        line_start = lowered.rfind("\n", 0, position) + 1
        line_end = lowered.find("\n", position)
        line = lowered[line_start : line_end if line_end != -1 else len(lowered)]
        if "".join(line.split(";", 1)[0].split()) == "[moleculetype]":
            offsets.append(line_start)
        position = lowered.find("moleculetype", position + 12)
    return offsets


def _write_file(path: Path, content: str) -> None:
    with open(path, "w", **_ENCODING) as handle:
        handle.write(content)


def _has_content(text: str) -> bool:
    return any(line.split(";", 1)[0].strip() for line in text.splitlines())


def _molecule_name(block: str) -> str | None:
    position = block.find("\n") + 1
    while position:
        end = block.find("\n", position)
        line = block[position : end if end != -1 else len(block)]
        position = end + 1
        line = line.split(";", 1)[0].strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("["):
            return None
        return line.split()[0]
    return None


def _relative_path(name: str, used: set[str]) -> str:
    shard = hashlib.sha1(name.encode(**_ENCODING)).hexdigest()[:2]
    stem = _UNSAFE_CHARS.sub("_", name)
    relative = f"{shard}/{stem}.itp"
    counter = 1
    while relative.lower() in used:  # case-insensitive file systems
        counter += 1
        relative = f"{shard}/{stem}__{counter}.itp"
    used.add(relative.lower())
    return relative


def _write_index(library_dir: Path, index: dict[str, str], source: Path) -> None:
    lines = [f"# source: {source.resolve()}", f"# molecules: {len(index)}"]
    lines += [f"{name}\t{relative}" for name, relative in index.items()]
    tmp = library_dir / (INDEX_FILE_NAME + ".tmp")
    tmp.write_text("\n".join(lines) + "\n", **_ENCODING)
    os.replace(tmp, library_dir / INDEX_FILE_NAME)


class MoleculeLibrary:
    """Read access to a library created by :func:`split_itp`."""

    def __init__(self, library_dir: str | Path):
        self.library_dir = Path(library_dir)
        index_file = self.library_dir / INDEX_FILE_NAME
        if not index_file.is_file():
            raise LibraryError(f"no molecule library at {self.library_dir} ({INDEX_FILE_NAME} missing)")
        self._paths: dict[str, str] = {}
        for line in index_file.read_text(**_ENCODING).splitlines():
            if line and not line.startswith("#"):
                name, _, relative = line.partition("\t")
                self._paths[name] = relative

    def __contains__(self, name: str) -> bool:
        return name in self._paths

    def __len__(self) -> int:
        return len(self._paths)

    @property
    def names(self) -> list[str]:
        return list(self._paths)

    def path(self, name: str) -> Path:
        return self.library_dir / self._paths[name]


def sample_molecules(library: MoleculeLibrary, n: int, seed: int) -> list[str]:
    """Pick ``n`` distinct molecule names at random, in library order."""
    if not 0 < n <= len(library):
        raise LibraryError(f"cannot sample {n} of {len(library)} molecules")
    names = library.names
    chosen = set(random.Random(seed).sample(names, n))
    return [name for name in names if name in chosen]


GROUP_SEPARATOR = "_"


def group_of(name: str) -> str:
    """The group of a molecule: the part of its name before the first ``_``.

    For names formatted ``<length>_<molid>`` this is the length.
    """
    group, separator, _ = name.partition(GROUP_SEPARATOR)
    if not separator:
        raise LibraryError(f"molecule {name!r} has no {GROUP_SEPARATOR!r} to take its group from")
    return group


def sample_per_group(
    library: MoleculeLibrary, per_group: int, seed: int
) -> tuple[list[str], dict[str, int]]:
    """Pick ``per_group`` molecules at random from every group (see :func:`group_of`).

    Returns:
        The chosen names in library order, and the size of every group.

    Raises:
        LibraryError: If a name has no group or a group is smaller than ``per_group``.
    """
    if per_group < 1:
        raise LibraryError("per_group must be at least 1")
    groups: dict[str, list[str]] = {}
    for name in library.names:
        groups.setdefault(group_of(name), []).append(name)

    sizes = {group: len(groups[group]) for group in sorted(groups, key=_group_sort_key)}
    too_small = {group: size for group, size in sizes.items() if size < per_group}
    if too_small:
        listing = ", ".join(f"{group} ({size})" for group, size in too_small.items())
        raise LibraryError(f"groups with fewer than {per_group} molecules: {listing}")

    rng = random.Random(seed)
    chosen: set[str] = set()
    for group in sizes:  # fixed order keeps the result reproducible
        chosen.update(rng.sample(groups[group], per_group))
    return [name for name in library.names if name in chosen], sizes


def _group_sort_key(group: str) -> tuple[bool, int, str]:
    return (not group.isdigit(), int(group) if group.isdigit() else 0, group)
