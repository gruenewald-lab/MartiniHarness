"""Helpers for GROMACS topology (.top/.itp) files."""

from pathlib import Path


def read_moleculetypes(itp_file: str | Path) -> list[str]:
    """Return the names of all ``[ moleculetype ]`` entries in an ITP file.

    This is a lightweight parser: it does not follow ``#include`` statements
    or evaluate ``#ifdef`` blocks.
    """
    names = []
    in_moleculetype = False
    for raw_line in Path(itp_file).read_text().splitlines():
        line = raw_line.split(";", 1)[0].strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("["):
            in_moleculetype = line.strip("[] ").lower() == "moleculetype"
            continue
        if in_moleculetype:
            names.append(line.split()[0])
            in_moleculetype = False
    return names


def render_topology(
    includes: list[Path],
    system_name: str,
    molecules: list[tuple[str, int]],
) -> str:
    """Render the text of a system topology file.

    Args:
        includes: ITP files to ``#include``, in order (force field first).
        system_name: Title written to the ``[ system ]`` section.
        molecules: ``(name, count)`` pairs for the ``[ molecules ]`` section.

    Returns:
        The topology file content.
    """
    lines = [f'#include "{path}"' for path in includes]
    lines += ["", "[ system ]", system_name, "", "[ molecules ]"]
    lines += [f"{name:<10} {count}" for name, count in molecules]
    return "\n".join(lines) + "\n"
