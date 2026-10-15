"""Build-matrix expansion with GitHub Actions semantics."""

import itertools
from typing import Any

Combination = dict[str, Any]


def expand_matrix(
    axes: dict[str, list[Any]],
    include: list[Combination] | None = None,
    exclude: list[Combination] | None = None,
) -> list[Combination]:
    """Expand a build matrix into a list of combinations.

    Follows the rules of GitHub Actions' ``strategy.matrix``:

    1. The Cartesian product of all ``axes`` forms the original combinations.
    2. Every ``exclude`` entry removes the combinations it (partially) matches.
    3. Every ``include`` entry is merged into each original combination whose
       axis values it does not contradict. Values that were added by an earlier
       ``include`` may be overwritten, original axis values never are. If an
       entry matches no original combination it is appended as a new one.

    Args:
        axes: Mapping of axis name to the list of values it takes.
        include: Entries that extend existing combinations or add new ones.
        exclude: Entries that remove matching combinations.

    Returns:
        The expanded combinations, in a deterministic order.
    """
    include = include or []
    exclude = exclude or []

    if axes:
        combinations = [
            dict(zip(axes, values)) for values in itertools.product(*axes.values())
        ]
    else:
        combinations = []

    combinations = [
        combo
        for combo in combinations
        if not any(_matches(combo, entry) for entry in exclude)
    ]

    originals = list(combinations)
    for entry in include:
        matched = False
        for combo in originals:
            if all(combo[key] == value for key, value in entry.items() if key in axes):
                combo.update(entry)
                matched = True
        if not matched:
            combinations.append(dict(entry))

    return combinations


def _matches(combo: Combination, entry: Combination) -> bool:
    return all(key in combo and combo[key] == value for key, value in entry.items())
