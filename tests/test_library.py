import textwrap

import pytest

from martini_harness import cli
from martini_harness.config import ConfigError, load_campaign
from martini_harness.library import (
    INDEX_FILE_NAME,
    PREAMBLE_FILE_NAME,
    LibraryError,
    MoleculeLibrary,
    group_of,
    sample_molecules,
    sample_per_group,
    split_itp,
)
from martini_harness.simulation import build_simulations, run_campaign
from martini_harness.topology import read_moleculetypes

from .conftest import EXAMPLES, FakeBackend

MULTI_ITP = """\
; generated library
#define SOMETHING

[ moleculetype ]
; name nrexcl
alpha 1
[ atoms ]
1 TC3 1 ALP C1 1 0.0

[moleculetype]
beta 1
[ atoms ]
1 TC3 1 BET C1 1 0.0
#ifdef FLEXIBLE
[ bonds ]
#endif

  [ MOLECULETYPE ] ; odd spacing and case
weird/name 1
[ atoms ]
1 TC3 1 WRD C1 1 0.0

[ moleculetype ]
alpha 1
[ atoms ]
1 TC4 1 ALP C1 1 0.0

[ system ]
test
[ molecules ]
alpha 1
"""


@pytest.fixture
def library(tmp_path):
    source = tmp_path / "multi.itp"
    source.write_text(MULTI_ITP)
    result = split_itp(source, tmp_path / "lib")
    return result, MoleculeLibrary(tmp_path / "lib")


def test_split_writes_one_file_per_molecule(library):
    result, lib = library
    assert result.written == 3
    assert lib.names == ["alpha", "beta", "weird/name"]
    for name in lib.names:
        assert read_moleculetypes(lib.path(name)) == [name]
    assert "#ifdef FLEXIBLE" in lib.path("beta").read_text()
    assert lib.path("weird/name").name == "weird_name.itp"
    assert all(len(lib.path(name).parent.name) == 2 for name in lib.names)  # hash shards


def test_split_keeps_first_duplicate_and_drops_topology_sections(library):
    result, lib = library
    assert result.duplicates == ["alpha"]
    assert "TC3" in lib.path("alpha").read_text()
    assert all("[ system ]" not in lib.path(name).read_text() for name in lib.names)
    assert any("[ system ]" in warning for warning in result.warnings)


def test_split_saves_preamble(library, tmp_path):
    result, _ = library
    assert result.preamble == tmp_path / "lib" / PREAMBLE_FILE_NAME
    assert "#define SOMETHING" in result.preamble.read_text()


def test_split_refuses_to_overwrite_without_force(library, tmp_path):
    with pytest.raises(LibraryError, match="already contains a library"):
        split_itp(tmp_path / "multi.itp", tmp_path / "lib")
    assert split_itp(tmp_path / "multi.itp", tmp_path / "lib", force=True).written == 3


def test_split_reports_unbalanced_conditionals(tmp_path):
    source = tmp_path / "bad.itp"
    source.write_text("[ moleculetype ]\na 1\n#ifdef X\n[ moleculetype ]\nb 1\n#endif\n")
    result = split_itp(source, tmp_path / "lib")
    assert sum("unbalanced" in warning for warning in result.warnings) == 2


def test_split_handles_crlf_line_endings(tmp_path):
    source = tmp_path / "crlf.itp"
    source.write_bytes(b"[ moleculetype ]\r\nmol1 1\r\n[ moleculetype ]\r\nmol2 1\r\n")
    split_itp(source, tmp_path / "lib")
    assert MoleculeLibrary(tmp_path / "lib").names == ["mol1", "mol2"]


def test_split_without_molecules_fails(tmp_path):
    source = tmp_path / "empty.itp"
    source.write_text("[ atomtypes ]\n")
    with pytest.raises(LibraryError, match="no \\[ moleculetype \\]"):
        split_itp(source, tmp_path / "lib")


def test_missing_library(tmp_path):
    with pytest.raises(LibraryError, match=INDEX_FILE_NAME):
        MoleculeLibrary(tmp_path)


def test_sampling_is_reproducible(tmp_path):
    source = tmp_path / "many.itp"
    source.write_text("".join(f"[ moleculetype ]\nm{i:03d} 1\n" for i in range(100)))
    split_itp(source, tmp_path / "lib")
    lib = MoleculeLibrary(tmp_path / "lib")
    first = sample_molecules(lib, 10, seed=3)
    assert first == sample_molecules(lib, 10, seed=3)
    assert first != sample_molecules(lib, 10, seed=4)
    assert len(set(first)) == 10 and first == sorted(first)
    with pytest.raises(LibraryError):
        sample_molecules(lib, 101, seed=0)


def test_split_sample_and_run_campaign_from_cli(tmp_path):
    big = tmp_path / "big.itp"
    row42 = (EXAMPLES / "row42.itp").read_text()
    big.write_text("\n".join(row42.replace("row42", f"mol{i}") for i in range(20)))

    assert cli.main(["split-itp", str(big), "-o", str(tmp_path / "library")]) == 0
    targets = tmp_path / "campaign" / "targets.yaml"
    assert cli.main(["sample", str(tmp_path / "library"), "-n", "3", "--seed", "1", "-o", str(targets)]) == 0
    assert "seed 1" in targets.read_text()

    campaign = tmp_path / "campaign" / "campaign.yaml"
    campaign.write_text(
        textwrap.dedent(
            """
            imports: [targets.yaml]
            output_dir: out
            simulation: {n_solvent: 50}
            matrix:
              solvent: [W, OCO]
            """
        )
    )
    simulations = build_simulations(load_campaign(campaign))
    assert len(simulations) == 6
    spec = simulations[0].spec
    assert spec.includes[-1].parent.parent == tmp_path / "library"
    assert read_moleculetypes(spec.includes[-1]) == [spec.target]

    backend = FakeBackend()
    assert all(o.success for o in run_campaign(simulations, backend))
    top = (tmp_path / "campaign" / "out" / spec.name / "system.top").read_text()
    assert str(big) not in top and str(spec.includes[-1]) in top


def test_unknown_target_mentions_library(tmp_path):
    source = tmp_path / "lib.itp"
    source.write_text("[ moleculetype ]\nmol1 1\n")
    split_itp(source, tmp_path / "lib")
    campaign = tmp_path / "campaign.yaml"
    campaign.write_text("molecule_library: lib\nmatrix: {target: [mol2], solvent: [W]}\n")
    with pytest.raises(ConfigError, match="molecule library"):
        build_simulations(load_campaign(campaign))


def _length_library(tmp_path, sizes):
    """Library with names <length>_<molid>; sizes maps length -> number of molecules."""
    source = tmp_path / "lengths.itp"
    source.write_text(
        "".join(
            f"[ moleculetype ]\n{length}_{molid:04d} 1\n"
            for length, count in sizes.items()
            for molid in range(count)
        )
    )
    split_itp(source, tmp_path / "lib")
    return MoleculeLibrary(tmp_path / "lib")


def test_group_of():
    assert group_of("12_0345") == "12"
    assert group_of("12_03_45") == "12"
    with pytest.raises(LibraryError, match="no '_'"):
        group_of("row42")


def test_sample_per_group_draws_equal_counts(tmp_path):
    lib = _length_library(tmp_path, {5: 50, 12: 200, 8: 20})
    names, sizes = sample_per_group(lib, 10, seed=1)

    assert sizes == {"5": 50, "8": 20, "12": 200}  # numeric order
    counts = {}
    for name in names:
        counts[group_of(name)] = counts.get(group_of(name), 0) + 1
    assert counts == {"5": 10, "8": 10, "12": 10}
    assert names == [name for name in lib.names if name in set(names)]  # library order
    assert names == sample_per_group(lib, 10, seed=1)[0]
    assert names != sample_per_group(lib, 10, seed=2)[0]


def test_sample_per_group_reports_small_groups(tmp_path):
    lib = _length_library(tmp_path, {5: 50, 8: 3, 9: 1})
    with pytest.raises(LibraryError, match=r"fewer than 4 molecules: 8 \(3\), 9 \(1\)"):
        sample_per_group(lib, 4, seed=0)


def test_sample_cli_per_group(tmp_path):
    _length_library(tmp_path, {5: 30, 6: 30})
    out = tmp_path / "targets.yaml"
    assert cli.main(["sample", str(tmp_path / "lib"), "--per-group", "7", "--seed", "3", "-o", str(out)]) == 0
    text = out.read_text()
    assert "7 from each of 2 groups" in text
    targets = [line.split('"')[1] for line in text.splitlines() if line.startswith("    - ")]
    assert sorted(group_of(name) for name in targets) == ["5"] * 7 + ["6"] * 7
