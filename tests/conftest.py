import shutil
import textwrap
from pathlib import Path

import pytest

EXAMPLES = Path(__file__).parent.parent / "examples"


class FakeBackend:
    """Stands in for polyply/GROMACS: records calls and writes the expected outputs."""

    def __init__(self, fail_on: set[str] | None = None):
        self.calls: list[tuple[str, dict]] = []
        self.fail_on = fail_on or set()

    def generate_coordinates(self, topology, output, density, name, log):
        self.calls.append(("coords", {"topology": topology, "output": output, "density": density}))
        assert topology.exists()
        output.write_text(f"fake coordinates for {name}\n")

    def grompp(self, flags, log):
        self.calls.append(("grompp", flags))
        for key in ("f", "c", "p"):
            assert Path(flags[key]).exists(), f"grompp input -{key} missing: {flags[key]}"
        if Path(flags["o"]).stem in self.fail_on:
            raise RuntimeError(f"grompp failed for {flags['o']}")
        Path(flags["o"]).write_text("tpr")

    def mdrun(self, flags, log):
        self.calls.append(("mdrun", flags))
        deffnm = Path(flags["deffnm"])
        for suffix in (".gro", ".log", ".edr"):
            deffnm.with_name(deffnm.name + suffix).write_text("output")

    def stages_run(self) -> list[tuple[str, str]]:
        """(simulation directory, stage) for every mdrun call."""
        return [
            (Path(flags["deffnm"]).parent.name, Path(flags["deffnm"]).name)
            for kind, flags in self.calls
            if kind == "mdrun"
        ]


@pytest.fixture
def fake_backend():
    return FakeBackend()


@pytest.fixture
def campaign_file(tmp_path):
    """Write a campaign YAML (plus row42.itp) into tmp_path and return its path."""
    shutil.copy(EXAMPLES / "row42.itp", tmp_path / "row42.itp")

    def write(body: str) -> Path:
        path = tmp_path / "campaign.yaml"
        path.write_text(textwrap.dedent(body))
        return path

    return write


BASIC_CAMPAIGN = """
    title: test
    output_dir: out
    molecules:
      row42: row42.itp
    simulation:
      n_solvent: 100
    matrix:
      target: [row42]
      solvent: [W, CLF]
"""
