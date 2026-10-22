import os
from pathlib import Path

import pytest

from martini_harness.backend import BackendError, GromacsWrapperBackend


class FakeTool:
    command_name = "fake"

    def __init__(self, returncode=0):
        self.returncode = returncode
        self.calls = []

    def __call__(self, stdout=None, stderr=None, **flags):
        self.calls.append((os.getcwd(), flags))
        Path("step1b.pdb").write_text("dump")  # GROMACS writes these to the cwd
        stdout.write("output\n")
        return self.returncode, None, None


def test_tools_run_inside_the_simulation_directory(tmp_path, monkeypatch):
    workdir = tmp_path / "sim"
    workdir.mkdir()
    monkeypatch.chdir(tmp_path)
    tool = FakeTool()

    GromacsWrapperBackend._run(tool, {"s": workdir / "a.tpr", "nt": 2}, workdir / "a.mdrun.out")

    cwd, flags = tool.calls[0]
    assert Path(cwd).resolve() == workdir.resolve()
    assert flags == {"s": str(workdir / "a.tpr"), "nt": 2}
    assert (workdir / "step1b.pdb").exists() and not (tmp_path / "step1b.pdb").exists()
    assert Path.cwd().resolve() == tmp_path.resolve()
    assert (workdir / "a.mdrun.out").read_text() == "output\n"


def test_nonzero_exit_raises_and_restores_cwd(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "sim").mkdir()
    with pytest.raises(BackendError, match="exited with code 1"):
        GromacsWrapperBackend._run(FakeTool(returncode=1), {}, tmp_path / "sim" / "x.out")
    assert Path.cwd().resolve() == tmp_path.resolve()
