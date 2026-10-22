import json
from pathlib import Path

import pytest

from martini_harness import cli
from martini_harness.check import (
    COMPLETE,
    FAILED,
    INCOMPLETE,
    NOT_STARTED,
    check_production,
    expected_end_time,
    last_xvg_time,
)
from martini_harness.config import load_campaign
from martini_harness.simulation import build_simulations, run_campaign
from martini_harness.state import STATE_FILE_NAME

from .conftest import BASIC_CAMPAIGN, FakeBackend

HEADER = '# dhdl\n@    title "dH/dl"\n@ s0 legend "Thermodynamic state"\n'

SHORT_AWH = BASIC_CAMPAIGN.replace("solvent: [W, CLF]", "solvent: [W, CLF, OCO, TCM]") + """
    stages:
      - {name: min, mdp: min.mdp}
      - name: awh
        mdp: awh.mdp
        substitute: {TOCOUPLE: "{target}"}
        mdp_params: {nsteps: 1000, dt: 0.02, nstdhdl: 100}
"""
END_TIME = 20.0  # 1000 steps * 0.02 ps


def _xvg(path: Path, last_time: float, step: float = 2.0) -> None:
    times = [i * step for i in range(int(last_time / step) + 1)]
    path.write_text(HEADER + "".join(f"{t:.4f} 0 -1.5 2.5\n" for t in times))


class XvgBackend(FakeBackend):
    """Writes <stage>.xvg up to a fraction of the run, per simulation directory."""

    def __init__(self, written_fraction: dict[str, float], **kwargs):
        super().__init__(**kwargs)
        self.written_fraction = written_fraction

    def mdrun(self, flags, log):
        super().mdrun(flags, log)
        deffnm = Path(flags["deffnm"])
        fraction = self.written_fraction.get(deffnm.parent.name, 1.0)
        _xvg(deffnm.with_name(deffnm.name + ".xvg"), END_TIME * fraction)


def test_last_xvg_time_reads_last_complete_line(tmp_path):
    xvg = tmp_path / "a.xvg"
    _xvg(xvg, 10.0)
    assert last_xvg_time(xvg) == 10.0

    with open(xvg, "a") as handle:
        handle.write("12.0000 0 -1.")  # partially written frame
    assert last_xvg_time(xvg) == 10.0

    (tmp_path / "header.xvg").write_text(HEADER)
    assert last_xvg_time(tmp_path / "header.xvg") is None
    assert last_xvg_time(tmp_path / "missing.xvg") is None


def test_last_xvg_time_on_large_file(tmp_path):
    xvg = tmp_path / "big.xvg"
    _xvg(xvg, 50_000.0, step=1.0)
    assert xvg.stat().st_size > 1_000_000
    assert last_xvg_time(xvg) == 50_000.0


@pytest.mark.parametrize(
    "mdp, override, expected",
    [
        ("nsteps = 500000\ndt = 0.02\n", None, 10_000.0),  # nstdhdl defaults to 50
        ("nsteps = 1050\ndt = 0.02\nnstdhdl = 100\n", None, 20.0),  # last multiple of nstdhdl
        ("nsteps = 1000\ndt = 0.02\ntinit = 5\nnstdhdl = 100\n", None, 25.0),
        ("nsteps = 1000\ndt = 0.02\nnstdhdl = 100\n", 500, 10.0),  # mdrun -nsteps
        ("nsteps = -1\ndt = 0.02\n", None, None),
        ("nsteps = 1000\nnstdhdl = 0\n", None, None),
    ],
)
def test_expected_end_time(mdp, override, expected):
    assert expected_end_time(mdp, override) == expected


def test_check_classifies_simulations(campaign_file):
    path = campaign_file(SHORT_AWH)
    backend = XvgBackend({"row42_in_CLF": 0.5})
    simulations = build_simulations(load_campaign(path))
    # row42_in_TCM is never run; row42_in_OCO fails at awh.
    to_run = [sim for sim in simulations if sim.name != "row42_in_TCM"]
    run_campaign([s for s in to_run if s.name != "row42_in_OCO"], backend)
    run_campaign([s for s in to_run if s.name == "row42_in_OCO"], FakeBackend(fail_on={"awh"}))

    checks = {c.name: c for c in (check_production(sim, "awh") for sim in simulations)}
    assert checks["row42_in_W"].category == COMPLETE
    assert checks["row42_in_W"].end_time == END_TIME
    assert checks["row42_in_CLF"].category == INCOMPLETE
    assert checks["row42_in_CLF"].progress == 0.5
    assert (checks["row42_in_OCO"].category, checks["row42_in_OCO"].failed_step) == (FAILED, "awh")
    assert (checks["row42_in_TCM"].category, checks["row42_in_TCM"].next_step) == (NOT_STARTED, "topology")


def test_check_uses_mdout_written_by_grompp(campaign_file, tmp_path):
    path = campaign_file(SHORT_AWH)
    simulations = build_simulations(load_campaign(path))
    run_campaign(simulations[:1], XvgBackend({}))
    workdir = tmp_path / "out" / "row42_in_W"
    # grompp's mdout says the run was twice as long as the xvg reaches.
    (workdir / "awh_mdout.mdp").write_text("nsteps = 2000\ndt = 0.02\nnstdhdl = 100\n")
    check = check_production(simulations[0], "awh")
    assert (check.category, check.progress) == (INCOMPLETE, 0.5)


def test_check_cli(campaign_file, capsys):
    path = campaign_file(SHORT_AWH)
    run_campaign(build_simulations(load_campaign(path))[:2], XvgBackend({"row42_in_CLF": 0.5}))

    assert cli.main(["check", str(path)]) == 0
    out = capsys.readouterr().out
    assert "Production stage 'awh': 4 simulation(s)" in out
    assert "complete          1  ( 25.0%)" in out
    assert "median progress 50%" in out
    assert "not started       2" in out

    assert cli.main(["check", str(path), "--list", "complete"]) == 0
    assert capsys.readouterr().out == "row42_in_W\n"

    assert cli.main(["check", str(path), "--details"]) == 0
    assert "row42_in_CLF  incomplete    50.0%" in capsys.readouterr().out

    assert cli.main(["check", str(path), "--stage", "min"]) == 0  # no min.xvg: never complete
    assert cli.main(["check", str(path), "--stage", "nope"]) == 2


ERRORS = Path(__file__).parent / "data" / "gromacs_errors"


class CrashingBackend(XvgBackend):
    """mdrun of the awh stage fails like real GROMACS, per simulation directory."""

    def __init__(self, crashes: dict[str, str]):
        super().__init__({})
        self.crashes = crashes

    def mdrun(self, flags, log):
        deffnm = Path(flags["deffnm"])
        sample = self.crashes.get(deffnm.parent.name)
        if sample and deffnm.name == "awh":
            Path(log).write_text((ERRORS / sample).read_text())
            raise RuntimeError(f"mdrun failed; see {log}")
        super().mdrun(flags, log)


def test_failure_reasons_in_state_status_and_check(campaign_file, tmp_path, capsys):
    path = campaign_file(SHORT_AWH)
    crashes = {"row42_in_CLF": "lpf.mdrun.out", "row42_in_OCO": "bgpu.mdrun.out", "row42_in_TCM": "lpf.mdrun.out"}
    outcomes = run_campaign(build_simulations(load_campaign(path)), CrashingBackend(crashes))
    reasons = {o.name: o.reason for o in outcomes}
    assert reasons["row42_in_W"] is None
    assert reasons["row42_in_CLF"].startswith("non-finite forces (pforce)")
    assert reasons["row42_in_OCO"].startswith("GPU: Bonded interactions on GPU are not supported")

    status = {s.name: s for s in build_simulations(load_campaign(path))[1].status()}
    assert status["awh"].reason.startswith("non-finite forces (pforce)")

    # An older batch has no reason in its state files: check reads the logs instead.
    state_file = tmp_path / "out" / "row42_in_OCO" / STATE_FILE_NAME
    state = json.loads(state_file.read_text())
    del state["steps"]["awh"]["reason"]
    state_file.write_text(json.dumps(state))

    assert cli.main(["check", str(path)]) == 0
    out = capsys.readouterr().out
    assert "failed            3  ( 75.0%)  at awh 3" in out
    assert "reasons: non-finite forces (pforce) 2, GPU 1" in out
    assert cli.main(["check", str(path), "--details"]) == 0
    assert "row42_in_OCO  failed           -  at awh  GPU: Bonded interactions" in capsys.readouterr().out


def test_crashed_incomplete_run_is_reported(campaign_file, tmp_path, capsys):
    path = campaign_file(SHORT_AWH)
    run_campaign(build_simulations(load_campaign(path))[:1], XvgBackend({"row42_in_W": 0.5}))
    # The harness itself was killed, so no failure was recorded, but mdrun's output has the error.
    workdir = tmp_path / "out" / "row42_in_W"
    (workdir / "awh.mdrun.out").write_text((ERRORS / "lincs.mdrun.out").read_text())

    assert cli.main(["check", str(path)]) == 0
    assert "1 with an error in the log: LINCS, then non-finite energy 1" in capsys.readouterr().out
