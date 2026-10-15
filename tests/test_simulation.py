import json

from martini_harness.config import load_campaign
from martini_harness.simulation import build_simulations, run_campaign
from martini_harness.state import STATE_FILE_NAME

from .conftest import BASIC_CAMPAIGN, FakeBackend

ALL_STEPS = ["topology", "coords", "min", "eq", "awh"]


def _simulations(path):
    return build_simulations(load_campaign(path))


def test_run_creates_files_in_each_workdir(campaign_file, tmp_path, fake_backend):
    outcomes = run_campaign(_simulations(campaign_file(BASIC_CAMPAIGN)), fake_backend)

    assert [o.success for o in outcomes] == [True, True]
    assert outcomes[0].ran == ALL_STEPS
    workdir = tmp_path / "out" / "row42_in_W"
    top = (workdir / "system.top").read_text()
    assert "row42.itp" in top and "W          100" in top
    awh_mdp = (workdir / "awh.mdp").read_text()
    assert "TOCOUPLE" not in awh_mdp
    assert any(line.startswith("couple-moltype") and "row42" in line for line in awh_mdp.splitlines())
    spec = json.loads((workdir / "simulation.json").read_text())
    assert spec["solvent"] == "W"
    assert fake_backend.stages_run() == [
        ("row42_in_W", "min"), ("row42_in_W", "eq"), ("row42_in_W", "awh"),
        ("row42_in_CLF", "min"), ("row42_in_CLF", "eq"), ("row42_in_CLF", "awh"),
    ]


def test_stages_chain_structures_and_respect_input(campaign_file, tmp_path, fake_backend):
    body = BASIC_CAMPAIGN + """
    stages:
      - {name: min, mdp: min.mdp, maxwarn: 1, mdrun: {nsteps: 10}}
      - {name: eq, mdp: eq.mdp}
      - {name: eq2, mdp: eq.mdp, input: min, grompp: {r: "{name}.gro"}}
    """
    run_campaign(_simulations(campaign_file(body))[:1], fake_backend)

    grompps = [flags for kind, flags in fake_backend.calls if kind == "grompp"]
    mdruns = [flags for kind, flags in fake_backend.calls if kind == "mdrun"]
    assert [g["c"].name for g in grompps] == ["init.gro", "min.gro", "min.gro"]
    assert grompps[0]["maxwarn"] == 1 and grompps[1]["maxwarn"] == 0
    assert grompps[2]["r"] == "row42_in_W.gro"
    assert mdruns[0]["nsteps"] == 10 and "nsteps" not in mdruns[1]


def test_second_run_skips_finished_steps(campaign_file, fake_backend):
    path = campaign_file(BASIC_CAMPAIGN)
    run_campaign(_simulations(path), fake_backend)

    rerun = FakeBackend()
    outcomes = run_campaign(_simulations(path), rerun)
    assert rerun.calls == []
    assert all(o.ran == [] and o.skipped == ALL_STEPS for o in outcomes)
    assert all(s.status == "done" for sim in _simulations(path) for s in sim.status())


def test_changed_stage_reruns_it_and_everything_after(campaign_file, fake_backend):
    path = campaign_file(BASIC_CAMPAIGN)
    run_campaign(_simulations(path), fake_backend)

    changed = BASIC_CAMPAIGN + """
    stages:
      - {name: min, mdp: min.mdp}
      - {name: eq, mdp: eq.mdp, maxwarn: 2, mdp_params: {nsteps: 100}}
      - {name: awh, mdp: awh.mdp, substitute: {TOCOUPLE: "{target}"}}
    """
    path = campaign_file(changed)
    statuses = [s.status for s in _simulations(path)[0].status()]
    assert statuses == ["done", "done", "done", "stale", "pending"]

    rerun = FakeBackend()
    outcome = run_campaign(_simulations(path), rerun)[0]
    assert outcome.skipped == ["topology", "coords", "min"]
    assert outcome.ran == ["eq", "awh"]


def test_changed_matrix_parameter_reruns_from_coordinates(campaign_file, fake_backend):
    path = campaign_file(BASIC_CAMPAIGN)
    run_campaign(_simulations(path), fake_backend)

    path = campaign_file(BASIC_CAMPAIGN.replace("n_solvent: 100", "n_solvent: 100\n      density: 950"))
    rerun = FakeBackend()
    outcome = run_campaign(_simulations(path), rerun)[0]
    assert outcome.skipped == ["topology"]
    assert outcome.ran == ["coords", "min", "eq", "awh"]


def test_missing_output_reruns_from_that_step(campaign_file, tmp_path, fake_backend):
    path = campaign_file(BASIC_CAMPAIGN)
    run_campaign(_simulations(path), fake_backend)
    (tmp_path / "out" / "row42_in_W" / "eq.gro").unlink()

    outcomes = run_campaign(_simulations(path), FakeBackend())
    assert outcomes[0].ran == ["eq", "awh"]
    assert outcomes[1].ran == []


def test_failure_stops_simulation_but_not_campaign(campaign_file, tmp_path):
    path = campaign_file(BASIC_CAMPAIGN)
    backend = FakeBackend(fail_on={"eq"})
    outcomes = run_campaign(_simulations(path), backend)

    assert [(o.success, o.failed) for o in outcomes] == [(False, "eq"), (False, "eq")]
    assert outcomes[0].ran == ["topology", "coords", "min"]
    assert [s.status for s in _simulations(path)[0].status()] == [
        "done", "done", "done", "failed", "pending",
    ]
    state = json.loads((tmp_path / "out" / "row42_in_W" / STATE_FILE_NAME).read_text())
    assert "grompp failed" in state["steps"]["eq"]["error"]
    assert "awh" not in state["steps"]

    # Resuming after the problem is fixed only reruns the failed stage onwards.
    retry = FakeBackend()
    outcomes = run_campaign(_simulations(path), retry)
    assert all(o.success for o in outcomes)
    assert outcomes[0].ran == ["eq", "awh"]


def test_rerunning_a_step_invalidates_later_state(campaign_file, tmp_path, fake_backend):
    """If a run is interrupted after an upstream rerun, stale downstream results are not reused."""
    path = campaign_file(BASIC_CAMPAIGN)
    run_campaign(_simulations(path), fake_backend)
    workdir = tmp_path / "out" / "row42_in_W"
    (workdir / "min.gro").unlink()

    class Interrupt(FakeBackend):
        def mdrun(self, flags, log):
            super().mdrun(flags, log)
            if flags["deffnm"].name == "min":
                raise KeyboardInterrupt

    try:
        run_campaign(_simulations(path)[:1], Interrupt())
    except KeyboardInterrupt:
        pass
    statuses = [s.status for s in _simulations(path)[0].status()]
    assert statuses == ["done", "done", "pending", "pending", "pending"]


def test_force_reruns_everything(campaign_file, fake_backend):
    path = campaign_file(BASIC_CAMPAIGN)
    run_campaign(_simulations(path), fake_backend)
    outcomes = run_campaign(_simulations(path), FakeBackend(), force=True)
    assert all(o.ran == ALL_STEPS for o in outcomes)


def test_global_mdrun_flags_skip_minimization(campaign_file, fake_backend):
    body = BASIC_CAMPAIGN + """
    mdrun: {bonded: gpu, ntomp: 8, ntmpi: 1}
    stages:
      - {name: min, mdp: min.mdp}
      - {name: eq, mdp: eq.mdp, maxwarn: 2, mdrun: {ntomp: 4}}
      - {name: awh, mdp: awh.mdp, substitute: {TOCOUPLE: "{target}"}}
    """
    run_campaign(_simulations(campaign_file(body))[:1], fake_backend)

    mdruns = {flags["deffnm"].name: flags for kind, flags in fake_backend.calls if kind == "mdrun"}
    assert not {"bonded", "ntomp", "ntmpi"} & set(mdruns["min"])
    assert (mdruns["eq"]["bonded"], mdruns["eq"]["ntomp"], mdruns["eq"]["ntmpi"]) == ("gpu", 4, 1)
    assert (mdruns["awh"]["bonded"], mdruns["awh"]["ntomp"], mdruns["awh"]["ntmpi"]) == ("gpu", 8, 1)


def test_changing_global_mdrun_flags_does_not_rerun(campaign_file, fake_backend):
    path = campaign_file(BASIC_CAMPAIGN)
    run_campaign(_simulations(path), fake_backend)

    path = campaign_file(BASIC_CAMPAIGN + "    mdrun: {ntomp: 16}\n")
    outcomes = run_campaign(_simulations(path), FakeBackend())
    assert all(o.ran == [] for o in outcomes)
