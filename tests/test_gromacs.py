"""End-to-end run with real polyply and GROMACS on a tiny system."""

import shutil

import pytest

from martini_harness.config import load_campaign
from martini_harness.simulation import build_simulations, run_campaign

from .conftest import BASIC_CAMPAIGN

pytestmark = [
    pytest.mark.gromacs,
    pytest.mark.skipif(shutil.which("gmx") is None, reason="gmx not on PATH"),
]

SHORT_STAGES = """
    stages:
      - {name: min, mdp: min.mdp, mdp_params: {nsteps: 50}}
      - {name: eq, mdp: eq.mdp, maxwarn: 2, mdp_params: {nsteps: 50, dt: 0.005}}
      - name: awh
        mdp: awh.mdp
        substitute: {TOCOUPLE: "{target}"}
        mdp_params: {define: "", nsteps: 50, nstcalcenergy: 50, nstenergy: 50, awh-nstout: 50}
"""


def test_short_campaign_runs_and_resumes(campaign_file, tmp_path):
    body = BASIC_CAMPAIGN.replace("solvent: [W, CLF]", "solvent: [W]") + SHORT_STAGES
    path = campaign_file(body.replace("n_solvent: 100", "n_solvent: 300\n      density: 997"))

    outcome = run_campaign(build_simulations(load_campaign(path)))[0]
    assert outcome.success, outcome.error
    workdir = tmp_path / "out" / "row42_in_W"
    for name in ("init.gro", "min.gro", "eq.gro", "awh.gro", "awh.edr", "eq_mdout.mdp"):
        assert (workdir / name).exists(), name

    outcome = run_campaign(build_simulations(load_campaign(path)))[0]
    assert outcome.ran == []
