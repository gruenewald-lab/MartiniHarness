"""End-to-end run with real polyply and GROMACS on a tiny system."""

import shutil
import subprocess

import pytest

from martini_harness.analysis import OK, run_analysis
from martini_harness.config import AnalysisConfig, load_campaign
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

    # 50 steps of 0.02 ps: gmx awh writes awh_t0.xvg and awh_t1.xvg; the last
    # line of the latter is the PMF at the last lambda state.
    analysis = AnalysisConfig(
        name="dG",
        stage="awh",
        command="gmx awh -s awh.tpr -f awh.edr -o {tmpdir}/awh",
        output="{tmpdir}/awh_t{end_time}.xvg",
    )
    result = run_analysis(analysis, build_simulations(load_campaign(path))[0])
    assert result.status == OK, result.error
    assert float(result.values["(nm, deg or lambda state)"]) == 20
    assert "PMF" in result.values

    # The free energy is PMF(first state, coupled) - PMF(last state), not the last value.
    delta = analysis.model_copy(update={"result": "first-last", "columns": ["dG_solv"]})
    simulation = build_simulations(load_campaign(path))[0]
    result = run_analysis(delta, simulation)
    assert result.status == OK, result.error
    xvg = simulation.spec.workdir / "rerun"
    xvg.mkdir()
    subprocess.run(
        ["gmx", "awh", "-quiet", "-s", "awh.tpr", "-f", "awh.edr", "-o", "rerun/awh"],
        cwd=simulation.spec.workdir, check=True, capture_output=True,
    )
    rows = [
        line.split() for line in (xvg / "awh_t1.xvg").read_text().splitlines()
        if line.strip() and line[0] not in "#@"
    ]
    expected = float(rows[0][1]) - float(rows[-1][1])
    assert float(result.values["dG_solv"]) == pytest.approx(expected, abs=1e-4)
