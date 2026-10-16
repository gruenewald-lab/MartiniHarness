from pathlib import Path

import pytest

from martini_harness.config import DEFAULT_FORCEFIELD_DIR, DEFAULT_MDP_DIR
from martini_harness.mdp import (
    TemplateError,
    get_mdp_param,
    is_minimization,
    render_mdp,
    render_value,
    set_mdp_params,
)
from martini_harness.topology import read_moleculetypes, render_topology

TEMPLATE = """\
integrator   = sd   ; stochastic dynamics
nsteps       = 1000
ref_t        = 298.15
couple-moltype = TOCOUPLE
"""


def test_substitute_placeholder():
    text = render_mdp(TEMPLATE, substitute={"TOCOUPLE": "row42"})
    assert "couple-moltype = row42" in text


def test_missing_placeholder_is_an_error():
    with pytest.raises(TemplateError, match="NOPE"):
        render_mdp(TEMPLATE, substitute={"NOPE": "x"})


def test_params_replace_existing_and_append_new():
    text = set_mdp_params(TEMPLATE, {"nsteps": 10, "ref-t": 310, "gen_vel": True})
    lines = text.splitlines()
    assert any(line.split("=")[0].strip() == "nsteps" and line.endswith("= 10") for line in lines)
    assert any(line.startswith("ref-t") and line.endswith("= 310") for line in lines)
    assert "ref_t        = 298.15" not in text
    assert lines[-1].startswith("gen_vel") and lines[-1].endswith("= yes")
    assert "integrator   = sd   ; stochastic dynamics" in text


def test_render_value_recurses_and_reports_unknown_variables():
    variables = {"target": "row42", "n": 3}
    assert render_value({"a": ["{target}", 5], "b": "{n}x"}, variables) == {"a": ["row42", 5], "b": "3x"}
    with pytest.raises(TemplateError, match="solvent"):
        render_value("{solvent}", variables)


def test_read_moleculetypes():
    names = read_moleculetypes(DEFAULT_FORCEFIELD_DIR / "martini_v3.0.0_solvents_v1.itp")
    assert {"W", "OCO", "CLF"} <= set(names)
    assert read_moleculetypes(Path(__file__).parent.parent / "examples" / "row42.itp") == ["row42"]


def test_render_topology():
    text = render_topology([Path("/ff/a.itp")], "sys", [("row42", 1), ("W", 10)])
    assert text.splitlines()[0] == '#include "/ff/a.itp"'
    assert "[ molecules ]\nrow42      1\nW          10\n" in text


def test_is_minimization():
    assert get_mdp_param(TEMPLATE, "Couple_Moltype") == "TOCOUPLE"
    assert get_mdp_param(TEMPLATE, "tcoupl") is None
    assert is_minimization("integrator = steep ; comment\n")
    assert not is_minimization("nsteps = 10\n")  # GROMACS default integrator is md
    assert is_minimization((DEFAULT_MDP_DIR / "min.mdp").read_text())
    assert not any(is_minimization((DEFAULT_MDP_DIR / f).read_text()) for f in ("eq.mdp", "awh.mdp"))
