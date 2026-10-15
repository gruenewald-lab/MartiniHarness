import pytest

from martini_harness.config import DEFAULT_FORCEFIELD_DIR, ConfigError, load_campaign
from martini_harness.simulation import build_simulations

from .conftest import BASIC_CAMPAIGN, EXAMPLES


def test_sample_campaign_loads():
    config = load_campaign(EXAMPLES / "sample_campaign.yaml")
    simulations = build_simulations(config)
    assert [sim.name for sim in simulations] == ["row42_in_W", "row42_in_OCO", "row42_in_CLF"]
    clf = simulations[2].spec
    assert (clf.n_solvent, clf.density) == (6000, 1480.0)


def test_paths_resolve_relative_to_campaign_file(campaign_file, tmp_path):
    config = load_campaign(campaign_file(BASIC_CAMPAIGN))
    assert config.output_dir == tmp_path / "out"
    assert config.molecules["row42"] == tmp_path / "row42.itp"
    assert config.topology.includes[0] == DEFAULT_FORCEFIELD_DIR / "martini_v3.0.0.itp"


def test_unknown_keys_are_rejected(campaign_file):
    with pytest.raises(ConfigError, match="n_solvnet"):
        load_campaign(campaign_file(BASIC_CAMPAIGN.replace("n_solvent", "n_solvnet")))


def test_stage_input_must_reference_earlier_stage(campaign_file):
    body = BASIC_CAMPAIGN + """
    stages:
      - {name: min, mdp: min.mdp, input: eq}
      - {name: eq, mdp: eq.mdp}
    """
    with pytest.raises(ConfigError, match="earlier stage"):
        load_campaign(campaign_file(body))


@pytest.mark.parametrize("name", ["coords", "has space"])
def test_invalid_stage_names(campaign_file, name):
    body = BASIC_CAMPAIGN + f"""
    stages:
      - {{name: "{name}", mdp: min.mdp}}
    """
    with pytest.raises(ConfigError):
        load_campaign(campaign_file(body))


def test_matrix_overrides_simulation_parameters(campaign_file):
    body = BASIC_CAMPAIGN + """
      include:
        - {solvent: CLF, n_solvent: 42, density: 1480}
    """
    sims = {sim.name: sim.spec for sim in build_simulations(load_campaign(campaign_file(body)))}
    assert sims["row42_in_CLF"].n_solvent == 42
    assert sims["row42_in_W"].n_solvent == 100


@pytest.mark.parametrize(
    "replace, message",
    [
        (("solvent: [W, CLF]", "solvent: [W, NOPE]"), "solvent 'NOPE' is not defined"),
        (("target: [row42]", "target: [row99]"), "target 'row99' is not defined"),
        (("n_solvent: 100", 'name: "{target}"'), "duplicate simulation name"),
        (("n_solvent: 100", 'name: "{temperature}"'), "unknown variable"),
        (("row42: row42.itp", "row42: missing.itp"), "ITP file not found"),
    ],
)
def test_build_simulations_reports_errors(campaign_file, replace, message):
    config = load_campaign(campaign_file(BASIC_CAMPAIGN.replace(*replace)))
    with pytest.raises(ConfigError, match=message):
        build_simulations(config)


def test_missing_mdp_placeholder_is_reported_before_running(campaign_file):
    body = BASIC_CAMPAIGN + """
    stages:
      - name: min
        mdp: min.mdp
        substitute: {TOCOUPLE: "{target}"}
    """
    with pytest.raises(ConfigError, match="placeholder 'TOCOUPLE' not found"):
        build_simulations(load_campaign(campaign_file(body)))
