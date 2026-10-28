import shutil
import textwrap
from pathlib import Path

import pytest

from martini_harness.config import DEFAULT_FORCEFIELD_DIR, ConfigError, load_campaign
from martini_harness.simulation import build_simulations, run_campaign, split_count

from .conftest import BASIC_CAMPAIGN, EXAMPLES


def test_sample_campaign_loads():
    config = load_campaign(EXAMPLES / "sample_campaign.yaml")
    simulations = build_simulations(config)
    assert [sim.name for sim in simulations] == ["row42_in_W", "row42_in_wetOCO", "row42_in_CLF"]
    clf = simulations[2].spec
    assert (clf.n_solvent, clf.density) == (6000, 1480.0)
    assert simulations[1].spec.solvent_molecules == [("OCO", 7360), ("W", 640)]


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


MIXTURE_CAMPAIGN = BASIC_CAMPAIGN.replace("solvent: [W, CLF]", "solvent: [wetOCO]") + """
    solvents:
      wetOCO: {OCO: 920, W: 80}
"""


def test_solvent_mixture_divides_n_solvent(campaign_file, tmp_path, fake_backend):
    (sim,) = build_simulations(load_campaign(campaign_file(MIXTURE_CAMPAIGN)))
    assert sim.name == "row42_in_wetOCO"
    assert sim.spec.solvent_molecules == [("OCO", 92), ("W", 8)]
    run_campaign([sim], fake_backend)
    molecules = (tmp_path / "out" / sim.name / "system.top").read_text().split("[ molecules ]")[1]
    assert molecules.split() == ["row42", "1", "OCO", "92", "W", "8"]


@pytest.mark.parametrize(
    "amounts, n_solvent, expected",
    [
        ({"A": 1, "B": 1, "C": 1}, 100, [("A", 34), ("B", 33), ("C", 33)]),
        ({"A": 0.92, "B": 0.08}, 8000, [("A", 7360), ("B", 640)]),
        ({"A": 2, "B": 1}, 7, [("A", 5), ("B", 2)]),
    ],
)
def test_split_count_adds_up(amounts, n_solvent, expected):
    assert split_count(amounts, n_solvent) == expected


@pytest.mark.parametrize(
    "replace, message",
    [
        (("W: 80", "NOPE: 80"), "solvent 'NOPE' is not defined"),
        (("n_solvent: 100", "n_solvent: 5"), "n_solvent 5 leaves no W molecules"),
    ],
)
def test_solvent_mixture_build_errors(campaign_file, replace, message):
    config = load_campaign(campaign_file(MIXTURE_CAMPAIGN.replace(*replace)))
    with pytest.raises(ConfigError, match=message):
        build_simulations(config)


@pytest.mark.parametrize(
    "mixture, message",
    [("{OCO: 920, W: 0}", "amount of 'W' must be positive"), ("{}", "at least one molecule")],
)
def test_solvent_mixture_validation(campaign_file, mixture, message):
    body = MIXTURE_CAMPAIGN.replace("{OCO: 920, W: 80}", mixture)
    with pytest.raises(ConfigError, match=message):
        load_campaign(campaign_file(body))


def test_missing_mdp_placeholder_is_reported_before_running(campaign_file):
    body = BASIC_CAMPAIGN + """
    stages:
      - name: min
        mdp: min.mdp
        substitute: {TOCOUPLE: "{target}"}
    """
    with pytest.raises(ConfigError, match="placeholder 'TOCOUPLE' not found"):
        build_simulations(load_campaign(campaign_file(body)))


# --- imports and stage_overrides -------------------------------------------


def with_imports(imports: str) -> str:
    return f"    imports: {imports}\n" + BASIC_CAMPAIGN


@pytest.fixture
def write_file(tmp_path):
    """Write a (dedented) file relative to tmp_path and return its path."""

    def write(relative: str, body: str) -> Path:
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(textwrap.dedent(body))
        return path

    return write


def test_imports_example_matches_sample_campaign():
    sample = load_campaign(EXAMPLES / "sample_campaign.yaml")
    imported = load_campaign(EXAMPLES / "imports_campaign.yaml")
    assert imported.stages == sample.stages
    assert imported.molecules == sample.molecules
    assert imported.matrix == sample.matrix


def test_imports_deep_merge_mappings_and_replace_lists(campaign_file, write_file):
    write_file(
        "base.yaml",
        """
        title: base
        simulation: {n_target: 3, density: 500}
        topology:
          includes: [martini_v3.0.0.itp, martini_v3.0.0_solvents_v1.itp]
        stages:
          - {name: min, mdp: min.mdp}
          - {name: eq, mdp: eq.mdp}
        """,
    )
    write_file("later.yaml", "simulation: {density: 600}\ntitle: later\n")
    config = load_campaign(
        campaign_file(
            with_imports("[base.yaml, later.yaml]")
            + """
    topology:
      includes: [martini_v3.0.0.itp]
    stages:
      - {name: only, mdp: min.mdp}
    """
        )
    )
    assert config.title == "test"  # importing file wins
    # simulation merges key by key: base, later import, then campaign file.
    assert (config.simulation.n_target, config.simulation.density) == (3, 600)
    assert config.simulation.n_solvent == 100
    assert [stage.name for stage in config.stages] == ["only"]
    assert [path.name for path in config.topology.includes] == ["martini_v3.0.0.itp"]


def test_imported_matrix_axes_merge_by_key(campaign_file, write_file):
    write_file("matrix.yaml", "matrix:\n  solvent: [W, OCO, CLF]\n  temperature: [300, 310]\n")
    config = load_campaign(campaign_file(with_imports("matrix.yaml")))
    assert config.matrix.axes == {
        "solvent": ["W", "CLF"],
        "temperature": [300, 310],
        "target": ["row42"],
    }


def test_imported_paths_resolve_relative_to_imported_file(tmp_path, write_file):
    write_file("lib/mdps/min.mdp", "integrator = steep\n")
    shutil.copy(EXAMPLES / "row42.itp", write_file("lib/deep/row42.itp", ""))
    write_file(
        "lib/deep/molecules.yaml",
        """
        molecules:
          row42: row42.itp
        """,
    )
    write_file(
        "lib/pipeline.yaml",
        """
        imports: [deep/molecules.yaml]
        mdp_dir: mdps
        output_dir: out
        stages:
          - {name: min, mdp: min.mdp}
        """,
    )
    campaign = write_file(
        "campaigns/campaign.yaml",
        """
        imports: [../lib/pipeline.yaml]
        matrix:
          target: [row42]
          solvent: [W]
        """,
    )
    config = load_campaign(campaign)
    assert config.molecules == {"row42": tmp_path / "lib/deep/row42.itp"}
    assert config.mdp_dir == tmp_path / "lib/mdps"
    assert config.output_dir == tmp_path / "lib/out"
    assert config.mdp_path(config.stages[0]) == tmp_path / "lib/mdps/min.mdp"
    assert [sim.name for sim in build_simulations(config)] == ["row42_in_W"]


def test_importing_file_paths_override_imported_paths(campaign_file, write_file, tmp_path):
    write_file("lib/settings.yaml", "output_dir: lib_out\nmdp_dir: mdps\n")
    config = load_campaign(campaign_file(with_imports("[lib/settings.yaml]")))
    assert config.output_dir == tmp_path / "out"
    assert config.mdp_dir == tmp_path / "lib/mdps"


def test_diamond_imports_are_allowed(campaign_file, write_file):
    write_file("common.yaml", "simulation: {density: 750}\n")
    write_file("a.yaml", "imports: [common.yaml]\n")
    write_file("b.yaml", "imports: [common.yaml]\n")
    config = load_campaign(campaign_file(with_imports("[a.yaml, b.yaml]")))
    assert config.simulation.density == 750


def test_import_cycle_is_reported(campaign_file, write_file):
    write_file("a.yaml", "imports: [b.yaml]\n")
    write_file("b.yaml", "imports: [a.yaml]\n")
    with pytest.raises(ConfigError, match=r"import cycle: .*a\.yaml -> .*b\.yaml -> .*a\.yaml"):
        load_campaign(campaign_file(with_imports("[a.yaml]")))


def test_self_import_is_a_cycle(campaign_file):
    with pytest.raises(ConfigError, match="import cycle"):
        load_campaign(campaign_file(with_imports("[campaign.yaml]")))


def test_missing_import_is_reported(campaign_file, write_file):
    write_file("a.yaml", "imports: [nope.yaml]\n")
    with pytest.raises(ConfigError, match=r"a\.yaml: imported file not found: nope\.yaml"):
        load_campaign(campaign_file(with_imports("[a.yaml]")))


@pytest.mark.parametrize(
    "body, message",
    [
        ("- just\n- a list\n", "top level must be a mapping"),
        ("stages: [unclosed\n", "cannot read"),
    ],
)
def test_bad_imported_file_is_named(campaign_file, write_file, body, message):
    write_file("bad.yaml", body)
    with pytest.raises(ConfigError, match=rf"bad\.yaml.*{message}|{message}.*bad\.yaml"):
        load_campaign(campaign_file(with_imports("[bad.yaml]")))


def test_imports_must_be_file_names(campaign_file):
    with pytest.raises(ConfigError, match="imports must be a file name or a list"):
        load_campaign(campaign_file(with_imports("{a: b}")))


def test_stage_overrides_merge_into_imported_stage(campaign_file, write_file):
    write_file(
        "pipeline.yaml",
        """
        stages:
          - name: min
            mdp: min.mdp
          - name: awh
            mdp: awh.mdp
            maxwarn: 1
            substitute: {TOCOUPLE: "{target}"}
            mdp_params: {nsteps: 1000, dt: 0.02}
        """,
    )
    config = load_campaign(
        campaign_file(
            with_imports("[pipeline.yaml]")
            + """
    stage_overrides:
      awh:
        mdp_params: {nsteps: 500}
        mdrun: {nt: 8}
        maxwarn: 3
    """
        )
    )
    awh = config.stages[1]
    assert awh.mdp_params == {"nsteps": 500, "dt": 0.02}
    assert awh.mdrun == {"nt": 8}
    assert awh.maxwarn == 3
    assert awh.substitute == {"TOCOUPLE": "{target}"}
    assert config.stages[0].mdp_params == {}


def test_stage_overrides_apply_to_default_pipeline(campaign_file):
    body = BASIC_CAMPAIGN + "    stage_overrides:\n      eq: {maxwarn: 5}\n"
    config = load_campaign(campaign_file(body))
    assert [(stage.name, stage.maxwarn) for stage in config.stages] == [
        ("min", 0),
        ("eq", 5),
        ("awh", 0),
    ]


@pytest.mark.parametrize(
    "override, message",
    [
        ("nope: {maxwarn: 1}", "no stage named 'nope'"),
        ("eq: {name: other}", "cannot be renamed"),
        ("eq: {maxwran: 1}", r"stage_overrides\.eq\.maxwran: Extra inputs"),
    ],
)
def test_stage_overrides_errors(campaign_file, override, message):
    body = BASIC_CAMPAIGN + f"    stage_overrides:\n      {override}\n"
    with pytest.raises(ConfigError, match=message):
        load_campaign(campaign_file(body))
