import logging
import sys
import textwrap

import pytest

from martini_harness import cli
from martini_harness.analysis import FAILED, NOT_FINISHED, OK, run_analysis
from martini_harness.config import AnalysisConfig, ConfigError, load_campaign
from martini_harness.simulation import build_simulations, run_campaign

from .conftest import BASIC_CAMPAIGN, FakeBackend

# Stands in for `gmx awh`: writes <prefix>_t<time>.xvg for every --times value,
# each ending in the line "20.0000 <time / 10>".
FAKE_AWH = '''
import argparse, pathlib, sys
parser = argparse.ArgumentParser()
parser.add_argument("-o")
parser.add_argument("--times", type=float, nargs="+", default=[0.0, 10.0, 20.0])
parser.add_argument("--fail", action="store_true")
args = parser.parse_args()
print("selection:", sys.stdin.read().strip())
if args.fail:
    sys.exit(3)
for t in args.times:
    path = pathlib.Path(f"{args.o}_t{t:g}.xvg")
    path.write_text(
        '@    xaxis  label "(nm, deg or lambda state)"\\n'
        '@ s0 legend "PMF"\\n'
        + "".join(f"{x:8.4f}  {t / 10 * x / 20:g}\\n" for x in range(21))
    )
'''

AWH_CAMPAIGN = BASIC_CAMPAIGN + """
    stages:
      - {name: min, mdp: min.mdp}
      - name: awh
        mdp: awh.mdp
        substitute: {TOCOUPLE: "{target}"}
        mdp_params: {nsteps: 1000, dt: 0.02}
"""


@pytest.fixture
def fake_awh(tmp_path):
    script = tmp_path / "fake_awh.py"
    script.write_text(FAKE_AWH)
    return f"{sys.executable} {script}"


def _analysis_yaml(command: str, extra: str = "") -> str:
    return textwrap.indent(
        textwrap.dedent(f"""\
            analyses:
              - name: dG
                command: {command} -o {{tmpdir}}/awh{extra}
                output: "{{tmpdir}}/awh_t{{end_time}}.xvg"
        """),
        "    ",
    )


def test_awh_analysis_reads_last_line_at_end_time(campaign_file, fake_awh, tmp_path):
    path = campaign_file(AWH_CAMPAIGN + _analysis_yaml(fake_awh))
    config = load_campaign(path)
    simulations = build_simulations(config)
    run_campaign(simulations[:1], FakeBackend())

    analysis = config.analyses[0]
    assert analysis.stage == "awh"  # defaults to the last stage
    result = run_analysis(analysis, simulations[0])
    assert result.status == OK, result.error
    # end time 1000 * 0.02 = 20 ps -> awh_t20.xvg, not the earlier awh_t10.xvg
    assert result.values == {"(nm, deg or lambda state)": "20.0000", "PMF": "2"}

    workdir = tmp_path / "out" / "row42_in_W"
    assert not list(workdir.glob("analysis_tmp_*"))  # scratch directory removed
    assert (workdir / "dG.analysis.out").is_file()

    assert run_analysis(analysis, simulations[1]).status == NOT_FINISHED


def test_awh_analysis_failures(campaign_file, fake_awh):
    path = campaign_file(AWH_CAMPAIGN + _analysis_yaml(fake_awh, " --times 0 10"))
    config = load_campaign(path)
    simulations = build_simulations(config)
    run_campaign(simulations[:1], FakeBackend())

    result = run_analysis(config.analyses[0], simulations[0])
    assert result.status == FAILED
    assert "awh_t20.xvg not written" in result.error

    failing = config.analyses[0].model_copy(update={"command": f"{fake_awh} --fail"})
    result = run_analysis(failing, simulations[0])
    assert (result.status, "exited with code 3" in result.error) == (FAILED, True)

    missing = config.analyses[0].model_copy(update={"command": "no-such-program-xyz"})
    assert "cannot run no-such-program-xyz" in run_analysis(missing, simulations[0]).error


def test_analysis_stdin_columns_and_variables(campaign_file, fake_awh, tmp_path):
    path = campaign_file(AWH_CAMPAIGN)
    simulations = build_simulations(load_campaign(path))
    run_campaign(simulations[:1], FakeBackend())

    analysis = AnalysisConfig(
        name="custom",
        stage="awh",
        command=fake_awh.split() + ["-o", "{tmpdir}/{target}", "--times", "20"],
        output="{tmpdir}/{target}_t20.xvg",
        stdin="{solvent}\n",
        columns=["lambda_state", "dG"],
    )
    result = run_analysis(analysis, simulations[0])
    assert result.values == {"lambda_state": "20.0000", "dG": "2"}
    assert "selection: W" in (tmp_path / "out" / "row42_in_W" / "custom.analysis.out").read_text()

    unknown = analysis.model_copy(update={"output": "{nope}.xvg"})
    assert "unknown variable(s) nope" in run_analysis(unknown, simulations[0]).error


def test_analyze_cli_writes_table(campaign_file, fake_awh, tmp_path, caplog):
    path = campaign_file(AWH_CAMPAIGN + _analysis_yaml(fake_awh))
    run_campaign(build_simulations(load_campaign(path))[:1], FakeBackend())

    caplog.set_level(logging.INFO)
    assert cli.main(["analyze", str(path)]) == 0
    lines = (tmp_path / "out" / "dG.tsv").read_text().splitlines()
    assert lines[0].split("\t") == [
        "name", "target", "solvent", "n_target", "n_solvent", "density",
        "status", "(nm, deg or lambda state)", "PMF", "error",
    ]
    assert lines[1].split("\t") == ["row42_in_W", "row42", "W", "1", "100", "900.0", "ok", "20.0000", "2", ""]
    assert lines[2].split("\t")[6:] == ["not-finished", "", "", ""]
    assert "dG: 1 ok, 1 with awh not finished, 0 failed" in caplog.text

    assert cli.main(["analyze", str(path), "--analysis", "nope"]) == 2
    assert cli.main(["analyze", str(campaign_file(AWH_CAMPAIGN))]) == 2  # no analyses defined


def test_analyze_cli_exit_code_on_failure(campaign_file, fake_awh, tmp_path):
    path = campaign_file(AWH_CAMPAIGN + _analysis_yaml(fake_awh, " --fail"))
    run_campaign(build_simulations(load_campaign(path))[:1], FakeBackend())
    assert cli.main(["analyze", str(path), "--only", "row42_in_W"]) == 1
    row = (tmp_path / "out" / "dG.tsv").read_text().splitlines()[1].split("\t")
    assert row[6] == "failed" and "exited with code 3" in row[-1]


@pytest.mark.parametrize(
    "analyses, message",
    [
        ("[{name: a, stage: eq, command: x, output: y}]", "no stage named 'eq'"),
        ("[{name: a, command: x, output: y}, {name: a, command: z, output: y}]", "duplicate analysis"),
        ("[{name: a/b, command: x, output: y}]", "analysis name"),
    ],
)
def test_invalid_analyses(campaign_file, analyses, message):
    path = campaign_file(AWH_CAMPAIGN + f"    analyses: {analyses}\n")
    with pytest.raises(ConfigError, match=message):
        load_campaign(path)


def test_end_time_uses_gromacs_number_format(campaign_file):
    body = AWH_CAMPAIGN.replace("nsteps: 1000, dt: 0.02", "nsteps: 500000, dt: 0.02")
    path = campaign_file(body)
    simulations = build_simulations(load_campaign(path))
    run_campaign(simulations[:1], FakeBackend())
    analysis = AnalysisConfig(
        name="t", command=["echo", "{end_time}"], output="t.analysis.out", stage="awh"
    )
    assert run_analysis(analysis, simulations[0]).values == {"x": "10000"}


@pytest.mark.parametrize(
    "result, column, expected",
    [("first-last", 1, "-2"), ("last-first", 1, "2"), ("first-last", 0, "-20")],
)
def test_signed_difference_of_first_and_last_line(campaign_file, fake_awh, result, column, expected):
    path = campaign_file(AWH_CAMPAIGN)
    simulations = build_simulations(load_campaign(path))
    run_campaign(simulations[:1], FakeBackend())
    analysis = AnalysisConfig(
        name="dG",
        stage="awh",
        command=fake_awh.split() + ["-o", "{tmpdir}/awh"],
        output="{tmpdir}/awh_t{end_time}.xvg",  # PMF 0 ... 2 over lambda state 0 ... 20
        result=result,
        column=column,
    )
    assert run_analysis(analysis, simulations[0]).values == {"dG": expected}

    named = analysis.model_copy(update={"columns": ["dG_solv"]})
    assert run_analysis(named, simulations[0]).values == {"dG_solv": expected}

    missing = analysis.model_copy(update={"column": 9})
    assert "no numeric column 9" in run_analysis(missing, simulations[0]).error


def test_difference_needs_one_column_name():
    with pytest.raises(ValueError, match="single column name"):
        AnalysisConfig(name="a", command="x", output="y", result="first-last", columns=["a", "b"])
