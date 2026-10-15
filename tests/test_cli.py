from martini_harness import cli

from .conftest import BASIC_CAMPAIGN, EXAMPLES, FakeBackend


def test_dry_run_lists_simulations(campaign_file, capsys):
    assert cli.main(["run", str(campaign_file(BASIC_CAMPAIGN)), "--dry-run"]) == 0
    out = capsys.readouterr().out
    assert "row42_in_W" in out and "row42_in_CLF" in out
    assert "2/2 simulation(s) have steps to run" in out


def test_only_restricts_simulations(campaign_file, capsys):
    path = str(campaign_file(BASIC_CAMPAIGN))
    assert cli.main(["status", path, "--only", "row42_in_CLF"]) == 0
    out = capsys.readouterr().out
    assert "row42_in_CLF" in out and "row42_in_W" not in out
    assert cli.main(["status", path, "--only", "nope"]) == 2


def test_invalid_config_exits_with_2(campaign_file):
    assert cli.main(["run", str(campaign_file("matrix: {}\n")), "--dry-run"]) == 2


def test_run_exit_code_reflects_failures(campaign_file, monkeypatch, tmp_path):
    path = str(campaign_file(BASIC_CAMPAIGN))
    monkeypatch.setattr(cli, "run_campaign", _with_backend(FakeBackend(fail_on={"awh"})))
    assert cli.main(["run", path]) == 1
    monkeypatch.setattr(cli, "run_campaign", _with_backend(FakeBackend()))
    assert cli.main(["run", path]) == 0
    assert cli.main(["run", path, "-o", str(tmp_path / "elsewhere")]) == 0
    assert (tmp_path / "elsewhere" / "row42_in_W" / "awh.gro").exists()


def test_list_molecules(capsys):
    assert cli.main(["list-molecules", str(EXAMPLES / "row42.itp")]) == 0
    assert capsys.readouterr().out.startswith("row42\t")


def _with_backend(backend):
    from martini_harness.simulation import run_campaign

    return lambda simulations, force=False: run_campaign(simulations, backend, force=force)
