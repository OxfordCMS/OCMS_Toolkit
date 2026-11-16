"""Unit tests for the ocms_toolkit command-line interface.

These tests need neither Slurm nor the pipeline tools.
"""

import sys
from pathlib import Path
from unittest import mock

import pytest

from ocmstoolkit import cli


@pytest.fixture
def pipelines():
    return cli.discover_pipelines()


def test_discovers_all_pipelines(pipelines):
    assert set(pipelines) == {"fastqc", "sra", "subsample"}


def test_every_pipeline_has_a_config_template(pipelines):
    for name in pipelines:
        assert cli.config_template(name).is_file(), name


def test_every_template_has_the_account_placeholder(pipelines):
    for name in pipelines:
        config = cli.load_config(cli.config_template(name))
        assert config["resources"]["slurm_account"] == cli.ACCOUNT_PLACEHOLDER


def test_discovers_scripts():
    assert set(cli.discover_scripts()) == {
        "fastq2filteredfastq", "new_project", "rename_and_link"}


def test_pipeline_summary_is_first_comment(pipelines):
    assert cli.pipeline_summary(pipelines["fastqc"]).startswith("Run FastQC")


def test_write_config_fills_account(tmp_path):
    dest = tmp_path / "pipeline_fastqc.yaml"
    assert cli.write_config("fastqc", dest, account="kir.prj")
    config = cli.load_config(dest)
    assert config["resources"]["slurm_account"] == "kir.prj"
    # comments are kept
    assert dest.read_text().startswith("# Configuration")


def test_write_config_keeps_placeholder_without_account(tmp_path):
    dest = tmp_path / "pipeline_fastqc.yaml"
    assert not cli.write_config("fastqc", dest, account=None)
    assert cli.load_config(dest)["resources"]["slurm_account"] == \
        cli.ACCOUNT_PLACEHOLDER


def test_write_config_does_not_overwrite(tmp_path):
    dest = tmp_path / "pipeline_fastqc.yaml"
    cli.write_config("fastqc", dest)
    with pytest.raises(cli.CLIError):
        cli.write_config("fastqc", dest)
    cli.write_config("fastqc", dest, force=True)


def test_write_config_unknown_pipeline(tmp_path):
    with pytest.raises(cli.CLIError):
        cli.write_config("nonexistent", tmp_path / "x.yaml")


def test_validate_slurm_rejects_placeholder():
    config = {"resources": {"slurm_account": cli.ACCOUNT_PLACEHOLDER}}
    with pytest.raises(cli.CLIError):
        cli.validate_slurm(config, "user")


def test_validate_slurm_rejects_unknown_account():
    config = {"resources": {"slurm_account": "other.prj"}}
    with mock.patch.object(cli, "slurm_accounts", return_value={"kir.prj"}):
        with pytest.raises(cli.CLIError):
            cli.validate_slurm(config, "user")


def test_validate_slurm_accepts_known_account():
    config = {"resources": {"slurm_account": "kir.prj"}}
    with mock.patch.object(cli, "slurm_accounts", return_value={"kir.prj"}):
        assert cli.validate_slurm(config, "user") == []


def test_validate_slurm_warns_if_unverifiable():
    config = {"resources": {"slurm_account": "kir.prj"}}
    with mock.patch.object(cli, "slurm_accounts", return_value=None):
        assert len(cli.validate_slurm(config, "user")) == 1


def test_format_default_resources_quotes_text():
    assert cli.format_default_resources(
        {"slurm_account": "kir.prj", "x": 3, "empty": None}) == \
        ["slurm_account='kir.prj'", "x=3"]


def test_format_default_resources_rejects_lists():
    with pytest.raises(cli.CLIError):
        cli.format_default_resources({"slurm_account": ["a"]})


CONFIG = {"resources": {"slurm_account": "kir.prj", "slurm_partition": "short"}}


def test_slurm_command_does_not_pass_cores(tmp_path):
    # --cores would cap the threads of every Slurm job
    cmd = cli.build_snakemake_command(
        Path("p.smk"), tmp_path / "p.yaml", CONFIG, jobs=20)
    assert "--cores" not in cmd
    assert cmd[cmd.index("--executor") + 1] == "slurm"
    assert cmd[cmd.index("--jobs") + 1] == "20"
    assert "slurm_account='kir.prj'" in cmd
    assert "--local-cores" not in cmd


def test_slurm_command_cores_limit_local_rules(tmp_path):
    cmd = cli.build_snakemake_command(
        Path("p.smk"), tmp_path / "p.yaml", CONFIG, jobs=20, cores=4)
    assert cmd[cmd.index("--local-cores") + 1] == "4"
    assert "--cores" not in cmd


def test_local_command(tmp_path):
    cmd = cli.build_snakemake_command(
        Path("p.smk"), tmp_path / "p.yaml", CONFIG, jobs=20, local=True)
    assert "--executor" not in cmd
    assert cmd[cmd.index("--cores") + 1] == "all"


def test_extra_arguments_are_appended(tmp_path):
    cmd = cli.build_snakemake_command(
        Path("p.smk"), tmp_path / "p.yaml", CONFIG, jobs=1, dry_run=True,
        extra=["--unlock"])
    assert cmd[-2:] == ["--dry-run", "--unlock"]


def test_config_task_writes_file(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    with mock.patch.object(cli, "default_slurm_account", return_value=None):
        assert cli.main(["subsample", "config"]) == 0
    assert (tmp_path / "pipeline_subsample.yaml").is_file()


def test_full_without_config_fails(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    assert cli.main(["fastqc", "full", "--local"]) == 1
    assert "pipeline_fastqc.yaml not found" in capsys.readouterr().err


def test_old_cgat_syntax_gets_a_hint(capsys):
    assert cli.main(["fastqc", "make", "full"]) == 2
    assert "full --jobs" in capsys.readouterr().err


def test_extra_arguments_only_for_full(capsys):
    with pytest.raises(SystemExit):
        cli.main(["fastqc", "config", "--", "--unlock"])


def test_script_dispatch_passes_arguments(monkeypatch):
    received = {}

    def fake_main(argv):
        received["argv"] = list(argv)
        received["sys.argv"] = list(sys.argv)
        return 0

    module = mock.Mock(main=fake_main)
    with mock.patch.object(cli.importlib, "import_module", return_value=module):
        assert cli.main(["rename_and_link", "-i", "raw"]) == 0
    assert received["argv"] == ["ocms_toolkit rename_and_link", "-i", "raw"]
    assert received["sys.argv"] == received["argv"]
