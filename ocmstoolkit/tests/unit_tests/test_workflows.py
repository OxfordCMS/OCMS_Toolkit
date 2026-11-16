"""Dry runs of every pipeline on the test FASTQ files.

These check that each Snakefile parses and that its DAG resolves. They need
Snakemake, but not Slurm or the tools the pipelines run.
"""

import shutil
from pathlib import Path
from unittest import mock

import pytest

from ocmstoolkit import cli

FASTQ_DIR = Path(__file__).resolve().parent.parent / "fastqs"

ENA_SCRIPT = (
    "wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/001/SRR1230001/"
    "SRR1230001_1.fastq.gz\n"
)


@pytest.fixture
def rundir(tmp_path, monkeypatch):
    shutil.copytree(FASTQ_DIR, tmp_path / "input.dir")
    (tmp_path / "ena.sh").write_text(ENA_SCRIPT)
    monkeypatch.chdir(tmp_path)
    return tmp_path


@pytest.mark.parametrize("pipeline", ["fastqc", "subsample", "sra"])
def test_dry_run(rundir, pipeline, capfd):
    with mock.patch.object(cli, "default_slurm_account", return_value=None):
        assert cli.main([pipeline, "config"]) == 0
    configfile = rundir / f"pipeline_{pipeline}.yaml"
    if pipeline == "sra":
        configfile.write_text(configfile.read_text().replace(
            'ena_script: ""', 'ena_script: "ena.sh"'))

    assert cli.main([pipeline, "full", "--local", "--dry"]) == 0
    out, err = capfd.readouterr()
    assert "rule full" in out + err
