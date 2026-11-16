import re

import pytest

from ocmstoolkit.modules import Workflow as W
from ocmstoolkit.modules import Sra


@pytest.mark.parametrize("value, expected", [
    ("4G", 4000), ("15G", 15000), ("500M", 500), ("1.5G", 1500),
    ("2GB", 2000), ("512k", 1), (2000, 2000),
])
def test_memory_to_mb(value, expected):
    assert W.memory_to_mb(value) == expected


@pytest.mark.parametrize("value", ["lots", "", "0G", "-1G", True])
def test_memory_to_mb_rejects(value):
    with pytest.raises(ValueError):
        W.memory_to_mb(value)


def test_job_resources_defaults():
    assert W.job_resources({}) == {"mem_mb_per_cpu": 4000}
    assert W.job_threads(None) == 1


def test_job_resources_runtime():
    assert W.job_resources({"job_memory": "15G", "job_runtime": 90}) == \
        {"mem_mb_per_cpu": 15000, "runtime": 90}


def test_job_threads_rejects_zero():
    with pytest.raises(ValueError):
        W.job_threads({"job_threads": "0"})


@pytest.mark.parametrize("name, expected", [
    ("x.fastq.gz", "x"), ("x.fastq.1.gz", "x.fastq.1"), ("dir/y.fq.gz", "y"),
])
def test_fastqc_basename(name, expected):
    assert W.fastqc_basename(name) == expected


def test_list_input_files(tmp_path):
    for name in ["a.fastq.1.gz", "a.fastq.2.gz", "notes.txt"]:
        (tmp_path / name).write_text("x")
    found = W.list_input_files(str(tmp_path), re.compile(r"^(\S+)\.fastq.*gz$"))
    assert sorted(m.group(1) for m in found) == ["a", "a"]
    assert len(found) == 2


def test_list_input_files_errors(tmp_path):
    with pytest.raises(ValueError):
        W.list_input_files(str(tmp_path / "missing"), re.compile("x"))
    with pytest.raises(ValueError):
        W.list_input_files(str(tmp_path), re.compile("x"))


def test_wildcard_choices():
    pattern = W.wildcard_choices(["a", "a.1", "b+c"])
    assert re.fullmatch(pattern, "b+c")
    assert re.fullmatch(pattern, "a.1")
    assert not re.fullmatch(pattern, "a.2")
    assert not re.fullmatch(W.wildcard_choices([]), "")


ENA_SCRIPT = """#!/bin/bash

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/001/SRR1230001/SRR1230001_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/001/SRR1230001/SRR1230001_2.fastq.gz
"""


def test_parse_ena_script(tmp_path):
    script = tmp_path / "ena.sh"
    script.write_text(ENA_SCRIPT)
    parsed = Sra.parse_ena_script(str(script))
    assert list(parsed) == [".sentinel.SRR1230001_1.fastq.gz",
                            ".sentinel.SRR1230001_2.fastq.gz"]
    command, fq, srr = parsed[".sentinel.SRR1230001_1.fastq.gz"]
    assert command.startswith("wget -nc ftp://")
    assert fq == "SRR1230001_1.fastq.gz"
    assert srr == "SRR1230001"
    assert Sra.downloads(str(script))["SRR1230001_2.fastq.gz"].endswith(
        "SRR1230001_2.fastq.gz")


def test_parse_ena_script_rejects_duplicates_and_empty(tmp_path):
    script = tmp_path / "ena.sh"
    line = ENA_SCRIPT.splitlines()[2]
    script.write_text(f"{line}\n{line}\n")
    with pytest.raises(ValueError):
        Sra.parse_ena_script(str(script))
    script.write_text("#!/bin/bash\n\n")
    with pytest.raises(ValueError):
        Sra.parse_ena_script(str(script))
