# Run FastQC on each FASTQ file and summarise the results with MultiQC
"""
pipeline_fastqc.smk
===================

Runs FastQC on every FASTQ file in the input directory, then MultiQC to
summarise the results in ``multiqc_report.html``.

Usage
-----

From the directory in which the results should be written::

    ocms_toolkit fastqc config            # write pipeline_fastqc.yaml
    ocms_toolkit fastqc full --dry        # show the jobs that would run
    ocms_toolkit fastqc full --jobs 20    # submit the jobs to Slurm

Input
-----

FASTQ files in ``input_dir`` (default ``input.dir``), named:

``<sample>.fastq.gz``
    single-end reads;
``<sample>.fastq.1.gz``, ``<sample>.fastq.2.gz``
    paired-end reads.

Singleton files (``<sample>.fastq.3.gz``) are not checked.

Requirements
------------

FastQC and MultiQC must be on the ``PATH`` of the shell from which the
pipeline is launched (Slurm jobs inherit it), for example by loading the
FastQC and MultiQC modules before running the pipeline.

Output
------

``fastqc.dir/<track>.fastqc``
    FastQC results, extracted; ``<track>`` is ``<sample>``, ``<sample>.1``
    or ``<sample>.2``;
``fastqc.dir/<track>.fastqc.html``, ``fastqc.dir/<track>.fastqc.zip``
    FastQC report and archive;
``multiqc_report.html``, ``multiqc_data/``
    MultiQC summary of all FastQC results.
"""

import re

from snakemake.exceptions import WorkflowError

from ocmstoolkit.modules import Workflow as W

INPUT_DIR = str(config.get("input_dir") or "input.dir")
FASTQC = config.get("fastqc") or {}
MULTIQC = config.get("multiqc") or {}

# <sample>.fastq.gz, <sample>.fastq.1.gz and <sample>.fastq.2.gz
FASTQ_RE = re.compile(r"^(?P<sample>\S+)\.fastq(?P<read>\.[12])?\.gz$")

try:
    # track (<sample> or <sample>.<read>) -> input FASTQ file
    TRACKS = {
        match["sample"] + (match["read"] or ""): path
        for match, path in W.list_input_files(INPUT_DIR, FASTQ_RE).items()
    }
    FASTQC_THREADS = W.job_threads(FASTQC)
    FASTQC_RESOURCES = W.job_resources(FASTQC)
    MULTIQC_THREADS = W.job_threads(MULTIQC)
    MULTIQC_RESOURCES = W.job_resources(MULTIQC)
except ValueError as err:
    raise WorkflowError(f"pipeline_fastqc: {err}")

wildcard_constraints:
    track=W.wildcard_choices(TRACKS),


localrules: full


rule full:
    input:
        "multiqc_report.html",


rule fastqc:
    """Run FastQC on one FASTQ file."""
    input:
        lambda wildcards: TRACKS[wildcards.track],
    output:
        report=directory("fastqc.dir/{track}.fastqc"),
        html="fastqc.dir/{track}.fastqc.html",
        zip="fastqc.dir/{track}.fastqc.zip",
    params:
        # FastQC writes to a job-specific directory first; its output is then
        # renamed after the track.
        tmpdir="fastqc.dir/{track}.tmp",
        name=lambda wildcards: W.fastqc_basename(TRACKS[wildcards.track]),
    log:
        "logs/fastqc/{track}.log",
    threads: FASTQC_THREADS
    resources:
        **FASTQC_RESOURCES,
    shell:
        """
        rm -rf {params.tmpdir}
        mkdir -p {params.tmpdir}
        fastqc --extract --threads {threads} --outdir {params.tmpdir} \
            {input} > {log} 2>&1
        mv {params.tmpdir}/{params.name}_fastqc {output.report}
        mv {params.tmpdir}/{params.name}_fastqc.html {output.html}
        mv {params.tmpdir}/{params.name}_fastqc.zip {output.zip}
        rmdir {params.tmpdir}
        """


rule multiqc:
    """Summarise the FastQC results of all tracks."""
    input:
        expand("fastqc.dir/{track}.fastqc", track=TRACKS),
        expand("fastqc.dir/{track}.fastqc.zip", track=TRACKS),
    output:
        report="multiqc_report.html",
        data=directory("multiqc_data"),
    log:
        "multiqc.log",
    threads: MULTIQC_THREADS
    resources:
        **MULTIQC_RESOURCES,
    shell:
        # MultiQC's default output names are multiqc_report.html and
        # multiqc_data/; --force replaces the output of an earlier run.
        "multiqc --fullnames --force fastqc.dir > {log} 2>&1"
