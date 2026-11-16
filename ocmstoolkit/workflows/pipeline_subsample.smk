# Randomly subsample FASTQ files to a fixed number of reads with seqtk
"""
pipeline_subsample.smk
======================

Uses ``seqtk sample`` to randomly subsample every FASTQ file in the input
directory to a fixed number of reads (``seqtk.depth``).

The same random seed is used for every file. Each file is subsampled
independently, but the two files of a read pair are subsampled identically,
so read pairs are preserved.

Usage
-----

From the directory in which the results should be written::

    ocms_toolkit subsample config            # write pipeline_subsample.yaml
    ocms_toolkit subsample full --dry        # show the jobs that would run
    ocms_toolkit subsample full --jobs 20    # submit the jobs to Slurm

Input
-----

Files matching ``*.fastq*.gz`` in ``input_dir`` (default ``input.dir``):
single-end ``<sample>.fastq.gz`` or paired-end ``<sample>.fastq.1.gz`` and
``<sample>.fastq.2.gz``.

Requirements
------------

seqtk must be on the ``PATH`` of the shell from which the pipeline is
launched (Slurm jobs inherit it), for example after
``module load seqtk/1.4-GCC-12.2.0``.

Output
------

``subsampled.dir/<sample>_subsampled.fastq[.N].gz``
    one subsampled file for each input file.
"""

import re

from snakemake.exceptions import WorkflowError

from ocmstoolkit.modules import Workflow as W

INPUT_DIR = str(config.get("input_dir") or "input.dir")
SEQTK = config.get("seqtk") or {}

# <prefix>.<suffix>, where the suffix starts with "fastq" and ends with "gz"
FASTQ_RE = re.compile(r"^(?P<prefix>\S+)\.(?P<suffix>fastq\S*gz)$")

try:
    # output file name -> input FASTQ file
    SUBSAMPLED = {
        f"{match['prefix']}_subsampled.{match['suffix']}": path
        for match, path in W.list_input_files(INPUT_DIR, FASTQ_RE).items()
    }
    DEPTH = SEQTK.get("depth")
    if DEPTH in (None, ""):
        raise ValueError("set seqtk.depth, the number of reads to keep.")
    SEED = int(SEQTK.get("seed", 100))
    SEQTK_THREADS = W.job_threads(SEQTK)
    SEQTK_RESOURCES = W.job_resources(SEQTK)
except ValueError as err:
    raise WorkflowError(f"pipeline_subsample: {err}")

wildcard_constraints:
    fastq=W.wildcard_choices(SUBSAMPLED),


localrules: full


rule full:
    input:
        expand("subsampled.dir/{fastq}", fastq=SUBSAMPLED),


rule subsample_fastq:
    """Subsample one FASTQ file with a fixed seed."""
    input:
        lambda wildcards: SUBSAMPLED[wildcards.fastq],
    output:
        "subsampled.dir/{fastq}",
    params:
        depth=DEPTH,
        seed=SEED,
    log:
        # Kept out of subsampled.dir, which is often the input of the next
        # pipeline.
        "logs/subsample_fastq/{fastq}.log",
    threads: SEQTK_THREADS
    resources:
        **SEQTK_RESOURCES,
    shell:
        """
        seqtk sample -s{params.seed} {input} {params.depth} 2> {log} \
            | gzip > {output}
        """
