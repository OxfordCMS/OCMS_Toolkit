'''
Workflow.py - helper functions shared by the ocms_toolkit Snakemake workflows
=============================================================================

The functions in this module do not depend on Snakemake, so that they can be
unit tested on their own. The Snakefiles in ``ocmstoolkit/workflows`` convert
the ``ValueError`` raised here into a Snakemake ``WorkflowError``.

Job resources
-------------

Each tool section of a pipeline configuration file may contain:

``job_threads``
    Number of CPUs requested for each job (default 1).

``job_memory``
    Memory *per CPU*, as in cgatcore, e.g. ``4G`` or ``500M``. A job receives
    ``job_memory`` x ``job_threads`` in total. It is passed to Slurm as
    ``--mem-per-cpu``.

``job_runtime``
    Optional wall time in minutes. If it is not set, the default wall time
    of the Slurm partition applies.
'''

import os
import re

#: Memory per CPU used when a tool section does not set job_memory. This is
#: the cgatcore default.
DEFAULT_JOB_MEMORY = "4G"

_MEMORY_RE = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*([KMGT]?)B?\s*$", re.IGNORECASE)
_MEMORY_FACTOR_MB = {"K": 1 / 1000, "": 1, "M": 1, "G": 1000, "T": 1000 ** 2}


def memory_to_mb(value):
    '''Convert a memory specification such as ``15G`` to megabytes.

    The conversion follows cgatcore (1G = 1000M). A bare number is taken to
    be megabytes.
    '''
    if isinstance(value, bool):
        raise ValueError(f"invalid memory specification: {value!r}")
    if isinstance(value, (int, float)):
        megabytes = value
    else:
        match = _MEMORY_RE.match(str(value))
        if not match:
            raise ValueError(
                f"invalid memory specification: {value!r}"
                " (expected a value such as 4G or 500M)")
        number, unit = match.groups()
        megabytes = float(number) * _MEMORY_FACTOR_MB[unit.upper()]
    megabytes = int(-(-megabytes // 1))  # round up
    if megabytes <= 0:
        raise ValueError(f"memory must be greater than zero: {value!r}")
    return megabytes


def job_threads(section):
    '''Number of threads for a job, from ``job_threads`` (default 1).'''
    section = section or {}
    threads = section.get("job_threads") or 1
    try:
        threads = int(threads)
    except (TypeError, ValueError):
        raise ValueError(f"job_threads must be a whole number: {threads!r}")
    if threads < 1:
        raise ValueError(f"job_threads must be at least 1: {threads!r}")
    return threads


def job_resources(section):
    '''Snakemake resources for a tool section of the configuration.

    Returns a dictionary with ``mem_mb_per_cpu`` and, if ``job_runtime`` is
    set, ``runtime`` (minutes). The Slurm account and partition are not
    included: the ocms_toolkit command applies them to every job with
    ``--default-resources``.
    '''
    section = section or {}
    resources = {
        "mem_mb_per_cpu": memory_to_mb(
            section.get("job_memory") or DEFAULT_JOB_MEMORY),
    }
    runtime = section.get("job_runtime")
    if runtime not in (None, ""):
        try:
            runtime = int(runtime)
        except (TypeError, ValueError):
            raise ValueError(
                f"job_runtime must be a whole number of minutes: {runtime!r}")
        if runtime < 1:
            raise ValueError(f"job_runtime must be at least 1: {runtime!r}")
        resources["runtime"] = runtime
    return resources


def list_input_files(indir, pattern):
    '''Return ``{match: path}`` for files in ``indir`` matching ``pattern``.

    ``pattern`` is a compiled regular expression applied to the file name.
    Files (or symbolic links to files) are returned in sorted order. A
    ``ValueError`` is raised if ``indir`` does not exist or contains no
    matching file.
    '''
    if not os.path.isdir(indir):
        raise ValueError(
            f"input directory '{indir}' not found. Set input_dir in the"
            " configuration file, relative to the directory the pipeline"
            " is run from, or as an absolute path.")
    found = {}
    for name in sorted(os.listdir(indir)):
        match = pattern.match(name)
        if match and os.path.isfile(os.path.join(indir, name)):
            found[match] = os.path.join(indir, name)
    if not found:
        raise ValueError(
            f"no input files found in '{indir}'. FASTQ files should be named"
            " <sample>.fastq.gz (single end) or <sample>.fastq.1.gz and"
            " <sample>.fastq.2.gz (paired end).")
    return found


def fastqc_basename(path):
    '''Name that FastQC gives its output for ``path``, without ``_fastqc``.

    FastQC removes these extensions, in this order, from the file name.
    '''
    name = os.path.basename(path)
    for extension in (".gz", ".bz2", ".txt", ".fastq", ".fq", ".csfastq",
                      ".sam", ".bam"):
        if name.endswith(extension):
            name = name[:-len(extension)]
    return name


def wildcard_choices(values):
    '''Regular expression accepting exactly the strings in ``values``.

    Used as a wildcard constraint so that wildcards only match the files the
    workflow found.
    '''
    values = sorted(set(values), key=len, reverse=True)
    if not values:
        return "(?!)"
    return "|".join(re.escape(value) for value in values)
