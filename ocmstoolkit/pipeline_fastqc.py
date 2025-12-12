"""====================
FastQC pipeline
====================


The readqc pipeline imports unmapped reads from one or more input
files and performs basic quality control steps.

Quality metrics are based on the FastQC tools, see
http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/ for further
details.

Usage
=====

Configuration
-------------

See :file:`pipeline.yml` for setting configuration values affecting
the workflow (pre-processing or no pre-processing) and options for
various pre-processing tools.

Input
-----

Reads are imported by placing files or linking to files in the :term:
`working directory`.

The default file format assumes the following convention:

   <sample>.<suffix>

fastq.gz
   Single-end reads in fastq format.

fastq.1.gz, fastq.2.gz
   Paired-end reads in fastq format.
   The two fastq files must be sorted by read-pair.

.. note::

   Quality scores need to be of the same scale for all input files.
   Thus it might be difficult to mix different formats.

Pipeline output
----------------

The major output is a set of HTML pages and plots reporting on the quality of
the sequence archive.

Example
=======

"""

# import ruffus
from ruffus import transform, merge, follows, mkdir, regex, suffix, \
    jobs_limit, subdivide, collate, active_if, originate, split, formatter

# import useful standard python modules
import sys
import os
import re
import shutil
import sqlite3
import glob
import pandas as pd

# import modules from the cgat code collection
import cgatcore.experiment as E
#import cgatpipelines.tasks.mapping as mapping
from cgatcore import pipeline as P
#import cgatpipelines.tasks.readqc as readqc
#import cgatpipelines.tasks.preprocess as preprocess
import cgatcore.iotools as iotools


# Initialize the pipeline
P.initialize()

# Define input files and preprocessing steps list of acceptable input
# formats

PARAMS = P.get_parameters(['pipeline.yml'])
indir = PARAMS.get("general_input.dir", "input.dir")
SEQUENCEFILES = (f"{indir}/*fastq.*gz")

# optional group agnostic to paired end or single end reads
SEQUENCEFILES_REGEX = fr"{indir}/(\S+)(\.fastq)(\.[1-2])?.gz"

@follows(mkdir("fastqc.dir"))
@transform(SEQUENCEFILES,
           regex(SEQUENCEFILES_REGEX),
           r"fastqc.dir/\1\3.fastqc")
def runFastQC(infile, outfile, seq_regex=SEQUENCEFILES_REGEX):
    '''run FastQC on each input file.
    '''
    
    outdir = os.path.join("fastqc.dir/", os.path.basename(outfile))
    m = re.search(seq_regex, infile)
    if m.group(3) is None:
        fastqc_out = [m.group(1), '_fastqc']
    else:
        fastqc_out = [m.group(1), m.group(2), m.group(3), '_fastqc']
    fastqc_out = ''.join(filter(None, fastqc_out))
    
    to_move = os.path.join("fastqc.dir", fastqc_out)
    to_move_html = os.path.join("fastqc.dir", fastqc_out+".html")
    to_move_zip = os.path.join("fastqc.dir", fastqc_out+".zip")
        
    prefix = ''.join(filter(None, [m.group(1),m.group(3)]))
    out1 = os.path.join("fastqc.dir", prefix + ".fastqc.html")
    out2 = os.path.join("fastqc.dir", prefix + ".fastqc.zip")

    statement = ("fastqc --extract --outdir=fastqc.dir %(infile)s;"
                 " mv %(to_move)s %(outfile)s;"
                 " mv %(to_move_html)s %(out1)s;"
                 " mv %(to_move_zip)s %(out2)s")
                 #" rm -rf %(to_move)s")

    P.run(statement)

@merge(runFastQC, "multiqc_report.html")
def build_report(infiles, outfile):
    '''build report'''
    statement = '''multiqc . -s'''
    P.run(statement,
          job_threads = PARAMS['job_threads'],
          job_memory = PARAMS['job_memory'])

@follows(build_report)
@split("multiqc_data/multiqc_fastqc.txt",
       "multiqc_data/failed_qc_metrics/*.txt")
def find_failed_samples(infile, outfiles):
    """
    Identify samples that failed QC metrics.
    Output: one file per metric containing failed sample names.
    """

    df = pd.read_csv(infile, sep="\t")

    # Explicit list of QC metric columns based on your MultiQC file
    qc_metrics = [
        "basic_statistics",
        "per_base_sequence_quality",
        "per_tile_sequence_quality",
        "per_sequence_quality_scores",
        "per_base_sequence_content",
        "per_sequence_gc_content",
        "per_base_n_content",
        "sequence_length_distribution",
        "sequence_duplication_levels",
        "overrepresented_sequences",
        "adapter_content",
    ]

    output_dir = "multiqc_data/failed_qc_metrics"
    os.makedirs(output_dir, exist_ok=True)

    output_files = []

    for metric in qc_metrics:
        failed_samples = df[df[metric].astype(str).str.lower() == "fail"]["Sample"]

        if failed_samples.empty:
            continue

        outpath = os.path.join(output_dir, f"failed_{metric}.txt")
        failed_samples.to_csv(outpath, index=False, header=False)
        output_files.append(outpath)

    return output_files

@merge(find_failed_samples,
       "multiqc_data/failed_qc_metrics/failed_combined_samples.txt")
def combine_failed_samples(infiles, outfile):
    """
    Combine all failed-sample lists across QC metrics into one unique list.
    """

    failed_samples = set()

    # Read each metric-level failure file
    for infile in infiles:
        print(f"Processing file: {infile}")

        if not os.path.isfile(infile):
            print(f"Warning: file not found (skipping): {infile}")
            continue

        with open(infile, "r") as f:
            for line in f:
                failed_samples.add(line.strip())

    # Write the combined unique samples
    with open(outfile, "w") as out:
        for sample in sorted(failed_samples):
            out.write(sample + "\n")

@transform(combine_failed_samples,
           suffix(".txt"),
           "_report.html")
def rerun_failed_multiqc(infile, outfile):
    """
    Rerun MultiQC on the FastQC directories corresponding to failed samples.
    """

    # Read combined failed samples
    with open(infile, "r") as f:
        failed_samples = [line.strip() for line in f]

    fastqc_dirs = []

    for sample in failed_samples:
        #remove .gz
        name = sample.replace(".gz", "")
        #remove .fastq
        name = name.replace(".fastq", "")
        # append .fastqc to match directory names
        fastqc_dir = os.path.join("fastqc.dir", name + ".fastqc")

        print(f"Checking: {fastqc_dir}")

        if os.path.isdir(fastqc_dir):
            fastqc_dirs.append(fastqc_dir)
        else:
            print(f"WARNING: FastQC directory not found: {fastqc_dir}")

    if not fastqc_dirs:
        raise RuntimeError("No FastQC directories were found for rerunning MultiQC.")

    # Output directory for failed MultiQC
    output_dir = os.path.join("multiqc_data", "failed_qc_metrics", "multiqc_failed.dir")
    os.makedirs(output_dir, exist_ok=True)

    # Run MultiQC
    statement = f"multiqc {' '.join(fastqc_dirs)} -s -o {output_dir}"
    print(f"Running: {statement}")
    P.run(statement)

@follows(rerun_failed_multiqc)
def full():
    pass

def main(argv=None):
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
