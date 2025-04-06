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

# optional group agnostic to paired end or single end reads
SEQUENCEFILES = ["*.fastq.*gz"]
SEQUENCEFILES_REGEX = r"(\S+)(\.fastq)(\.[1-2])?.gz"

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
    P.run(statement)
@follows(build_report)
@split("multiqc_data/multiqc_fastqc.txt", "multiqc_data/failed_qc_metrics/*.txt"
)
def find_failed_samples(infile, outfiles):
    '''Identify samples that failed specific QC metrics and generate separate fi
les for each metric.'''

    # Read the MultiQC FastQC summary file
    df = pd.read_csv(infile, sep="\t")

    # List of QC metric columns — adjust index range if your MultiQC version cha
nges
    qc_metrics = df.columns[12:23]

    output_dir = "multiqc_data/failed_qc_metrics"
    os.makedirs(output_dir, exist_ok=True)

    output_files=[]

    for metric in qc_metrics:
        failed_samples = df[df[metric].str.contains("fail", case=False, na=False
)]["Sample"]
        if not failed_samples.empty:
            output_file = os.path.join(output_dir, f"failed_{metric}.txt")
            failed_samples.to_csv(output_file, index=False, header=False)
            output_files.append(output_file)
    return output_files

@merge(find_failed_samples, "multiqc_data/failed_qc_metrics/failed_combined_samp
les.txt")
def combine_failed_samples(infiles, outfile):
    '''Combine all failed samples across QC metrics and remove duplicates.'''

    all_failed_samples = set()

    # Read each file and add samples to the set
    for infile in infiles:
        print(f"Processing file: {infile}")  # Debugging statement
        if not os.path.isfile(infile):
            print(f"File not found: {infile}")  # Debugging statement
            continue
        with open(infile, "r") as f:
            for line in f:
                all_failed_samples.add(line.strip())  # Add to set to avoid dupl
icates

    # Write the unique samples to the output file
    with open(outfile, "w") as f_out:
        for sample in sorted(all_failed_samples):
            f_out.write(sample + "\n")


@transform(combine_failed_samples,
           regex(r"(\S+).txt"),
           r"\1_report.html")
def rerun_failed_multiqc(infile, outfile):
    '''Rerun MultiQC on the combined list of failed samples.'''

    # Read the combined list of failed samples
    with open(infile, "r") as f:
        failed_samples = [line.strip() for line in f]

    # Extract unique sample names without the `.fastq.` and `.gz` extensions
    sample_names = set()
    for sample in failed_samples:
        # Remove `.fastq.` and replace `.1.gz` or `.2.gz` with `_1` or `_2`
        transformed_sample = sample.replace(".fastq.", "_").replace(".gz", "")
        sample_names.add(transformed_sample)

    # Locate the FastQC directories for these samples in fastqc.dir
    failed_dirs = []
    for sample in sample_names:
        dir_path = os.path.join("fastqc.dir", sample + ".fastqc")  # Fix: append
 `.fastqc` to directory name
        print(f"Checking directory: {dir_path}")  # Debugging statement
        if os.path.exists(dir_path):
            failed_dirs.append(dir_path)
        else:
            print(f"Directory does not exist: {dir_path}")  # Debugging statemen
t

    # Ensure there are directories to run MultiQC
    if failed_dirs:
        # Define the output directory inside `failed_qc_metrics`
        output_dir = os.path.join("multiqc_data", "failed_qc_metrics", "multiqc_
failed.dir")
        os.makedirs(output_dir, exist_ok=True)

        # Run MultiQC on the failed sample directories
        statement = f"multiqc {' '.join(failed_dirs)} -s -o {output_dir}"
        print(f"Running command: {statement}")  # Debugging statement
        P.run(statement)
    else:
        raise RuntimeError("No valid FastQC directories found for the failed samples.")
@follows(rerun_failed_multiqc)
def full():
    '''Final target.'''
    pass

def main(argv=None):
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
