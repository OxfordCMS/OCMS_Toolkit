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
def full():
    pass

def main(argv=None):
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
