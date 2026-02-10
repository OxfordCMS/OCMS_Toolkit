'''
subsample_fastq.py
====================

:Author: Sandi Yen
:Tags: Python

Overview
========

This script uses seqtk to randomly subsample fastq files (with seed)

Usage
=====

Script takes in all fastq.*gz files in input.dir and subsamples to a specified read depth.

Example::

    ocms_toolkit subsample_fastq make full


Configuration
-------------
ocms_toolkit subsample_fastq config

Input files
-----------
Input files should be fastq files. Can be single (fastq.gz) or paired-end reads (fastq.1.gz, fastq.2.gz). 


Requirements
------------
module load seqtk/1.4-GCC-12.2.0

Pipeline output
===============
subsampled.dir containing subsampled fastq files.


Glossary
========

..glossary::


Code
====

'''

import sys
import os
import re
import glob
from pathlib import Path
from ruffus import *
from cgatcore import pipeline as P 

# get all sequence files within directory to process
PARAMS = P.get_parameters(['pipeline.yml'])
indir = PARAMS.get("general_input.dir", "input.dir")
SEQUENCEFILES = (f"{indir}/*fastq.*gz")
SEQUENCEFILES_REGEX = regex(fr"{indir}/(\S+)\.(fastq.*gz)")


######################################################
######################################################
######################################################
# subsamples fastq files to specified depth

# produces a subsampled.dir which contains
# subsampled fastq files
######################################################
######################################################
######################################################

@follows(mkdir("subsampled.dir"))
@transform(SEQUENCEFILES, 
         SEQUENCEFILES_REGEX,
         r"subsampled.dir/\1_subsampled.\2")

def subsample_fastq(infile, outfile):
    """Expects single end to be fastq.gz and paired end fastq files to be 
       in format fastq.1.gz, fastq.2.gz
    """

    depth = PARAMS['seqtk_depth']
    # subsample file with seed
    statement = '''seqtk sample -s100 %(infile)s %(depth)s | gzip %(outfile)s'''
    P.run(statement,
          job_threads = PARAMS['seqtk_job_threads'],
          job_memory = PARAMS['seqtk_job_memory'])

@follows(subsample_fastq)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
