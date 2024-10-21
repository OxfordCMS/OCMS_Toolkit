'''
subsample_fastq.py
====================

:Author: Sandi Yen
:Tags: Python

Overview
========

This script concatenates fastq files of paired-end reads into one fastq file. It's written as a pipeline so paired-end fastqs can be processed as a job.

Usage
=====

Script takes in all fastq.*.gz files in current directory and concatenates paired-end reads for each sample into one fastq file. Concatenated fastqs written to  concat_fastq.dir

Example::

    ocms subsample_fastq make full


Configuration
-------------
No configuration required

Input files
-----------
Input files should be fastq.1.gz, fastq.2.gz


Requirements
------------
module load seqtk/1.4-GCC-12.2.0

Pipeline output
===============
reads from fastq.1.gz and fastq.2.gz concatenated into a fastq.gz file


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
SEQUENCEFILES = ("input.dir/*fastq.*gz")
SEQUENCEFILES_REGEX = regex(r"input.dir/(\S+)\.(fastq.*gz)")

PARAMS = P.get_parameters(['pipeline.yml'])

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

    depth = PARAMS['depth']
    # subsample file with seed
    fq = re.sub("\\.gz$", "", outfile)
    statement = '''seqtk sample -s100 %(infile)s %(depth)s > %(fq)s &&
                   gzip %(fq)s'''
    P.run(statement,
          job_threads = PARAMS['job_threads'],
          job_memory = PARAMS['job_memory'])

@follows(subsample_fastq)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
