"""
zstd_compression.py
====================

:Author: Holly Roach
:Tags: Python

Overview
========

This script uses zstd to compress gunzipped files into zstandard compressed files.
zstd has been shown to achive higher fastq compression than gzip:
http://www.bioinformaticszen.com/post/use-zstd-for-raw-fastq/
documentation: https://manpages.debian.org/testing/zstd/zstd.1.en.html


Usage
=====

Script takes in all *.gz files in input.dir, and uncompresses and re-compresses 
the file using zstd to obtain a better compression ratio.

Example::

    ocms_toolkit zstd_compression make full


Configuration
-------------
ocms_toolkit zstd_compression config

Input files
-----------
Input files should be gunzipped (.gz) files located in input.dir


Requirements
------------
module load zstd/1.5.5-GCCcore-12.3.0

Pipeline output
===============
01_input_md5sum.dir contains record of gunzipped prior to extraction.
02_compressed.dir contains zstd compressed files.
03_output_md5sum.dir contains record of gunzipped file after extraction from zstd


Glossary
========

..glossary::


Code
====

"""

import sys
import re
from pathlib import Path
from ruffus import follows, mkdir, transform, regex
from cgatcore import pipeline as P 

# get all gunzip files within directory to process
FILES = ("input.dir/*.gz")
FILES_REGEX = regex(r".*input\.dir\/(\S+)\.gz")

PARAMS = P.get_parameters(["pipeline.yml"])

###############################################################################
# Create md5sums for each input file
###############################################################################
@follows(mkdir("01_input_md5sum.dir"))
@transform(
    FILES, 
    FILES_REGEX,
    r"01_input_md5sum.dir/\1.md5sum.txt"
)

def input_md5sum(infile, outfile):
    """Return md5sum for uncompressed input files"""

    # create statment for running md5sum
    statement = (
        "gzip"
        " --decompress"
        " --keep"
        f" {infile}"
        f" | md5sum > {outfile}"
    )

    # create script for slurm job submission
    P.run(
        statement,
        job_threads=PARAMS["md5sum_job_threads"],
        job_memory=PARAMS["md5sum_job_memory"],
    )



@follows(input_md5sum)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

