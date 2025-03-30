'''
zstd_compression.py
====================

:Author: Holly Roach
:Tags: Python

Overview
========

This script uses zstd to compress gunzipped files into zstandard compressed files.
zstd has been shown to achive higher fastq compression than gzip: http://www.bioinformaticszen.com/post/use-zstd-for-raw-fastq/
documentation: https://manpages.debian.org/testing/zstd/zstd.1.en.html


Usage
=====

Script takes in all *.gz files in input.dir, and uncompresses and re-compresses the file using zstd to obtain a better compression ratio

Example::

    ocms_toolkit zstd_compression make full


Configuration
-------------
ocms_toolkit zstd_compression config

Input files
-----------
Input files should be gunzipped (.gz) files 


Requirements
------------
module load zstd/1.5.5-GCCcore-12.3.0

Pipeline output
===============
compressed.dir containing zstd compressed files.


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

# get all gunzip files within directory to process
FILES = ("input.dir/.*gz")
FILES_REGEX = regex(r"input.dir/(\S+)\.*gz")

PARAMS = P.get_parameters(['pipeline.yml'])

######################################################
######################################################
######################################################
# compress .gz files using a specified level of compression

# produces a compressed.dir which contains
# zstd compressed files
######################################################
######################################################
######################################################


###############################################################################
# Create md5sums for each input file
###############################################################################
@follows(mkdir("01_input_md5sum.dir"))
@transform(FILES, 
         FILES_REGEX,
         r"01_input_md5sum.dir/\1_md5sum.txt")

def input_md5sum(infile, outfile):
    """Return md5sum for input files"""

    print(f"infile: {infile}")      # infile: input.dir/samples_pooled_corrected.megahit.contigs.fasta
    print(f"outfile: {outfile}")    # outfile: prokka_output.dir/samples_pooled_corrected.megahit/

    # capture sample id from output dir name
    sample_id = re.sub("01_input_md5sum.dir/", "", outfile)   # Sample ID: samples_pooled_corrected.megahit/
    print(f"Sample ID: {sample_id}")

    statement = f"gzip -dk {infile} | md5sum > {outfile}"

    P.run(statement,
          # md5sum can't be multi-threaded
          job_threads = 1,
          job_memory = 15G)

###############################################################################
# Re-compress input using zstd
###############################################################################
@follows(mkdir("02_compressed.dir"))
@transform(FILES, 
         FILES_REGEX,
         r"02_compressed.dir/\1.gz")

def zstd_compress(infile, outfile):
    """Uncompress .gz files using gunzip"""

    # define level of compression
    compression_lvl = PARAMS['compression_lvl']

    # subsample file with seed
    fq = re.sub("\\.gz$", "", outfile)
    statement = '''seqtk sample -s100 %(infile)s %(depth)s > %(fq)s &&
                   gzip %(fq)s'''
    P.run(statement,
          job_threads = PARAMS['job_threads'],
          job_memory = PARAMS['job_memory'])
    
###############################################################################
# Check md5sums for new zstd compressed files
###############################################################################
@follows(zstd_compress)
def check_md5sum(infile, outfile):
    """Expects input files to be gunzipped (.gz)
    """

    # define level of compression
    compression_lvl = PARAMS['compression_lvl']

    # subsample file with seed
    fq = re.sub("\\.gz$", "", outfile)
    statement = '''seqtk sample -s100 %(infile)s %(depth)s > %(fq)s &&
                   gzip %(fq)s'''
    P.run(statement,
          job_threads = PARAMS['job_threads'],
          job_memory = PARAMS['job_memory'])

@follows(check_md5sum)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
