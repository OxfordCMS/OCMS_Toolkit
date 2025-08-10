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
01_input_md5sum.dir contains md5sum hash of uncompressed .gz files
02_compressed.dir contains zstd compressed (.zst) files 
03_check_md5sum.dir contains stdout of md5sum check using uncompressed .zst files


Glossary
========

..glossary::


Code
====

"""

import sys
import os
import re
from pathlib import Path
from ruffus import follows, mkdir, transform, regex, collate
from cgatcore import pipeline as P 

# get all gunzip files within directory to process
FILES = ("input.dir/*.gz")
FILES_REGEX = regex(r"input\.dir\/(\S+)\.gz")

PARAMS = P.get_parameters(["pipeline.yml"])

###############################################################################
# Create md5sums for uncompressed input files
###############################################################################
@follows(mkdir("01_input_md5sum.dir"))
@transform(
    FILES, 
    FILES_REGEX,
    r"01_input_md5sum.dir/\1.md5"
)

def input_md5sum(infile, outfile):
    """Return md5sum for uncompressed input files"""

    # capture uncompressed filename from infile
    # file_name = re.search(r"input\.dir\/(\S+)\.gz", infile).group(1)

    # correct format for md5sum: <md5sum_checksum><space><space><file_name>
    # file_name = f"  {file_name}"

    # create statment for running md5sum
    statement = (
        "zstd"
        " --force"
        " --decompress"
        " --keep"
        " --stdout"
        f" {infile}"
        " | md5sum"
        # f" | sed 's/[[:blank:]]*-/{file_name}/g'"
        # creates error when doing md5sum check as file does not exit
        f" > {outfile}"
    )

    # create script for slurm job submission
    P.run(
        statement,
        job_threads=PARAMS["md5sum_job_threads"],
        job_memory=PARAMS["md5sum_job_memory"],
    )


###############################################################################
# Extract and re-compress input using zstd
###############################################################################
@follows(input_md5sum, mkdir("02_compressed.dir"))
@transform(
    FILES, 
    FILES_REGEX,
    r"02_compressed.dir/\1.zst"
)

def zstd_compress(infile, outfile):
    """Uncompress .gz files using gzip and then re-compress files using zstd"""

    # define level of compression
    compression_lvl = PARAMS['zstd_compression_lvl']

    threads = PARAMS['zstd_job_threads']

    if compression_lvl >= 20 and compression_lvl <= 22 :
        # create statment for running zstd using ultra compression
        statement = (
            "zcat"
            f" {infile}"
            " | zstd"
            " --compress"
            f" -{compression_lvl}"
            " --ultra"
            f" --threads={threads}"
            f" -o {outfile}"
        )
    elif compression_lvl > 0 and compression_lvl <= 19 : 
         # create statment for running zstd not using ultra compression
        statement = (
            "zcat"
            f" {infile}"
            " | zstd"
            " --compress"
            f" -{compression_lvl}"
            f" --threads={threads}"
            f" -o {outfile}"
        )
    else :
        # raise an error if invalid compression level entred
        print(f"ERROR: Invalid compression_lvl defined in pipeline.yml: {compression_lvl}")
        assert PARAMS['zstd_compression_lvl'] == 0 or PARAMS['zstd_compression_lvl'] > 22, \
            "Exception occurred: Invalid compression_lvl entered in pipeline.yml, must be between 1 - 22"        
    
    # create script for slurm job submission
    P.run(statement,
          job_threads = PARAMS['zstd_job_threads'],
          job_memory = PARAMS['zstd_job_memory'])
    

###############################################################################
# Check md5sums for uncompressed output files
###############################################################################
@follows(zstd_compress, mkdir("03_check_md5sum.dir"))
@transform(
    "02_compressed.dir/*.zst", 
    regex(r"02_compressed\.dir\/(\S+)\.zst"),
    r"03_check_md5sum.dir/\1.md5sum.out"
)

def check_md5sum(infile, outfile):
    """Uncompress .zstd file and perform md5sum checksum to check if it is the 
    same as the original uncompressed input files"""

    # capture filename from infile
    md5sum_file = re.search(r"02_compressed\.dir\/(\S+)\.zst", infile).group(1)

    # create md5sum filename
    md5sum_file = f"01_input_md5sum.dir/{md5sum_file}.md5"

    # create statment for extracting md5sum
    statement = (
        "zstd"
        " --decompress"
        " --keep"
        " --stdout"
        f" {infile}"
        " | md5sum"
        f" -c {md5sum_file}"
        f" > {outfile}"
        " 2>/dev/null"
    )

    # create script for slurm job submission
    P.run(
        statement,
        job_threads=PARAMS["md5sum_job_threads"],
        job_memory=PARAMS["md5sum_job_memory"],
    )

###############################################################################
# Creates a list of all samples with successfull or failed md5sum checks
###############################################################################
@follows(check_md5sum)
@collate(
    "03_check_md5sum.dir/*.md5sum.out",
    regex("03_check_md5sum.dir\/(\S+)\.md5sum.out"),
    "03_check_md5sum.dir/successfull.md5sum.out",
)

def list_md5_checks(infiles, outfile):
    """Creates list of samples where decompressed fastq.gz and decompressed 
    fastq.zstd files have the same md5sum"""

    # create empty list to store list of samples with successful md5sum checks
    success = []

    # create empty list to store list of samples with failed md5sum checks
    fail = []

    for infile in infiles : 

        # capture sample id from infile
        sample_id = re.search(r"03_check_md5sum\.dir\/(\S+)\.md5sum.out", infile).group(1)

        # read in outcome of md5sum check
        md5sum_check = open(os.path.abspath(infile),"r")

        md5sum_check = md5sum_check.read()

        if "OK" in md5sum_check :
            # add sample id to a list containing all successful md5sum checkss
            success.append(sample_id)

        elif "FAILED" in md5sum_check :
            # add sample id to a list containing all failed md5sum checkss
            fail.append(sample_id)

    # save file containing samples with successful md5sum checks
    with open(outfile, 'w') as output:
        # join the list elements into a single string with a newline character
        success = '\n'.join(success)
        # Write the data to the file
        output.write(success)

    # if any samples failed the md5sum check, save them as a list
    if len(fail) > 0 :
        # create file name for list of failed md5sum checks
        outfile2 = re.sub("successfull", "failed", outfile)
        # save file containing samples with failed md5sum checks
        with open(outfile2, 'w') as output:
            # join the list elements into a single string with a newline character
            fail = '\n'.join(fail)
            # Write the data to the file
            output.write(fail)    

@follows(list_md5_checks)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

