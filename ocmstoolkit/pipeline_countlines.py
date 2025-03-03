'''
pipeline_countlines.py
====================

:Author: Sandi Yen
:Tags: Python

Overview
========

This script counts number of lines in all files in input.dir

Usage
=====

Script takes in all files in input.dir and counts then number of lines 
(or reads)

Example::

    ocms_toolkit countlines make full


Configuration
-------------
ocms_toolkit countlines config

Input files
-----------
text files, fastq files, fasta files

Requirements
------------

Pipeline output
===============


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
indir = PARAMS.get('input.dir','input.dir')
@mkdir('countlines.dir')
@transform(fr"{indir}/*",
           regex(f"{indir}/(.*)"),
           r"countlines.dir/\1.tmp")
def countLines(infile, outfile):

    with open(outfile, "w") as tmpfile:
        entry = os.path.basename(infile)
        tmpfile.write(f"{entry}\n")

    if infile.endswith(".gz"):
        statement = f"zcat {infile} | wc -l"
    else:
        statement = f"cat {infile} | wc -l"

    type = PARAMS['type']

    if type == 'fasta':
        statement = f"echo $({statement})/2 | bc"
    elif type == 'fastq':
        statement = f"echo $({statement})/4 | bc"

    statement = statement + f" >> {outfile}"

    P.run(statement, 
          job_threads=PARAMS['job_threads'],
          job_memory=PARAMS['job_memory'])


def clearTmp():
    statement = "rm countlines.dir/*.tmp"
    P.run(statement)

@posttask(clearTmp)
@merge(countLines,
       "countlines.dir/merged_countlines.tsv")
def mergeCountLines(infiles, outfile):
    count_dict = {}
    for infile in infiles:
        with open(infile, "r") as tmpfile:
            lines = tmpfile.readlines()
            lines = [x.rstrip("\n") for x in lines]
            # filename as key count as value
            fname = os.path.basename(lines[0])
            count_dict[fname] = lines[1]

    header = ['file', 'count']
    with open(outfile, "w") as merged:
        merged.write("\t".join(header) + "\n")
        for k,v in count_dict.items():
            entry = [k, v]
            merged.write("\t".join(entry) + "\n")
    
@follows(mergeCountLines)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv))