'''
================================================================================
pipeline_sra.py
================================================================================
:Author: Sandi Yen
:Tags: Python

Overview
========

This pipeline downloads genomes/metagenomes 
using SRR accession numbers. Depending on where the sequences were deposited, 
download method varies slightly. This pipeline provides a means of parallelizing 
downloads from various repositories. Currently supported repositories are: ENA

Usage
=====

Pipeline downloads genomes/metagenomes based on 
SRA accession numbers (SRRXXXXXXXX)

Example::

    ocms_toolkit sra make ena


Configuration
-------------
ocms_toolkit sratoolkit config

Input files
-----------

ENA Task:
ENA browser generates a download script for all the SRR accessions in a 
bioproject. Point to this file in pipeline.yml and the pipeline will
download each accessing into their respective subdirectories


Requirements
------------

Pipeline output
===============
Each accession number is downloaded into its own respective subdirectory.


Glossary
========

..glossary::


Code
====
'''

from ruffus import *
import os, sys
from cgatcore import pipeline as P
from ocmstoolkit.modules import Sra as SRA

PARAMS = P.get_parameters(['pipeline.yml'])

try:
    # parse ena script
    sra_ftp, sra_fqs = SRA.parse_ena_script(PARAMS['ena_script'])
    sentinels = []
    for accession, fqs in sra_fqs.items():
        # make subdirectories
        if not os.path.exists(accession):
            os.mkdir(accession)
        for fq in fqs:
            sentinel = f"{accession}/{fq}.sentinel"
            sentinels.append(sentinel)
except (IOError, FileNotFoundError):
    pass

@originate(sentinels, 'sra_download.log')
def enaDownload(sentinel, logfile):
    '''
    parse ENA download script and download each in a seperate job.
    Jobs are not sent to cluster as download requires internet access.
    '''
    
    # initialize log file
    open(logfile, 'w').close()
    accession = os.path.dirname(sentinel)

    # get fastq file names and ftp links
    fqs = sra_fqs[accession]
    ftp = sra_ftp[accession]

    # go into accession subdirectory
    statements = [f"cd {accession}"]
    # wget ftp commands
    statements.extend(ftp)
    # sentinel file for each fastq file
    statements.append("cd ../")
    entry = [f"echo 'downloaded {x}' >> {logfile}" for x in fqs]
    statements.extend(entry)
    statements.append(f"touch {sentinel}")
    
    statement = " && ".join(statements)

    P.run(statement, without_cluster = True)

@follows(enaDownload)
@transform(r"(\S+)/*.gz",
           regex("\S+\/(.*gz$)"),
           r"check_sums/\1.md5")
def generateMD5(infile, outfile):
    statement = f"md5sum {infile} > {outfile}"

    P.run(statement)

def check_report(infile):
     '''
     parses the md5 check sums report to see if all downloads passed md5 checks
     '''

     P.run(f"cat {infile} | grep -v OK > failed_check_sums.txt")

@posttask(check_report)
@transform(generateMD5,
           regex("check_sums/(\S+).md5"),
           "check_sums/report.txt")
def verifyMD5(infile, outfile):
    statement = f"md5sum -c {infile} > {outfile}"

    P.run(statement)

@follows(verifyMD5)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv))