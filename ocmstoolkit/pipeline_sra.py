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
from datetime import datetime
from pathlib import Path
import subprocess
import glob
from cgatcore import pipeline as P
from ocmstoolkit.modules import Sra as SRA

PARAMS = P.get_parameters(['pipeline.yml'])

try:
    # parse ena script
    srr_dict = SRA.parse_ena_script(PARAMS['ena_script'])
    sentinels = [os.path.join("downloaded.dir", x) for x in srr_dict.keys()]
        
except (IOError, FileNotFoundError):
    pass

@mkdir("downloaded.dir")
@originate(sentinels, 'README.md')
def enaDownload(sentinel, readme):
    '''
    parse ENA download script and download each in a seperate job.
    Jobs are not sent to cluster as download requires internet access.
    At the same time, produces a README.md to document the data download
    '''
    
    # initialize readme
    with open(readme, 'w') as f_readme:
        date = datetime.today().strftime('%Y-%m-%d')
        msg = [f"BioProject: {PARAMS['bioproject']}",
               f"Download pipeline run on: {date}",
               "Download log:\n"]
        f_readme.write('\n'.join(msg))


    # get fastq file names and ftp links
    key = os.path.basename(sentinel)
    ftp = srr_dict[key][0]
    fq = srr_dict[key][1]
    timestamp = datetime.today().strftime("%Y-%m-%d %X")
    # go into accession subdirectory
    statements = [f"cd downloaded.dir"]
    # wget ftp commands
    statements.append(ftp)
    # sentinel file for each fastq file
    statements.append("cd ../")
    entry = f"echo '{timestamp}    downloaded {fq}' > {sentinel}"
    statements.append(entry)
    
    statement = " && ".join(statements)

    P.run(statement, without_cluster = True)

@merge(enaDownload,
       ".sentinel.updateReadme")
def updateReadme(sentinels, outfile):
    '''
    Update README.md with download message in the sentinel
    '''
    for sentinel in sentinels:
        with (open(sentinel, 'r') as curr_sentinel,
              open("README.md", 'a') as readme):
            for line in curr_sentinel:
                readme.write(line)
    Path(outfile).touch()

@mkdir('checksums.dir')
@transform(enaDownload,
           regex("(\S+)\/.sentinel.(.*gz)"),
           r"checksums.dir/\2.md5")
def generateMD5(sentinel, md5):
    '''
    Generate md5 file for each download
    '''
    key = os.path.basename(sentinel)
    fq = srr_dict[key][1]
    statement = f"md5sum downloaded.dir/{fq} > {md5}"

    P.run(statement)

def merge_sentinels():
    '''
    combines all check sums sentinel file into one report and 
    parses the md5 check sums report to see if all downloads passed md5 checks
    Any downloads that filed the md5 check sum written to 
    checksums.dir/failed_check_sums.txt
    '''
    sentinels = glob.glob("checksums.dir/.sentinel*.md5")
    with open("checksums.dir/report.txt", 'a') as report:
        for sentinel in sentinels:
            with open(sentinel, 'r') as checksum:
                entry = checksum.readline()
            report.write(entry)

@posttask(merge_sentinels)
@transform(generateMD5,
           regex("checksums.dir/(\S+).md5"),
           r"checksums.dir/.sentinel.\1.md5")
def verifyMD5(md5, sentinel):
    '''
    Verify md5 sums and write to a checksums sentinel file
    '''
    statement = f"md5sum -c {md5} > {sentinel}"

    P.run(statement)

def update_readme_md5():
    '''
    Updates readme with checksums results. lists out files that failed
    check sums if applicable.
    '''
    # update readme
    with (open("README.md", 'a') as readme,
          open("checksums.dir/failed_check_sums.txt", "r") as f_failed):
        failed = f_failed.readlines()

        if failed:
            readme_entry = ["failed check sums:\n"]
            for line in failed:
                entry = line.split(" ")[0]
                readme_entry.extend(entry+"\n")
            readme.write('\n'.join(readme_entry))
        else:
            readme.write("all downloads passed md5 check sums")

@follows(verifyMD5)
@posttask(update_readme_md5)
@transform("checksums.dir/report.txt",
           regex(r"checksums\.dir/report\.txt"),
           "checksums.dir/failed_check_sums.txt")
def reportCheckSums(report, outfile):
    
    with open(outfile, 'w') as failed_check_sums:
        statement = f"grep -v OK {report}"
        subprocess.run(statement.split(" "))

@merge(enaDownload, "README.md")
def cleanUp(sentinels, readme):
    '''
    Delete downloaded files and update readme. Does this serially as is a 
    quick command line operation so not bothering to parallelizing
    '''
    fastqs = []
    for sentinel in sentinels:
        key = os.path.basename(sentinel)
        fq = srr_dict[key][1]
        fastqs.append(fq)
    timestamp = datetime.today().strftime("%Y-%m-%d %X")
    
    statements = []
    for fq in fastqs:
        statements.append(f"rm {fq}")
        statements.append(f"echo {timestamp}    deleting {fq} >> {readme}")
    
    statement = " && ".join(statements)
    
    P.run(statement)

@follows(reportCheckSums)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv))