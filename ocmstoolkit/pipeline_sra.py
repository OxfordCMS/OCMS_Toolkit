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

import ruffus
import os
from cgatcore import pipeline as P
from ocmstoolkit.modules import Sra as SRA

PARAMS = P.get_parameters(['pipeline.yml'])

# parse ena script
sra_ftp, sra_fqs = SRA.parse_ena_script(PARAMS['ena_script'])
accessions = list(sra_fqs.keys())
fastqs = list(sra_fqs.values())

@mkdir(accessions)
@originate(fastqs)
def ENA(srr_accession):
    statement = sra_ftp[srr_accession]

    P.run(statement, to_cluster = False)