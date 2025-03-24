'''
countlines.py
==============

:Author: Sandi Yen
:Tags: Python

Overview
========

Use this script to count the number of lines, or sequences in a file.

Usage
-----

The script takes 4 parameters. 
--infile file to be counted. Can be piped into ocmstoolkit countlines or can \
provide file name to the --infile flag
--outfile optional. name of output file that records name of file (see --fname) \
and count in tab separated format (fname    count)
--type determines type of counting. Set to 'lines' to count number of lines, \
'fasta' to count number of fasta sequences (lines/2), or 'fastq' to count number \
of fastq sequences (lines/4). defaults to 'lines'. 
--fname only used for logging purposes. the name of the file being counted. \
If fname is not set and file is piped in, output file will only contain number \
of lines.

Example::
    ocms_toolkit countlines --infile=<INFILE> --outfile=<OUTFILE> --type=<TYPE> --fname=<FILENAME>
    ocms_toolkit countlines sample1.fastq.1.gz -t fastq
    ocms_toolkit countlines sample1.fasta.1.gz -t fasta
    ocms_toolkit countlines example.tsv -t lines
    cat example.tsv | ocms_toolkit countlines -t fastq -f example
    

Type::
  
    ocms_toolkit countlines --help

for comand line help.

'''

import sys
import os
import subprocess
import argparse
import gzip

def main(argv=None):
    """
    Count the number of lines, or sequences in a file and returns count in stdout.
    parses command line options in sys.argv, unless *argv* is given.
    """
    
    if argv is None:
        argv = sys.argv
    
    # set up command line parser
    parser = argparse.ArgumentParser(description=__doc__)
    
    parser.add_argument("infile",
                        default=sys.stdin, nargs="?",
                        help="file to count")
    parser.add_argument("-t", "--type", dest="type", type=str, const='lines',
                        nargs='?', choices=['lines','fastq','fasta'], 
                        default='lines', help="type of counting. \
                            'lines' to count lines. \
                            'fastq' to count fastq sequences \
                            'fasta' to count fasta sequences")
    
    # unpack commandline arguments
    args = parser.parse_args()

    # stream in file or stdin
    if args.infile is sys.stdin: 
        infile = sys.stdin.buffer
    else:
        infile = open(args.infile, "rb")

    # Streaming data to subprocess
    with subprocess.Popen(
        ["wc", "-l"],
        stdin=infile,
        stdout=subprocess.PIPE,
        text=True  # Ensures output is treated as text
    ) as proc:
        nlines = proc.communicate()[0].rstrip()
    
    if args.infile is not sys.stdin:
        infile.close()
    
    # if type is fasta or fastq, make sure nlines is an even number
    if args.type == "lines":
        assert (nlines % 2) != 0, (
            f"File {args.infile} contains an odd number of lines but \
                --type is not set to 'lines'.")

    if args.type == 'fasta':
        nlines = int(nlines)/2
    elif args.type == 'fastq':
        nlines = int(nlines)/4
   
    # print to stdout
    out = str(int(nlines))
    sys.stdout.write(out + "\n")