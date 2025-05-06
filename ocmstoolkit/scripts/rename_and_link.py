'''
rename_and_link.py
==============

:Author: Sandi Yen
:Tags: Python

Overview
========

Use this script to sym link files and renames the sym links based on a id_mapping file. This is helpful for sym linking files that have very long/cumbersome barcodes produced by the sequencer. This is a stripped back version of combine lanes.

Usage
-----

The script takes four parameters. Specifying input directory --indir containing files to be symlinked. --suffix specifies the file extension of files to be symlinked (i.e. .fastq.1.gz, .fastq.gz etc.) Output directory --outdir is location where symlinks will be created. --mapping specifyies the file that maps the original barcodes and the new ids to be used when renaming. --log name of logfile, default = read.map

Example::
    ocms rename_and_link --indir=<INDIR> --src_suffix=<SRC_SUFFIX> --<TARGET_SUFFIX> --outdir=<OUTDIR> --mapping<ID-MAPPING>
    ocms rename_and_link -i raw -s _1.fastq.gz -t .fastq.1.gz -o renamed -m id_mapping.tsv -l read1.map
    
    # run directory
    raw/
        /raw/long_barcode1.fastq.1.gz
        /raw/long_barcode2.fastq.1.gz
        id_mapping.tsv
    renamed/
        /renamed/clean_id1.fastq.1.gz
        /renamed/clean_id2.fastq.1.gz

    # id_mapping.tsv
    long_barcode1    clean_id1
    long_barcode2    clean_id2

Type::
  
    ocms rename_and_link --help

for comand line help.

'''

import sys
import os
import re
import glob
import argparse

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """
    
    if argv is None:
        argv = sys.argv
    
    # set up command line parser
    parser = argparse.ArgumentParser(description=__doc__)
    
    parser.add_argument("-i", "--indir", dest="indir", type=str,
                        required=True, help="input directory")
    parser.add_argument("-s", "--src_suffix", dest="src_suffix", type=str,
                        required=True,
                        help=("extension that encompasses all files to be symlinked"
                              " (i.e. .fastq.1.gz, _1.fastq.gz). glob expressions"
                              " currently not supported."))
    parser.add_argument("-t", "--targ_suffix", dest="target_suffix", type=str,
                        default=None,
                        help=("extension of target files."
                              " (i.e. .fastq.1.gz). glob expressions"
                              " currently not supported."))
    parser.add_argument("-o", "--outdir", dest="outdir", type=str,
                        required=True, help="output directory")
    parser.add_argument("-m", "--mapping", dest="mapping", type=str,
                        required=True,
                        help=("tab-seperated file mapping file names to new names."
                              " First column is barcode (no file extension)."
                              " Second column is sample ID that files will be renamed to"
                              " (no file extension)"))
    parser.add_argument("-l", "--log", dest="logfile", type=str,
                        default='read.map',
                        help="name of log file, default=read.map")

    # unpack commandline arguments
    args = parser.parse_args()

    # run directory
    rundir = os.path.abspath(".")
    
    # make output directories
    os.makedirs(args.outdir, exist_ok=True)
    
    # get input files
    readfiles = glob.glob(os.path.join(rundir, args.indir, "*"+args.src_suffix))
    readfiles.sort()

    # checks for existing logfile in output directory
    assert not os.path.isfile(os.path.join(rundir, args.outdir,args.logfile)), (
        f"{args.logfile} already exists in the output directory."
                         " Specify different name for log file"
                         " with --log")

    # check for target extension
    if args.target_suffix is None:
        args.target_suffix = args.src_suffix

    #read identifier map into dictionary
    mapping_dict = {}
    inf = open(os.path.join(rundir, args.indir, args.mapping), encoding="utf_8_sig")
    for line in inf.readlines():
        data = line[:-1].split("\t")
        mapping_dict[data[0]] = data[1]
    
    found_barcode = {}
    out_map = open(os.path.join(args.outdir, args.logfile), "w")
    for readfile in readfiles:
        prefix = os.path.basename(readfile).replace(args.src_suffix, "")
        
        # make sure that the file has a mapping id
        assert prefix in mapping_dict, (f"ERROR: File {readfile}s exists but there is no "
                                        f"sample ID for {prefix} in --mapping file")

        if prefix in found_barcode:
            continue
        else:
            new_name = mapping_dict[prefix] + args.target_suffix
            new_name = os.path.join(rundir, args.outdir, new_name)
            found_barcode[prefix] = [readfile, new_name]
            out_map.write("\t".join(found_barcode[prefix]) + "\n")

        original = found_barcode[prefix][0]
        newname = found_barcode[prefix][1]

        # check symlink doesn't already exist
        assert not os.path.exists(newname), (
            "Remove previously created files from directory."
            " Will not force overwriting")
        statement = f"ln -s {original} {newname}"
        os.system(statement)
        
        
    def main(argv=None):
        main(argv)
        
if __name__ == "__main__":
    sys.exit(main(sys.argv))
