# module for utility functions useful in handling shotgun data

import os
import re
import glob
import errno
import cgatcore.iotools as IOTools

# Check that the input files correspond
def get_fastns(datadir='.', out='fastn1'):
    '''
    Gets all fastn files from specified directory. If detects paired ends
    will check all samples are paired.
    out = 'fastn1' - returns list of fastn1s (works for both se and pe cases)
    out = 'fastn2' - returns list of fastn2s
    out = 'all' - returns list of fastn1s and None in case of single end; 
                  returns fastn1 and fastn2 in case of paired end
    '''
    assert out in ['fastn1','fastn2','all'], (
        "get_fastns argument out accepts 'fastn1','fastn2', or 'all'"
    )

    # look for all fastn.1.gz, fastn.2.gz, fastn.gz
    fn_regex = re.compile(r"(\S+)(\.fast[a,q].*gz)")
    # search object
    fastn_search = [re.search(fn_regex, x) for x in os.listdir(datadir)]

    # build dicionary of paired files
    fastn_dict = {}
    for i in fastn_search:
        curr_fastn = i.groups() # sample name, suffix, end
        # if already have matched, append file name
        if curr_fastn[0] in fastn_dict:
            fastn_dict[curr_fastn[0]].append("".join(curr_fastn))
        # add file name to dictionary
        else:
            fastn_dict[curr_fastn[0]] = ["".join(curr_fastn)]

    # getting fastn.1.gz
    fastn1s = [fastn_dict[x][0] for x in fastn_dict.keys()]
    fastn2s = []
    # check if any fastn.1.gz found
    paired = sum([len(fastn_dict[x]) != 1 for x in fastn_dict.keys()])
    
    # if paired, check if all files are paired
    if paired:
        # temporary dictionary to store unmatched fastns
        unmatched = {}
        for key, fn in fastn_dict.items():
            for curr_fn in fn:
                # if paired file found removed from unmatched dictionary
                if key in unmatched.keys():
                    unmatched.pop(key)
                    fastn2s.append(curr_fn)
                # if haven't found paired file yet, add to in dictionary
                else: 
                    unmatched[key] = curr_fn
        assert not unmatched, "Paired files detected but not all input files \
            have pairs"
        
    assert fastn_dict, (
        f"No input files detected in {datadir} directory."
        " Check the file suffixes follow the notation"
        " fastq.1.gz, fastq.2.gz for paired end reads,"
        " and fast.gz for single end reads")

    # add path back onto file names
    fastn1s = [os.path.join(datadir, x) for x in fastn1s]
    fastn2s = [os.path.join(datadir, x) for x in fastn2s]
    
    if out == 'fastn1':
        return fastn1s
    elif out == 'fastn2':
        return fastn2s
    elif out == 'all' and paired:
        return fastn1s, fastn2s
    elif out == 'all' and not paired:
        return fastn1s, None

def symlnk(inf, outf):
    try:
        os.symlink(os.path.abspath(inf), outf)
    except OSError as e:
        if e.errno == errno.EEXIST:
            os.remove(outf)
            os.symlink(inf, outf)

class BaseTool:
    """
    Base class for all tool classes for handling pipeline infiles, outfiles, 
    and PARAMS. All tool classes should inherit this base class.
    """
    def __init__(self, infile, outfile, **PARAMS):
        self.infile = infile
        self.outfile = outfile
        self.indir = os.path.dirname(os.path.abspath(infile))
        self.outdir = os.path.dirname(os.path.abspath(outfile))
        self.PARAMS = PARAMS

class MetaFastn:
    """
    A data class for handling fastn files. This data class is intented to be 
    used when handling fastn files (not iterating through fastn files). 
    Intended to work with single, paired, or paired + singleton fastq or 
    fasta files. 

    Some elements pulled form CGATMetaSequencing by Matt Jackson.

    Some options are  assumed to be passed via kwargs, as this and 
    inherited classes are written to work with a PARAMS dict 
    generated from a pipeline.yml config file.
    """

    def __init__(self, fastn1):
        self.fastn1 = fastn1
        self.fastn2 = None
        self.fastn3 = None
        self.fn1_suffix = None
        self.fn2_suffix = None
        self.fn3_suffix = None
        self.prefixstrip = None
        self.prefix = None
        self.fileformat = None
        self.head = []

        '''check file can be opened on init & capture header 
           (first 5 lines used for interleave and format checks)
        '''
        try:
            self.openfile = IOTools.open_file(self.fastn1)
        except FileNotFoundError as e:
            msg = f"cannot open file {self.fastn1}"
            raise Exception(msg) from e
        self.head = [self.openfile.readline().rstrip("\n") for x in range(5)]

        '''autocheck file format and pairedness, 
           read count must be specified seperately
        '''
        self.isPaired()
        self.hasSingleton()
        self.getFormat()
        self.getSuffix()
        if self.prefixstrip is None:
            self.prefix = os.path.basename(self.fastn1.rstrip(self.fn1_suffix))
        else:
            self.prefix = os.path.basename(self.fastn1.rstrip(self.prefixstrip))

    '''check if paired and if containts interleaved pairs or matching files'''
    def isPaired(self):
        if self.fastn1.endswith(".1.gz"):
            paired_name = self.fastn1.replace(".1",".2")
            assert len(glob.glob(paired_name)) > 0, (
                f"cannot find read 2 file at location {paired_name}"
                f" associated with read 1 file {self.fastn1}")
            self.fastn2 = paired_name

    '''check for singletons'''
    def hasSingleton(self):
        fq3_name = self.fastn1.replace(".1",".3")
        if os.path.exists(fq3_name):
            # check file is not empty
            if os.stat(fq3_name).st_size != 0:
                self.fastn3 = fq3_name

    '''check it is fasta or fastq and if compressed'''    
    def getFormat(self):
        extensions=("fasta","fastq")
        for i in extensions:    
            if self.fastn1.endswith((i+".1.gz",i+".gz")):
                if i == "fastq":
                    self.fileformat="fastq"
                else:
                    self.fileformat="fasta"
           
        msg = f"file {self.fastn1} is not of the correct format (fasta or fastq)."
        assert self.fileformat, msg
        if self.fileformat == "fasta":
            assert self.head[0][0] == ">", (
                "invalid header on first line for fasta format")
        else:
            assert self.head[0][0] == "@", (
                "invalid header on first line for fastq format")
            
    '''get fastq1 file suffix '''
    def getSuffix(self):
        # set suffix
        # if self.fastq1.endswith(".fastq.1.gz"):
        if self.fastn2 is not None:
            assert self.fastn1.endswith(f".{self.fileformat}.1.gz"), (
                f"Paired-end {self.fileformat} files must be in notation"
                f" '.{self.fileformat}.1.gz'")
            self.fn1_suffix = f".{self.fileformat}.1.gz"
            self.fn2_suffix = f'.{self.fileformat}.2.gz'
            self.fn3_suffix = f'.{self.fileformat}.3.gz'
        else:
            assert self.fastn1.endswith(f".{self.fileformat}.gz"), (
                "Single-end fastq files must be in notation "
                f"'{self.fileformat}.gz'"
            )
            self.fn1_suffix = f".{self.fileformat}.gz"

class MetaBam:
    '''
    Data class for bam and sam files
    '''
    def __init__(self, infile):
        if infile.endswith(".bam"):
            self.bamfile = infile
            self.samfile = re.sub("bam$", "sam", infile)
        elif infile.endswith(".sam"):
            self.samfile = infile
            self.bamfile = re.sub("sam$", "bam", infile)

        assert infile.endswith(('.bam','.sam')), (
            "Class MetaBam is only handles .bam and .sam files")