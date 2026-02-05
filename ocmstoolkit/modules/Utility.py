# module for utility functions useful in handling shotgun data

import os
import re
import glob
import errno
import cgatcore.iotools as IOTools
import cgatcore.experiment as E

# Check that the input files correspond
def get_fastns(datadir='.', *args: int):
    '''
    Gets all fastn files from specified directory. does not perform any checks
    for pairedness. Can return fastn1, fastn2, fastn3 with 1, 2, 3 respectively.
    Defaults to returning fastn1s.
    Examples: 
    get_fastns() # returns fastn1s in current directory
    get_fastns(., 1,2) # returns fastn1s and fastn2s in current directory
    get_fastns(., 1,2,3) # returns fastn1s, fastn2s and fastn3s
    '''
    # defaults to returning fastn1
    if not args: args = (1,)
    
    # convert args to set
    args_set = set(args)

    # valid sets of fastns
    valid_sets = [ {0}, {1}, {1,2}, {1,2,3}]
    assert args_set in valid_sets, (
        "specify fastns to return as 0; 1; 1,2; or 1,2,3."
        " i.e. get_fastns(datadir='.', 1,2) for paired end reads;"
        " i.e. get_fastns(datadir='.', 0) for single end reads"
    )

    # look for all fastn.1.gz, fastn.2.gz, fastn.gz
    fn_regex = re.compile(r"(\S+)(\.fast[a,q])(.*gz)")
    # search object
    fastn_search = [re.search(fn_regex, x) for x in os.listdir(datadir)]

    assert fastn_search, (
        f"No fastq or fasta files detected in {datadir} directory."
        " Check the file suffixes follow the notation"
        " fastq.1.gz, fastq.2.gz for paired end reads,"
        " and fastq.gz for single end reads")

    # build list of fastn1, 2, 3 files
    fastn0s = []
    fastn1s = []
    fastn2s = []
    fastn3s = []
    for i in fastn_search:
        prefix, format, suffix = i.groups() # sample name, suffix, end
        # detect end
        if suffix == ".gz":
            fastn0s.append(i.group(0))
        elif suffix == ".1.gz":
            fastn1s.append(i.group(0))
        elif suffix == '.2.gz':
            fastn2s.append(i.group(0))
        elif suffix == '.3.gz':
            fastn3s.append(i.group(0))
    
    # add path back onto file names
    fastn0s = [os.path.join(datadir, x) for x in fastn0s]
    fastn1s = [os.path.join(datadir, x) for x in fastn1s]
    fastn2s = [os.path.join(datadir, x) for x in fastn2s]
    fastn3s = [os.path.join(datadir, x) for x in fastn3s]

    if args_set == {0}:
        return fastn0s
    elif args_set == {1}:
        return fastn1s
    elif args_set == {1,2}:
        return fastn1s, fastn2s
    elif args_set == {1,2,3}:
        return fastn1s, fastn2s, fastn3s

def relink(inf, outf):
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
    A data class for handling fastn files. This data class is intented to
    be used when handling fastn files (not iterating through fastn files). 
    Intended to work with single, paired, or paired + singleton fastq or 
    fasta files. 

    For compatibility with pipelines, file names are expected to conform
    to naming conventions:
      single-end read files end with *.fastn.gz 
      read one of paired-end read files end with *.fastn.1.gz
      read two of paired-end read files end with *.fastn.2.gz
      singletons end with *.fastn.3.gz
    
    """

    def __init__(self, fastn1):
        self.fastn1 = fastn1
        self.fastn2 = None
        self.fastn3 = None
        self.fn1_suffix = None
        self.fn2_suffix = None
        self.fn3_suffix = None
        self.prefix = None
        self.file_format = None
        self.head = []

        # First check file can be opened on init & capture header 
        # (first 5 lines used for interleave and format checks)
        try:
            self.openfile = IOTools.open_file(self.fastn1)
        except FileNotFoundError as e:
            msg = f"cannot open file {self.fastn1}"
            raise Exception(msg) from e
        self.head = [self.openfile.readline().rstrip("\n") for x in range(5)]
        try:
             self.openfile.close()
        except Exception:
            pass
        finally:
            if hasattr(self, "openfile"):
                del self.openfile

        # Set the file format to either be fastq or fasta        
        self.get_format()

        # Check for the presence of (non-empty) mate files
        self.check_for_mates()
        
        # Set the file suffix attributes and the file name
        self.get_suffix()
        if self.prefixstrip is None:
            self.prefix = os.path.basename(self.fastn1.rstrip(self.fn1_suffix))
        else:
            self.prefix = os.path.basename(self.fastn1.rstrip(self.prefixstrip))

            
    def get_format(self):
        '''Check file is fasta or fastq and if compressed'''    
        extensions=("fasta","fastq")
        for i in extensions:    
            if self.fastn1.endswith((i+".1.gz",i+".gz")):
                self.file_format = i
                break
            else:
                continue
                       
        msg = f"File {self.fastn1} is not of the correct format (fasta or fastq)."
        assert self.file_format, msg

        if self.file_format == "fasta":
            assert self.head[0][0] == ">", (
                "Invalid header on first line for fasta format")
        else:
            assert self.head[0][0] == "@", (
                "Invalid header on first line for fastq format")

            
    def check_for_mates(self):
        '''
        Check if paired read and singleton read files exist.
        
        If singleton read file exists check that paired read file also
        exists, as singletons without pairs are not permissible.

        If singleton read file exists, but is empty, then treat the sample
        like no singletons exist as some tools/pipeline steps will create
        empty files by defult.'''
        if self.fastn1.endswith(".1.gz"):
            # Start by checking for mate pairs
            pair_name = self.fastn.rstrip(".1.gz") + ".2.gz"
            if os.path.exists(pair_name):
                self.fastn2 = pair_name

            # Subsequently check for singletons
            singleton_name = self.fastn.rstrip(".1.gz") + ".3.gz"
            if.os.path.exists(singleton_name):
                assert self.fastn2, f"Singleton file {singleton_name} exists without mate pair"
                # Only record presence of singleton file if it is not empty
                if len(IOTools.open_file(singleton_name).read(1)) > 0:              
                    self.fastn3 = singleton_name
                else:
                    E.warn(f"File {singleton_name} exists, but is empty"
                           " and is therefore being ignored")


    def get_suffix(self):
        '''Set file suffix attributes and prefix attribute'''

        if self.fastn1.endswith(self.file_format + ".1.gz"):
            self.fn1_suffix = '.' + self.file_format + '.1.gz'
            if self.fastn2:
                self.fn2_suffix = '.' + self.file_format + '.2.gz'
            if self.fastn3:
                self.fn3_suffix = '.' + self.file_format + '.3.gz'
        else:
            msg = f"Files expected to end with either {self.file_format}.gz or {self.file_format}.1.gz"
            assert self.fastn1.endswith('.' + self.file_format + '.gz'), msg
            self.fn1_suffix = '.' + self.file_format + '.gz'

        self.prefix = os.path.basename(self.fastn1).rstrip(self.fastn1_suffix)

        
class MetaBam:
    '''
    Data class for bam and sam files. 
    Handles filenames for when converting between cram, bam, and sam files.
    **Does not perform the actual conversion itself**
    '''
    def __init__(self, infile):
        if infile.endswith(".bam"):
            self.bamfile = infile
            self.samfile = re.sub("bam$", "sam", infile)
            self.cramfile = re.sub("bam$", "cram", infile)
        elif infile.endswith(".sam"):
            self.samfile = infile
            self.bamfile = re.sub("sam$", "bam", infile)
            self.cramfile = re.sub("sam$", "cram", infile)
        elif infile.endswith(".cram"):
            self.cramfile = infile
            self.bamfile = re.sub("cram$", "bam", infile)
            self.samfile = re.sub("cram$", "sam", infile)

        assert infile.endswith(('.bam','.sam','.cram')), (
            "Class MetaBam is only handles .bam, .sam, and .cram files")
