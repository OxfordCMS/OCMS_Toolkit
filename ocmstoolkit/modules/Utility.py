# module for utility functions useful in handling shotgun data

import os
import re
import glob
import errno
import cgatcore.iotools as IOTools
import ast
import csv
from collections import defaultdict

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
        self.has_singleton = None
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
        self.is_paired()
        self.find_singleton()
        self.get_format()
        self.get_suffix()
        if self.prefixstrip is None:
            self.prefix = os.path.basename(self.fastn1.rstrip(self.fn1_suffix))
        else:
            self.prefix = os.path.basename(self.fastn1.rstrip(self.prefixstrip))

    def is_paired(self):
        '''check if paired
        '''
        if self.fastn1.endswith(".1.gz"):
            paired_name = self.fastn1.replace(".1",".2")
            assert len(glob.glob(paired_name)) > 0, (
                f"cannot find read 2 file at location {paired_name}"
                f" associated with read 1 file {self.fastn1}")
            self.fastn2 = paired_name
            # fastn3 is defined if data is paired, regardless of whether or
            # or not singletons actually exist - this is because a fastn3 file
            # is created automatically in pipelines even if there are no singletons
            self.fastn3 = self.fastn1.replace(".1",".3")

    '''check for singletons'''
    def find_singleton(self):
        fq3_name = self.fastn1.replace(".1",".3")
        if os.path.exists(fq3_name):
            # check file is not empty
            if os.stat(fq3_name).st_size != 0:
                self.has_singleton = True
            else:
                self.has_singleton = False


    '''check it is fasta or fastq and if compressed'''    
    def get_format(self):
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
    def get_suffix(self):
        # set suffix
        # if self.fastq1.endswith(".fastq.1.gz"):
        if self.fastn2 is not None:
            assert self.fastn1.endswith(f".{self.fileformat}.1.gz"), (
                f"Paired-end {self.fileformat} files must be in notation"
                f" '.{self.fileformat}.1.gz'")
            self.fn1_suffix = f".{self.fileformat}.1.gz"
            self.fn2_suffix = f'.{self.fileformat}.2.gz'
            # define fn3_suffix regardless if file exists
            self.fn3_suffix = f'.{self.fileformat}.3.gz'
        else:
            assert self.fastn1.endswith(f".{self.fileformat}.gz"), (
                "Single-end fastq files must be in notation "
                f"'{self.fileformat}.gz'"
            )
            self.fn1_suffix = f".{self.fileformat}.gz"

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

class CheckMParser:
    """
    Minimal parser for CheckM storage/bin_stats_ext.tsv.

    Usage:
        p = CheckMParser("/path/to/03_checkm.dir/sample/tool")
        p.parse_and_write(
            hq_comp=90.0, hq_cont=5.0,
            mq_comp=50.0, mq_cont=10.0,
            score_mult=5.0, top_n=10,
            sample_name="sample", tool_name="metabat2/maxbin2/etc"
        )

    Outputs (written into checkm_out_dir):
      - bin_stats_parsed.tsv
      - top_bins.tsv
      - hq_bin_ids.txt
      - quality_summary.txt
    """
    def __init__(self, checkm_out_dir):
        self.outdir = checkm_out_dir
        self.bin_stats_ext = os.path.join(checkm_out_dir, "storage", "bin_stats_ext.tsv")
        self.parsed_out = os.path.join(checkm_out_dir, "bin_stats_parsed.tsv")
        self.top_out = os.path.join(checkm_out_dir, "top_bins.tsv")
        self.hq_out = os.path.join(checkm_out_dir, "hq_bin_ids.txt")
        self.summary_out = os.path.join(checkm_out_dir, "quality_summary.txt")

    @staticmethod
    def _tofloat(x, default=0.0):
        try:
            return float(x)
        except Exception:
            return default

    @staticmethod
    def _toint(x, default=0):
        try:
            return int(x)
        except Exception:
            try:
                return int(float(x))
            except Exception:
                return default

    def parse_bin_stats_ext(self):
        """Return list of rows (dict) parsed from bin_stats_ext.tsv (dict-per-line)."""
        rows = []
        if not os.path.exists(self.bin_stats_ext):
            return rows

        with open(self.bin_stats_ext, "r") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    bid, dict_text = line.split("\t", 1)
                    stats = ast.literal_eval(dict_text)
                    if not isinstance(stats, dict):
                        continue
                except Exception:
                    continue

                comp = self._tofloat(stats.get("Completeness", stats.get("completeness", 0)))
                cont = self._tofloat(stats.get("Contamination", stats.get("contamination", 0)))
                gsize = self._toint(stats.get("Genome size", stats.get("Genome_size", 0)))
                gc = self._tofloat(stats.get("GC", stats.get("gc", 0)))
                nsc = self._toint(stats.get("# scaffolds", stats.get("# contigs", 0)))
                n50 = self._toint(stats.get("N50 (scaffolds)", stats.get("N50 (contigs)", 0)))
                genes = self._toint(stats.get("# predicted genes", 0))

                rows.append({
                    "bin_id": bid,
                    "Completeness": comp,
                    "Contamination": cont,
                    "Genome_size": gsize,
                    "GC": gc,
                    "N_scaffolds": nsc,
                    "N50": n50,
                    "Predicted_genes": genes
                })
        return rows

    def parse_and_write(self, *, hq_comp=90.0, hq_cont=5.0,
                        mq_comp=50.0, mq_cont=10.0,
                        score_mult=5.0, top_n=10,
                        sample_name=None, tool_name=None):
        """
        Parse bin_stats_ext and write outputs.
        score_mult: float or None. If None -> DO NOT compute score; sort by completeness desc, contamination asc.
        """
        rows = self.parse_bin_stats_ext()

        # write an informative empty marker if nothing parsed
        if not rows:
            with open(self.parsed_out, "w") as f:
                f.write("# no parsed rows (bin_stats_ext.tsv missing or empty)\n")
            with open(self.summary_out, "w") as f:
                f.write(f"Sample: {sample_name or 'unknown'}\nTool: {tool_name or 'unknown'}\nNo bins parsed.\n")
            return

        # classify
        for r in rows:
            comp = r["Completeness"]
            cont = r["Contamination"]
            if comp >= hq_comp and cont <= hq_cont:
                r["quality"] = "HQ"
            elif comp >= mq_comp and cont <= mq_cont:
                r["quality"] = "MQ"
            else:
                r["quality"] = "Low"

            if score_mult is not None:
                r["score"] = r["Completeness"] - score_mult * r["Contamination"]

        # sort
        if score_mult is None:
            # sort by completeness desc, contamination asc
            rows.sort(key=lambda x: (-x["Completeness"], x["Contamination"]))
        else:
            # sort by score desc, completeness desc, contamination asc
            rows.sort(key=lambda x: (x.get("score", 0.0), x["Completeness"], -x["Contamination"]), reverse=True)

        # write parsed TSV
        headers = ["bin_id","Completeness","Contamination","Genome_size","GC","N_scaffolds","N50","Predicted_genes","quality"]
        if score_mult is not None:
            headers.insert(headers.index("quality"), "score")
        with open(self.parsed_out, "w", newline="") as outf:
            w = csv.DictWriter(outf, fieldnames=headers, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for r in rows:
                w.writerow(r)

        # write top N
        with open(self.top_out, "w", newline="") as tf:
            tw = csv.writer(tf, delimiter="\t")
            if score_mult is None:
                tw.writerow(["rank","bin_id","Completeness","Contamination","quality"])
                for i, r in enumerate(rows[:top_n], start=1):
                    tw.writerow([i, r["bin_id"], r["Completeness"], r["Contamination"], r["quality"]])
            else:
                tw.writerow(["rank","bin_id","score","Completeness","Contamination","quality"])
                for i, r in enumerate(rows[:top_n], start=1):
                    tw.writerow([i, r["bin_id"], r["score"], r["Completeness"], r["Contamination"], r["quality"]])

        # write HQ ids
        with open(self.hq_out, "w") as hf:
            for r in rows:
                if r["quality"] == "HQ":
                    hf.write(r["bin_id"] + "\n")

        # write summary
        with open(self.summary_out, "w") as sf:
            sf.write(f"Sample: {sample_name or 'unknown'}\n")
            sf.write(f"Tool: {tool_name or 'unknown'}\n")
            sf.write(f"Thresholds: HQ: completeness>={hq_comp}, contamination<={hq_cont}; MQ: completeness>={mq_comp}, contamination<={mq_cont}\n")
            if score_mult is None:
                sf.write("Score: disabled (sorting by completeness/contamination)\n\n")
            else:
                sf.write(f"Score formula: completeness - {score_mult}*contamination\n\n")
            sf.write(f"Total bins parsed: {len(rows)}\n")
            sf.write(f"HQ: {sum(1 for r in rows if r['quality']=='HQ')}\n")
            sf.write(f"MQ: {sum(1 for r in rows if r['quality']=='MQ')}\n")
            sf.write(f"Low: {sum(1 for r in rows if r['quality']=='Low')}\n")
            if any(r["quality"] == "HQ" for r in rows):
                sf.write("\nHigh-quality bins (top):\n")
                for r in rows:
                    if r["quality"] == "HQ":
                        if score_mult is None:
                            sf.write(f"{r['bin_id']}\tComp:{r['Completeness']}\tCont:{r['Contamination']}\n")
                        else:
                            sf.write(f"{r['bin_id']}\tComp:{r['Completeness']}\tCont:{r['Contamination']}\tScore:{r['score']}\n")

def aggregate_checkm_summaries(root_dir="03_checkm.dir"):
    """
    Simple aggregator for per-tool parsed CheckM files.

    Writes:
      - 03_checkm.dir/<sample>/combined_summary.tsv
      - 03_checkm.dir/<sample>/combined_summary.txt
      - 03_checkm.dir/global_binning_summary.tsv
    """
    # find all parsed files
    pattern = os.path.join(root_dir, "*", "*", "bin_stats_parsed.tsv")
    parsed_files = sorted(glob.glob(pattern))
    if not parsed_files:
        return  # nothing to do

    # keep a global list of rows (dicts) and per-sample stats
    global_rows = []
    sample_stats = defaultdict(lambda: defaultdict(int))

    for pf in parsed_files:
        # derive sample/tool from relative path: <root_dir>/<sample>/<tool>/bin_stats_parsed.tsv
        rel = os.path.relpath(pf, root_dir)
        parts = rel.split(os.sep)
        if len(parts) < 3:
            # unexpected layout: skip
            continue
        sample, tool = parts[0], parts[1]

        # read parsed TSV for this tool
        with open(pf, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = list(reader)

        # append to per-sample combined_summary.tsv (create header on first write)
        combined_tsv = os.path.join(root_dir, sample, "combined_summary.tsv")
        os.makedirs(os.path.dirname(combined_tsv), exist_ok=True)
        write_header = not os.path.exists(combined_tsv)
        with open(combined_tsv, "a", newline="") as outfh:
            fieldnames = ["sample", "tool"] + (reader.fieldnames or [])
            w = csv.DictWriter(outfh, fieldnames=fieldnames, delimiter="\t")
            if write_header:
                w.writeheader()
            for r in rows:
                r_out = {"sample": sample, "tool": tool}
                r_out.update(r)
                w.writerow(r_out)

        # update per-sample counts and collect for global table
        for r in rows:
            q = r.get("quality", "").strip() or "Low"  # default to Low if missing
            sample_stats[sample][f"{tool}_total"] += 1
            sample_stats[sample][f"{tool}_{q}"] += 1

            r2 = {"sample": sample, "tool": tool}
            r2.update(r)
            global_rows.append(r2)

    # write per-sample human-readable summary .txt files
    for sample, stats in sample_stats.items():
        summary_txt = os.path.join(root_dir, sample, "combined_summary.txt")
        with open(summary_txt, "w") as outfh:
            outfh.write(f"Sample: {sample}\n\n")
            # find tools seen for this sample
            tools = set(k.split("_")[0] for k in stats.keys())
            for tool in sorted(tools):
                total = stats.get(f"{tool}_total", 0)
                hq = stats.get(f"{tool}_HQ", 0)
                mq = stats.get(f"{tool}_MQ", 0)
                low = stats.get(f"{tool}_Low", 0)
                outfh.write(f"{tool}:\n")
                outfh.write(f"  Total bins: {total}\n")
                outfh.write(f"  HQ: {hq}\n")
                outfh.write(f"  MQ: {mq}\n")
                outfh.write(f"  Low: {low}\n\n")

    # write a single global TSV with all rows
    global_tsv = os.path.join(root_dir, "global_binning_summary.tsv")
    if global_rows:
        # determine fieldnames: start with sample,tool then all keys from first row
        first = global_rows[0]
        other_keys = [k for k in first.keys() if k not in ("sample", "tool")]
        fieldnames = ["sample", "tool"] + other_keys
        with open(global_tsv, "w", newline="") as outfh:
            w = csv.DictWriter(outfh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for r in global_rows:
                w.writerow(r)
