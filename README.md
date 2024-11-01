# About OCMS_Toolkit
OCMS_Toolkit is a repository of scripts useful for bioinformatics work at OCMS.

# Requirements
python 3.8

# Install
OCMS_Toolkit can be installed by running the following:
```
git clone https://github.com/OxfordCMS/OCMS_Toolkit
cd OCMS_Toolkit
python setup.py install
```
This will place relevant modules in your path and enable the use of a command line interface.

# OCMS_Toolkit
## scripts
### new_project
Use this script to initiate a new project on BMRC. `new_project` sets up the relevant directories and symlinks on BMRC in group and/or user directories that conforms to how OCMS works on the BMRC.

```
ocms_toolkit new_project --project_name=NEW_PROJECT --level=both
```

`-p` or `--project_name` is the new project name
`-l` or `--level` is the level at which new projects should be made. Takes `group`, `user`, or `both`. `--level=group` creates the directories in `projects` and `archive`. `l-level=user` creates directories in `devel` and `work`. `--level=both` makes all directories. You may want to set `--level=user` if the project has already been created in `project` and `archive` and you just need the directories in your own user space.

### rename_and_link
User this script to sym link files and rename the sym links based ona mapping file. This is helpful for symlinking files that have very long/cumbersome barcodes produced by the sequencer. This is a stripped back version of combine_lanes.py

This script takes five parameters.
`-i` or `--indir` Specifying input directory containing files to be symlinked.
`-s` or `--suffix` specifies the file extension of files to be symlinked (i.e. `.fastq.1.gz`, `.fastq.gz` etc.)
`-o` or `--outdir` Output directory is the location where symlinks will be created.
`-m` or `--mapping` specifies the file that maps the original barcodes and the new IDs to be used when renaming
`-l` or `--log` specifies name of logfile produced. default=read.map
```
    ocms rename_and_link --indir=<INDIR> --suffix=<SUFFIX> --outdir=<OUTDIR> --mapping<ID-MAPPING>
    ocms rename_and_link -i raw -s .fastq.1.gz -o renamed -m id_mapping.tsv -l read1.map
    raw/
	/raw/long_barcode1.fastq.1.gz
	/raw/long_barcode2.fastq.1.gz
    renamed/
	/renamed/clean_id1.fastq.1.gz
	/renamed/clean_id2.fastq.1.gz

    mapping.tsv
    long_barcode1    clean_id1
    long_barcode2    clean_id2
```

## pipelines
### pipeline_subsample.py
This script uses seqtk to randomly subsample fastq files (with seed). Script takes in all fastq.*gz in input.dir and subsamples to a specified read depth.

#### Configuration
`ocms_toolkit subsample_fastq config`

#### Requirements
`module load seqtk/1.4-GCC-12.2.0`