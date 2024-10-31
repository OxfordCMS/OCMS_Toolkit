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

## pipelines
### pipeline_subsample.py
This script uses seqtk to randomly subsample fastq files (with seed). Script takes in all fastq.*gz in input.dir and subsamples to a specified read depth.

#### Configuration
`ocms_toolkit subsample_fastq config`

#### Requirements
`module load seqtk/1.4-GCC-12.2.0`