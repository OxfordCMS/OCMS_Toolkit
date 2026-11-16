[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15350987.svg)](https://doi.org/10.5281/zenodo.15350987)

# About OCMS_Toolkit
OCMS_Toolkit is a repository of scripts and Snakemake pipelines useful for bioinformatics work at OCMS.

# Requirements
- Python >= 3.11 (required by Snakemake 9)
- Snakemake >= 9 and the Snakemake Slurm executor plugin (installed with the package)
- The tools run by each pipeline (FastQC, MultiQC, seqtk), on the `PATH` of the shell from which the pipeline is launched, for example as environment modules

# Installation
OCMS_Toolkit can be installed by running the following:
```
git clone https://github.com/OxfordCMS/OCMS_Toolkit
cd OCMS_Toolkit
pip install .
```
This will place relevant modules in your path and enable the use of a command line interface, `ocms_toolkit`. List the available pipelines and scripts with:

```
ocms_toolkit --help
```

# OCMS_Toolkit
## scripts
<details>
<summary>new project</summary>

### new_project
Use this script to initiate a new project on BMRC. `new_project` sets up the relevant directories and symlinks on BMRC in group and/or user directories that conforms to how OCMS works on the BMRC.

```
ocms_toolkit new_project --project_name=NEW_PROJECT --level=both
```

`-p` or `--project_name` is the new project name
`-l` or `--level` is the level at which new projects should be made. Takes `group`, `user`, or `both`. `--level=group` creates the directories in `projects` and `archive`. `--level=user` creates directories in `devel` and `work`. `--level=both` makes all directories. You may want to set `--level=user` if the project has already been created in `project` and `archive` and you just need the directories in your own user space.
</details>

<details>
<summary>rename and link</summary>

### rename_and_link
Use this script to symlink files and rename the symlinks based on a mapping file. This is helpful for symlinking files that have very long/cumbersome barcodes produced by the sequencer. This is a stripped back version of combine_lanes.py

This script takes the following parameters.
`-i` or `--indir` Specifying input directory containing files to be symlinked.
`-s` or `--src_suffix` specifies the file extension of files to be symlinked (i.e. `.fastq.1.gz`, `_1.fastq.gz` etc.)
`-t` or `--targ_suffix` specifies the file extension of the symlinks (defaults to `--src_suffix`)
`-o` or `--outdir` Output directory is the location where symlinks will be created.
`-m` or `--mapping` specifies the file (in the input directory) that maps the original barcodes and the new IDs to be used when renaming
`-l` or `--log` specifies name of logfile produced. default=read.map
```
    ocms_toolkit rename_and_link --indir=<INDIR> --src_suffix=<SUFFIX> --outdir=<OUTDIR> --mapping=<ID-MAPPING>
    ocms_toolkit rename_and_link -i raw -s .fastq.1.gz -o renamed -m id_mapping.tsv -l read1.map
    raw/
	/raw/long_barcode1.fastq.1.gz
	/raw/long_barcode2.fastq.1.gz
	/raw/id_mapping.tsv
    renamed/
	/renamed/clean_id1.fastq.1.gz
	/renamed/clean_id2.fastq.1.gz

    id_mapping.tsv
    long_barcode1    clean_id1
    long_barcode2    clean_id2
```
</details>

<details>
<summary>fastq2filteredfastq</summary>

### fastq2filteredfastq
Filters a FASTQ file to keep the reads that are mapped (or unmapped) in a BAM file.

```
ocms_toolkit fastq2filteredfastq --in_fastq=x.fastq.gz --outdir=filtered.dir \
    --bamfile=x.bam --filter_mapping=unmapped
ocms_toolkit fastq2filteredfastq -i x.fastq.1.gz -o filtered.dir -b x.bam \
    -f mapped --paired
```
</details>

## pipelines
The pipelines are Snakemake workflows. Each is run from the directory in which its results should be written, in three steps:

1. Create the configuration file `pipeline_<name>.yaml` in the current directory. Your default Slurm account is filled in where it can be found.

   ```
   ocms_toolkit <name> config
   ```

2. Edit `pipeline_<name>.yaml`, then check which jobs would run:

   ```
   ocms_toolkit <name> full --dry
   ```

3. Run the pipeline, submitting at most 20 concurrent jobs to Slurm:

   ```
   ocms_toolkit <name> full --jobs 20
   ```

`ocms_toolkit <name> full --local --cores N` runs the pipeline on the current node instead of submitting jobs to Slurm. Options after a bare `--` are passed unchanged to Snakemake, for example to remove a stale lock after an interrupted run:

```
ocms_toolkit <name> full -- --unlock
```

Each configuration file has a `resources` section with the Slurm account and partition, which apply to every job, and a section for each tool with:

| Setting | Meaning |
|---------------|---------|
| `job_threads` | CPUs per job |
| `job_memory`  | memory per CPU, e.g. `4G`; a job receives `job_memory` x `job_threads` in total (as in the earlier cgatcore pipelines) |
| `job_runtime` | wall time in minutes; if empty, the default of the Slurm partition applies |

The `make` syntax of the earlier cgatcore pipelines (`ocms_toolkit fastqc make full -p 10 -v 6`) is replaced by `ocms_toolkit fastqc full --jobs 10`. Configuration files from the cgatcore pipelines (`pipeline.yml`) are not read; create a new one with `ocms_toolkit <name> config`.

<details>
<summary>pipeline fastqc</summary>

### fastqc
Runs FastQC on each FASTQ file, followed by MultiQC to produce a summary report, `multiqc_report.html`, in the working directory.

#### Input files
FASTQ files in `input_dir` (default `input.dir`): `<sample>.fastq.gz` for single-end reads, or `<sample>.fastq.1.gz` and `<sample>.fastq.2.gz` for paired-end reads. Singleton files (`.fastq.3.gz`) are not checked.

#### Requirements
FastQC and MultiQC, for example from the BMRC modules.

#### Usage
```
ocms_toolkit fastqc config
ocms_toolkit fastqc full --jobs 24
```

#### Outputs
```
fastqc.dir/
|- <sample>[.1|.2].fastqc/        extracted FastQC results
|- <sample>[.1|.2].fastqc.html
|- <sample>[.1|.2].fastqc.zip
multiqc_report.html
multiqc_data/
```
</details>

<details>
<summary>pipeline subsample</summary>

### subsample
This pipeline uses `seqtk` to randomly subsample each FASTQ file to a fixed number of reads. Subsampling uses the same seed (default 100) for every file, so while each file is subsampled independently, the subsampling is the same for paired-end reads.

#### Configuration
```
ocms_toolkit subsample config
```

`seqtk.depth` is the number of reads kept from each file, and `seqtk.seed` the random seed. Input files should use the notation `.fastq.gz`, or `.fastq.1.gz` and `.fastq.2.gz`, and be located in `input_dir` (default `input.dir`).

#### Requirements
`module load seqtk/1.4-GCC-12.2.0`

#### Usage
```
ocms_toolkit subsample full --jobs 10
```

#### Outputs
```
subsampled.dir/
|- <sample>_subsampled.fastq.1.gz
|- <sample>_subsampled.fastq.2.gz
```
</details>

<details>
<summary>pipeline sra</summary>

### sra
This pipeline parses the download script generated by the ENA Browser ("Download all"), runs each download as a separate job, and performs MD5 checksums on the downloads.

Downloads need internet access, so they run on the node from which the pipeline is launched rather than as Slurm jobs; `download.max_parallel` (default 4) sets how many run at the same time. The checksum jobs are submitted to Slurm.

#### Configuration
```
ocms_toolkit sra config
```

`ena_script` is the download script generated by the ENA Browser; `bioproject` is recorded in `README.md`.

#### Requirements
No modules required. Internet access required.

#### Usage
```
ocms_toolkit sra full --jobs 20
```

#### Outputs
- `downloaded.dir/` contains the downloaded files, and a sentinel file recording when each was downloaded.
- `checksums.dir/` contains an MD5 file for each download and the result of checking it, all results in `report.txt`, and any files that failed in `failed_check_sums.txt`.
- `README.md` records the BioProject, the date, the files downloaded, and the checksum results.

The optional `clean_up` step deletes the downloaded files and records the deletion in `README.md`:

```
ocms_toolkit sra full -- clean_up
```
</details>

# Development
Install the package in editable mode with the test dependencies, and run the tests:

```
pip install -e ".[test]"
pytest ocmstoolkit/tests/unit_tests
```

A pipeline consists of two files with the same name:

- `ocmstoolkit/workflows/pipeline_<name>.smk`, whose first line is a `#` comment describing the pipeline (shown in `ocms_toolkit --help`);
- `ocmstoolkit/config/pipeline_<name>.yaml`, the configuration template, with a `resources` section containing `slurm_account: "your_slurm_account"` and `slurm_partition`.

The pipeline is then available as `ocms_toolkit <name>` without changes to the command-line interface. The configuration is supplied to Snakemake with `--configfile`, so Snakefiles should not contain a `configfile:` directive. Scripts placed in `ocmstoolkit/scripts/` with a `main(argv)` function are available as `ocms_toolkit <script>`.
