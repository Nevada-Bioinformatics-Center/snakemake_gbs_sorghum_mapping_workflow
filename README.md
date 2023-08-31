# snakemake_gbs_sorghum_mapping_workflow

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.2.1-brightgreen.svg)](https://snakemake.bitbucket.io)

This workflow performs fastqc on an input PROJECT directory downloaded from basespace. It will merge the 
FASTQ files between lanes, then run fastqc on all merged data and compile an aggregate report with multiqc.

## Authors

* Hans Vasquez-Gross (@hans-vg)

## Usage

### Simple

#### Step 1: Install workflow

clone this workflow to your local computer


#### Step 2: Configure workflow

Configure the workflow according to your needs by editing the config.yaml to configure your input basespace PROJECT directory.

#### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --cluster qsub --jobs 100

or

    snakemake --use-conda --drmaa --jobs 100

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.
