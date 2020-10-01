# RNA-Seq Workflow

This is an RNA-Seq workflow written in `snakemake`, and tested on the University of Adelaide HPC (phoenix).
This has been written for single end reads only and will not behave correctly for paired-end reads.

## Outline

The steps currently implemented are

1. Download the DNA reference, GTF file
2. Generate the STAR index for the reference
3. Trim files using `AdapterRemoval`
4. Align Trimmed files using STAR
5. Count reads using `featureCounts`

Additional tasks performed are

1. Initialising a `workflowr` directory structure
2. FastQC reports are also generated for both raw and trimmed reads

Snakemake rules rely heavily on those provided at https://github.com/snakemake/snakemake-wrappers.
Version numbers for these are also available at https://snakemake-wrappers.readthedocs.io/

## Essential Files

In order to run this workflow, please ensure that you have

1. Placed unprocessed fastq files in the directory `data/raw/fastq`
2. Placed a `tsv` file (usually called `samples.tsv`) in the `config` folder
    + This file **must** contain a column called `sample`
3. Edited `config.yml` in the `config folder to ensure all parameters are correct
