# RNA-Seq Workflow

This is an RNA-Seq workflow written in `snakemake`, and tested on the University of Adelaide HPC (phoenix).
This has been written for single end reads only and will not behave correctly for paired-end reads.

It is assumed that the reference will be obtained from [Ensembl](www.ensembl.org) and this workflow is not compatible with any other source, without further modification.

__The compilation of some qc*Rmd files currently depends on bioc-devel and these rules may misbehave until the next Bioconductor release.__
Manual compilation using `workflowr::wflow_publish('path/to/rmd', 'Commit Message')` will be effective if using `ngsReports` > 1.5.4.

## Outline

The steps currently implemented are

1. Download the DNA reference, GTF file
2. Generate the STAR index for the reference
3. Trim files using `AdapterRemoval`
4. Align Trimmed files using STAR
5. Count reads using `featureCounts`
6. Compile initial QC Rmarkdown documents using `workflowr`

Additional tasks performed are

1. Initialising a `workflowr` directory structure
2. FastQC reports are also generated for both raw and trimmed reads
3. The `rulegraph` will be generated as a dot file and a pdf

Snakemake rules rely heavily on those provided at https://github.com/snakemake/snakemake-wrappers.
Version numbers for these are also available at https://snakemake-wrappers.readthedocs.io/

## Starting A New Analysis

This is a Template repository so to begin a new analysis:

1. Click `Use this template` at the top of the page
2. Decide on your new repository name
3. Clone to your local machine
4. Add an R-project file by opening a new project in the root directory
5. Use `git push` to update the remote repository
6. Clone to `phoenix`

This order ensures that we have a `git` repository ready to go and an `Rproj` file.
The report-writing steps from `workflowr` require both of these things to have been performed beforehand.

## Essential Files

In order to run this workflow, please ensure that you have

1. Placed unprocessed fastq files in the directory `data/raw/fastq`
2. Placed a `tsv` file (usually called `samples.tsv`) in the `config` folder
    + This file **must** contain a column called `sample`
3. Edited `config.yml` in the `config folder to ensure all parameters are correct

## Running the workflow

Once you have placed your data in `data/faw/fastq`, edited the file `samples.tsv` and checked `config/config.yml`, please run the workflow as follows:

First, activate your `snakemake` conda/mamba environment

```
micromamba activate snakemake
```

Now check you have everything correct

```
snakemake -n
```

Building the conda environments can also be helpful in advance as you can be sure every environment builds successfully.
Unfortunately, problems with internet connections whilst building environments can be terminal for a slurm job.
This only needs to be performed once, and is optional but advised.

```
snakemake --use-conda --create-envs-only
```

### Running on `phoenix`

To run the complete workflow on phoenix, edit the script `scripts/run_snakemake.sh`:

1. Ensure you have suitable resources
2. Ensure you have modified the `slurm` output directories
3. Ensure you have modified the project root directory

Now submit the script to the queueing system

```
sbatch scripts/run_snakemake.sh
```

### Running locally

If you are running locally, there is no need to use the job submission script, so just run the following with the most appropriate number of threads.

```
snakemake --cores 8 --use-conda
```
