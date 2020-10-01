#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 6
#SBATCH --time=1:00:00
#SBATCH --mem=48GB
#SBATCH -o /home/a1018048/slurm/snakemake_rnaseq/%x_%j.out
#SBATCH -e /home/a1018048/slurm/snakemake_rnaseq/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stephen.pederson@adelaide.edu.au

## Cores
CORES=12
if [ -d "/hpcfs" ]; then
	module load arch/arch/haswell
	module load arch/haswell
	module load modulefiles/arch/haswell
	HPC="/hpcfs"
else
    if [ -d "/fast" ]; then
        HPC=/fast
    else
        exit 1
    fi
fi

## Project Root
PROJ=${HPC}/users/a1018048/snakemake_rnaseq

## The environment containing snakemake
micromamba activate snakemake
cd ${PROJ}

## Create dot and pdf files for visualisation
snakemake --dag > output/dag.dot
dot -Tpdf output/dag.dot > output/dag.pdf
snakemake --rulegraph > output/rulegraph.dot
dot -Tpdf output/rulegraph.dot > output/rulegraph.pdf

## Run snakemake
snakemake \
  --cores ${CORES} \
  --use-conda \
  --wrapper-prefix 'https://raw.githubusercontent.com/snakemake/snakemake-wrappers/'

## Add files to git
bash ${PROJ}/scripts/update_git.sh
