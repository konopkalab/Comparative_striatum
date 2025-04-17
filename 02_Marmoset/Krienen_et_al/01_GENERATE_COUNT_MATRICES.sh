#!/bin/bash
#
# CREATED USING THE BIOHPC PORTAL on Fri Jun 21 2019 22:56:40 GMT-0500 (Central Daylight Time)
#
# This file is batch script used to run commands on the BioHPC cluster.
# The script is submitted to the cluster using the SLURM `sbatch` command.
# Lines starting with # are comments, and will not be run.
# Lines starting with #SBATCH specify options for the scheduler.
# Lines that do not start with # or #SBATCH are commands that will run.

# Name for the job that will be visible in the job queue and accounting tools.
#SBATCH --job-name Sample_M1

# Name of the SLURM partition that this job should run on.
#SBATCH -p 256GB       # partition (queue)
# Number of nodes required to run this job
#SBATCH -N 1

# Memory (RAM) requirement/limit in MB.
#SBATCH --mem 254976      # Memory Requirement (MB)

# Time limit for the job in the format Days-H:M:S
# A job that reaches its time limit will be cancelled.
# Specify an accurate time limit for efficient scheduling so your job runs promptly.
#SBATCH -t 2-0:0:00

# The standard output and errors from commands will be written to these files.
# %j in the filename will be replace with the job number when it is submitted.
#SBATCH -o count_Sample_M1_%j.out
#SBATCH -e count_Sample_M1_%j.err

module load cellranger/6.0.0

cd /home2/gkonop/project/01_CELLRANGER_OUT/MARMOSET_KRIENEN_FASTQ/PUTAMEN

refdir="/home2/gkonop/workdir/REFERENCES/cellranger_reference_RNA/caljac3_human_liftoff/reference_genome"
fastqdir="/home2/gkonop/project/00_FASTQ/MARMOSET_KRIENEN/PUTAMEN/marm027_rxn1"

# Run cellranger
cellranger count --id=celllranger_count_marm027_RXN1 --transcriptome=$refdir --fastqs=$fastqdir --sample=marm027_RXN1 --include-introns


fastqdir="/home2/gkonop/project/00_FASTQ/MARMOSET_KRIENEN/PUTAMEN/marm027_rxn2"

# Run cellranger
cellranger count --id=celllranger_count_marm027_RXN2 --transcriptome=$refdir --fastqs=$fastqdir --sample=marm027_RXN2 --include-introns


# Run cellranger
fastqdir="/home2/gkonop/project/00_FASTQ/MARMOSET_KRIENEN/CAUDATE/marm027_rxn1"

cellranger count --id=celllranger_count_marm027_RXN1 --transcriptome=$refdir --fastqs=$fastqdir --sample=marm027_RXN1 --include-introns

# Run cellranger
fastqdir="/home2/gkonop/project/00_FASTQ/MARMOSET_KRIENEN/CAUDATE/marm027_rxn2"

cellranger count --id=celllranger_count_marm027_RXN2 --transcriptome=$refdir --fastqs=$fastqdir --sample=marm027_RXN2 --include-introns


# Run cellranger
fastqdir="/home2/gkonop/project/00_FASTQ/MARMOSET_KRIENEN/CAUDATE/marm028_rxn1"

cellranger count --id=celllranger_count_marm028_RXN1 --transcriptome=$refdir --fastqs=$fastqdir --sample=marm028_RXN1 --include-introns

# Run cellranger
fastqdir="/home2/gkonop/project/00_FASTQ/MARMOSET_KRIENEN/CAUDATE/marm028_rxn2"

cellranger count --id=celllranger_count_marm028_RXN2 --transcriptome=$refdir --fastqs=$fastqdir --sample=marm028_RXN2 --include-introns

