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
#SBATCH --job-name pr3_sc2

# Name of the SLURM partition that this job should run on.
#SBATCH -p super       # partition (queue)
# Number of nodes required to run this job
#SBATCH -N 1

# Memory (RAM) requirement/limit in MB.
#SBATCH --mem 254976      # Memory Requirement (MB)

# Time limit for the job in the format Days-H:M:S
# A job that reaches its time limit will be cancelled.
# Specify an accurate time limit for efficient scheduling so your job runs promptly.
#SBATCH -t 1-0:0:00

# The standard output and errors from commands will be written to these files.
# %j in the filename will be replace with the job number when it is submitted.
#SBATCH -o pr3_sc2_%j.out
#SBATCH -e pr3_sc2_%j.err

module load cellranger/6.0.0

cd /home2/gkonop/workdir/REFERENCES/cellranger_reference_RNA
mkdir panTro5_human_liftoff
cd panTro5_human_liftoff

# Create cellranger reference
cellranger mkref --genome=reference_genome --fasta=~/workdir/REFERENCES/FASTA_ONLY/panTro5.fa --genes=~/workdir/REFERENCES/GTF_ONLY/hg38_liftoverTo_panTro5.gtf









