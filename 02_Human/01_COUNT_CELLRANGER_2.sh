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
#SBATCH --job-name count_put

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
#SBATCH -o count_put_%j.out
#SBATCH -e count_put_%j.err

module load cellranger/3.1.0

cd /home2/gkonop/project/00_FASTQ/HUMAN/COUNT

refdir="/home2/gkonop/workdir/REFERENCES/cellranger_reference_RNA/GRCh38-1.2.0_premrna"
fastqdir="/home2/gkonop/project/00_FASTQ/HUMAN/FASTQ"

# Run all caudates
cellranger count --id=Sample_242973 --transcriptome=$refdir --fastqs=$fastqdir --sample=Sample_242973

cellranger count --id=Sample_242975 --transcriptome=$refdir --fastqs=$fastqdir --sample=Sample_242975

cellranger count --id=Sample_242976 --transcriptome=$refdir --fastqs=$fastqdir --sample=Sample_242976

cellranger count --id=Sample_242981 --transcriptome=$refdir --fastqs=$fastqdir --sample=Sample_242981

cellranger count --id=Sample_242983 --transcriptome=$refdir --fastqs=$fastqdir --sample=Sample_242983

cellranger count --id=Sample_242989 --transcriptome=$refdir --fastqs=$fastqdir --sample=Sample_242989

cellranger count --id=Sample_242990 --transcriptome=$refdir --fastqs=$fastqdir --sample=Sample_242990

cellranger count --id=Sample_242993 --transcriptome=$refdir --fastqs=$fastqdir --sample=Sample_242993







