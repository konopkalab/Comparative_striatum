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
#SBATCH --job-name Samp_M1

# Name of the SLURM partition that this job should run on.
#SBATCH -p 256GB       # partition (queue)
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
#SBATCH -o count_Samp_M1_%j.out
#SBATCH -e count_Samp_M1_%j.err

module load cellranger/6.0.0

cd /home2/gkonop/project/01_CELLRANGER_OUT/MACAQUE_OUR

refdir="/home2/gkonop/workdir/REFERENCES/cellranger_reference_RNA/rheMac10_human_liftoff/reference_genome"
fastqdir="/home2/gkonop/project/00_BCL/NHP_STRIATUM_050323/2023_05_03_NOV6_12830_0/FASTQ_NOV6"

# Run cellranger
cellranger count --id=celllranger_count_Sample_30 --transcriptome=$refdir --fastqs=$fastqdir --sample=Sample_30 --include-introns

cellranger count --id=celllranger_count_Sample_26 --transcriptome=$refdir --fastqs=$fastqdir --sample=Sample_26 --include-introns

cellranger count --id=celllranger_count_Sample_29 --transcriptome=$refdir --fastqs=$fastqdir --sample=Sample_29 --include-introns

cellranger count --id=celllranger_count_Sample_33 --transcriptome=$refdir --fastqs=$fastqdir --sample=Sample_33 --include-introns

cellranger count --id=celllranger_count_Sample_38 --transcriptome=$refdir --fastqs=$fastqdir --sample=Sample_38 --include-introns

cellranger count --id=celllranger_count_Sample_48 --transcriptome=$refdir --fastqs=$fastqdir --sample=Sample_48 --include-introns

cellranger count --id=celllranger_count_Sample_53 --transcriptome=$refdir --fastqs=$fastqdir --sample=Sample_53 --include-introns

cellranger count --id=celllranger_count_Sample_25 --transcriptome=$refdir --fastqs=$fastqdir --sample=Sample_25 --include-introns

cellranger count --id=celllranger_count_Sample_42 --transcriptome=$refdir --fastqs=$fastqdir --sample=Sample_42 --include-introns

cellranger count --id=celllranger_count_Sample_34 --transcriptome=$refdir --fastqs=$fastqdir --sample=Sample_34 --include-introns

cellranger count --id=celllranger_count_Sample_43 --transcriptome=$refdir --fastqs=$fastqdir --sample=Sample_43 --include-introns

cellranger count --id=celllranger_count_Sample_37 --transcriptome=$refdir --fastqs=$fastqdir --sample=Sample_37 --include-introns

cellranger count --id=celllranger_count_Sample_47 --transcriptome=$refdir --fastqs=$fastqdir --sample=Sample_47 --include-introns

cellranger count --id=celllranger_count_Sample_52 --transcriptome=$refdir --fastqs=$fastqdir --sample=Sample_52 --include-introns
