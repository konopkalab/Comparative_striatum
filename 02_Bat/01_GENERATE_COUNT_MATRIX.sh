#!/bin/bash
#
# CREATED USING THE BIOHPC PORTAL on Sat Jan 22 2022 23:35:26 GMT-0600 (Central Standard Time)
#
# This file is batch script used to run commands on the BioHPC cluster.
# The script is submitted to the cluster using the SLURM `sbatch` command.
# Lines starting with # are comments, and will not be run.
# Lines starting with #SBATCH specify options for the scheduler.
# Lines that do not start with # or #SBATCH are commands that will run.

# Name for the job that will be visible in the job queue and accounting tools.
#SBATCH --job-name count

# Name of the SLURM partition that this job should run on.
#SBATCH -p 256GB       # partition (queue)
# Number of nodes required to run this job
#SBATCH -N 1

# Time limit for the job in the format Days-H:M:S
# A job that reaches its time limit will be cancelled.
# Specify an accurate time limit for efficient scheduling so your job runs promptly.
#SBATCH -t 2-0:0:00

# The standard output and errors from commands will be written to these files.
# %j in the filename will be replace with the job number when it is submitted.
#SBATCH -o count_%j.out
#SBATCH -e count_%j.err

module load cellranger/6.1.2
cd /home2/gkonop/project/01_CELLRANGER_OUT/BAT_OUT

refdir="/home2/gkonop/workdir/REFERENCES/FROM_DEVIN_BAT_REF/phyllostomus-discolor/phyllostomus-discolor"
fastqdir="/home2/gkonop/project/00_FASTQ/BAT_WT"

# Run all caudates
cellranger count --id=Put_3 --include-introns --transcriptome=$refdir --fastqs=$fastqdir --sample=Put_3
cellranger count --id=Put_4 --include-introns --transcriptome=$refdir --fastqs=$fastqdir --sample=Put_4
cellranger count --id=Caud_4 --include-introns --transcriptome=$refdir --fastqs=$fastqdir --sample=Caud_4
#####
cd /home2/gkonop/project/00_BCL/NOV31_13135_091923/18035-83

refdir="/home2/gkonop/workdir/REFERENCES/cellranger_reference_ARC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
libdir="/home2/gkonop/project/00_BCL/NOV31_13135_091923/18035-83/libraries_63C468.csv"

cellranger-arc count --id=63C468_Cerebellum_CR_OUT \
                       --reference=$refdir \
                       --libraries=$libdir \
                       --localcores=23 \
                       --localmem=256
