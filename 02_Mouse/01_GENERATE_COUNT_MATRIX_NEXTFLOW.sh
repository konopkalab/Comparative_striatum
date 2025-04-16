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
#SBATCH --job-name preprocess

# Name of the SLURM partition that this job should run on.
#SBATCH -p 256GB       # partition (queue)
# Number of nodes required to run this job
#SBATCH -N 1

# Memory (RAM) requirement/limit in MB.
#SBATCH --mem 252928      # Memory Requirement (MB)

# Time limit for the job in the format Days-H:M:S
# A job that reaches its time limit will be cancelled.
# Specify an accurate time limit for efficient scheduling so your job runs promptly.
#SBATCH -t 2-0:0:00

# The standard output and errors from commands will be written to these files.
# %j in the filename will be replace with the job number when it is submitted.
#SBATCH -o preprocess_%j.out
#SBATCH -e preprocess_%j.err


cd /home2/gkonop/project/BAM_DOWNLOADED/KRIENEN_2020/MOUSE
gtfdir='/home2/gkonop/workdir/REFERENCES/GTF_ONLY/Mus_musculus.GRCm38.102.chr.gtf.gz'
tmpdir='/home2/gkonop/project/TMP'
pubdir='/home2/gkonop/project/PREPROCESS'
tomtx='/home2/gkonop/SCRIPTS/TO_MTX_DROPSEQ.R'


# FOR SRR TO COUNT MATRIX
for i in `ls -d SRR*`
do

/home2/gkonop/workdir/PROGRAMS/./nextflow-22.04.4-all run ~/SCRIPTS/MOUSE_SCRIPTS/Mouse_Aligned_To_Count.nf \
		--bamfile /home2/gkonop/project/BAM_DOWNLOADED/KRIENEN_2020/MOUSE/${i}/*.bam \
		--tmpdir $tmpdir \
		--pubdir $pubdir \
		--tomtx $tomtx \
		--pref $i

done
