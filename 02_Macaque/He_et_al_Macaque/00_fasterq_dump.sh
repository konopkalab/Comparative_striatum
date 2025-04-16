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


# Echo GEO accession number
echo 'GEO accession number: GSE167920'

# Go to the NCBI GEO website and download the metadata associated with all fastq files (Run Selector)

# Load the modules
module load sra_toolkit/3.0.2
module load pigz

# Go to directory
cd /home2/gkonop/project/00_FASTQ/MACAQUE_HE2021

# Write all SRRs to be downloaded in a file
grep 'Putamen' He2021_GSE167920_Metadata.txt | cut -f1 -d, > todownload.txt

# Download
for i in `cat todownload.txt`
do
fasterq-dump -e 23 --progress --include-technical --split-files $i
pigz $i*
done

# Reshape for cellranger count
#mv SRR13551786_2.fastq.gz SRR13551786_S1_L001_R1_001.fastq.gz
#mv SRR13551786_3.fastq.gz SRR13551786_S1_L001_R2_001.fastq.gz

#mv SRR13551787_2.fastq.gz SRR13551787_S1_L001_R1_001.fastq.gz
#mv SRR13551787_3.fastq.gz SRR13551787_S1_L001_R2_001.fastq.gz

#mv SRR13551794_2.fastq.gz SRR13551794_S1_L001_R1_001.fastq.gz
#mv SRR13551794_3.fastq.gz SRR13551794_S1_L001_R2_001.fastq.gz

#mv SRR13551795_2.fastq.gz SRR13551795_S1_L001_R1_001.fastq.gz
#mv SRR13551795_3.fastq.gz SRR13551795_S1_L001_R2_001.fastq.gz


# Merge resequencing of the same library preps
#cat SRR13551786_S1_L001_R1_001.fastq.gz SRR13551787_S1_L001_R1_001.fastq.gz > v321_S1_L001_R1_001.fastq.gz
#cat SRR13551786_S1_L001_R2_001.fastq.gz SRR13551787_S1_L001_R2_001.fastq.gz > v321_S1_L001_R2_001.fastq.gz

#cat SRR13551794_S1_L001_R1_001.fastq.gz SRR13551795_S1_L001_R1_001.fastq.gz > v324_S1_L001_R1_001.fastq.gz
#cat SRR13551794_S1_L001_R2_001.fastq.gz SRR13551795_S1_L001_R2_001.fastq.gz > v324_S1_L001_R2_001.fastq.gz
