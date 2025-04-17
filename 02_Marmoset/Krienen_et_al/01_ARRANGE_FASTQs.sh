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

cd /home2/gkonop/project/00_FASTQ/MARMOSET_KRIENEN/PUTAMEN

#cd marm022
#cat */*read1* > marm022_S1_L001_R1_001.fastq.gz
#cat */*read2* > marm022_S1_L001_R2_001.fastq.gz

#cd marm023
#for i in *tar; do tar -xvf $i; done
#rm *.tar
#cat *read1* > marm023_S1_L001_R1_001.fastq.gz
#cat *read2* > marm023_S1_L001_R2_001.fastq.gz

# Concatanate and remove all small files
cd marm029_rxn1
for i in *tar; do tar -xvf $i; done
cat */*read1* > marm029_RXN1_S1_L001_R1_001.fastq.gz
cat */*read2* > marm029_RXN1_S1_L001_R2_001.fastq.gz
rm *.tar
ls -I "*fastq.gz" | xargs rm -rf


cd marm029_rxn2
for i in *tar; do tar -xvf $i; done
cat */*read1* > marm029_RXN2_S1_L001_R1_001.fastq.gz
cat */*read2* > marm029_RXN2_S1_L001_R2_001.fastq.gz
rm *.tar
ls -I "*fastq.gz" | xargs rm -rf


# Carry the files to separate folders
mkdir marm027_rxn1
mkdir marm027_rxn2
mv TMP/*marm027*rxn1* marm027_rxn1/
mv TMP/*marm027*rxn2* marm027_rxn2/

cd marm027_rxn1
for i in *tar; do tar -xvf $i; done
cat */*read1* > marm027_RXN1_S1_L001_R1_001.fastq.gz
cat */*read2* > marm027_RXN1_S1_L001_R2_001.fastq.gz

rm *.tar
ls -I "*fastq.gz" | xargs rm -rf

cd marm027_rxn2
for i in *tar; do tar -xvf $i; done
cat */*read1* > marm027_RXN2_S1_L001_R1_001.fastq.gz
cat */*read2* > marm027_RXN2_S1_L001_R2_001.fastq.gz

rm *.tar
ls -I "*fastq.gz" | xargs rm -rf
