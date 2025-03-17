#!/bin/bash
#
# CREATED USING THE BIOHPC PORTAL on Tue Apr 28 2020 00:37:21 GMT-0500 (Central Daylight Time)
#
# This file is batch script used to run commands on the BioHPC cluster.
# The script is submitted to the cluster using the SLURM `sbatch` command.
# Lines starting with # are comments, and will not be run.
# Lines starting with #SBATCH specify options for the scheduler.
# Lines that do not start with # or #SBATCH are commands that will run.

# Name for the job that will be visible in the job queue and accounting tools.
#SBATCH --job-name fastq

# Name of the SLURM partition that this job should run on.
#SBATCH -p super       # partition (queue)
# Number of nodes required to run this job
#SBATCH -N 1

# Time limit for the job in the format Days-H:M:S
# A job that reaches its time limit will be cancelled.
# Specify an accurate time limit for efficient scheduling so your job runs promptly.
#SBATCH -t 2-0:0:00

# The standard output and errors from commands will be written to these files.
# %j in the filename will be replace with the job number when it is submitted.
#SBATCH -o fastq_%j.out
#SBATCH -e fastq_%j.err

# Do not use an earlier version
module load cellranger/6.0.0

cd /home2/gkonop/project/00_BCL/NHP_STRIATUM_040723
tar -xvf 230322_A00330_0146_AHMVTJDSX5.tar

#mv scratch/s167561/illumina_hiseq/230123_A00672_0126_AHG3JNDSX5 NOV19_12545_HUMAN_BASAL_GANGLIA
#rm -rf scratch/
#cd NOV19_12545_HUMAN_BASAL_GANGLIA

cellranger mkfastq --run=230322_A00330_0146_AHMVTJDSX5 \
                        --csv=indexsheet.csv














