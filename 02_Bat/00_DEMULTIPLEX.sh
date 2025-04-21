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
#SBATCH --job-name nov19

# Name of the SLURM partition that this job should run on.
#SBATCH -p 256GB       # partition (queue)
# Number of nodes required to run this job
#SBATCH -N 1

# Time limit for the job in the format Days-H:M:S
# A job that reaches its time limit will be cancelled.
# Specify an accurate time limit for efficient scheduling so your job runs promptly.
#SBATCH -t 1-0:0:00

# The standard output and errors from commands will be written to these files.
# %j in the filename will be replace with the job number when it is submitted.
#SBATCH -o nov19_%j.out
#SBATCH -e nov19_%j.err

# Do not use an earlier version
module load bcl2fastq/2.20.0

FLOWCELL_DIR="/home2/gkonop/project/00_BCL/NOV31_13135_091923/18035-83/18035-83-run-081123-XP-fc2-Lane23-08142023_211619"
OUTPUT_DIR="/home2/gkonop/project/00_BCL/NOV31_13135_091923/18035-83/FASTQ_NOV31"

INPUT_DIR="/home2/gkonop/project/00_BCL/NOV31_13135_091923/18035-83/18035-83-run-081123-XP-fc2-Lane23-08142023_211619"

SAMPLE_SHEET_PATH="/home2/gkonop/project/00_BCL/NOV31_13135_091923/18035-83/indexsheet_NOV31_L3.csv"

BASES_MASK="Y150n*,I10n*,n*,Y150n*"

bcl2fastq --use-bases-mask=${BASES_MASK} \
          --create-fastq-for-index-reads \
          --minimum-trimmed-read-length=8 \
          --mask-short-adapter-reads=8 \
          --ignore-missing-positions \
          --ignore-missing-controls \
          --ignore-missing-filter \
          --ignore-missing-bcls \
          -R ${FLOWCELL_DIR} \
          -r 23 -w 6 \
          --output-dir=${OUTPUT_DIR} \
          --input-dir=${INPUT_DIR} \
          --sample-sheet=${SAMPLE_SHEET_PATH}

# For Another Lane
FLOWCELL_DIR="/home2/gkonop/project/00_BCL/NOV31_13135_091923/18035-83/18035-83-run-081023-XP-fc1-Lane78-08142023_211619"

INPUT_DIR="/home2/gkonop/project/00_BCL/NOV31_13135_091923/18035-83/18035-83-run-081023-XP-fc1-Lane78-08142023_211619"

SAMPLE_SHEET_PATH="/home2/gkonop/project/00_BCL/NOV31_13135_091923/18035-83/indexsheet_NOV31_L8.csv"

BASES_MASK="Y150n*,I10n*,n*,Y150n*"

bcl2fastq --use-bases-mask=${BASES_MASK} \
          --create-fastq-for-index-reads \
          --minimum-trimmed-read-length=8 \
          --mask-short-adapter-reads=8 \
          --ignore-missing-positions \
          --ignore-missing-controls \
          --ignore-missing-filter \
          --ignore-missing-bcls \
          -R ${FLOWCELL_DIR} \
          -r 23 -w 6 \
          --output-dir=${OUTPUT_DIR} \
          --input-dir=${INPUT_DIR} \
          --sample-sheet=${SAMPLE_SHEET_PATH}


# Move the files to the FASTQ path
rm *I*
rm Undetermined*
mv *R*gz ~/project/00_FASTQ/BAT_WT/

######

FLOWCELL_DIR="/home2/gkonop/project/00_BCL/NOV31_13135_091923/18035-83/18035-83-run-081123-XP-fc2-Lane23-08142023_211619"

OUTPUT_DIR="/home2/gkonop/project/00_BCL/NOV31_13135_091923/18035-83/FASTQ_NOV31"

INPUT_DIR="/home2/gkonop/project/00_BCL/NOV31_13135_091923/18035-83/18035-83-run-081123-XP-fc2-Lane23-08142023_211619"

SAMPLE_SHEET_PATH="/home2/gkonop/project/00_BCL/NOV31_13135_091923/18035-83/indexsheet_NOV31_L3_EMMA.csv"

BASES_MASK="Y150n*,I10n*,n*,Y150n*"

bcl2fastq --use-bases-mask=${BASES_MASK} \
          --create-fastq-for-index-reads \
          --minimum-trimmed-read-length=8 \
          --mask-short-adapter-reads=8 \
          --ignore-missing-positions \
          --ignore-missing-controls \
          --ignore-missing-filter \
          --ignore-missing-bcls \
          -R ${FLOWCELL_DIR} \
          -r 23 -w 6 \
          --output-dir=${OUTPUT_DIR} \
          --input-dir=${INPUT_DIR} \
          --sample-sheet=${SAMPLE_SHEET_PATH}



# For Another Lane
FLOWCELL_DIR="/home2/gkonop/project/00_BCL/NOV31_13135_091923/18035-83/18035-83-run-081023-XP-fc1-Lane78-08142023_211619"

OUTPUT_DIR="/home2/gkonop/project/00_BCL/NOV31_13135_091923/18035-83/FASTQ_NOV31"

INPUT_DIR="/home2/gkonop/project/00_BCL/NOV31_13135_091923/18035-83/18035-83-run-081023-XP-fc1-Lane78-08142023_211619"

SAMPLE_SHEET_PATH="/home2/gkonop/project/00_BCL/NOV31_13135_091923/18035-83/indexsheet_NOV31_L8_EMMA.csv"

BASES_MASK="Y150n*,I10n*,n*,Y150n*"

bcl2fastq --use-bases-mask=${BASES_MASK} \
          --create-fastq-for-index-reads \
          --minimum-trimmed-read-length=8 \
          --mask-short-adapter-reads=8 \
          --ignore-missing-positions \
          --ignore-missing-controls \
          --ignore-missing-filter \
          --ignore-missing-bcls \
          -R ${FLOWCELL_DIR} \
          -r 23 -w 6 \
          --output-dir=${OUTPUT_DIR} \
          --input-dir=${INPUT_DIR} \
          --sample-sheet=${SAMPLE_SHEET_PATH}


# Move the files to the FASTQ path
rm *I*
rm Undetermined*
mv *R*gz ~/project/00_FASTQ/BAT_WT/




