#!/bin/bash

# Name for the job that will be visible in the job queue and accounting tools.
#SBATCH --job-name count_ferret_rxn2

# Name of the SLURM partition that this job should run on.
#SBATCH -p 256GB       # partition (queue)
# Number of nodes required to run this job
#SBATCH -N 1



# Time limit for the job in the format Days-H:M:S
# A job that reaches its time limit will be cancelled.
# Specify an accurate time limit for efficient scheduling so your job runs promptly.
#SBATCH -t 7-00:0:00

# The standard output and errors from commands will be written to these files.
# %j in the filename will be replace with the job number when it is submitted.
#SBATCH -o count_ferret_rxn2_%j.out
#SBATCH -e count_ferret_rxn2_%j.err

# packages
module load cellranger/8.0.1

# change dir
#cd /project/Neuroinformatics_Core/Konopka_lab/gkonop/00_BAM_DOWNLOADED/KRIENEN_2020/FERRET/SRR11921037/cellranger_count/

# count
#cellranger count --id=repaired_krienen_ferret_rxn1 \
#           --transcriptome=/home2/gkonop/projectShared/For_Gozde/cellranger_reference_RNA/Mustela_putorius_furo_genome \
#           --fastqs=/home2/gkonop/project/00_BAM_DOWNLOADED/KRIENEN_2020/FERRET/SRR11921037/fastqs/ \
#           --sample=repaired_krienen_ferret_rxn1 \
#           --create-bam=true \
#           --localcores=8 \
#           --localmem=100  
     
# change dir
cd /project/Neuroinformatics_Core/Konopka_lab/gkonop/00_BAM_DOWNLOADED/KRIENEN_2020/FERRET/SRR11921038/cellranger_count/      

# count
cellranger count --id=repaired_krienen_ferret_rxn2 \
           --transcriptome=/home2/gkonop/projectShared/For_Gozde/cellranger_reference_RNA/Mustela_putorius_furo_genome \
           --fastqs=/home2/gkonop/project/00_BAM_DOWNLOADED/KRIENEN_2020/FERRET/SRR11921038/fastqs/ \
           --sample=repaired_krienen_ferret_rxn2 \
           --create-bam=true \
           --localcores=8 \
           --localmem=100  
