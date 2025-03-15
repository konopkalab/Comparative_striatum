#!/bin/bash


# Download the marmoset genome
wget http://hgdownload.soe.ucsc.edu/goldenPath/calJac3/bigZips/calJac3.fa.gz

# Subset the small chromosomes to test if the tool works
cd /home2/gkonop/workdir/PROGRAMS
git clone https://github.com/lh3/seqtk.git;
cd seqtk; make
seqtk subseq test.fa test.txt > hg38_chr21_22.fa

# Run liftoff on the subset
cd /endosome/work/Neuroinformatics_Core/gkonop/REFERENCES/TMP
liftoff -g hg38.gtf calJac3_chr21_22.fa hg38_chr21_22.fa -m /endosome/work/Neuroinformatics_Core/gkonop/PROGRAMS/minimap2-2.24_x64-linux/minimap2 -o hg38_to_calJac3.gtf

# Run liftoff on the whole genome
cd /endosome/work/Neuroinformatics_Core/gkonop/REFERENCES
liftoff -g '/endosome/work/Neuroinformatics_Core/gkonop/REFERENCES/cellranger_reference_ATAC/Hsa_GRCh38/HomSap_GRCh38/genes/genes.gtf' FASTA_ONLY/calJac3.fa /endosome/work/Neuroinformatics_Core/gkonop/REFERENCES/cellranger_reference_ATAC/Hsa_GRCh38/HomSap_GRCh38/fasta/genome.fa -m /endosome/work/Neuroinformatics_Core/gkonop/PROGRAMS/minimap2-2.24_x64-linux/minimap2 -o GTF_ONLY/hg38_liftoverTo_calJac3.gtf



