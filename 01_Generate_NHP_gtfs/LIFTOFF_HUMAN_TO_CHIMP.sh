#!/bin/bash


# Download the chimp genome
cd /home2/gkonop/workdir/REFERENCES/FASTA_ONLY
wget https://hgdownload.soe.ucsc.edu/goldenPath/panTro5/bigZips/panTro5.fa.gz
gunzip panTro5.fa.gz

# Subset the small chromosomes (to test if the tool works)
cd /home2/gkonop/workdir/REFERENCES/TMP
/home2/gkonop/workdir/PROGRAMS/seqtk/seqtk subseq /endosome/work/Neuroinformatics_Core/gkonop/REFERENCES/cellranger_reference_ATAC/Hsa_GRCh38/HomSap_GRCh38/fasta/genome.fa test.txt > hg38_chr20.fa
/home2/gkonop/workdir/PROGRAMS/seqtk/seqtk subseq /home2/gkonop/workdir/REFERENCES/FASTA_ONLY/panTro5.fa test.txt > panTro5_chr20.fa

# Subset the gtf file (to test if the tool works)
grep 'chr20' /endosome/work/Neuroinformatics_Core/gkonop/REFERENCES/cellranger_reference_ATAC/Hsa_GRCh38/HomSap_GRCh38/genes/genes.gtf > hg38_chr20.gtf

# Run liftoff on the subset
liftoff -g hg38_chr20.gtf panTro5_chr20.fa hg38_chr20.fa -m /endosome/work/Neuroinformatics_Core/gkonop/PROGRAMS/minimap2-2.24_x64-linux/minimap2 -o hg38_to_panTro5.gtf

# Run liftoff on the whole genome
cd /endosome/work/Neuroinformatics_Core/gkonop/REFERENCES
liftoff -g '/endosome/work/Neuroinformatics_Core/gkonop/REFERENCES/cellranger_reference_ATAC/Hsa_GRCh38/HomSap_GRCh38/genes/genes.gtf' FASTA_ONLY/panTro5.fa /endosome/work/Neuroinformatics_Core/gkonop/REFERENCES/cellranger_reference_ATAC/Hsa_GRCh38/HomSap_GRCh38/fasta/genome.fa -m /endosome/work/Neuroinformatics_Core/gkonop/PROGRAMS/minimap2-2.24_x64-linux/minimap2 -o GTF_ONLY/hg38_liftoverTo_panTro5.gtf
