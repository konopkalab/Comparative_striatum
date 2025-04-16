# load packages
module load samtools/gcc/1.10
module load cellranger/8.0.1
module load picard/2.10.3
module load seqtk/1.2-r94
module load bbmap/38.46

# convert bam file to fastq using awk

#!/bin/bash
# packages
module load samtools/gcc/1.10
 
# go to the correct directory
#cd /home2/gkonop/project/00_BAM_DOWNLOADED/KRIENEN_2020/FERRET/SRR11921037
cd /home2/gkonop/project/00_BAM_DOWNLOADED/KRIENEN_2020/FERRET/SRR11921038

# Convert BAM to SAM
#samtools view -h ferretP42_str_nova_rxn1.bam > input.sam
# Convert BAM to SAM
samtools view -h ferretP42_str_nova_rxn2.bam > input.sam

# Generate Read 1 and Read 2 FASTQ Files
awk '
BEGIN { OFS = "\t" }
/^@/ { next }
{
    read_id = "@"$1
    seq = $10
    qual = $11

    # Extract cell barcode and UMI from the XC and XM tags
    for (i = 12; i <= NF; i++) {
        if ($i ~ /^XC:Z:/) {
            cell_barcode = substr($i, 6)
        } else if ($i ~ /^XM:Z:/) {
            umi = substr($i, 6)
        }
    }

    # Read 1: Cell barcode and UMI
    read1_seq = cell_barcode umi
    read1_qual = substr(qual, 1, length(read1_seq))

    # Read 2: Original sequence and quality
    read2_seq = seq
    read2_qual = qual

    ##############################
    # Print Read 1 in FASTQ format
    print read_id "\n" read1_seq "\n+\n" read1_qual > "./fastqs/krienen_ferret_rxn2_S1_L001_R1_001.fastq"

    # Print Read 2 in FASTQ format
    print read_id "\n" read2_seq "\n+\n" read2_qual > "./fastqs/krienen_ferret_rxn2_S1_L001_R2_001.fastq"
}
' input.sam

# Compress the FASTQ Files
#gzip read1.fastq
#gzip read2.fastq

# change the names
#cd /home2/gkonop/project/00_BAM_DOWNLOADED/KRIENEN_2020/FERRET/SRR11921038/fastqs
#mv read1.fastq krienen_ferret_rxn2_S1_L001_R1_001.fastq
#mv read2.fastq krienen_ferret_rxn2_S1_L001_R2_001.fastq

# use bbmap suite to remove the unpaired reads
# rxn1
repair.sh in1=/project/Neuroinformatics_Core/Konopka_lab/gkonop/00_BAM_DOWNLOADED/KRIENEN_2020/FERRET/SRR11921037/fastqs/krienen_ferret_rxn1_filtered_S1_L001_R1_001.fastq in2=/project/Neuroinformatics_Core/Konopka_lab/gkonop/00_BAM_DOWNLOADED/KRIENEN_2020/FERRET/SRR11921037/fastqs/krienen_ferret_rxn1_filtered_S1_L001_R2_001.fastq out1=repaired_R1.fq.gz out2=repaired_R2.fq.gz

# rxn2
repair.sh in1=/project/Neuroinformatics_Core/Konopka_lab/gkonop/00_BAM_DOWNLOADED/KRIENEN_2020/FERRET/SRR11921038/fastqs/filtered_krienen_ferret_rxn2_S1_L001_R1_001.fastq in2=/project/Neuroinformatics_Core/Konopka_lab/gkonop/00_BAM_DOWNLOADED/KRIENEN_2020/FERRET/SRR11921038/fastqs/filtered_krienen_ferret_rxn2_S1_L001_R2_001.fastq out1=repaired_krienen_ferret_rxn2_S1_L001_R1_001.fq.gz out2=repaired_krienen_ferret_rxn2_S1_L001_R2_001.fq.gz
