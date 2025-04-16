rm(list = ls())
library(dplyr)
library(plyr)
library(Seurat)
library(ggpubr)
library(reshape2)
library(tidyr)
library(tidyverse)
library(DropletQC)
library(Matrix)
library(EnsDb.Mmusculus.v79)
source('~/SCRIPTS/utility_functions.R')

samps = paste0('SRR1192', 1005:1012)
folddir = '/home2/gkonop/project/01_DROPSEQ_PREPROCESS/KRIENEN_MOUSE'
outdir = '~/project/02_MATRIX_FOR_CELLBENDER/KRIENEN_MOUSE'

####
## Loop through the samples
####

for(i in 2:length(samps)){

	# Genes
	gns = read.table(paste0(folddir, '/', samps[i], '_count_gns.txt'), sep = '\t')
	gnsids = genes(EnsDb.Mmusculus.v79) %>% as.data.frame()
	gns$V2 = gns$V1
	gns$V1 = gnsids[match(gns[,1], gnsids$symbol), 'gene_id']
	gns$V3 = 'Gene Expression'

	# Barcodes
	barcs = read.table(paste0(folddir, '/', samps[i], '_count_barc.txt'), sep = '\t')
	barcs$V1 = paste0(barcs$V1, '-1')

	# Matrix
	mat = readMM(paste0(folddir, '/', samps[i], '_count.mtx'))

	# Remove genes with no stable gene id
	#whichKeep = which(!(is.na(gns$V1)))
	#gns = gns[whichKeep,]
	#mat = mat[whichKeep,]

	# Write to its own folder
	dir.create(paste0(outdir, '/', samps[i]))
	write.table(gns, paste0(outdir, '/', samps[i], '/features.tsv'), col.names = F, sep = '\t', quote = F, row.names = F)
	write.table(barcs, paste0(outdir, '/', samps[i], '/barcodes.tsv'), col.names = F, sep = '\t', quote = F, row.names = F)
	writeMM(mat, paste0(outdir, '/', samps[i], '/matrix.mtx'))

}
