rm(list = ls())
library(dplyr)
library(plyr)
library(Seurat)
library(ggpubr)
library(reshape2)
library(tidyr)
library(tidyverse)
library(DropletQC)
source('~/SCRIPTS/utility_functions.R')

#####
## LOAD COUNT MATRIX
#####

# Set the samples and read them with QC
samplePath = '/home2/gkonop/project/03_MATRIX_FROM_CELLBENDER/MARMOSET_KRIENEN/LIFTOFF'
samps = dir(samplePath)
samps = samps[grepl('marm', samps)]

allSeurL = list()
for(i in 1:length(samps)){

	#  Read CellBender output
	mat = Read10X_h5(paste0(samplePath, '/', samps[i], '/CellBender_out_filtered.h5'))

	# Soft filtering. Remove lower than 200 UMI
	mat = mat[, colSums(mat) > 200]
	cbBarcode = colnames(mat)
	colnames(mat) = gsub('-1', '', colnames(mat)) %>% paste0(., '_', samps[i])

	# Create seurat object
	tmpSeur = CreateSeuratObject(counts = mat, min.cells = 10, min.features = 1)

	# Calculate intronic read ratio. Skip this since we did not keep the BAMs after assigning to liftoffed genes
	nf2 = nuclear_fraction_annotation(
		annotation_path = '/home2/gkonop/workdir/REFERENCES/GTF_ONLY/hg38_liftoverTo_calJac3.gtf',
		bam = paste0('~/project/01_CELLRANGER_OUT/MARMOSET_KRIENEN_FASTQ/celllranger_count_', samps[i], '/outs/possorted_genome_bam.bam'),
		barcodes = cbBarcode,
		cores = 23, verbose = T)

	nf2$barc = colnames(tmpSeur)
	tmpSeur$intronRat = nf2$nuclear_fraction
	tmpSeur$orig.ident = samps[i]
	tmpSeur$Species = 'marmoset'
	allSeurL[[i]] = tmpSeur

	print(i)
}

allseur = Reduce(merge, allSeurL)
saveRDS(allseur, '~/workdir/PRIMATE_CAUDATE/MARMOSET_INTEGRATION/KRIENEN_SEURAT_RAW.RDS')
allseur = readRDS('~/workdir/PRIMATE_CAUDATE/MARMOSET_INTEGRATION/KRIENEN_SEURAT_RAW.RDS')

# Add meta data
meta = rio::import_list('~/project/00_BAM_DOWNLOADED/KRIENEN_2020/MARMOSET/krienen_marmoset_meta.xlsx')
meta = meta[[2]]
allseur$Individual = gsub('_RXN.*', '', allseur$orig.ident)

seurmeta = allseur@meta.data
seurmeta2 = cbind(seurmeta, meta[match(seurmeta$Individual, meta$BICCN_id),])
allseur@meta.data = seurmeta2

####
## CLUSTERING 1
####

pref = 'Marmoset_Caudate_Krienen'
seurM = allseur
seurM = subset(seurM, subset = nCount_RNA > 0)
seurM_corr = clusterBatchCorrect(seurM, npcs = 30, res = 0.1, pref = pref, sampleName = 'orig.ident')

# QC PLOTS

pdf(paste0(pref, "_Clusters.pdf"))
DimPlot(seurM_corr, group.by = 'seurat_clusters', label = T, raster = T, label.size = 7) +
NoLegend() +
ggtitle('')
dev.off()

pdf(paste0(pref, "_Depth_Gene.pdf"))
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_Depth_UMI.pdf"))
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'nCount_RNA', ylim = c(0,10000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_NuclearFraction.pdf"))
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# MARKER GENE PLOTS
markerGenes = c('NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA')

DefaultAssay(seurM_corr) = 'RNA'
seurM_corr = NormalizeData(seurM_corr)
pdf(paste0(pref, "_MarkersViolin.pdf"), width = 20, height = 15)
VlnPlot(seurM_corr, features = markerGenes, pt.size = 0, raster = T)
dev.off()

pdf(paste0(pref, "_MarkersDotPlot.pdf"), width = 20, height = 10)
DotPlot(seurM_corr, features = markerGenes) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

# Save (optional)
saveRDS(seurM_corr, paste0('~/workdir/PRIMATE_CAUDATE/MARMOSET_KRIENEN/03_CELLBENDERED_ANNOTATE/', pref, '_clustering1.RDS'))


####
## CLUSTERING 2
####

pref = 'Marmoset_Caudate_Krienen'
vasculature = c(4,8)
empty_drop = c(10)
doublet = c(11,12)
seurM = subset(seurM_corr, subset = seurat_clusters %in% c(vasculature, empty_drop, doublet), invert = T)
seurM_corr = clusterBatchCorrect(seurM, npcs = 30, res = 0.1, pref = pref, sampleName = 'orig.ident')

# QC PLOTS
pdf(paste0(pref, "_Clusters.pdf"))
DimPlot(seurM_corr, group.by = 'seurat_clusters', label = T, raster = T, label.size = 7) +
NoLegend() +
ggtitle('')
dev.off()

pdf(paste0(pref, "_Depth_Gene.pdf"))
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_Depth_UMI.pdf"))
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'nCount_RNA', ylim = c(0,10000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_NuclearFraction.pdf"))
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# MARKER GENE PLOTS
markerGenes = c('NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA')

DefaultAssay(seurM_corr) = 'RNA'
seurM_corr = NormalizeData(seurM_corr)
pdf(paste0(pref, "_MarkersViolin.pdf"), width = 20, height = 15)
VlnPlot(seurM_corr, features = markerGenes, pt.size = 0, raster = T)
dev.off()

pdf(paste0(pref, "_MarkersDotPlot.pdf"), width = 20, height = 10)
DotPlot(seurM_corr, features = markerGenes) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

# Save (optional)
saveRDS(seurM_corr, paste0('~/workdir/PRIMATE_CAUDATE/MARMOSET_KRIENEN/03_CELLBENDERED_ANNOTATE/', pref, '_clustering2.RDS'))


####
## CLUSTERING 3
####

pref = 'Marmoset_Caudate_Krienen'
vasculature = c()
empty_drop = c()
doublet = c(9)
seurM = subset(seurM_corr, subset = seurat_clusters %in% c(vasculature, empty_drop, doublet), invert = T)
seurM_corr = clusterBatchCorrect(seurM, npcs = 30, res = 0.1, pref = pref, sampleName = 'orig.ident')

# QC PLOTS
pdf(paste0(pref, "_Clusters.pdf"))
DimPlot(seurM_corr, group.by = 'seurat_clusters', label = T, raster = T, label.size = 7) +
NoLegend() +
ggtitle('')
dev.off()

pdf(paste0(pref, "_Depth_Gene.pdf"))
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_Depth_UMI.pdf"))
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'nCount_RNA', ylim = c(0,10000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_NuclearFraction.pdf"))
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# MARKER GENE PLOTS
markerGenes = c('NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA')

DefaultAssay(seurM_corr) = 'RNA'
seurM_corr = NormalizeData(seurM_corr)
pdf(paste0(pref, "_MarkersViolin.pdf"), width = 20, height = 15)
VlnPlot(seurM_corr, features = markerGenes, pt.size = 0, raster = T)
dev.off()

pdf(paste0(pref, "_MarkersDotPlot.pdf"), width = 20, height = 10)
DotPlot(seurM_corr, features = markerGenes) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

# Save (optional)
saveRDS(seurM_corr, paste0('~/workdir/PRIMATE_CAUDATE/MARMOSET_KRIENEN/03_CELLBENDERED_ANNOTATE/', pref, '_clustering3.RDS'))

### Putamen
library(dplyr)
library(plyr)
library(Seurat)
library(ggpubr)
library(reshape2)
library(tidyr)
library(tidyverse)
library(DropletQC)
source('~/SCRIPTS/utility_functions.R')

#####
## LOAD COUNT MATRIX
#####

# Set the samples and read them with QC
samplePath = '/home2/gkonop/project/03_MATRIX_FROM_CELLBENDER/MARMOSET_KRIENEN/LIFTOFF/PUTAMEN'
samps = dir(samplePath)
samps = samps[grepl('marm', samps)]

path_cr = '~/project/01_CELLRANGER_OUT/MARMOSET_KRIENEN_FASTQ/PUTAMEN/'

allSeurL = list()
for(i in 1:length(samps)){

	#  Read CellBender output
	mat = Read10X_h5(paste0(samplePath, '/', samps[i], '/CellBender_out_filtered.h5'))

	# Soft filtering. Remove lower than 200 UMI
	mat = mat[, colSums(mat) > 200]
	cbBarcode = colnames(mat)
	colnames(mat) = gsub('-1', '', colnames(mat)) %>% paste0(., '_', samps[i])

	# Create seurat object
	tmpSeur = CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0)

	# Load intronic read ratio
	fls = list.files( paste0(path_cr, 'cellranger_count_', samps[i], '/outs'), pattern = 'Intronic', full.names=T)

	if(any(grepl('RDS', fls)) == T){

		nf1 = readRDS(fls[grepl('RDS', fls)])
		print('Intronic_read_ratio_loaded')
		
	} else{

		# Calculate intronic read ratio
		nf1 = nuclear_fraction_tags(
			bam = paste0(path_cr, 'cellranger_count_', samps[i], '/outs/possorted_genome_bam.bam'),
			barcodes = cbBarcode,
			cores = 23, verbose = T)

		saveRDS(nf1, paste0(path_cr, 'cellranger_count_', samps[i], '/outs/Intronic_read_ratio.RDS'))
	}

	# Keep the final barcodes
	nf1 = nf1[match(cbBarcode, rownames(nf1)),]

	tmpSeur$intronRat = nf1
	tmpSeur$orig.ident = samps[i]
	tmpSeur$tissue = 'Putamen'
	tmpSeur$Species = 'Marmoset'
	allSeurL[[i]] = tmpSeur

	print(i)
}

allseur = Reduce(merge, allSeurL)
write_rds(allseur, '~/workdir/01_MARMOSET/MARMOSET_KRIENEN/03_CELLBENDERED_ANNOTATE/PUTAMEN/KRIENEN_PUTAMEN_RAW.RDS')
allseur = read_rds('~/workdir/01_MARMOSET/MARMOSET_KRIENEN/03_CELLBENDERED_ANNOTATE/PUTAMEN/KRIENEN_PUTAMEN_RAW.RDS')

# Add meta data
meta = rio::import_list('~/project/00_BAM_DOWNLOADED/KRIENEN_2020/MARMOSET/krienen_marmoset_meta.xlsx')
meta = meta[[2]]
allseur$Individual = gsub('_RXN.*', '', allseur$orig.ident)

seurmeta = allseur@meta.data
seurmeta2 = cbind(seurmeta, meta[match(seurmeta$Individual, meta$BICCN_id),])
allseur@meta.data = seurmeta2

####
## SET VARIABLES
####

pref = 'Marmoset_Putamen_Krienen'

markerGenes = c('SLC17A7', 'SATB2', 'RBFOX3','NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'OTOF', 'CACNG5', 'PCDH8', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA','IL29')

####
## CLUSTERING 1
####

seurM = allseur
seurM = subset(seurM, subset = nCount_RNA > 0)
seurM_corr = clusterBatchCorrect(seurM, npcs = 30, res = 1, pref = pref, sampleName = 'orig.ident')

# QC PLOTS
# change working directory
setwd("/home2/gkonop/workdir/01_MARMOSET/MARMOSET_KRIENEN/03_CELLBENDERED_ANNOTATE/PUTAMEN/CLUSTERING1")

# plots
pdf(paste0(pref, "_GB_Clusters.pdf"))
DimPlot(seurM_corr, group.by = 'seurat_clusters', label = T, raster = T, label.size = 7) +
NoLegend() +
ggtitle('')
dev.off()

pdf(paste0(pref, "_GB_Depth_Gene.pdf"))
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_GB_Depth_UMI.pdf"))
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'nCount_RNA', ylim = c(0,10000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_GB_NuclearFraction.pdf"))
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# MARKER GENE PLOTS
DefaultAssay(seurM_corr) = 'RNA'
seurM_corr = NormalizeData(seurM_corr)

pdf(paste0(pref, "_GB_MarkersDotPlot.pdf"), width = 20, height = 10)
DotPlot(seurM_corr, features = markerGenes) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

qw = FindMarkers(seurM_corr, ident.1 = 20, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
DotPlot(seurM_corr, features = rownames(qw)[1:30]) + rotate_x_text(90)

# Save (optional)
write_rds(seurM_corr, '~/workdir/01_MARMOSET/MARMOSET_KRIENEN/03_CELLBENDERED_ANNOTATE/PUTAMEN/GB_KRIENEN_PUTAMEN_CLUSTERING1.RDS')

####
## CLUSTERING 2
####

## Cluster_ID,CellType
#0,MOL
#1,Astrocyte
#2,iSPN
#3,dSPN
#4,MOL
#5,dSPN
#6*,Empty
#7,iSPN
#8,dSPN
#9,iSPN
#10,iSPN
#11,OPC
#12,Astrocyte
#13*,Endo
#14,Microglia
#15,dSPN
#16,MOL
#17,eSPN
#18*,Ast/Endo
#19*,non-SPN/MOL
#20*,dsPN/Ast
#21*,dSPN/iSPN
#22*,non-SPN/MOL
#23*,MOL/Ast
#24,non-SPN
#25,non-SPN
#26,non-SPN
#27*,MOL/OPC
#28*,Peri
#29*,Ast/non-SPN?
#30,non-SPN
#31,iSPN
#32,non-SPN
#33*,SPN/non-SPN


vasculature = c(13,28)
empty_drop = c(6)
doublet = c(18,19,20,21,22,23,27,29,33)
seurM = subset(seurM_corr, subset = seurat_clusters %in% c(vasculature, empty_drop, doublet), invert = T)
seurM_corr = clusterBatchCorrect(seurM, npcs = 30, res = 1, pref = pref, sampleName = 'orig.ident')

# set working dir
setwd("/home2/gkonop/workdir/01_MARMOSET/MARMOSET_KRIENEN/03_CELLBENDERED_ANNOTATE/PUTAMEN/CLUSTERING2")

# QC PLOTS
pdf(paste0(pref, "_GB_Clusters.pdf"))
DimPlot(seurM_corr, group.by = 'seurat_clusters', label = T, raster = T, label.size = 7) +
NoLegend() +
ggtitle('')
dev.off()

pdf(paste0(pref, "_GB_Depth_Gene.pdf"))
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_GB_Depth_UMI.pdf"))
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'nCount_RNA', ylim = c(0,10000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_GB_NuclearFraction.pdf"))
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# MARKER GENE PLOTS
DefaultAssay(seurM_corr) = 'RNA'
seurM_corr = NormalizeData(seurM_corr)

pdf(paste0(pref, "_GB_MarkersDotPlot.pdf"), width = 20, height = 10)
DotPlot(seurM_corr, features = markerGenes) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

# Save (optional)
write_rds(seurM_corr, '~/workdir/01_MARMOSET/MARMOSET_KRIENEN/03_CELLBENDERED_ANNOTATE/PUTAMEN/KRIENEN_PUTAMEN_CLUSTERING2.RDS')


qw = FindMarkers(seurM_corr, ident.1 = 11, ident.2 = 17, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
DotPlot(seurM_corr, features = rownames(qw)[1:30]) + rotate_x_text(90)

####
## BROADLY ANNOTATE
####

mapnames = setNames(c('MOL', 'Astrocyte', 'SPN', 'SPN', 'SPN', 'MOL',
			'SPN', 'SPN', 'SPN', 'OPC', 'SPN',
			'Microglia', 'Astrocyte', 'MOL', 'MOL', 'SPN',
			'Non_SPN', 'Non_SPN', 'Non_SPN', 'Astrocyte', 'Non_SPN', 'SPN', 'Non_SPN', 'SPN'),
		      c(0:23))

seurM_corr[["newannot"]] = mapnames[seurM_corr[["seurat_clusters"]][,1]]
Idents(seurM_corr) = seurM_corr$newannot

# set working dir
setwd("/home2/gkonop/workdir/01_MARMOSET/MARMOSET_KRIENEN/03_CELLBENDERED_ANNOTATE/PUTAMEN/ANNOTATED")

# PLOTS
pdf(paste0(pref, "_GB_Clusters.pdf"))
DimPlot(seurM_corr, group.by = 'newannot', label = T, raster = T, label.size = 7) +
NoLegend() +
ggtitle('')
dev.off()

pdf(paste0(pref, "_GB_Depth_Gene.pdf"))
ggboxplot(seurM_corr[[]], x = 'newannot', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_GB_Depth_UMI.pdf"))
ggboxplot(seurM_corr[[]], x = 'newannot', y = 'nCount_RNA', ylim = c(0,80000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_GB_NuclearFraction.pdf"))
ggboxplot(seurM_corr[[]], x = 'newannot', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()


# MARKER GENE PLOTS
DefaultAssay(seurM_corr) = 'RNA'
seurM_corr = NormalizeData(seurM_corr)
pdf(paste0(pref, "_GB_MarkersDotPlot.pdf"), width = 20, height = 10)
DotPlot(seurM_corr, features = markerGenes) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

# Save
stackedbarplot(seurM_corr[[]], groupfill = 'newannot', groupx = 'orig.ident', fn = 'Marmoset_Stacked_Putamen_Annotated')
write_rds(seurM_corr, '/home2/gkonop/workdir/01_MARMOSET/MARMOSET_KRIENEN/03_CELLBENDERED_ANNOTATE/PUTAMEN/KRIENEN_PUTAMEN_ANNOTATED.RDS')

