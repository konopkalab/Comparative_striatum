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

# Only the putamen samples
allsamps = rio::import('~/workdir/TABLES/Primate_Tables_From_EmilyOh/Chimp Macaque Language Summary.xlsx')
colnames(allsamps) = gsub(' ', '_', colnames(allsamps))
allsamps$Species = gsub(' ', '_', allsamps$Species)
putamenSamps = allsamps[allsamps$Brain_Region == 'Putamen' & allsamps$Species == 'Macaque',]
rownames(putamenSamps) = paste0('Sample_', putamenSamps[,1])

# Set the samples and read them with QC
fls = list.files('~/project/03_MATRIX_FROM_CELLBENDER/MACAQUE_OUR', pattern = "filtered.h5", recursive=T, full.names = T)
tmp = gsub('.*Macaque_|\\/CellBender.*', '', fls)
fls = fls[which(tmp %in% rownames(putamenSamps))]

allSeurL = list()
for(i in 1:length(fls)){

	# Sample name
	samp = gsub('.*Macaque_', '', fls[i]) %>% gsub('\\/.*', '', .)

	#  Read CellBender output
	mat = Read10X_h5(fls[i])

	# Soft filtering. Remove lower than 200 UMI
	mat = mat[, colSums(mat) > 200]
	cbBarcode = colnames(mat)
	colnames(mat) = gsub('-1', '', colnames(mat)) %>% paste0(., '-', samp)

	# Create seurat object
	tmpSeur = CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0)

	# Calculate intronic read ratio
	nf2 = nuclear_fraction_annotation(
		annotation_path = '/home2/gkonop/workdir/REFERENCES/GTF_ONLY/hg38_liftoverTo_rheMac10.gtf',
		bam = paste0('~/project/01_CELLRANGER_OUT/MACAQUE_OUR/celllranger_count_', samp, '/outs/possorted_genome_bam.bam'),
		barcodes = cbBarcode,
		cores = 23, verbose = T)

	nf2$barc = colnames(tmpSeur)
	tmpSeur$intronRat = nf2$nuclear_fraction
	tmpSeur$orig.ident = samp
	tmpSeur$tissue = 'putamen'
	tmpSeur$Species = 'macaque'
	allSeurL[[i]] = tmpSeur

	print(i)
}

allseur = Reduce(merge, allSeurL)

write_rds(allseur, '~/workdir/01_MACAQUE/OUR/02_CELLBENDERED_CLUSTER/Macaque_Putamen_Seurat_CB_Raw.RDS')
allseur = read_rds('~/workdir/MACAQUE/OUR/02_CELLBENDERED_CLUSTER/Macaque_Putamen_Seurat_CB_Raw.RDS')

# Add meta data
seurMeta = allseur@meta.data
meta = rio::import('~/workdir/TABLES/Primate_Tables_From_EmilyOh/Chimp Macaque Language Summary.xlsx')
rownames(meta) = paste0('Sample_', meta[,1])
colnames(meta) = gsub(' ', '_', colnames(meta)) %>% gsub('\\(|\\)', '', .)
seurMeta = cbind(seurMeta, meta[match(seurMeta$orig.ident, rownames(meta)), c('Brain_Region', 'Age', 'Sex')])
allseur@meta.data = seurMeta

####
## SET VARIABLES
####

pref = 'Macaque_Putamen_Our'

markerGenes = c('NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA', 'SLC17A7', 'SATB2', 'IL19', 'DNAH9', 'DNAH12')


####
## CLUSTERING 1
####

seurM = allseur
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
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'nCount_RNA', ylim = c(0,50000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_NuclearFraction.pdf"))
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()


# MARKER GENE PLOTS
DefaultAssay(seurM_corr) = 'RNA'
seurM_corr = NormalizeData(seurM_corr)

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
dir.create('~/workdir/01_MACAQUE/OUR/02_CLUSTER_ANNOTATE/PUTAMEN_ONLY/CLUSTERING_1')
write_rds(seurM_corr, paste0('~/workdir/01_MACAQUE/OUR/02_CELLBENDERED_CLUSTER/', pref, '_Clustering1.RDS'))
seurM_corr = read_rds(paste0('~/workdir/01_MACAQUE/OUR/02_CELLBENDERED_CLUSTER/', pref, '_Clustering1.RDS'))

qw = FindMarkers(seurM_corr, ident.1 = 13, only.pos = T, logfc.threshold = 0.5, min.pct =0.5)
DotPlot(seurM_corr, features = rownames(qw)[1:30]) + rotate_x_text(90)


####
## CLUSTERING 2
####

vasculature = c(7,13)
empty_drop = c(14)
doublet = c(3)
seurM = subset(seurM_corr, subset = seurat_clusters %in% c(vasculature, empty_drop, doublet), invert = T)
seurM_corr = clusterBatchCorrect(seurM, npcs = 30, res = 0.5, pref = pref, sampleName = 'orig.ident')

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
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'nCount_RNA', ylim = c(0,50000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_NuclearFraction.pdf"))
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()


# MARKER GENE PLOTS
DefaultAssay(seurM_corr) = 'RNA'
seurM_corr = NormalizeData(seurM_corr)
pdf(paste0(pref, "_MarkersDotPlot.pdf"), width = 25, height = 10)
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
dir.create('~/workdir/01_MACAQUE/OUR/02_CELLBENDERED_CLUSTER/PUTAMEN/CLUSTERING_2')
write_rds(seurM_corr, paste0('~/workdir/01_MACAQUE/OUR/02_CELLBENDERED_CLUSTER/', pref, '_Clustering2.RDS'))

qw = FindMarkers(seurM_corr, ident.1 = 3, only.pos = T, logfc.threshold = 0.5, min.pct =0.5)

####
## CLUSTERING 3
####

vasculature = c(19,20)
empty_drop = c()
doublet = c(5,10,12,13,14,17,24)
seurM = subset(seurM_corr, subset = seurat_clusters %in% c(vasculature, empty_drop, doublet), invert = T)
seurM_corr = clusterBatchCorrect(seurM, npcs = 30, res = 0.5, pref = pref, sampleName = 'orig.ident')

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
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'nCount_RNA', ylim = c(0,50000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_NuclearFraction.pdf"))
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()


# MARKER GENE PLOTS
DefaultAssay(seurM_corr) = 'RNA'
seurM_corr = NormalizeData(seurM_corr)
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
dir.create('~/workdir/01_MACAQUE/OUR/02_CELLBENDERED_CLUSTER/PUTAMEN/CLUSTERING_3')
write_rds(seurM_corr, paste0('~/workdir/01_MACAQUE/OUR/02_CELLBENDERED_CLUSTER/', pref, '_Clustering3.RDS'))

#qw = FindMarkers(seurM_corr, ident.1 = 12, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)

####
## CLUSTERING 4
####

vasculature = c()
empty_drop = c(4,10)
doublet = c()
neurogenesis = c()
seurM = subset(seurM_corr, subset = seurat_clusters %in% c(vasculature, empty_drop, doublet, neurogenesis), invert = T)
seurM_corr = clusterBatchCorrect(seurM, npcs = 30, res = 0.5, pref = pref, sampleName = 'orig.ident')

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
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'nCount_RNA', ylim = c(0,50000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_NuclearFraction.pdf"))
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()


# MARKER GENE PLOTS
DefaultAssay(seurM_corr) = 'RNA'
seurM_corr = NormalizeData(seurM_corr)
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
dir.create('~/workdir/01_MACAQUE/OUR/02_CELLBENDERED_CLUSTER/PUTAMEN/CLUSTERING_4')
write_rds(seurM_corr, paste0('~/workdir/01_MACAQUE/OUR/02_CELLBENDERED_CLUSTER/', pref, '_Clustering4.RDS'))
