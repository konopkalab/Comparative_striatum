library(dplyr)
library(plyr)
library(Seurat)
library(ggpubr)
library(reshape2)
library(tidyr)
library(tidyverse)
library(DropletQC)
source('~/SCRIPTS/utility_functions.R')

# Load datasets
ourseur = readRDS("/home2/gkonop/workdir/01_MACAQUE/OUR/02_CELLBENDERED_CLUSTER/CAUDATE/ANNOTATION/Macaque_Caudate_Our_Annotated_FINAL.RDS")
heseur = readRDS('/home2/gkonop/workdir/01_MACAQUE/HE/03_CLUSTER_ANNOTATE/CAUDATE_ONLY/Macaque_Caudate_He_Clustering3.RDS')

# Merge datasets
seurM = merge(ourseur, heseur)

# Set variables

markerGenes = c('SLC17A7', 'SATB2', 'RBFOX3', 'NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA')

pref = 'MACAQUE_OUR_HE_INTEGRATED'

####
## INTEGRATE - SEURAT
####

# Perform SCT
seurML = SplitObject(seurM, split.by = "orig.ident")
for (i in 1:length(seurML)) {
	seurML[[i]] = SCTransform(seurML[[i]], verbose = T, ncells = 1000)
	print(i)
}

# Feature selection
options(future.globals.maxSize= 99891289600)
allseur_feats = SelectIntegrationFeatures(object.list = seurML, nfeatures = 2000)
commonfeats = allseur_feats

# Integration
seurML = PrepSCTIntegration(object.list = seurML, anchor.features = commonfeats, verbose = FALSE)
allseur_anchors = FindIntegrationAnchors(object.list = seurML, normalization.method = "SCT", anchor.features = commonfeats, verbose = FALSE)
allseur_integrated = IntegrateData(anchorset = allseur_anchors, normalization.method = "SCT", verbose = FALSE)

# Clustering
allseur_integrated = RunPCA(allseur_integrated, verbose = FALSE)
allseur_integrated = RunUMAP(allseur_integrated, dims = 1:30)

# change the dir
setwd("/home2/gkonop/workdir/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/CAUDATE/CLUSTERING_1")
# Basic Plots
pdf(paste0(pref, "_GB_SAMPLE.pdf"))
DimPlot(allseur_integrated, group.by = 'orig.ident', raster = T)
dev.off()

allseur_integrated$type = ifelse(grepl('SRR', allseur_integrated$orig.ident), 'HE', 'OUR')
pdf(paste0(pref, "_GB_DATASET.pdf"))
DimPlot(allseur_integrated, group.by = 'type', raster = T)
dev.off()

allseur_integrated = FindNeighbors(allseur_integrated, dims = 1:30)
allseur_integrated = FindClusters(allseur_integrated, resolution = 1)

pdf(paste0(pref, "_GB_CLUSTER.pdf"))
DimPlot(allseur_integrated, label = T, raster = T) + NoLegend()
dev.off()

stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_GB_STACK_SAMPLES'))
stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'type', paste0(pref, '_GB_STACK_DATASETS'))


# MARKER GENE PLOTS
DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)
pdf(paste0(pref, "_GB_MarkersViolin.pdf"), width = 20, height = 15)
VlnPlot(allseur_integrated, features = markerGenes, pt.size = 0, raster = T)
dev.off()

pdf(paste0(pref, "_GB_MarkersDotPlot.pdf"), width = 20, height = 10)
DotPlot(allseur_integrated, features = markerGenes) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "_GB_NuclearFraction.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_GB_Depth_Gene.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

saveRDS(allseur_integrated, paste0('/home2/gkonop/workdir/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/CAUDATE/CLUSTERING_1/INTEGRATED_GB_Clustering1.RDS'))
allseur_integrated <- readRDS('/home2/gkonop/workdir/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/CAUDATE/CLUSTERING_1/INTEGRATED_GB_Clustering1.RDS')

####
## BROADLY ANNOTATE
####

mapnames = setNames(c('MOL', 'SPN', 'MOL', 'Astrocyte', 'SPN', 'MOL', 'SPN',
			 'Astrocyte','OPC', "MOL", 'Microglia', "SPN", "SPN", 'SPN', 'non-SPN', 'Astrocyte', 
			'non-SPN', 'Microglia','non-SPN', 'SPN', 'COP','non-SPN', 
			'non-SPN', 'non-SPN', 	'SPN'),
		      c(0:24))

allseur_integrated[["newannot"]] = mapnames[allseur_integrated[["seurat_clusters"]][,1]]
Idents(allseur_integrated) = allseur_integrated$newannot

# set the work directory
setwd("/home2/gkonop/workdir/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/CAUDATE/ANNOTATION")
# PLOTS
pdf(paste0(pref, "_GB_Clusters.pdf"))
DimPlot(allseur_integrated, group.by = 'newannot', label = T, raster = T, label.size = 7) +
NoLegend() +
ggtitle('')
dev.off()

pdf(paste0(pref, "_GB_Depth_Gene.pdf"))
ggboxplot(allseur_integrated[[]], x = 'newannot', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_GB_Depth_UMI.pdf"))
ggboxplot(allseur_integrated[[]], x = 'newannot', y = 'nCount_RNA', ylim = c(0,80000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_GB_NuclearFraction.pdf"))
ggboxplot(allseur_integrated[[]], x = 'newannot', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()


# MARKER GENE PLOTS
DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)
pdf(paste0(pref, "_GB_MarkersViolin.pdf"), width = 20, height = 25)
VlnPlot(allseur_integrated, features = markerGenes, pt.size = 0, raster = T)
dev.off()

pdf(paste0(pref, "_GB_MarkersDotPlot.pdf"), width = 20, height = 10)
DotPlot(allseur_integrated, features = markerGenes) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

stackedbarplot(allseur_integrated[[]], groupfill = 'newannot', groupx = 'orig.ident', fn = paste0(pref, '_GB_STACKED_ANNOTATED'))

# Save
saveRDS(allseur_integrated, paste0('~/workdir/01_MACAQUE/INTEGRATED/', pref, '_GB_Annotated.RDS'))

# Save only He et al
heseur = subset(allseur_integrated, subset = type == 'HE')
saveRDS(heseur, '~/workdir/01_MACAQUE/HE/03_CLUSTER_ANNOTATE/Macaque_Caudate_annotated_FINAL.RDS')



####
## Update matrix if the initial filtering was done incorrectly
####

allseur_orig = readRDS('~/workdir/01_MACAQUE/HE/03_CLUSTER_ANNOTATE/Macaque_Caudate_annotated_FINAL.RDS')
meta_orig = allseur_orig@meta.data

# Only the caudate samples
srrMeta = read.table('~/workdir/01_MACAQUE/HE/Macaque_He_SRA_Metadata.txt', header = T)
srrMeta$Age = as.numeric(gsub('.*,', '', srrMeta[,1]))
srrMeta$Sample =  srrMeta[, 1] %>% gsub(',.*', '', .)

samps = srrMeta[grepl('Caudate', srrMeta$Study.TISSUE), 'Sample']
allSeurL = list()
for(i in 1:length(samps)){

	#  Read CellBender output
	mat = Read10X_h5(paste0('~/project/03_MATRIX_FROM_CELLBENDER/MACAQUE_HE/', samps[i], '/CellBender_out_filtered.h5'))

	# Add sample name to the cell barcode name
	cbBarcode = colnames(mat)
	colnames(mat) = gsub('-1', '', colnames(mat)) %>% paste0(., '-', samps[i])

	# Create seurat object
	tmpSeur = CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0)
	tmpSeur$orig.ident = samps[i]
	tmpSeur$tissue = 'caudate'
	tmpSeur$Species = 'macaque'
	allSeurL[[i]] = tmpSeur

	print(i)
}

allseur = Reduce(merge, allSeurL)
allseur_filt = subset(allseur, cells = rownames(meta_orig))

mat = allseur_filt@assays$RNA@counts
newSeur = CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0, meta.data = meta_orig)

saveRDS(newSeur, '~/workdir/01_MACAQUE/HE/03_CLUSTER_ANNOTATE/Macaque_Caudate_annotated_FINAL_ALLGENES.RDS')
