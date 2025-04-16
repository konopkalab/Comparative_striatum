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

# Load datasets
ourseur = read_rds('~/workdir/01_MACAQUE/OUR/02_CELLBENDERED_CLUSTER/Macaque_Putamen_Our_Clustering4.RDS')
heseur = read_rds('~/workdir/01_MACAQUE/HE/03_CLUSTER_ANNOTATE/PUTAMEN_ONLY/Macaque_Putamen_He_Clustering5.RDS')

# Merge datasets
seurM = merge(ourseur, heseur)

# Set variables
markerGenes = c('NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA', 'SATB2', 'SLC17A7')
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

# Basic Plots
pdf(paste0(pref, "_SAMPLE.pdf"))
DimPlot(allseur_integrated, group.by = 'orig.ident', raster = T)
dev.off()

allseur_integrated$type = ifelse(grepl('SRR', allseur_integrated$orig.ident), 'HE', 'OUR')
pdf(paste0(pref, "_DATASET.pdf"))
DimPlot(allseur_integrated, group.by = 'type', raster = T)
dev.off()

allseur_integrated = FindNeighbors(allseur_integrated, dims = 1:30)
allseur_integrated = FindClusters(allseur_integrated, resolution = 0.5)

pdf(paste0(pref, "_CLUSTER.pdf"))
DimPlot(allseur_integrated, label = T, raster = T) + NoLegend()
dev.off()

stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_STACK_SAMPLES'))
stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'type', paste0(pref, '_STACK_DATASETS'))


# MARKER GENE PLOTS
DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)

pdf(paste0(pref, "_MarkersDotPlot.pdf"), width = 20, height = 10)
DotPlot(allseur_integrated, features = markerGenes) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "_NuclearFraction.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_Depth_Gene.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

write_rds(allseur_integrated, paste0('~/workdir/01_MACAQUE/INTEGRATED/PUTAMEN/', pref, '_Clustering1_Putamen.RDS'))

####
## CLUSTERING 2
####

allseur_integrated = subset(allseur_integrated, subset = seurat_clusters == 17, invert = T)
DefaultAssay(allseur_integrated) = 'integrated'

# Clustering
allseur_integrated = RunPCA(allseur_integrated, verbose = FALSE)
allseur_integrated = RunUMAP(allseur_integrated, dims = 1:30)

# Basic Plots
pdf(paste0(pref, "_SAMPLE.pdf"))
DimPlot(allseur_integrated, group.by = 'orig.ident', raster = T)
dev.off()

allseur_integrated$type = ifelse(grepl('SRR', allseur_integrated$orig.ident), 'HE', 'OUR')
pdf(paste0(pref, "_DATASET.pdf"))
DimPlot(allseur_integrated, group.by = 'type', raster = T)
dev.off()

allseur_integrated = FindNeighbors(allseur_integrated, dims = 1:30)
allseur_integrated = FindClusters(allseur_integrated, resolution = 0.5)

pdf(paste0(pref, "_CLUSTER.pdf"))
DimPlot(allseur_integrated, label = T, raster = T) + NoLegend()
dev.off()

stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_STACK_SAMPLES'))
stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'type', paste0(pref, '_STACK_DATASETS'))


# MARKER GENE PLOTS
DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)

pdf(paste0(pref, "_MarkersDotPlot.pdf"), width = 20, height = 10)
DotPlot(allseur_integrated, features = markerGenes) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "_NuclearFraction.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_Depth_Gene.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

# Remove cluster 17 (oligo-microglia doublet)
allseur_integrated = subset(allseur_integrated, subset = seurat_clusters == 17, invert = T)
write_rds(allseur_integrated, paste0('~/workdir/01_MACAQUE/INTEGRATED/PUTAMEN/', pref, '_Clustering2_Putamen.RDS'))

DefaultAssay(allseur_integrated) = 'integrated'
allseur_integrated <- BuildClusterTree(object = allseur_integrated)
PlotClusterTree(object = allseur_integrated)

####
## BROADLY ANNOTATE
####

mapnames = setNames(c('MOL', 'Astrocyte', 'SPN', 'MOL', 'SPN', 'OPC',
			'SPN', 'Microglia', 'SPN', 'Non_SPN', 'SPN',
			'Non_SPN', 'Non_SPN', 'Astrocyte', 'Non_SPN', 'Non_SPN',
			'OPC'),
		      c(0:16))

allseur_integrated[["newannot"]] = mapnames[allseur_integrated[["seurat_clusters"]][,1]]
Idents(allseur_integrated) = allseur_integrated$newannot

# PLOTS
pdf(paste0(pref, "_Clusters.pdf"))
DimPlot(allseur_integrated, group.by = 'newannot', label = T, raster = T, label.size = 7) +
NoLegend() +
ggtitle('')
dev.off()

pdf(paste0(pref, "_Depth_Gene.pdf"))
ggboxplot(allseur_integrated[[]], x = 'newannot', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_Depth_UMI.pdf"))
ggboxplot(allseur_integrated[[]], x = 'newannot', y = 'nCount_RNA', ylim = c(0,80000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_NuclearFraction.pdf"))
ggboxplot(allseur_integrated[[]], x = 'newannot', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()


# MARKER GENE PLOTS
DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)

pdf(paste0(pref, "_MarkersDotPlot.pdf"), width = 20, height = 10)
DotPlot(allseur_integrated, features = markerGenes) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

stackedbarplot(allseur_integrated[[]], groupfill = 'newannot', groupx = 'orig.ident', fn = paste0(pref, '_STACKED_ANNOTATED'))

# Save
write_rds(allseur_integrated, paste0('~/workdir/01_MACAQUE/INTEGRATED/PUTAMEN/', pref, '_Annotated_Putamen.RDS'))
