rm(list = ls())
library(dplyr)
library(plyr)
library(Seurat)
library(ggpubr)
library(reshape2)
library(tidyr)
library(tidyverse)
library(DropletQC)
library(readr)
source('~/SCRIPTS/utility_functions.R')

# Load datasets
linseur = readRDS('~/workdir/PRIMATE_CAUDATE/MARMOSET_LIN/03_CELLBENDERED_ANNOTATE/Marmoset_Caudate_Lin_clustering2.RDS')

kriseur = readRDS('~/workdir/PRIMATE_CAUDATE/MARMOSET_KRIENEN/03_CELLBENDERED_ANNOTATE/Marmoset_Caudate_Krienen_clustering3.RDS')

# Merge datasets
seurM = merge(linseur, kriseur)

# Set variables
markerGenes = c('NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA')
pref = 'MARMOSET_KRIENEN_LIN_INTEGRATED'

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

allseur_integrated$type = ifelse(allseur_integrated$orig.ident %in% c('v321', 'v324'), 'LIN', 'KRIENEN')
pdf(paste0(pref, "_DATASET.pdf"))
DimPlot(allseur_integrated, group.by = 'type', raster = T)
dev.off()

allseur_integrated = FindNeighbors(allseur_integrated, dims = 1:30)
allseur_integrated = FindClusters(allseur_integrated, resolution = 0.1)

pdf(paste0(pref, "_CLUSTER.pdf"))
DimPlot(allseur_integrated, label = T, raster = T) + NoLegend()
dev.off()

stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_STACK_SAMPLES'))
stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'type', paste0(pref, '_STACK_DATASETS'))


# MARKER GENE PLOTS
DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)
pdf(paste0(pref, "_MarkersViolin.pdf"), width = 20, height = 15)
VlnPlot(allseur_integrated, features = markerGenes, pt.size = 0, raster = T)
dev.off()

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

saveRDS(allseur_integrated, '~/workdir/PRIMATE_CAUDATE/MARMOSET_INTEGRATION/Marmoset_Integrated_Clustering_1.RDS')

qw = FindMarkers(allseur_integrated, ident.1 = 11, ident.2 = 6, only.pos = T, logfc.threshold = 0.5)


####
## CLUSTERING 2
####

# Filtering
vasculature = c()
empty_drop = c()
doublet = c(11,12)
tmp = allseur_integrated
allseur_integrated = subset(allseur_integrated, subset = seurat_clusters %in% c(vasculature, empty_drop, doublet), invert = T)

# Clustering
allseur_integrated = RunPCA(allseur_integrated, verbose = FALSE, assay = 'integrated')
allseur_integrated = RunUMAP(allseur_integrated, dims = 1:30, assay = 'integrated')

# Basic Plots
pdf(paste0(pref, "_SAMPLE.pdf"))
DimPlot(allseur_integrated, group.by = 'orig.ident', raster = T)
dev.off()

allseur_integrated$type = ifelse(allseur_integrated$orig.ident %in% c('v321', 'v324'), 'LIN', 'KRIENEN')
pdf(paste0(pref, "_DATASET.pdf"))
DimPlot(allseur_integrated, group.by = 'type', raster = T)
dev.off()

DefaultAssay(allseur_integrated) = 'integrated'
allseur_integrated = FindNeighbors(allseur_integrated, dims = 1:30, assay = 'integrated')
allseur_integrated = FindClusters(allseur_integrated, resolution = 0.1)

pdf(paste0(pref, "_CLUSTER.pdf"))
DimPlot(allseur_integrated, label = T, raster = T) + NoLegend()
dev.off()

stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_STACK_SAMPLES'))
stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'type', paste0(pref, '_STACK_DATASETS'))


# MARKER GENE PLOTS
DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)
pdf(paste0(pref, "_MarkersViolin.pdf"), width = 20, height = 15)
VlnPlot(allseur_integrated, features = markerGenes, pt.size = 0, raster = T)
dev.off()

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

saveRDS(allseur_integrated, '~/workdir/PRIMATE_CAUDATE/MARMOSET_INTEGRATION/Marmoset_Integrated_Clustering_2.RDS')


####
## BROADLY ANNOTATE
####

allseur_integrated = readRDS('~/workdir/PRIMATE_CAUDATE/MARMOSET_INTEGRATION/Marmoset_Integrated_Clustering_2.RDS')

mapnames = setNames(c('MOL', 'SPN', 'Astrocyte', 'SPN', 'SPN', 'OPC',
			'Microglia', 'SPN', 'Other_GABA', 'SPN', 'Other_GABA',
			'Other_GABA', 'Other_GABA'),
		      c(0:12))

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
pdf(paste0(pref, "_MarkersViolin.pdf"), width = 20, height = 25)
VlnPlot(allseur_integrated, features = markerGenes, pt.size = 0, raster = T)
dev.off()

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
saveRDS(allseur_integrated, paste0('~/workdir/PRIMATE_CAUDATE/MARMOSET_INTEGRATION/', pref, '_FINAL_ANNOTATED.RDS'))
