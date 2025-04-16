rm(list = ls())
library(Seurat)
library(Matrix)
library(ggplot2)
library(plyr)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggpubr)
library(reshape2)
library(rio)
library(data.table)
library(harmony)
library(GeneOverlap)
library(readr)
source('~/SCRIPTS/utility_functions.R')
set.seed(1234)

# Subset non-spns
seurM = read_rds('~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/Human_CaudatePutamen_NEURONS_Annotated.RDS')
to_exclude = c('SPN_IEG', 'dSPN', 'iSPN', 'eSPN_DRD2', 'eSPN_FOXP2', 'Excitatory', 'SPN')
seurM = subset(seurM, subset = newannot %in% to_exclude, invert = T)

####
## SET VARIABLES FOR CLUSTERING
####

pref = 'Human_CaudatePutamen_NONSPN'

markerGenes = c('NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CFTR', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA', 'SATB2', 'SLC17A7', 'EYA2', 'PDGFD', 'PTHLH', 'PVALB')

####
## CLUSTER-1
####

seurM_corr = clusterBatchCorrect(seurM, npcs = 10, res = 0.5, pref = pref, sampleName = 'orig.ident')
stackedbarplot(seurM_corr[[]], 'seurat_clusters', 'Tissue', paste0(pref, '_Tissue'), wd = 10, hg = 10)

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

# Move the plots
outfls = list.files('~/', pattern = pref, full.names = F)
todir = '~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN_NONSPN/CLUSTERING_1/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)

# Save (optional)
write_rds(seurM_corr, '~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/Human_CaudatePutamen_NONSPN_Clustering1.RDS')
seurM_corr = read_rds('~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/Human_CaudatePutamen_NONSPN_Clustering1.RDS')

qw = FindMarkers(seurM_corr, ident.1 = 4, only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
pdf('Cluster_4_Markers.pdf', width = 12)
DotPlot(seurM_corr, features = rownames(qw)[1:30]) + rotate_x_text(45)
dev.off()

####
## CLUSTER-2
####

# Remove doublets or spns and recluster
toremove = c(2,10,12,13,14)

seurM = subset(seurM_corr, subset = seurat_clusters %in% toremove, invert = T)
seurM_corr = clusterBatchCorrect(seurM, npcs = 10, res = 0.5, pref = pref, sampleName = 'orig.ident')
stackedbarplot(seurM_corr[[]], 'seurat_clusters', 'Tissue', paste0(pref, '_Tissue'), wd = 10, hg = 10)

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

# Move the plots
outfls = list.files('~/', pattern = pref, full.names = F)
todir = '~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN_NONSPN/CLUSTERING_2/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)

# Save (optional)
write_rds(seurM_corr, '~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/Human_CaudatePutamen_NONSPN_Clustering2.RDS')
seurM_corr = read_rds('~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/Human_CaudatePutamen_NONSPN_Clustering2.RDS')

qw = FindMarkers(seurM_corr, ident.1 = 4, only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
pdf('Cluster_4_Markers.pdf', width = 12)
DotPlot(seurM_corr, features = rownames(qw)[1:30]) + rotate_x_text(45)
dev.off()


####
## ANNOTATE -- TEMPORARY
####

mapnames = setNames(c('TAC3', 'PVALB_PTHLH', 'UNK_FOXP2', 'CCK_VIP+', 'UNK_SPARCL1', 'SST_NPY',
			'PVALB_PTHLH', 'TAC3', 'CCK_VIP-', 'CHAT', 'PVALB_PTHLH',
			'PVALB_PTHLH', 'CCK_VIP-'),
		      c(0:12))

seurM_corr[["newannot"]] = mapnames[seurM_corr[["seurat_clusters"]][,1]]
Idents(seurM_corr) = seurM_corr$newannot

# PLOTS
pdf(paste0(pref, "_Clusters_Annotated.pdf"))
DimPlot(seurM_corr, group.by = 'newannot', label = T, raster = T, repel = T, label.size = 7) +
NoLegend() +
ggtitle('')
dev.off()

pdf(paste0(pref, "_Depth_Gene_Annotated.pdf"))
ggboxplot(seurM_corr[[]], x = 'newannot', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_Depth_UMI_Annotated.pdf"))
ggboxplot(seurM_corr[[]], x = 'newannot', y = 'nCount_RNA', ylim = c(0,80000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_NuclearFraction_Annotated.pdf"))
ggboxplot(seurM_corr[[]], x = 'newannot', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()


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

stackedbarplot(seurM_corr[[]], groupfill = 'newannot', groupx = 'orig.ident', fn = paste0(pref, '_Stacked_Annotated'))

stackedbarplot(seurM_corr[[]], groupfill = 'orig.ident', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_2'))

stackedbarplot(seurM_corr[[]], groupfill = 'Tissue', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Tissue'))

# Move the plots
outfls = list.files('~/', pattern = pref, full.names = F)
todir = '~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN_NONSPN/ANNOTATED/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)

# Save (optional)
write_rds(seurM_corr, '~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/Human_CaudatePutamen_NONSPN_Annotated.RDS')
seurM_corr = read_rds('~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/Human_CaudatePutamen_NONSPN_Annotated.RDS')
