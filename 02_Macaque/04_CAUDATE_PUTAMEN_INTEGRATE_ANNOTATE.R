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
library(DropletQC)
library(xlsx)
library(biomaRt)
source('~/SCRIPTS/utility_functions.R')
set.seed(1234)

### Load macaque data
macaque_caud <- readRDS("/home2/gkonop/workdir/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/MACAQUE_OUR_HE_INTEGRATED_Annotated_Caudate.RDS")
table(macaque_caud$Tissue)
table(macaque_caud$Species)
table(macaque_caud$newannot)
macaque_caud_SPN <- subset(macaque_caud, subset = newannot == "SPN")

macaque_put <- readRDS("/home2/gkonop/workdir/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/PUTAMEN/MACAQUE_OUR_HE_INTEGRATED_Annotated_Putamen.RDS")
macaque_put$Species <- rep("Macaque", nrow(macaque_put[[]]))
table(macaque_put$newannot)
table(macaque_put$Tissue)
table(macaque_put$Species)
macaque_put$tTissue <- NULL
macaque_put_SPN <- subset(macaque_put, subset = newannot == "SPN")

## integrate putamen and caudate
seurM <- merge(macaque_caud, macaque_put)
# Split data
seurM$id = paste0(seurM$orig.ident, '_', seurM$Tissue)
seurML = SplitObject(seurM, split.by = "id")

# normalize and identify variable features for each dataset independently
seurML = lapply(X = seurML, FUN = function(x) {
    x = NormalizeData(x)
    x = FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features = SelectIntegrationFeatures(object.list = seurML)
print(paste0('Number of genes to use for integration: ', length(features)))
seurML = lapply(X = seurML, FUN = function(x) {
    x = ScaleData(x, features = features, verbose = FALSE)
    x = RunPCA(x, features = features, verbose = FALSE)
})


## cca gives a mmerory problem try integrating with rpca
anchors = FindIntegrationAnchors(object.list = seurML, anchor.features = features, reduction = "rpca")
saveRDS(anchors, '/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/GB_All_cells_anchors_rpca_macaque_caud_put.RDS')
anchors = readRDS('/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/GB_SPN_anchors_rpca_macaque_caud_put.RDS')


# integrate 
macaque_integrated = IntegrateData(anchorset = anchors)
saveRDS(macaque_integrated, '/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/GB_Integrated_rpca_macaque_caud_put.RDS')
macaque_integrated <- readRDS( '/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/GB_Integrated_rpca_macaque_caud_put.RDS')

# set variables
marksToPlot = c('NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CFTR', 'CHAT', 'ADARB2', 'LHX6', 'MEIS2', 'TAC3', 'VIP', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA', 'SATB2', 'SLC17A7', 'EYA2', 'PDGFD', 'PTHLH', 'PVALB', 'SPARCL1', 'RSPO2', 'RMST')

# Clustering
DefaultAssay(macaque_integrated) = 'integrated'
macaque_integrated = ScaleData(macaque_integrated, verbose = FALSE)
macaque_integrated = RunPCA(macaque_integrated, verbose = FALSE)
macaque_integrated = RunUMAP(macaque_integrated, dims = 1:20, reduction = 'pca')

# Basic Plots
pref = 'GB_macaque_caudate_putamen_integrated_SPN_corrected'
setwd("/endosome/work/Neuroinformatics_Core/gkonop/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/Caudate_Putamen/CLUSTERING_1/")

macaque_integrated = FindNeighbors(macaque_integrated, dims = 1:20, reduction = 'pca')
macaque_integrated = FindClusters(macaque_integrated, resolution = 1)


pdf(paste0(pref, "_INTEGRATED_SAMPLES.pdf"))
DimPlot(macaque_integrated, group.by = "orig.ident", label = T, raster = T) + NoLegend()
dev.off()

pdf(paste0(pref, "_INTEGRATED_CLUSTERS.pdf"))
DimPlot(macaque_integrated, label = T, raster = T) + NoLegend()
dev.off()

stackedbarplot(macaque_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_INTEGRATED_STACK_SAMPLES'), wd = 20)
stackedbarplot(macaque_integrated[[]], groupx = 'seurat_clusters', groupfill = 'Species', paste0(pref, '_INTEGRATED_STACK_SPECIES'))

# Plot previously identified markers
DefaultAssay(macaque_integrated) = 'RNA'
macaque_integrated = NormalizeData(macaque_integrated)

pdf(paste0(pref, "_INTEGRATED_DOTPLOT.pdf"), width = 20, height = 10)
DotPlot(macaque_integrated, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "_INTEGRATED_DEPTH.pdf"))
ggboxplot(macaque_integrated[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'black') + NoLegend() +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_INTEGRATED_FEATUREPLOT.pdf"), width = 15, height = 35)
FeaturePlot(macaque_integrated, features = marksToPlot, sort = T, raster = T, pt.size = 2)
dev.off()

pdf(paste0(pref, "_INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(macaque_integrated[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

saveRDS(macaque_integrated, file= "/endosome/work/Neuroinformatics_Core/gkonop/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/Caudate_Putamen/CLUSTERING_1/macaque_integrated_cuadate_putamen_CLUSTERING_1.RDS")
macaque_integrated <- readRDS(file= "/endosome/work/Neuroinformatics_Core/gkonop/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/Caudate_Putamen/CLUSTERING_1/macaque_integrated_cuadate_putamen_CLUSTERING_1.RDS")

##########################################
###
## FILTER AND REPLOT 1
####
#Cluster#,Cell_type
#0,dSPN
#1,MOL
#2,iSPN
#3,MOL
#4,Astrocyte
#5,MOL
#6,MOL
#7,OPC
#8,MOL
#9,Microglia
#10,MOL(NRGN expressing)
#11,Astrocyte
#12,iSPN
#13,iSPN
#14,iSPN
#15,Astrocyte
#16,dSPN
#17,Non_SPN
#18,eSPN
#19,Non_SPN(TAC3)
#20*,MOL+Ast
#21,Astrocyte
#22,iSPN
#23,dSPN
#24,iSPN
#25,eSPN
#26*,Doublet
#27,Non_SPN
#28*,iSPN+eSPN+MOL
#29,dSPN
#30,eSPN
#31,Non_SPN
#32*,Exc_neu
#33*,iSPN+MOL
#34,Microglia
#35,MOL
#36,COP

## remove 
vasculature = c()
empty_drop = c()
doublet = c(20,26,28,32,33)
seurM = subset(macaque_integrated, subset = seurat_clusters %in% c(vasculature, empty_drop, doublet), invert = T)
seurM_corr = clusterBatchCorrect(seurM, npcs = 30, res = 1, pref = pref, sampleName = 'orig.ident')

setwd("/endosome/work/Neuroinformatics_Core/gkonop/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/Caudate_Putamen/CLUSTERING_2")

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
DefaultAssay(seurM_corr) = 'RNA'
seurM_corr = NormalizeData(seurM_corr)

pdf(paste0(pref, "_MarkersDotPlot.pdf"), width = 20, height = 10)
DotPlot(seurM_corr, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

# Save (optional)
saveRDS(seurM_corr, 'macaque_integrated_cuadate_putamen_CLUSTERING_2.RDS')
seurM_corr <- readRDS('macaque_integrated_cuadate_putamen_CLUSTERING_2.RDS')
######################################################################################################
###
## FILTER AND REPLOT 2
####
#Cluster#,Cell_type
#0,MOL
#1,iSPN
#2,MOL
#3,Astrocyte
#4,MOL
#5,dSPN
#6,iSPN
#7,OPC
#8,iSPN
#9,Astrocyte
#10,Microglia
#11,eSPN
#12,iSPN
#13,Non_SPN
#14,Non_SPN
#15,Astrocyte
#16,Non_SPN
#17,iSPN
#18,dSPN
#19,Non_SPN
#20,dSPN
#21,Non_SPN
#22*,MOL+Astrocyte
#23*,MOL+Astrocyte
#24,Astrocyte
#25,iSPN
#26,COP
#27*,Neu+MOL+Ast

## remove 
vasculature = c()
empty_drop = c()
doublet = c(22,23,27)
seurM = subset(seurM_corr, subset = seurat_clusters %in% c(vasculature, empty_drop, doublet), invert = T)
seurM_corr = clusterBatchCorrect(seurM, npcs = 30, res = 1, pref = pref, sampleName = 'orig.ident')

setwd("/home2/gkonop/workdir/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/Caudate_Putamen/CLUSTERING_3")


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
DefaultAssay(seurM_corr) = 'RNA'
seurM_corr = NormalizeData(seurM_corr)

pdf(paste0(pref, "_MarkersDotPlot.pdf"), width = 20, height = 10)
DotPlot(seurM_corr, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()


# Save
saveRDS(seurM_corr, 'macaque_integrated_cuadate_putamen_CLUSTERING_3.RDS')
seurM_corr <- readRDS('macaque_integrated_cuadate_putamen_CLUSTERING_3.RDS')
###################################################################################################
###
## ANNOTATE
####
#Cluster#,Cell_type
#0,MOL
#1,MOL
#2,MOL
#3,iSPN
#4,Astrocyte
#5,dSPN
#6,iSPN
#7,Astrocyte
#8,OPC
#9,iSPN
#10,Microglia
#11,eSPN
#12,MOL
#13,iSPN
#14,MOL
#15,Non_SPN
#16,Non_SPN
#17,Astrocyte
#18,iSPN
#19,Non_SPN
#20,dSPN
#21,Non_SPN
#22,dSPN
#23,Non_SPN
#24,Astrocyte
#25,iSPN
#26,COP
#27,Non_SPN

setwd("/home2/gkonop/workdir/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/Caudate_Putamen/ANNOTATED")

mapnames <- setNames(c("MOL", "MOL", "MOL", "SPN", "Astrocyte", "SPN", "SPN", 
"Astrocyte", "OPC", "SPN","Microglia","SPN","MOL","SPN",
"MOL", "Non_SPN", "Non_SPN", "Astrocyte","SPN","Non_SPN","SPN",
"Non_SPN","SPN","Non_SPN","Astrocyte", "SPN","OPC","Non_SPN"), c(0:27))

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

pdf(paste0(pref, "_MarkersDotPlot.pdf"), width = 25, height = 10)
DotPlot(seurM_corr, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

stackedbarplot(seurM_corr[[]], groupfill = 'newannot', groupx = 'Species', fn = paste0(pref, '_Stacked_Annotated_Species'))

stackedbarplot(seurM_corr[[]], groupfill = 'Species', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Species2'))

stackedbarplot(seurM_corr[[]], groupfill = 'newannot', groupx = 'orig.ident', fn = paste0(pref, '_Stacked_Annotated_Samples'))

stackedbarplot(seurM_corr[[]], groupfill = 'orig.ident', groupx = 'newannot', fn = paste0(pref, '_Stacked_AnnotatedSamples2'), wd = 20)

stackedbarplot(seurM_corr[[]], groupfill = 'Tissue', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Tissue'))


# Save
saveRDS(seurM_corr, 'macaque_integrated_cuadate_putamen_ANNOTATED.RDS')
seurM_corr <- readRDS('macaque_integrated_cuadate_putamen_ANNOTATED.RDS')
