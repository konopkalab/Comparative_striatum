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
###############
## integrate all
######################################################################################################
### Load marmoset data
marmoset_caud = readRDS('/home2/gkonop/workdir/01_MARMOSET/MARMOSET_INTEGRATION/GB_CAUDATE/ANNOTATION/LIN_KRIENEN_INTEGRATED_CAUDATE_FINAL_ANNOTATION.RDS')
table(marmoset_caud$Species)
table(marmoset_caud$newannot)
table(marmoset_caud$Tissue)
marmoset_caud_SPN <- subset(marmoset_caud, subset = newannot == "SPN")
dim(marmoset_caud_SPN) #[1] 48990 26880

marmoset_put = readRDS('/home2/gkonop/workdir/01_MARMOSET/MARMOSET_KRIENEN/03_CELLBENDERED_ANNOTATE/PUTAMEN/KRIENEN_PUTAMEN_ANNOTATED.RDS')
table(marmoset_put$Species)
table(marmoset_put$newannot)
table(marmoset_put$Tissue)
marmoset_put_SPN <-subset(marmoset_put, subset = newannot == "SPN")


## integrate putamen and caudate
seurM <- merge(marmoset_caud, marmoset_put)
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
saveRDS(anchors, '/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/GB_All_cells_anchors_rpca_marmoset_caud_put.RDS')
anchors = readRDS('/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/GB_SPN_anchors_rpca_marmoset_caud_put.RDS')


# integrate 
marmoset_integrated = IntegrateData(anchorset = anchors)
saveRDS(marmoset_integrated, '/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/GB_Integrated_rpca_marmoset_caud_put.RDS')
marmoset_integrated <- readRDS( '/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/GB_Integrated_rpca_marmoset_caud_put.RDS')


######
## CLUSTER AND FILTER 
######
# Clustering
DefaultAssay(marmoset_integrated) = 'integrated'
marmoset_integrated = ScaleData(marmoset_integrated, verbose = FALSE)
marmoset_integrated = RunPCA(marmoset_integrated, verbose = FALSE)
marmoset_integrated = RunUMAP(marmoset_integrated, dims = 1:20, reduction = 'pca')

# Basic Plots
pref = 'GB_marmoset_caudate_putamen_integrated'
setwd("/endosome/work/Neuroinformatics_Core/gkonop/01_MARMOSET/MARMOSET_INTEGRATION/Caudate_putamen/CLUSTERING_1")

marmoset_integrated = FindNeighbors(marmoset_integrated, dims = 1:20, reduction = 'pca')
marmoset_integrated = FindClusters(marmoset_integrated, resolution = 1)


pdf(paste0(pref, "_INTEGRATED_SAMPLES.pdf"))
DimPlot(marmoset_integrated, group.by = "orig.ident", label = T, raster = T) + NoLegend()
dev.off()

pdf(paste0(pref, "_INTEGRATED_CLUSTERS.pdf"))
DimPlot(marmoset_integrated, label = T, raster = T) + NoLegend()
dev.off()

stackedbarplot(marmoset_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_INTEGRATED_STACK_SAMPLES'), wd = 20)
stackedbarplot(marmoset_integrated[[]], groupx = 'seurat_clusters', groupfill = 'Species', paste0(pref, '_INTEGRATED_STACK_SPECIES'))

# Plot previously identified markers
DefaultAssay(marmoset_integrated) = 'RNA'
marmoset_integrated = NormalizeData(marmoset_integrated)

pdf(paste0(pref, "_INTEGRATED_DOTPLOT.pdf"), width = 20, height = 10)
DotPlot(marmoset_integrated, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "_INTEGRATED_DEPTH.pdf"))
ggboxplot(marmoset_integrated[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'black') + NoLegend() +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_INTEGRATED_FEATUREPLOT.pdf"), width = 15, height = 35)
FeaturePlot(marmoset_integrated, features = marksToPlot, sort = T, raster = T, pt.size = 2)
dev.off()

pdf(paste0(pref, "_INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(marmoset_integrated[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()


saveRDS(marmoset_integrated, file= "/endosome/work/Neuroinformatics_Core/gkonop/01_MARMOSET/MARMOSET_INTEGRATION/Caudate_putamen/CLUSTERING_1/marmoset_integrated_cuadate_putamen_CLUSTERING_1.RDS")
marmoset_integrated <- readRDS(file= "/endosome/work/Neuroinformatics_Core/gkonop/01_MARMOSET/MARMOSET_INTEGRATION/Caudate_putamen/CLUSTERING_1/marmoset_integrated_cuadate_putamen_CLUSTERING_1.RDS")
###
## CLUSTER AND FILTER
####
#Cluster#,Cell_type
#0,dSPN
#1,MOL
#2,Astrocyte
#3,iSPN
#4,iSPN
#5,MOL
#6,MOL
#7,dSPN
#8,dSPN
#9,OPC
#10,Astrocyte
#11,Microglia
#12,dSPN
#13,iSPN
#14,MOL
#15,eSPN
#16,Astrocyte
#17,MOL
#18,iSPN (TAC3 and VIP expressing)
#19,Non_SPN
#20,iSPN
#21,Astrocyte
#22,dSPN
#23,iSPN
#24*,Empty
#25,Non_SPN
#26*,Doublet
#27,iSPN(TAC1 expressing)
#28,Non_SPN
#29,MOL
#30,Non_SPN (CHAT expressing)
#31,COP


setwd("/endosome/work/Neuroinformatics_Core/gkonop/01_MARMOSET/MARMOSET_INTEGRATION/Caudate_putamen/CLUSTERING_2")

## remove 
vasculature = c()
empty_drop = c(24)
doublet = c(26)
seurM = subset(marmoset_integrated, subset = seurat_clusters %in% c(vasculature, empty_drop, doublet), invert = T)
seurM_corr = clusterBatchCorrect(seurM, npcs = 30, res = 1, pref = pref, sampleName = 'orig.ident')


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
saveRDS(seurM_corr, 'marmoset_integrated_cuadate_putamen_CLUSTERING_2.RDS')
seurM_corr <- readRDS('marmoset_integrated_cuadate_putamen_CLUSTERING_2.RDS')
###
## CLUSTER AND FILTER 2
####
#Cluster#,Cell_type
#0,dSPN
#1,iSPN
#2,MOL
#3,Astrocyte
#4,OPC
#5,MOL
#6,Astrocyte
#7,Microglia
#8,iSPN
#9,MOL
#10,dSPN
#11,MOL
#12,iSPN
#13,dSPN
#14,eSPN
#15,dSPN
#16,Astrocyte
#17,MOL
#18,iSPN
#19,Non_SPN
#20,Non_SPN
#21,dSPN
#22,Non_SPN
#23,iSPN
#24,eSPN
#25,iSPN
#26,Non_SPN
#27,Non_SPN



setwd("/endosome/work/Neuroinformatics_Core/gkonop/01_MARMOSET/MARMOSET_INTEGRATION/Caudate_putamen/ANNOTATED")

# annotate cell types
mapnames <- setNames(c("SPN", "SPN","MOL","Astrocyte","OPC","MOL","Astrocyte",
"Microglia", "SPN", "MOL", "SPN", "MOL", "SPN","SPN",
"SPN","SPN","Astrocyte","MOL","SPN","Non_SPN","Non_SPN",
"SPN","Non_SPN","SPN","SPN","SPN","Non_SPN","Non_SPN"), c(0:27))

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
saveRDS(seurM_corr, 'marmoset_integrated_cuadate_putamen_ANNOTATED.RDS')
