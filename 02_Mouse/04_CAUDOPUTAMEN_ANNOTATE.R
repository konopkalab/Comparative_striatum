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
samps = paste0('SRR1192', 1005:1012)

allSeurL = list()
for(i in 1:length(samps)){

	#  Read CellBender output
	mat = Read10X_h5(paste0('~/project/03_MATRIX_FROM_CELLBENDER/MOUSE_KRIENEN/', samps[i], '/CellBender_out_filtered.h5'))

	# Soft filtering. Remove lower than 200 UMI
	mat = mat[, colSums(mat) > 200]
	cbBarcode = colnames(mat) %>% gsub('-1', '', .)
	colnames(mat) = gsub('-1', '', colnames(mat)) %>% paste0(., '_', samps[i])

	# Remove barcodes that are likely false (not necessary for 10x)
	#mat = mat[,!(grepl('AAAAA|CCCCC|TTTTT|GGGGG', colnames(mat)))]

	# Create seurat object
	tmpSeur = CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0)

	print(i)
}
allseur = Reduce(merge, allSeurL)

write_rds(allseur, '~/workdir/01_MOUSE/03_CLUSTER_ANNOTATE/GB_Mouse_CaudoPutamen_Seurat_CB_Raw.RDS')
allseur = read_rds('~/workdir/01_MOUSE/03_CLUSTER_ANNOTATE/GB_Mouse_CaudoPutamen_Seurat_CB_Raw.RDS')

####
## SET VARIABLES
####

pref = 'Mouse_Caudate'

markerGenes = c('FOSL2', 'NPAS4', 'NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA', 'DCN', 'CP', 'BFSP2', 'TNR')
markerGenes = tolower(markerGenes)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

markerGenes = firstup(markerGenes)

####
## CLUSTERING 1
####

seurM = allseur
seurM = subset(seurM, subset = nCount_RNA > 0)
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
todir = '~/workdir/01_MOUSE/03_CLUSTER_ANNOTATE/CLUSTERING_1/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)

# Save (optional)
write_rds(seurM_corr, paste0('~/workdir/01_MOUSE/03_CLUSTER_ANNOTATE/', pref, '_Clustering1.RDS'))

qw = FindMarkers(seurM_corr, ident.1 = 22, only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
DotPlot(seurM_corr, features = rownames(qw)[1:30]) + rotate_x_text(90)

####
## CLUSTERING 2
####

vasculature = c(10,13,17,21,22)
neurogenesis = c(12)
empty_drop = c(15)
doublet = c()
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
write_rds(seurM_corr, paste0('~/workdir/01_MOUSE/03_CLUSTER_ANNOTATE/', pref, '_Clustering2.RDS'))

# Move the plots
outfls = list.files('~/', pattern = pref, full.names = F)
todir = '~/workdir/01_MOUSE/03_CLUSTER_ANNOTATE/CLUSTERING_2/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)


qw = FindMarkers(seurM_corr, ident.1 = 8, only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)

pdf('Markers_cluster_19.pdf', width = 20)
DotPlot(seurM_corr, features = rownames(qw)[1:20])
dev.off()



####
## CLUSTERING 3
####

vasculature = c()
neurogenesis = c()
empty_drop = c()
doublet = c(8)
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
write_rds(seurM_corr, paste0('~/workdir/01_MOUSE/03_CLUSTER_ANNOTATE/', pref, '_Clustering3.RDS'))

# Continute from here Oct 21, 2023
seurM_corr = read_rds(paste0('~/workdir/01_MOUSE/03_CLUSTER_ANNOTATE/', pref, '_Clustering3.RDS'))

qw = FindMarkers(seurM_corr, ident.1 = 17, only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)

pdf('Markers_cluster_19.pdf', width = 20)
DotPlot(seurM_corr, features = rownames(qw)[1:20])
dev.off()


####
## BROADLY ANNOTATE
####

seurM_corr = read_rds(paste0('~/workdir/01_MOUSE/03_CLUSTER_ANNOTATE/', pref, '_Clustering3.RDS'))
seurM_corr = subset(seurM_corr, subset = seurat_clusters == 17, invert = T)
seurM_corr$seurat_clusters = droplevels(seurM_corr$seurat_clusters)

mapnames = setNames(c('SPN', 'SPN', 'SPN', 'SPN', 'Astrocyte', 'SPN',
			'SPN', 'MOL', 'SPN', 'Microglia', 'SPN',
			'Non_SPN', 'OPC', 'Astrocyte', 'Non_SPN', 'Astrocyte',
			'Non_SPN', 'COP'),
		      c(0:16,18))

seurM_corr[["newannot"]] = mapnames[seurM_corr[["seurat_clusters"]][,1]]
Idents(seurM_corr) = seurM_corr$newannot

# PLOTS
pdf(paste0(pref, "_Clusters.pdf"))
DimPlot(seurM_corr, group.by = 'newannot', label = T, raster = T, label.size = 7, repel = T) +
NoLegend() +
ggtitle('')
dev.off()

pdf(paste0(pref, "_Depth_Gene.pdf"))
ggboxplot(seurM_corr[[]], x = 'newannot', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_Depth_UMI.pdf"))
ggboxplot(seurM_corr[[]], x = 'newannot', y = 'nCount_RNA', ylim = c(0,80000)) +
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

# Save
write_rds(seurM_corr, paste0('~/workdir/01_MOUSE/03_CLUSTER_ANNOTATE/', pref, '_Annotated_FINAL.RDS'))
