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

# Only the caudate samples
allsamps = rio::import('~/workdir/TABLES/Primate_Tables_From_EmilyOh/Chimp Macaque Language Summary.xlsx')
colnames(allsamps) = gsub(' ', '_', colnames(allsamps))
allsamps$Species = gsub(' ', '_', allsamps$Species)
caudSamps = allsamps[allsamps$Brain_Region == 'Caudate' & allsamps$Species == 'Macaque',]
rownames(caudSamps) = paste0('Sample_', caudSamps[,1])

# Set the samples and read them with QC
fls = list.files('~/project/03_MATRIX_FROM_CELLBENDER/MACAQUE_OUR', pattern = "filtered.h5", recursive=T, full.names = T)
tmp = gsub('.*Macaque_|\\/CellBender.*', '', fls)
fls = fls[which(tmp %in% rownames(caudSamps))]

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
	tmpSeur$tissue = 'caudate'
	tmpSeur$Species = 'macaque'
	allSeurL[[i]] = tmpSeur

	print(i)
}

allseur = Reduce(merge, allSeurL)

saveRDS(allseur, '~/workdir/01_MACAQUE/OUR/02_CELLBENDERED_CLUSTER/Macaque_Caudate_Seurat_CB_Raw.RDS')
allseur = readRDS('~/workdir/01_MACAQUE/OUR/02_CELLBENDERED_CLUSTER/Macaque_Caudate_Seurat_CB_Raw.RDS')

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

pref = 'Macaque_Caudate_Our'

markerGenes = c('SLC17A7', 'SATB2', 'RBFOX3', 'NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA')


####
## CLUSTERING 1
####

seurM = allseur
seurM_corr = clusterBatchCorrect(seurM, npcs = 30, res = 1, pref = pref, sampleName = 'orig.ident')

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
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'nCount_RNA', ylim = c(0,50000)) +
rotate_x_text(90) 
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
dir.create('~/workdir/01_MACAQUE/OUR/02_CLUSTER_ANNOTATE/CAUDATE_ONLY/CLUSTERING_1')
saveRDS(seurM_corr, paste0('~/workdir/01_MACAQUE/OUR/02_CELLBENDERED_CLUSTER/', pref, '_GB_Clustering1.RDS'))
seurM_corr = readRDS(paste0('~/workdir/01_MACAQUE/OUR/02_CELLBENDERED_CLUSTER/', pref, '_GB_Clustering1.RDS'))

############################################################################3
# using https://www.sciencedirect.com/science/article/pii/S2211124720301911?via%3Dihub for SPN markers
# using article: Innovations present in the primate  interneuron repertoire by Fenna M. Krienen et al for non-SPN markers
# * represents cluster to be removed

#Cluster,CellType
#0,Oli
#1,Oli
#2,dSPN
#3,iSPN
#4*,Inh_neu/Oli
#5,Oli
#6,non-SPN
#7,Ast
#8,OPC
#9,Ast
#10,Microglia
#11,non-SPN+dSPN+eSPN
#12,Oli
#13*,non-SPN+dSPN+eSPN+Ast
#14,eSPN
#15,non-SPN?
#16,non-SPN
#17,non-SPN
#18*,Oli+Ast
#19*,Endo
#20*,Oli+Microglia
#21*,Oli+non-SPN
#22*,Endo+Peri
#23,non-SPN
#24*,non-SPN+Microglia
#25*,Ependymal
#26*,Endo
#27*,OPC+Oli
#28,iSPN
#29,non-SPN
#30*,Microglia+Ast+Peri?
#31*,OPC+Ast
#32,non-SPN
#33*,Microglia+OPC+non-SPN
#34,OPC
#35*,Oli+microglia
#36,COP

############################################################################

####
## CLUSTERING 2
####

setwd("/home2/gkonop/workdir/01_MACAQUE/OUR/02_CELLBENDERED_CLUSTER/CAUDATE/CLUSTERING_2")
vasculature = c(19,22,25,26)
empty_drop = c()
doublet = c(4,13,18,20,21,24,27,30,31,33,35)
seurM = subset(seurM_corr, subset = seurat_clusters %in% c(vasculature, empty_drop, doublet), invert = T)
seurM_corr = clusterBatchCorrect(seurM, npcs = 30, res = 1, pref = pref, sampleName = 'orig.ident')

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
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'nCount_RNA', ylim = c(0,50000)) +
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


############################################################################3
# using https://www.sciencedirect.com/science/article/pii/S2211124720301911?via%3Dihub for SPN markers
# using article: Innovations present in the primate  interneuron repertoire by Fenna M. Krienen et al for non-SPN markers

#Cluster,CellType
#0,Oli
#1,Oli
#2,iSPN
#3,dSPN
#4,Oli
#5,Ast
#6*,non-SPN(low genes/umi)
#7,OPC
#8,Microglia
#9,Oli
#10*,non-SPN(low genes/umi)
#11,eSPN
#12,non-SPN
#13,dSPN
#14,Microglia
#15,iSPN
#16*,non-SPN(low genes/umi)
#17,non-SPN
#18,non-SPN
#19,non-SPN
#20,Oli
#21,iSPN
#22*,non-SPN+Oli
#23,non-SPN
#24,Oli
#25,non-SPN
#26*,Micro-Ast-Oli-OPC
#27*,eSPN+Micro
#28,COP
#29,iSPN
############################################################################


# Save (optional)
dir.create('~/workdir/01_MACAQUE/OUR/02_CLUSTER_ANNOTATE/CAUDATE_ONLY/CLUSTERING_2')
saveRDS(seurM_corr, paste0('~/workdir/01_MACAQUE/OUR/02_CELLBENDERED_CLUSTER/', pref, '_Clustering2.RDS'))

qw = FindMarkers(seurM_corr, ident.1 = 17, only.pos = T, logfc.threshold = 0.5, min.pct =0.5)

####
## CLUSTERING 3
####
setwd("/home2/gkonop/workdir/01_MACAQUE/OUR/02_CELLBENDERED_CLUSTER/CAUDATE/CLUSTERING_3")
vasculature = c()
empty_drop = c(6,10,16)
doublet = c(22,26,27)
seurM = subset(seurM_corr, subset = seurat_clusters %in% c(vasculature, empty_drop, doublet), invert = T)
seurM_corr = clusterBatchCorrect(seurM, npcs = 30, res = 1, pref = pref, sampleName = 'orig.ident')

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
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'nCount_RNA', ylim = c(0,50000)) +
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

############################################################################3
# using https://www.sciencedirect.com/science/article/pii/S2211124720301911?via%3Dihub for SPN markers
# using article: Innovations present in the primate  interneuron repertoire by Fenna M. Krienen et al for non-SPN markers

#Cluster,CellType
#0,Oli
#1,Oli
#2,Oli
#3,Oli
#4,iSPN
#5,dSPN
#6,MOL(Oli)
#7,Ast
#8,OPC
#9,iSPN
#10,dSPN
#11,Microglia
#12,Ast
#13,eSPN
#14,dSPN
#15,non-SPN
#16,iSPN
#17,non-SPN
#18,iSPN
#19,Ast
#20*,doublet
#21,non-SPN
#22,iSPN
#23,non-SPN
#24*,Exc_neu
#25,MOL
#26,non-SPN
#27,iSPN
#28,COP
#29,MOL
#30,Microglia
############################################################################


# Save (optional)
dir.create('~/workdir/01_MACAQUE/OUR/02_CLUSTER_ANNOTATE/CLUSTERING_3')
saveRDS(seurM_corr, paste0('~/workdir/01_MACAQUE/OUR/02_CELLBENDERED_CLUSTER/', pref, '_GB_Clustering3.RDS'))
seurM_corr <- readRDS(paste0('~/workdir/01_MACAQUE/OUR/02_CELLBENDERED_CLUSTER/', pref, '_GB_Clustering3.RDS'))

####
## CLUSTERING 4
####
setwd("/home2/gkonop/workdir/01_MACAQUE/OUR/02_CELLBENDERED_CLUSTER/CAUDATE/CLUSTERING_4")
vasculature = c()
empty_drop = c()
doublet = c(20,24)
neurogenesis = c()
seurM = subset(seurM_corr, subset = seurat_clusters %in% c(vasculature, empty_drop, doublet, neurogenesis), invert = T)
seurM_corr = clusterBatchCorrect(seurM, npcs = 30, res = 1, pref = pref, sampleName = 'orig.ident')

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
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'nCount_RNA', ylim = c(0,50000)) +
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
dir.create('~/workdir/MACAQUE/OUR/02_CELLBENDERED_CLUSTER/CLUSTERING_4')
saveRDS(seurM_corr, paste0('~/workdir/MACAQUE/OUR/02_CELLBENDERED_CLUSTER/', pref, '_Clustering4.RDS'))
 

####
## BROADLY ANNOTATE
####

############################################################################3
# using https://www.sciencedirect.com/science/article/pii/S2211124720301911?via%3Dihub for SPN markers
# using article: Innovations present in the primate  interneuron repertoire by Fenna M. Krienen et al for non-SPN markers

#Cluster,CellType
#0,MOL
#1,MOL
#2,MOL
#3,Ast
#4,dSPN
#5,iSPN
#6,iSPN
#7,OPC
#8,MOL
#9,Microglia
#10,dSPN
#11,eSPN
#12,Ast
#13,dSPN
#14,non-SPN
#15,iSPN
#16,non-SPN
#17,iSPN
#18,MOL
#19,Ast
#20,non-SPN
#21,iSPN
#22,non-SPN
#23,MOL
#24,non-SPN
#25,COP
#26,iSPN
#27,MOL
#28,Microglia

############################################################################


seurM_corr = readRDS(paste0('~/workdir/MACAQUE/OUR/02_CELLBENDERED_CLUSTER/', pref, '_GB_Clustering4.RDS'))

mapnames = setNames(c('MOL', 'MOL', 'MOL','Astrocyte', 'SPN', 'SPN', 'SPN','OPC', 'MOL', 'Microglia', 'SPN', 'SPN', 'Astrocyte','SPN','non-SPN','SPN','non-SPN','SPN','MOL',
			'Astrocyte', 'non-SPN', 'SPN','non-SPN','MOL', 'non-SPN','COP','SPN','MOL','Microglia'),
		      c(0:28))

seurM_corr[["newannot"]] = mapnames[seurM_corr[["seurat_clusters"]][,1]]
Idents(seurM_corr) = seurM_corr$newannot

setwd("/home2/gkonop/workdir/01_MACAQUE/OUR/02_CELLBENDERED_CLUSTER/CAUDATE/ANNOTATION")
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

stackedbarplot(seurM_corr[[]], groupfill = 'newannot', groupx = 'orig.ident', fn = paste0(pref, '_Caudate_Annotated'))
stackedbarplot(seurM_corr[[]], groupfill = 'newannot', groupx = 'Age', fn = paste0(pref, '_Age_Caudate_Annotated'))

# Save
saveRDS(seurM_corr, "/home2/gkonop/workdir/01_MACAQUE/OUR/02_CELLBENDERED_CLUSTER/CAUDATE/ANNOTATION/Macaque_Caudate_Our_Annotated_FINAL.RDS")
