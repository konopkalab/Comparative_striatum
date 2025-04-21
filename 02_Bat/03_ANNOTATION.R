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

# Set directories
path_cr = '~/project/01_CELLRANGER_OUT/BAT_OUT/'
path_cb = '~/project/03_MATRIX_FROM_CELLBENDER/BAT_OUR/'
samps = paste0('Caud_', 2:5) # 1st run is bad

#####
## LOAD COUNT MATRIX
#####

allSeurL = list()
for(i in 1:length(samps)){

	#  Read CellBender output
	mat = Read10X_h5(paste0(path_cb, samps[i], '/CellBender_out_filtered.h5'))

	# Soft filtering. Remove lower than 200 UMI
	mat = mat[, colSums(mat) > 200]
	cbBarcode = colnames(mat)
	colnames(mat) = gsub('-1', '', colnames(mat)) %>% paste0(., '-', samps[i])

	# Create seurat object
	tmpSeur = CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0)

	# Load intronic read ratio
	fls = list.files(paste0(path_cr, samps[i], '/outs'), pattern = 'Intronic', full.names=T)

	if(any(grepl('RDS', fls)) == T){

		nf1 = readRDS(fls[grepl('RDS', fls)])
		print('Intronic_read_ratio_loaded')
		
	} else{

		# Calculate intronic read ratio
		nf1 = nuclear_fraction_tags(
			bam = paste0(path_cr, samps[i], '/outs/possorted_genome_bam.bam'),
			barcodes = cbBarcode,
			cores = 23, verbose = T)

		saveRDS(nf1, paste0(path_cr, samps[i], '/outs/Intronic_read_ratio.RDS'))
	}

	# Keep the final barcodes
	nf1 = nf1[match(cbBarcode, rownames(nf1)),]

	tmpSeur$intronRat = nf1
	tmpSeur$orig.ident = samps[i]
	tmpSeur$tissue = 'Caudate'
	tmpSeur$Species = 'Bat'
	allSeurL[[i]] = tmpSeur

	print(i)
}

allseur = Reduce(merge, allSeurL)

dir.create('~/workdir/01_BAT/02_CLUSTER_ANNOTATE/CAUDATE')
write_rds(allseur, '~/workdir/01_BAT/02_CLUSTER_ANNOTATE/CAUDATE/Bat_Caudate_Seurat_CB_Raw.RDS')
allseur = read_rds('~/workdir/01_BAT/02_CLUSTER_ANNOTATE/CAUDATE/Bat_Caudate_Seurat_CB_Raw.RDS')

####
## SET VARIABLES
####

pref = 'Bat_Caudate'

markerGenes = c('NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'CAUDATER17', 'FLT1', 'DUSP1', 'COBLL1', 'CFTR', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA', 'SATB2', 'SLC17A7')

####
## CLUSTERING 1
####

seurM = allseur
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
write_rds(seurM_corr, '~/workdir/01_BAT/02_CLUSTER_ANNOTATE/CAUDATE/Bat_Caudate_Seurat_CB_Clustering1.RDS')

qw = FindMarkers(seurM_corr, ident.1 = 30, only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
DotPlot(seurM_corr, features = rownames(qw)[1:30]) + rotate_x_text(90)

# Move the plots
outfls = list.files('~/', pattern = pref, full.names = F)
todir = '~/workdir/01_BAT/02_CLUSTER_ANNOTATE/CAUDATE/CLUSTERING_1/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)

####
## CLUSTERING 2
####

vasculature = c(20,21,31)
empty_drop = c(11)
doublet = c(5,8,9,12,24,30)
others = c()
seurM = subset(seurM_corr, subset = seurat_clusters %in% c(vasculature, empty_drop, doublet, others), invert = T)
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
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'nCount_RNA', ylim = c(0,80000)) +
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
write_rds(seurM_corr, '~/workdir/01_BAT/02_CLUSTER_ANNOTATE/CAUDATE/Bat_Caudate_Seurat_CB_Clustering2.RDS')

qw = FindMarkers(seurM_corr, ident.1 = 22, only.pos = T, logfc.threshold=0.5, min.pct = 0.5)
DotPlot(seurM_corr, features = rownames(qw)[1:30]) + rotate_x_text(90)

# Move the plots
outfls = list.files('~/', pattern = pref, full.names = F)
todir = '~/workdir/01_BAT/02_CLUSTER_ANNOTATE/CAUDATE/CLUSTERING_2/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)


####
## CLUSTERING 3
####

vasculature = c()
empty_drop = c(6)
doublet = c(21)
others = c(10,12,15,22)
seurM = subset(seurM_corr, subset = seurat_clusters %in% c(vasculature, empty_drop, doublet, others), invert = T)
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
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'nCount_RNA', ylim = c(0,80000)) +
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
write_rds(seurM_corr, '~/workdir/01_BAT/02_CLUSTER_ANNOTATE/CAUDATE/Bat_Caudate_Seurat_CB_Clustering3.RDS')

qw = FindMarkers(seurM_corr, ident.1 = 7, only.pos = T, logfc.threshold=0.5, min.pct = 0.5)
DotPlot(seurM_corr, features = rownames(qw)[1:30]) + rotate_x_text(90)

# Move the plots
outfls = list.files('~/', pattern = pref, full.names = F)
todir = '~/workdir/01_BAT/02_CLUSTER_ANNOTATE/CAUDATE/CLUSTERING_3/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)

####
## BROADLY ANNOTATE
####

mapnames = setNames(c('SPN', 'SPN', 'Astrocyte', 'SPN', 'MOL', 'SPN',
			'OPC', 'Non_SPN', 'SPN', 'Astrocyte', 'Non_SPN',
			'Non_SPN', 'Astrocyte', 'Non_SPN', 'Non_SPN', 'SPN',
			'Non_SPN'),
		      c(0:10))

seurM_corr[["newannot"]] = mapnames[seurM_corr[["seurat_clusters"]][,1]]
Idents(seurM_corr) = seurM_corr$newannot

# PLOTS
pdf(paste0(pref, "_Clusters.pdf"))
DimPlot(seurM_corr, group.by = 'newannot', label = T, raster = T, label.size = 7) +
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

pdf(paste0(pref, "_NuclearFraction.pdf"))
ggboxplot(seurM_corr[[]], x = 'newannot', y = 'intronRat') +
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

stackedbarplot(seurM_corr[[]], groupfill = 'newannot', groupx = 'orig.ident', fn = paste0(pref, '_Stacked_Annotated'))

stackedbarplot(seurM_corr[[]], groupfill = 'orig.ident', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_2'))

# Move the plots
outfls = list.files('~/', pattern = pref, full.names = F)
todir = '~/workdir/01_BAT/02_CLUSTER_ANNOTATE/CAUDATE/ANNOTATION/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)

# Save
write_rds(seurM_corr, '~/workdir/01_BAT/02_CLUSTER_ANNOTATE/CAUDATE/Bat_Caudate_Seurat_CB_Annotated.RDS')

#### Putamen
# Set directories
path_cr = '~/project/01_CELLRANGER_OUT/BAT_OUT/'
path_cb = '~/project/03_MATRIX_FROM_CELLBENDER/BAT_OUR/'
samps = paste0('Put_', c(1,3:5)) # 2nd run is bad

#####
## LOAD COUNT MATRIX
#####

allSeurL = list()
for(i in 1:length(samps)){

	#  Read CellBender output
	mat = Read10X_h5(paste0(path_cb, samps[i], '/CellBender_out_filtered.h5'))

	# Soft filtering. Remove lower than 200 UMI
	mat = mat[, colSums(mat) > 200]
	cbBarcode = colnames(mat)
	colnames(mat) = gsub('-1', '', colnames(mat)) %>% paste0(., '-', samps[i])

	# Create seurat object
	tmpSeur = CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0)

	# Load intronic read ratio
	fls = list.files(paste0(path_cr, samps[i], '/outs'), pattern = 'Intronic', full.names=T)

	if(any(grepl('RDS', fls)) == T){

		nf1 = readRDS(fls[grepl('RDS', fls)])
		print('Intronic_read_ratio_loaded')
		
	} else{

		# Calculate intronic read ratio
		nf1 = nuclear_fraction_tags(
			bam = paste0(path_cr, samps[i], '/outs/possorted_genome_bam.bam'),
			barcodes = cbBarcode,
			cores = 23, verbose = T)

		saveRDS(nf1, paste0(path_cr, samps[i], '/outs/Intronic_read_ratio.RDS'))
	}

	# Keep the final barcodes
	nf1 = nf1[match(cbBarcode, rownames(nf1)),]

	tmpSeur$intronRat = nf1
	tmpSeur$orig.ident = samps[i]
	tmpSeur$tissue = 'Putamen'
	tmpSeur$Species = 'Bat'
	allSeurL[[i]] = tmpSeur

	print(i)
}

allseur = Reduce(merge, allSeurL)

dir.create('~/workdir/01_BAT/02_CLUSTER_ANNOTATE/PUTAMEN')
write_rds(allseur, '~/workdir/01_BAT/02_CLUSTER_ANNOTATE/PUTAMEN/Bat_Putamen_Seurat_CB_Raw.RDS')
allseur = read_rds('~/workdir/01_BAT/02_CLUSTER_ANNOTATE/PUTAMEN/Bat_Putamen_Seurat_CB_Raw.RDS')

####
## SET VARIABLES
####

pref = 'Bat_Putamen'

markerGenes = c('NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'PUTAMENR17', 'FLT1', 'DUSP1', 'COBLL1', 'CFTR', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA', 'SATB2', 'SLC17A7')

####
## CLUSTERING 1
####

seurM = allseur
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
write_rds(seurM_corr, '~/workdir/01_BAT/02_CLUSTER_ANNOTATE/PUTAMEN/Bat_Putamen_Seurat_CB_Clustering1.RDS')

qw = FindMarkers(seurM_corr, ident.1 = 0, only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
DotPlot(seurM_corr, features = rownames(qw)[1:30]) + rotate_x_text(90)

# Move the plots
outfls = list.files('~/', pattern = pref, full.names = F)
todir = '~/workdir/01_BAT/02_CLUSTER_ANNOTATE/PUTAMEN/CLUSTERING_1/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)

####
## CLUSTERING 2
####

vasculature = c(16)
empty_drop = c(27)
doublet = c(6,9,17,26)
others = c()
seurM = subset(seurM_corr, subset = seurat_clusters %in% c(vasculature, empty_drop, doublet, others), invert = T)
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
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'nCount_RNA', ylim = c(0,80000)) +
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
write_rds(seurM_corr, '~/workdir/01_BAT/02_CLUSTER_ANNOTATE/PUTAMEN/Bat_Putamen_Seurat_CB_Clustering2.RDS')

qw = FindMarkers(seurM_corr, ident.1 = 17, only.pos = T, logfc.threshold=0.5, min.pct = 0.5)
DotPlot(seurM_corr, features = rownames(qw)[1:30]) + rotate_x_text(90)

# Move the plots
outfls = list.files('~/', pattern = pref, full.names = F)
todir = '~/workdir/01_BAT/02_CLUSTER_ANNOTATE/PUTAMEN/CLUSTERING_2/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)

####
## CLUSTERING 3
####

vasculature = c(23)
empty_drop = c()
doublet = c(17)
others = c()
seurM = subset(seurM_corr, subset = seurat_clusters %in% c(vasculature, empty_drop, doublet, others), invert = T)
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
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'nCount_RNA', ylim = c(0,80000)) +
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
write_rds(seurM_corr, '~/workdir/01_BAT/02_CLUSTER_ANNOTATE/PUTAMEN/Bat_Putamen_Seurat_CB_Clustering3.RDS')

qw = FindMarkers(seurM_corr, ident.1 = 22, only.pos = T, logfc.threshold=0.5, min.pct = 0.5)
DotPlot(seurM_corr, features = rownames(qw)[1:30]) + rotate_x_text(90)

# Move the plots
outfls = list.files('~/', pattern = pref, full.names = F)
todir = '~/workdir/01_BAT/02_CLUSTER_ANNOTATE/PUTAMEN/CLUSTERING_3/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)



####
## BROADLY ANNOTATE
####

mapnames = setNames(c('MOL', 'SPN', 'SPN', 'Astrocyte', 'MOL', 'Microglia',
			'OPC', 'MOL', 'SPN', 'MOL', 'Non_SPN',
			'Astrocyte', 'Non_SPN', 'Non_SPN', 'Microglia', 'SPN',
			'OPC', 'Non_SPN', 'Astrocyte', 'Other_1', 'Non_SPN',
			'Non_SPN', 'Other_2'),
		      c(0:22))

seurM_corr[["newannot"]] = mapnames[seurM_corr[["seurat_clusters"]][,1]]
Idents(seurM_corr) = seurM_corr$newannot

# PLOTS
pdf(paste0(pref, "_Clusters.pdf"))
DimPlot(seurM_corr, group.by = 'newannot', label = T, raster = T, label.size = 7) +
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

pdf(paste0(pref, "_NuclearFraction.pdf"))
ggboxplot(seurM_corr[[]], x = 'newannot', y = 'intronRat') +
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

stackedbarplot(seurM_corr[[]], groupfill = 'newannot', groupx = 'orig.ident', fn = paste0(pref, '_Stacked_Annotated'))

stackedbarplot(seurM_corr[[]], groupfill = 'orig.ident', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_2'))

# Move the plots
outfls = list.files('~/', pattern = pref, full.names = F)
todir = '~/workdir/01_BAT/02_CLUSTER_ANNOTATE/PUTAMEN/ANNOTATION/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)

# Save
write_rds(seurM_corr, '~/workdir/01_BAT/02_CLUSTER_ANNOTATE/PUTAMEN/Bat_Putamen_Seurat_CB_Annotated.RDS')
