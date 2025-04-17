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
caudSamps = allsamps[allsamps$Brain_Region == 'Caudate' & allsamps$Species == 'Chimp',]
rownames(caudSamps) = paste0('Sample_', caudSamps[,1])

# Set the samples and read them with QC
samps = rownames(caudSamps)
fls = list.files('~/project/03_MATRIX_FROM_CELLBENDER/CHIMP_OUR', pattern = "filtered.h5", recursive=T, full.names = T)
fls = fls[grepl('Chimp_Sample', fls)]

allSeurL = list()
for(i in 1:length(fls)){

	# Sample name
	samp = gsub('.*Chimp_', '', fls[i]) %>% gsub('\\/.*', '', .)

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
		annotation_path = '/home2/gkonop/workdir/REFERENCES/GTF_ONLY/hg38_liftoverTo_panTro5.gtf',
		bam = paste0('~/project/01_CELLRANGER_OUT/CHIMP_OUR/celllranger_count_', samp, '/outs/possorted_genome_bam.bam'),
		barcodes = cbBarcode,
		cores = 23, verbose = T)

	nf2$barc = colnames(tmpSeur)
	tmpSeur$intronRat = nf2$nuclear_fraction
	tmpSeur$orig.ident = samp
	tmpSeur$tissue = 'caudate'
	tmpSeur$Species = 'chimp'
	allSeurL[[i]] = tmpSeur

	print(i)
}

allseur = Reduce(merge, allSeurL)

save_rds(allseur, '~/workdir/01_CHIMP/02_CELLBENDERED_CLUSTER/Chimp_Caudate_Seurat_CB_Raw.RDS')
allseur = read_rds('~/workdir/01_CHIMP/02_CELLBENDERED_CLUSTER/Chimp_Caudate_Seurat_CB_Raw.RDS')

# Update matrix if the initial filtering was done incorrectly
#allseur_orig = read_rds('~/workdir/01_CHIMP/02_CELLBENDERED_CLUSTER/Chimp_Caudate_Seurat_CB_Raw.RDS')
#meta_orig = allseur_orig@meta.data

#mat = allseur@assays$RNA@counts
#newSeur = CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0, meta.data = meta_orig)
#save_rds(newSeur, '~/workdir/01_CHIMP/02_CELLBENDERED_CLUSTER/Chimp_Caudate_Seurat_CB_Raw.RDS')

allseur = subset(allseur, subset = orig.ident %in% c('Sample_2'), invert = T)


# Add meta data
seurMeta = allseur@meta.data
meta = rio::import('~/workdir/TABLES/Primate_Tables_From_EmilyOh/Chimp Macaque Language Summary.xlsx')
rownames(meta) = paste0('Sample_', meta[,1])
colnames(meta) = gsub(' ', '_', colnames(meta)) %>% gsub('\\(|\\)', '', .)
seurMeta = cbind(seurMeta, meta[match(seurMeta$orig.ident, rownames(meta)), c('Brain_Region', 'Age', 'Sex')])
allseur@meta.data = seurMeta

####
## KEEP THE CELLS RETAINED BY CELLRANGER
####

# I normally do not implement this filter, but CellBender did not filter anything in this dataset. Probably because the sequencing depth is very high. 

#samps = unique(allseur$orig.ident)
#keepcells = list()
#for(i in 1:length(samps)){

	#  Read CellBender output
#	mat = Read10X_h5(paste0('~/project/01_CELLRANGER_OUT/CHIMP_OUR/celllranger_count_', samps[i], '/outs/filtered_feature_bc_matrix.h5'))

#	colnames(mat) = gsub('-1', '', colnames(mat)) %>% paste0(., '-', samps[i])
#	keepcells[[i]] = intersect(colnames(mat), colnames(allseur))
#}

#keepcellsV = Reduce(union, keepcells)
#allseur_filt = subset(allseur, cells = keepcellsV)

####
## SET VARIABLES
####

pref = 'Chimp_Caudate_Our'

markerGenes = c('NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA')

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
save_rds(seurM_corr, paste0('~/workdir/01_CHIMP/02_CELLBENDERED_CLUSTER/', pref, '_Clustering1.RDS'))
seurM_corr = read_rds(paste0('~/workdir/01_CHIMP/02_CELLBENDERED_CLUSTER/', pref, '_Clustering1.RDS'))

####
## CLUSTERING 2
####

vasculature = c(11,13)
empty_drop = c()
doublet = c(9)
seurM = subset(seurM_corr, subset = seurat_clusters %in% c(vasculature, empty_drop, doublet), invert = T)
seurM_corr = clusterBatchCorrect(seurM, npcs = 30, res = 0.5, pref = 'Chimp_Caudate', sampleName = 'orig.ident')

# QC PLOTS
pref = 'Chimp_Caudate'

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
dir.create('~/workdir/01_CHIMP/02_CELLBENDERED_CLUSTER/CLUSTERING_2')
save_rds(seurM_corr, paste0('~/workdir/01_CHIMP/02_CELLBENDERED_CLUSTER/', pref, '_Clustering2.RDS'))
seurM_corr = read_rds('~/workdir/01_CHIMP/02_CELLBENDERED_CLUSTER/Chimp_Caudate_Clustering2.RDS')

qw = FindMarkers(seurM_corr, ident.1 = 20, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)


####
## CLUSTERING 3
####

# Also remove sample 7 since it differs from other samples. It has low number of SPNs for some reason.
vasculature = c()
empty_drop = c()
doublet = c(10,20,21)
seurM_corr = subset(seurM_corr, subset = seurat_clusters %in% c(vasculature, empty_drop, doublet), invert = T)
seurM = subset(seurM_corr, subset = orig.ident %in% c('Sample_7'), invert = T)

seurM_corr = clusterBatchCorrect(seurM, npcs = 30, res = 0.5, pref = 'Chimp_Caudate', sampleName = 'orig.ident')

DefaultAssay(seurM_corr) = 'SCT'
seurM_corr = FindClusters(seurM_corr, verbose = FALSE, resolution = 0.1)
stackedbarplot(seurM_corr[[]], 'seurat_clusters', 'orig.ident', paste0(pref, '_AfterBatchCorrect_Stack'), wd = 10, hg = 10)

# QC PLOTS
pref = 'Chimp_Caudate'

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
dir.create('~/workdir/01_CHIMP/02_CELLBENDERED_CLUSTER/CLUSTERING_3')
save_rds(seurM_corr, paste0('~/workdir/01_CHIMP/02_CELLBENDERED_CLUSTER/', pref, '_Clustering3.RDS'))

seurM_corr = read_rds('~/workdir/01_CHIMP/02_CELLBENDERED_CLUSTER/Chimp_Caudate_Clustering3.RDS')

####
## BROADLY ANNOTATE
####

seurM_corr = read_rds(paste0('~/workdir/CHIMP/02_CELLBENDERED_CLUSTER/', pref, '_Clustering3.RDS'))

mapnames = setNames(c('MOL', 'SPN', 'Astrocyte', 'Microglia', 'OPC', 'MOL',
			'SPN', 'Other_GABA', 'SPN', 'MOL', 'Other_GABA',
			'Microglia', 'Other_GABA', 'Other_GABA'),
		      c(0:13))

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

# Save
stackedbarplot(seurM_corr[[]], groupfill = 'newannot', groupx = 'orig.ident', fn = 'Chimp_Stacked_Caudate_Annotated')
save_rds(seurM_corr, paste0('~/workdir/01_CHIMP/02_CELLBENDERED_CLUSTER/', pref, '_Annotated_FINAL.RDS'))


####
## Update matrix if the initial filtering was done incorrectly
####

allseur_orig = read_rds('~/workdir/01_CHIMP/02_CELLBENDERED_CLUSTER/Chimp_Caudate_Annotated_FINAL.RDS')
meta_orig = allseur_orig@meta.data

# Only the caudate samples
allsamps = rio::import('~/workdir/TABLES/Primate_Tables_From_EmilyOh/Chimp Macaque Language Summary.xlsx')
colnames(allsamps) = gsub(' ', '_', colnames(allsamps))
allsamps$Species = gsub(' ', '_', allsamps$Species)
caudSamps = allsamps[allsamps$Brain_Region == 'Caudate' & allsamps$Species == 'Chimp',]
rownames(caudSamps) = paste0('Sample_', caudSamps[,1])

# Set the samples and read them with QC
samps = rownames(caudSamps)
fls = list.files('~/project/03_MATRIX_FROM_CELLBENDER/CHIMP_OUR', pattern = "filtered.h5", recursive=T, full.names = T)
fls = fls[grepl('Chimp_Sample', fls)]

allSeurL = list()
for(i in 1:length(fls)){

	# Sample name
	samp = gsub('.*Chimp_', '', fls[i]) %>% gsub('\\/.*', '', .)

	#  Read CellBender output
	mat = Read10X_h5(fls[i])

	# Soft filtering. Remove lower than 200 UMI
	mat = mat[, colSums(mat) > 200]
	cbBarcode = colnames(mat)
	colnames(mat) = gsub('-1', '', colnames(mat)) %>% paste0(., '-', samp)

	# Create seurat object
	tmpSeur = CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0)
	tmpSeur$orig.ident = samp
	tmpSeur$tissue = 'caudate'
	tmpSeur$Species = 'chimp'
	allSeurL[[i]] = tmpSeur

	print(i)
}

allseur = Reduce(merge, allSeurL)
allseur_filt = subset(allseur, cells = rownames(meta_orig))

mat = allseur_filt@assays$RNA@counts
newSeur = CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0, meta.data = meta_orig)

save_rds(newSeur, '~/workdir/01_CHIMP/02_CELLBENDERED_CLUSTER/Chimp_Caudate_Annotated_FINAL_ALLGENES.RDS')

