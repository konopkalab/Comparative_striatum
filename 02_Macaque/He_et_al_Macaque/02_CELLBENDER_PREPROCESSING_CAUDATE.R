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

#####
## LOAD COUNT MATRIX
#####

# Only the caudate samples
srrMeta = read.table('~/workdir/MACAQUE/HE/Macaque_He_SRA_Metadata.txt', header = T)
srrMeta$Age = as.numeric(gsub('.*,', '', srrMeta[,1]))
srrMeta$Sample =  srrMeta[, 1] %>% gsub(',.*', '', .)
saveRDS(srrMeta, '~/workdir/PRIMATE_CAUDATE/MACAQUE_HE/Macaque_He_SRA_Metadata_FINAL.RDS')

samps = srrMeta[grepl('Caudate', srrMeta$Study.TISSUE), 'Sample']
allSeurL = list()
for(i in 1:length(samps)){

	#  Read CellBender output
	mat = Read10X_h5(paste0('~/project/MATRIX_FROM_CELLBENDER/MACAQUE_HE/', samps[i], '/CellBender_out_filtered.h5'))

	# Add sample name to the cell barcode name
	cbBarcode = colnames(mat)
	colnames(mat) = gsub('-1', '', colnames(mat)) %>% paste0(., '-', samps[i])

	# Create seurat object
	tmpSeur = CreateSeuratObject(counts = mat, min.cells = 10, min.features = 0)

	# Calculate intronic read ratio
	# Set TMPDIR='~/project' before running R.
	nf2 = nuclear_fraction_annotation(
		annotation_path = '/home2/gkonop/workdir/REFERENCES/cellranger_reference_RNA/rheMac10_human_liftoff/reference_genome/genes/genes.gtf.gz',
		bam = paste0('~/project/CELLRANGER_OUT/MACAQUE_HE2021/celllranger_count_', samps[i], '/outs/possorted_genome_bam.bam'),
		barcodes = cbBarcode,
		cores = 23, verbose = T)

	nf2$barc = colnames(tmpSeur)
	tmpSeur$intronRat = nf2$nuclear_fraction
	tmpSeur$orig.ident = samps[i]
	tmpSeur$tissue = 'caudate'
	tmpSeur$Species = 'macaque'
	allSeurL[[i]] = tmpSeur

	print(i)
}

allseur = Reduce(merge, allSeurL)

saveRDS(allseur, '~/workdir/MACAQUE/HE/03_CLUSTER_ANNOTATE/Macaque_Caudate_Seurat_CB_Raw.RDS')
allseur = readRDS('~/workdir/MACAQUE/HE/03_CLUSTER_ANNOTATE/Macaque_Caudate_Seurat_CB_Raw.RDS')

# Add meta data
seurMeta = allseur@meta.data
meta = readRDS('~/workdir/MACAQUE/HE/Macaque_He_SRA_Metadata_FINAL.RDS')
seurMeta$Age = meta[match(seurMeta$orig.ident, meta$Sample), 'Age']
allseur@meta.data = seurMeta

####
## KEEP THE CELLS RETAINED BY CELLRANGER
####

# I normally do not implement this filter, but CellBender did not filter anything in this dataset. Probably because the sequencing depth is very high. 

samps = unique(allseur$orig.ident)
keepcells = list()
for(i in 1:length(samps)){

	#  Read CellBender output
	mat = Read10X_h5(paste0('~/project/01_CELLRANGER_OUT/MACAQUE_HE2021/celllranger_count_', samps[i], '/outs/filtered_feature_bc_matrix.h5'))

	colnames(mat) = gsub('-1', '', colnames(mat)) %>% paste0(., '-', samps[i])
	keepcells[[i]] = intersect(colnames(mat), colnames(allseur))
}

keepcellsV = Reduce(union, keepcells)
allseur_filt = subset(allseur, cells = keepcellsV)

####
## SET VARIABLES
####

pref = 'Macaque_Caudate_He'

markerGenes = c('NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA')


####
## CLUSTERING 1
####

seurM = allseur_filt
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

seurM_corr$log10UMI = log10(seurM_corr$nCount_RNA)
pdf(paste0(pref, "_Depth_UMI_log10.pdf"))
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'log10UMI') +
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
pdf(paste0(pref, "_MarkersViolin.pdf"), width = 20, height = 15)
VlnPlot(seurM_corr, features = markerGenes, pt.size = 0, raster = T)
dev.off()

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

#qw = FindMarkers(seurM_corr, ident.1 = 12, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)

# Save (optional)
#dir.create('~/workdir/MACAQUE/HE/03_CLUSTER_ANNOTATE/CAUDATE_ONLY')
#dir.create('~/workdir/MACAQUE/HE/03_CLUSTER_ANNOTATE/CAUDATE_ONLY/CLUSTERING_1')
saveRDS(seurM_corr, paste0('~/workdir/MACAQUE/HE/03_CLUSTER_ANNOTATE/CAUDATE_ONLY/CLUSTERING_1/', pref, '_Clustering1.RDS'))


####
## CLUSTERING 2
####

vasculature = c(9,11)
empty_drop = c(6,14)
doublet = c(4,5,16)
seurM = subset(seurM_corr, subset = seurat_clusters %in% c(vasculature, empty_drop, doublet), invert = T)
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

seurM_corr$log10UMI = log10(seurM_corr$nCount_RNA)
pdf(paste0(pref, "_Depth_UMI_log10.pdf"))
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'log10UMI') +
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
pdf(paste0(pref, "_MarkersViolin.pdf"), width = 20, height = 15)
VlnPlot(seurM_corr, features = markerGenes, pt.size = 0, raster = T)
dev.off()

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

#qw = FindMarkers(seurM_corr, ident.1 = 12, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)

# Save (optional)
#dir.create('~/workdir/MACAQUE/HE/03_CLUSTER_ANNOTATE/CAUDATE_ONLY')
dir.create('~/workdir/MACAQUE/HE/03_CLUSTER_ANNOTATE/CAUDATE_ONLY/CLUSTERING_2')
saveRDS(seurM_corr, paste0('~/workdir/MACAQUE/HE/03_CLUSTER_ANNOTATE/CAUDATE_ONLY/', pref, '_Clustering2.RDS'))


####
## CLUSTERING 3
####

vasculature = c()
empty_drop = c(10,11)
doublet = c()
seurM = subset(seurM_corr, subset = seurat_clusters %in% c(vasculature, empty_drop, doublet), invert = T)
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

seurM_corr$log10UMI = log10(seurM_corr$nCount_RNA)
pdf(paste0(pref, "_Depth_UMI_log10.pdf"))
ggboxplot(seurM_corr[[]], x = 'seurat_clusters', y = 'log10UMI') +
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
pdf(paste0(pref, "_MarkersViolin.pdf"), width = 20, height = 15)
VlnPlot(seurM_corr, features = markerGenes, pt.size = 0, raster = T)
dev.off()

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

#qw = FindMarkers(seurM_corr, ident.1 = 12, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)

# Save (optional)
#dir.create('~/workdir/MACAQUE/HE/03_CLUSTER_ANNOTATE/CAUDATE_ONLY')
dir.create('~/workdir/MACAQUE/HE/03_CLUSTER_ANNOTATE/CAUDATE_ONLY/CLUSTERING_3')
saveRDS(seurM_corr, paste0('~/workdir/MACAQUE/HE/03_CLUSTER_ANNOTATE/CAUDATE_ONLY/', pref, '_Clustering3.RDS'))



### UNUSED ###

####
## BROADLY ANNOTATE
####

seurM_corr = readRDS(paste0('~/workdir/MACAQUE/HE/03_CLUSTER_ANNOTATE/CAUDATE_ONLY/SEURAT/CLUSTERING_5/', pref, '_Clustering5.RDS'))

mark13 = FindMarkers(seurM_corr, ident.1 = 13, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
mark15 = FindMarkers(seurM_corr, ident.1 = 15, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)

mapnames = setNames(c('Astrocyte', 'SPN', 'SPN', 'MOL', 'MOL', 'OPC',
			'Microglia', 'SPN', 'SPN', 'SPN', 'Astrocyte',
			'Astrocyte', 'Endothelia', 'Ependymal', 'Other_GABA', 'Mural',
			'Other_GABA', 'Other_GABA', 'Other_GABA', 'Other_GABA'),
		      c(0:19))

seurM_corr[["newannot"]] = mapnames[seurM_corr[["seurat_clusters"]][,1]]
Idents(seurM_corr) = seurM_corr$newannot

# PLOTS
pref = 'Macaque_Caudate'

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
pdf(paste0(pref, "_MarkersViolin.pdf"), width = 20, height = 25)
VlnPlot(seurM_corr, features = markerGenes, pt.size = 0, raster = T)
dev.off()

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
dir.create('~/workdir/MACAQUE/HE/03_CLUSTER_ANNOTATE/CAUDATE_ONLY/SEURAT/FINAL')
saveRDS(seurM_corr, paste0('~/workdir/MACAQUE/HE/03_CLUSTER_ANNOTATE/CAUDATE_ONLY/SEURAT/FINAL/', pref, '_Annotated_FINAL.RDS'))
