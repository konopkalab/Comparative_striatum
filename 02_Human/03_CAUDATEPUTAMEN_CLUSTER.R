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

#####
## LOAD COUNT MATRIX
#####

# Only the putamen samples
allsamps = rio::import('~/workdir/TABLES/Human_Tables/Language_Areas_Human_Tissue_Bank.xlsx')
colnames(allsamps) = gsub(' ', '_', colnames(allsamps))
allsamps$Brain_Region = gsub(' ', '_', allsamps$Brain_Region)
allsamps$Brain_Region = gsub('Caudate_Nucleus', 'Caudate', allsamps$Brain_Region)
cauputSamps = allsamps[allsamps$Brain_Region %in% c('Caudate', 'Putamen'),]
rownames(cauputSamps) = paste0('Sample_', cauputSamps[,1])

# Remove low quality samples
cauputSamps = cauputSamps[cauputSamps$Barcode != '242943',]

# Set the samples and read them with QC
samps = rownames(cauputSamps)

# Set variables
path_cr = '/home2/gkonop/project/00_FASTQ/HUMAN/COUNT/'
path_cb = '/home2/gkonop/project/03_MATRIX_FROM_CELLBENDER/HUMAN_OUR/'
path_gtf = '/home2/gkonop/workdir/REFERENCES/cellranger_reference_ATAC/Hsa_GRCh38/HomSap_GRCh38/genes/genes.gtf'
ncores = 23

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
		nf1 = nuclear_fraction_annotation(
			bam = paste0(path_cr, samps[i], '/outs/possorted_genome_bam.bam'),
			annotation_path = path_gtf,
			barcodes = cbBarcode,
			cores = ncores, verbose = T)

		saveRDS(nf1, paste0(path_cr, samps[i], '/outs/Intronic_read_ratio.RDS'))
	}

	# Keep the final barcodes
	nf1 = nf1[match(cbBarcode, rownames(nf1)),]

	tmpSeur$intronRat = nf1
	tmpSeur$orig.ident = samps[i]
	tmpSeur$Tissue = cauputSamps[samps[i], 'Brain_Region']
	tmpSeur$Species = 'Human'
	allSeurL[[i]] = tmpSeur

	print(i)
}

allseur = Reduce(merge, allSeurL)
write_rds(allseur, '~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/Human_CaudatePutamen_Raw.RDS')

# Add meta data
seurMeta = allseur@meta.data
meta = rio::import('~/workdir/TABLES/Human_Tables/Language_Areas_Human_Tissue_Bank.xlsx')
rownames(meta) = paste0('Sample_', meta[,1])
colnames(meta) = gsub(' ', '_', colnames(meta)) %>% gsub('_Label|_\\(min\\)', '', .)
meta$Brain_Region = gsub('Caudate Nucleus', 'Caudate', meta$Brain_Region)
seurMeta = cbind(seurMeta, meta[match(seurMeta$orig.ident, rownames(meta)), c('Brain_Region', 'Age', 'PMI', 'Sex', 'Race')])
seurMeta$Brain_Region = gsub(' ', '_', seurMeta$Brain_Region) %>% gsub('\\(|\\)', '', .)
allseur@meta.data = seurMeta

####
## SET VARIABLES FOR CLUSTERING
####

pref = 'Human_CaudatePutamen'

markerGenes = c('NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CFTR', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA', 'SATB2', 'SLC17A7')

####
## CLUSTER, REMOVE DOUBLETS, SPLIT NEURON-GLIA
####

seurM = allseur
seurM_corr = clusterBatchCorrect(seurM, npcs = 30, res = 0.5, pref = pref, sampleName = 'orig.ident')

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
todir = '~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN/CLUSTERING_1/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)

# Save (optional)
write_rds(seurM_corr, '~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/Human_CaudatePutamen_Clustering1.RDS')
seurM_corr = read_rds('~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/Human_CaudatePutamen_Clustering1.RDS')

qw = FindMarkers(seurM_corr, ident.1 = 21, only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
pdf('Cluster_4_Markers.pdf', width = 12)
DotPlot(seurM_corr, features = rownames(qw)[1:30]) + rotate_x_text(45)
dev.off()


####
## CLUSTERING 2
####

vasculature = c(16,20,21,22)
empty_drop = c(11,14)
doublet = c(18,21,23,26)
others = c()
seurM = subset(seurM_corr, subset = seurat_clusters %in% c(vasculature, empty_drop, doublet, others), invert = T)
seurM_corr = clusterBatchCorrect(seurM, npcs = 30, res = 0.5, pref = pref, sampleName = 'orig.ident')

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
todir = '~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN/CLUSTERING_2/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)

# Save (optional)
write_rds(seurM_corr, '~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/Human_CaudatePutamen_Clustering2.RDS')
seurM_corr = read_rds('~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/Human_CaudatePutamen_Clustering2.RDS')

####
## SUBSET AND BROADLY ANNOTATE
####

vasculature = c()
empty_drop = c()
doublet = c(19,20,27)
others = c()
seurM = subset(seurM_corr, subset = seurat_clusters %in% c(vasculature, empty_drop, doublet, others), invert = T)

# Find strong markers for major cell types
qw = FindMarkers(seurM, ident.1 = 6, only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
pdf('Cluster_6_Markers_Microglia.pdf', width = 12)
DotPlot(seurM, features = rownames(qw)[1:30]) + rotate_x_text(45)
dev.off()

qw = FindMarkers(seurM, ident.1 = 9, only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
pdf('Cluster_9_Markers_OPC.pdf', width = 12)
DotPlot(seurM, features = rownames(qw)[1:50]) + rotate_x_text(45)
dev.off()

qw = FindMarkers(seurM, ident.1 = c(0,2:5), only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
pdf('Cluster_02345_Markers_MOL.pdf', width = 14)
DotPlot(seurM, features = rownames(qw)[1:50]) + rotate_x_text(45)
dev.off()

qw = FindMarkers(seurM, ident.1 = c(1,8,12,13), only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
pdf('Cluster_Markers_Neurons.pdf', width = 14)
DotPlot(seurM, features = rownames(qw)[1:50]) + rotate_x_text(45)
dev.off()


# STRIATUM BROAD ANNOTATER
meta = seurM@meta.data

mic = c('APBB1IP', 'INPP5D', 'SP100', 'PTPRC', 'ARHGAP15')
opc = c('PCDH15', 'PTPRZ1', 'SOX6', 'COL11A1', 'MIR3681HG')
mol = c('MOG', 'PLP1', 'MBP', 'MOBP', 'RNF220')
ast = c('SLC1A2', 'SLC1A3', 'AQP4')
neu = c('SYT1', 'SNAP25', 'CACNA1E', 'DLGAP3', 'FNBP1L')

# Statistics on markers
marks = FindAllMarkers(seurM, features = c(mic, opc, mol, ast, neu), only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)

markL = list(mic, opc, mol, ast, neu)
names(markL) = c('MIC', 'OPC', 'MOL', 'AST', 'NEU')

# Cell type score of each cluster
dfL = list()
for(i in 1:length(markL)){

	df = marks[marks$gene %in% markL[[i]], 'cluster'] %>% table %>% as.data.frame
	colnames(df) = c('Cluster', 'Score')
	df$Score = df$Score / length(markL[[i]])

	lfc_score = marks[marks$gene %in% markL[[i]], ] %>% group_by(cluster, .drop = F) %>% dplyr::summarize(meanLFC = mean(avg_log2FC)) %>% as.data.frame
	lfc_score[is.nan(lfc_score$meanLFC), 'meanLFC'] = 0
	df$meanLFC = lfc_score$meanLFC
	df$Score = df$Score * df$meanLF
	df$CellType = names(markL)[i]
	dfL[[i]] = df
}

finaldf = do.call(rbind, dfL)

# Match cluster with annotation
cls = finaldf$Cluster %>% unique
annotL = list()
for(i in 1:length(cls)){
	subdf = finaldf[finaldf$Cluster == cls[i],]
	max_ind = subdf$Score %>% which.max
	ctype = subdf[max_ind, 'CellType']
	df2 = data.frame(Cluster = cls[i], CellType = ctype)
	annotL[[i]] = df2
}

annots = do.call(rbind, annotL)

# Some clusters are misclassified
annots[annots$Cluster == 22, 'CellType'] = 'NEU'
annots[annots$Cluster == 17, 'CellType'] = 'NEU'
annots[annots$Cluster == 21, 'CellType'] = 'MIC'

# Annotate
seurM$seurat_clusters = droplevels(seurM$seurat_clusters)
mapnames = setNames(annots$CellType, annots$Cluster)
seurM[["newannot"]] = mapnames[seurM[["seurat_clusters"]][,1]]
Idents(seurM) = seurM$newannot

# PLOTS
seurM = RunUMAP(seurM, dims = 1:30, reduction = 'harmony')

pdf(paste0(pref, "_Clusters_Annotated.pdf"))
DimPlot(seurM, group.by = 'newannot', label = T, raster = T, label.size = 7) +
NoLegend() +
ggtitle('')
dev.off()

pdf(paste0(pref, "_Clusters.pdf"))
DimPlot(seurM, group.by = 'seurat_clusters', label = T, raster = T, label.size = 7) +
NoLegend() +
ggtitle('')
dev.off()

# Save
write_rds(seurM, '~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/Human_CaudatePutamen_Annotated.RDS')

# Move the plots
outfls = list.files('~/', pattern = pref, full.names = F)
todir = '~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN/ANNOTATED/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)






