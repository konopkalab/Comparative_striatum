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

# Subset neurons
seurM = read_rds('~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/Human_CaudatePutamen_Annotated.RDS')
seurM = subset(seurM, subset = newannot == 'NEU')

####
## SET VARIABLES FOR CLUSTERING
####

pref = 'Human_CaudatePutamen_NEURONS'

markerGenes = c('NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CFTR', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA', 'SATB2', 'SLC17A7')

####
## CLUSTER-1
####

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
todir = '~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN_NEURONS/CLUSTERING_1/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)

# Save (optional)
write_rds(seurM_corr, '~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/Human_CaudatePutamen_NEURONS_Clustering1.RDS')
seurM_corr = read_rds('~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/Human_CaudatePutamen_NEURONS_Clustering1.RDS')

qw = FindMarkers(seurM_corr, ident.1 = 21, only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
pdf('Cluster_4_Markers.pdf', width = 12)
DotPlot(seurM_corr, features = rownames(qw)[1:30]) + rotate_x_text(45)
dev.off()



####
## SRIATUM NEURON-GLIA DOUBLET SUSPECTOR IN NEURONS
####

mic = c('APBB1IP', 'INPP5D', 'SP100', 'PTPRC', 'ARHGAP15')
opc = c('PCDH15', 'PTPRZ1', 'SOX6', 'COL11A1', 'MIR3681HG')
mol = c('MOG', 'PLP1', 'MBP', 'MOBP', 'RNF220')
ast = c('SLC1A2', 'SLC1A3', 'AQP4')

# Expressions of markers
Idents(seurM_corr) = seurM_corr$seurat_clusters
dotout = DotPlot(seurM_corr, features = c(mic, opc, mol, ast))
dotout = dotout$data

markL = list(mic, opc, mol, ast)
names(markL) = c('MIC', 'OPC', 'MOL', 'AST')

# Doublet score of each cluster
dfL = list()
for(i in 1:length(markL)){

	df = dotout[dotout$features.plot %in% markL[[i]], ] %>% group_by(id, .drop = F) %>% dplyr::summarize(meanPCT = mean(pct.exp)) %>% as.data.frame
	colnames(df) = c('Cluster', 'Mean_PCT')
	df$CellType = names(markL)[i]
	dfL[[i]] = df
}

finaldf = do.call(rbind, dfL)

# Predict doublets. This is a very crude way of doing this but it is very fast.
cls = finaldf$Cluster %>% unique
annotL = list()
for(i in 1:length(cls)){
	subdf = finaldf[finaldf$Cluster == cls[i],]

	if(sum(subdf$Mean_PCT > 50) == 0){
		ctype = 'Not_Doublet'
		df2 = data.frame(Cluster = cls[i], CellType = ctype)
		annotL[[i]] = df2
	} else{
		max_ind = subdf$Mean_PCT %>% which.max
		ctype = subdf[max_ind, 'CellType']
		df2 = data.frame(Cluster = cls[i], CellType = ctype)
		annotL[[i]] = df2
	}
}

annots = do.call(rbind, annotL)

# Adjust the results if they have distinct markers
annots[annots$Cluster == 22, 'CellType'] = 'Not_Doublet'

# Move the plots
outfls = list.files('~/', pattern = pref, full.names = F)
todir = '~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN_NEURONS/CLUSTERING_1/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)

#qw = FindMarkers(seurM_corr, ident.1 = 7, only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
#pdf('Cluster_4_Markers.pdf', width = 12)
#DotPlot(seurM_corr, features = rownames(qw)[1:30]) + rotate_x_text(45)
#dev.off()


####
## CLUSTER-2
####

# Remove doublets and recluster
doublets = annots[annots$CellType != 'Not_Doublet', 'Cluster'] %>% as.character %>% as.numeric
empty_droplets = c(7)

seurM = subset(seurM_corr, subset = seurat_clusters %in% c(doublets, empty_droplets), invert = T)
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
todir = '~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN_NEURONS/CLUSTERING_2/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)

# Save (optional)
write_rds(seurM_corr, '~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/Human_CaudatePutamen_NEURONS_Clustering2.RDS')
seurM_corr = read_rds('~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/Human_CaudatePutamen_NEURONS_Clustering2.RDS')


qw = FindMarkers(seurM_corr, ident.1 = 14, only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
pdf('Cluster_4_Markers.pdf', width = 12)
DotPlot(seurM_corr, features = rownames(qw)[1:30]) + rotate_x_text(45)
dev.off()


####
## CLUSTER-3
####

# Remove doublets and recluster
doublets = c(2,10)
empty_droplets = c(18)

seurM = subset(seurM_corr, subset = seurat_clusters %in% c(doublets, empty_droplets), invert = T)
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
todir = '~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN_NEURONS/CLUSTERING_3/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)

# Save (optional)
write_rds(seurM_corr, '~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/Human_CaudatePutamen_NEURONS_Clustering3.RDS')
seurM_corr = read_rds('~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/Human_CaudatePutamen_NEURONS_Clustering3.RDS')


####
## CLUSTER-4
####

# Remove doublets and recluster
doublets = c(13)
empty_droplets = c()

seurM = subset(seurM_corr, subset = seurat_clusters %in% c(doublets, empty_droplets), invert = T)
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
todir = '~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN_NEURONS/CLUSTERING_4/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)

# Save (optional)
write_rds(seurM_corr, '~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/Human_CaudatePutamen_NEURONS_Clustering4.RDS')
seurM_corr = read_rds('~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/Human_CaudatePutamen_NEURONS_Clustering4.RDS')



# PLOT NEURONAL CELL TYPE MARKERS (KRIENEN 2020)
marksToPlot = c('TAC3', 'TH', 'TRH', 'TRHDE', 'SOX6', 'PTK2B', 'PRDM1', 'EYA4', 'FRMD7', 'SHD', 'CCK', 'PVALB', 'PTHLH', 'PDGFD', 'EYA2', 'NTNG1', 'KLHL1', 'NDNF', 'SST', 'NPY', 'VIP', 'CHAT', 'SULF1', 'ANGPT1', 'RBFOX3', 'SYT1', 'SNAP25', 'CASZ1', 'OTOF', 'CACNG5', 'PCDH8', 'FOXP2', 'DRD1', 'DRD2', 'TAC1', 'PENK', 'MOG', 'MBP', 'PCDH15', 'NRGN', 'GAD1', 'SLC1A2', 'AQP4', 'APBB1IP', 'FLT1', 'SLC17A7', 'SATB2', 'MALAT1')

pdf(paste0(pref, "_NEURON_MarkersDotPlot.pdf"), width = 20, height = 10)
DotPlot(seurM_corr, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()





# Find top distinct markers
pct_dif_cutoff = 30
the_cls = 11

subseur = subset(seurM_corr, subset = seurat_clusters == the_cls)
Idents(subseur) = subseur$seurat_clusters
mat = subseur@assays$RNA@counts
feats = which( apply(mat, 1, function(x){sum(x>0) / length(x)}) > 0.25 ) %>% names
dotout = DotPlot(seurM_corr, features = feats)
dotout = dotout$data

cls = dotout$id %>% unique %>% as.character %>% as.numeric
cls = cls[cls != the_cls]
pctDiffL = list()
for(i in 1:length(cls)){

	pctdiff = dotout[dotout$id == the_cls, 'pct.exp'] - dotout[dotout$id == cls[i], 'pct.exp']
	pctDiffL[[i]] = data.frame(gns = dotout[dotout$id == the_cls, 'features.plot'], pctdiff = pctdiff)
}
finaldf = do.call(rbind, pctDiffL)

qw = finaldf %>% group_by(gns) %>% summarize( score = sum(pctdiff > pct_dif_cutoff) ) %>% as.data.frame
plotGns = qw[order(qw$score, decreasing = T)[1:30], 'gns'] %>% as.character

pdf(paste0(pref, "_Cluster_", the_cls, "_Top_Distinct_Markers.pdf"), width = 20, height = 10)
DotPlot(seurM_corr, features = plotGns) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()



####
## ANNOTATE
####

mapnames = setNames(c('iSPN', 'dSPN', 'dSPN', 'iSPN', 'eSPN_FOXP2', 'dSPN',
			'iSPN', 'PVALB_PDGFD', 'TAC3', 'Excitatory', 'Excitatory',
			'PVALB_PRDM1', 'eSPN_DRD2', 'SST_NPY', 'CCK_VIP', 'Excitatory',
			'CCK_EYA4', 'CHAT', 'dSPN'),
		      c(0:18))

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


# MARKER GENE PLOTS
marksToPlot = c('CCK', 'EYA4', 'VIP', 'PVALB', 'PDGFD', 'PTHLH', 'PRDM1', 'TAC3', 'TH', 'TRH', 'TRHDE', 'SOX6', 'PTK2B', 'EYA2', 'FRMD7', 'SHD', 'NTNG1', 'KLHL1', 'NDNF', 'SST', 'NPY', 'CHAT', 'SULF1', 'ANGPT1', 'RBFOX3', 'SYT1', 'SNAP25', 'CASZ1', 'OTOF', 'CACNG5', 'PCDH8', 'FOXP2', 'DRD1', 'DRD2', 'TAC1', 'PENK', 'MOG', 'MBP', 'PCDH15', 'NRGN', 'GAD1', 'SLC1A2', 'AQP4', 'APBB1IP', 'FLT1', 'SLC17A7', 'SATB2', 'MALAT1')

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

stackedbarplot(seurM_corr[[]], groupfill = 'newannot', groupx = 'orig.ident', fn = paste0(pref, '_Stacked_Annotated'))

stackedbarplot(seurM_corr[[]], groupfill = 'orig.ident', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_2'))

stackedbarplot(seurM_corr[[]], groupfill = 'Tissue', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Tissue'))

# Move the plots
outfls = list.files('~/', pattern = pref, full.names = F)
todir = '~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN_NEURONS/ANNOTATED/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)

# Save (optional)
write_rds(seurM_corr, '~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/Human_CaudatePutamen_NEURONS_Annotated.RDS')
seurM_corr = read_rds('~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/Human_CaudatePutamen_NEURONS_Annotated.RDS')






