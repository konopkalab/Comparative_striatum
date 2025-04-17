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

# Read annotated caudate and putamen data
caudSeur = read_rds('~/workdir/01_CHIMP/02_CELLBENDERED_CLUSTER/Chimp_Caudate_Annotated_FINAL.RDS')
putSeur = read_rds('~/workdir/01_CHIMP/02_CELLBENDERED_CLUSTER/Chimp_Putamen_Our_Annotated_FINAL.RDS')

# Subset neurons
caudSeur = subset(caudSeur, subset = newannot %in% c('SPN', 'Other_GABA'))
putSeur = subset(putSeur, subset = newannot %in% c('SPN', 'Other_GABA'))
seurM = merge(caudSeur, putSeur)
seurM$Tissue = seurM$tissue
seurM$tissue = NULL
rm(caudSeur, putSeur)
gc()

####
## SET VARIABLES FOR CLUSTERING
####

pref = 'Chimp_CaudatePutamen_NEURONS'

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
todir = '~/workdir/01_CHIMP/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN_NEURONS/CLUSTERING_1/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)

# Save (optional)
write_rds(seurM_corr, '~/workdir/01_CHIMP/03_CLUSTER_ANNOTATE/Chimp_CaudatePutamen_NEURONS_Clustering1.RDS')
seurM_corr = read_rds('~/workdir/01_CHIMP/03_CLUSTER_ANNOTATE/Chimp_CaudatePutamen_NEURONS_Clustering1.RDS')

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
#annots[annots$Cluster == 22, 'CellType'] = 'Not_Doublet'


####
## CLUSTER-2
####

# Remove doublets and recluster
doublets = annots[annots$CellType != 'Not_Doublet', 'Cluster'] %>% as.character %>% as.numeric
empty_droplets = c(19)

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
todir = '~/workdir/01_CHIMP/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN_NEURONS/CLUSTERING_2/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)

# Save (optional)
write_rds(seurM_corr, '~/workdir/01_CHIMP/03_CLUSTER_ANNOTATE/Chimp_CaudatePutamen_NEURONS_Clustering2.RDS')
seurM_corr = read_rds('~/workdir/01_CHIMP/03_CLUSTER_ANNOTATE/Chimp_CaudatePutamen_NEURONS_Clustering2.RDS')


qw = FindMarkers(seurM_corr, ident.1 = 2, only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
pdf('Cluster_4_Markers.pdf', width = 12)
DotPlot(seurM_corr, features = rownames(qw)[1:30]) + rotate_x_text(45)
dev.off()

####
## CLUSTER-3
####

# Remove doublets and recluster
doublets = c(2,9,14)
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
todir = '~/workdir/01_CHIMP/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN_NEURONS/CLUSTERING_3/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)

# Save (optional)
write_rds(seurM_corr, '~/workdir/01_CHIMP/03_CLUSTER_ANNOTATE/Chimp_CaudatePutamen_NEURONS_Clustering3.RDS')
seurM_corr = read_rds('~/workdir/01_CHIMP/03_CLUSTER_ANNOTATE/Chimp_CaudatePutamen_NEURONS_Clustering3.RDS')

qw = FindMarkers(seurM_corr, ident.1 = 17, only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
pdf('Cluster_4_Markers.pdf', width = 12)
DotPlot(seurM_corr, features = rownames(qw)[1:30]) + rotate_x_text(45)
dev.off()

####
## CLUSTER-4
####

# Remove doublets and recluster
doublets = c(12,17)
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
todir = '~/workdir/01_CHIMP/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN_NEURONS/CLUSTERING_4/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)

# Save (optional)
write_rds(seurM_corr, '~/workdir/01_CHIMP/03_CLUSTER_ANNOTATE/Chimp_CaudatePutamen_NEURONS_Clustering4.RDS')
seurM_corr = read_rds('~/workdir/01_CHIMP/03_CLUSTER_ANNOTATE/Chimp_CaudatePutamen_NEURONS_Clustering4.RDS')

#qw = FindMarkers(seurM_corr, ident.1 = 9, only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
#pdf('Cluster_4_Markers.pdf', width = 12)
#DotPlot(seurM_corr, features = rownames(qw)[1:30]) + rotate_x_text(45)
#dev.off()

####
## CLUSTER-5
####

# Remove doublets and recluster. Empty droplet here is just somewhat lower quality with high MT content
doublets = c(14)
empty_droplets = c(12)

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
todir = '~/workdir/01_CHIMP/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN_NEURONS/CLUSTERING_5/'
dir.create(todir)
sapply(outfls, function(x){file.copy(from = paste0('~/', x), to = paste0(todir, x))})
file.remove(outfls)

# Save (optional)
write_rds(seurM_corr, '~/workdir/01_CHIMP/03_CLUSTER_ANNOTATE/Chimp_CaudatePutamen_NEURONS_Clustering5.RDS')
seurM_corr = read_rds('~/workdir/01_CHIMP/03_CLUSTER_ANNOTATE/Chimp_CaudatePutamen_NEURONS_Clustering5.RDS')

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

####
## LABEL TRANSFER USING HUMAN AS REFERENCE
####

# Load datasets
humanSeur = read_rds('~/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/Human_CaudatePutamen_NEURONS_Annotated.RDS')
chimpSeur = read_rds('~/workdir/01_CHIMP/03_CLUSTER_ANNOTATE/Chimp_CaudatePutamen_NEURONS_Clustering5.RDS')

# Remove a doublet cluster
chimpSeur = subset(chimpSeur, subset = seurat_clusters == 12, invert = T)

# Merge clusters from the same cell type
chimpSeur$seurat_clusters = chimpSeur$seurat_clusters %>% as.character
chimpSeur$seurat_clusters = gsub('^2$', '0', chimpSeur$seurat_clusters)
chimpSeur$seurat_clusters = gsub('^3$', '1', chimpSeur$seurat_clusters)

# Update some annotations in human
#humanSeur$newannot = gsub('Exc_.*', 'Excitatory', humanSeur$newannot)

# Check how many genes will be used for label transfer
DefaultAssay(humanSeur) = 'SCT'
DefaultAssay(chimpSeur) = 'SCT'
sum(VariableFeatures(chimpSeur) %in% rownames(humanSeur@assays$SCT))

# Find transfer anchors
transfer_anchors = FindTransferAnchors(reference = humanSeur, query = chimpSeur,
					dims = 1:20,
					reference.assay = "SCT",
					query.assay = "SCT",
					reduction = "cca")

write_rds(transfer_anchors, '~/workdir/01_CHIMP/03_CLUSTER_ANNOTATE/Chimp_CaudatePutamen_NEURONS_TransferAnchors_FromHuman.RDS')

# CONTINUE FROM HERE OCT 21, 2023
transfer_anchors = read_rds('~/workdir/01_CHIMP/03_CLUSTER_ANNOTATE/Chimp_CaudatePutamen_NEURONS_TransferAnchors_FromHuman.RDS')


# Make cell type predictions based on anchors
pred = TransferData(anchorset = transfer_anchors,
			refdata = humanSeur$newannot, 
			weight.reduction = chimpSeur[["harmony"]],
			dims = 1:20)

# Add predictions to meta data
chimpSeur$pred = pred$predicted.id
chimpSeur$pred_score = pred$prediction.score.max

pdf("MaxPredictionScores_ToCHIMP_FromHUMAN.pdf")
hist(chimpSeur$pred_score)
dev.off()

# Plot label transfer
# Make the colors match
#newrna_int$pred <- factor(newrna_int$pred, levels = levels(allen_mtg_seur$cell_type))
Idents(chimpSeur) = chimpSeur$pred
Idents(humanSeur) = humanSeur$newannot
transferPlot(chimpSeur, humanSeur, "LABEL_TRANSFER_ToCHIMP_FromHUMAN.pdf")

# Examine label transfer per cluster. Increase resolution if cluster has two distinct labels.
chimpSeur$predicted.id = chimpSeur$pred
ltpred(chimpSeur, vars = c('predicted.id', 'seurat_clusters'), fn = "LABEL_TRANSFER_CLARITY_ToCHIMP_FromHUMAN.pdf")

# Jaccard similarity between prediction and annotation
meta = chimpSeur@meta.data
meta$newannot = meta$seurat_clusters
ctypes_target = names(table(meta$newannot))
ctypes_ref = names(table(humanSeur$newannot))
jacmatL = list()
for(i in 1:length(ctypes_target)){

	annotCells = meta[meta$newannot %in% ctypes_target[i],] %>% rownames
	jacL = list()

	for(j in 1:length(ctypes_ref)){
		predCells = meta[meta$pred %in% ctypes_ref[j],] %>% rownames
		jacL[[j]] = sum(annotCells %in% predCells) / length(union(annotCells,predCells))
		#print(jacL[[j]])
	}
	jacmatL[[i]] = unlist(jacL)
	names(jacmatL[[i]]) = ctypes_ref
}

jacmatDF = do.call(rbind, jacmatL)
rownames(jacmatDF) = ctypes_target
toplot = melt(jacmatDF)
colnames(toplot) = c('Annotation', 'Prediction', 'Jaccard_Index')

# Order prediction labels
levs = colnames(jacmatDF)[apply(jacmatDF, 1, function(x){which.max(x)}) %>% as.numeric] %>% rev %>% unique
levs = c(levs, unique(setdiff(toplot$Prediction, levs)))
toplot$Prediction = factor(toplot$Prediction, levels = levs)

# Find the best matches per annotation group
toplot$Annotation = factor(toplot$Annotation, levels = c(0,1,4:11,13,14,15))
toplot$Prediction = factor(toplot$Prediction, levels = table(humanSeur$newannot) %>% sort %>% names %>% rev)

annot_ctypes = levels(toplot$Annotation)
bestL = list()
for(i in 1:length(annot_ctypes)){
	sub = toplot[toplot$Annotation == annot_ctypes[i], ]
	ind = sub$Jaccard_Index %>% which.max
	if(sum(sub[ind, 'Jaccard_Index'] > 0.25) > 0){
		bestL[[i]] = paste0(sub[ind, 'Annotation'], '_', sub[ind, 'Prediction'])
	}
}
bests = unlist(bestL)
toplot$id = paste0(toplot$Annotation, '_', toplot$Prediction)

# Assign best prediction
toplot$is_best = ifelse(toplot$id %in% bests, '*', '')

pdf('PredictionPlot_JACCARD_ToCHIMP_FromHUMAN.pdf', width = 10, height =10)
ggscatter(toplot, x = 'Prediction', y = 'Annotation', color = 'Jaccard_Index', size = 10) +
scale_color_gradient2(low = "white", high = "blue", mid = "white", midpoint = 0, limit = c(0,1)) +
xlab('Prediction') + ylab('Annotation') +
theme(axis.text.x = element_text(size=20),
	axis.text.y = element_text(size=20),
	axis.title = element_text(size=20)) +
	rotate_x_text(90) + grids(linetype = "dashed", color = 'grey') +
geom_text(aes(label = is_best), vjust = 0.8, colour = "red", size = 15 )
dev.off()


# Annotate and replot
mapnames = setNames(toplot[toplot$is_best == '*', 'Prediction'] %>% as.character,
		      toplot[toplot$is_best == '*', 'Annotation'] %>% as.character)

toplot$Chimpanzee = toplot$Annotation
toplot$Human = toplot$Prediction
toplot[["Chimpanzee"]] = mapnames[toplot[["Annotation"]] %>% as.character]
toplot[toplot$Annotation == 8, 'Chimpanzee'] = 'SPN_IEG'

toplot$Human = factor(toplot$Human, levels = c('dSPN', 'iSPN', 'Excitatory', 'eSPN_FOXP2', 'PVALB_PDGFD', 'TAC3', 'PVALB_FLT3', 'eSPN_DRD2', 'SST_NPY', 'VIP', 'CCK_EYA4', 'CHAT'))
toplot$Chimpanzee = factor(toplot$Chimpanzee, levels = c('dSPN', 'iSPN', 'SPN_IEG', 'Excitatory', 'eSPN_FOXP2', 'PVALB_PDGFD', 'TAC3', 'PVALB_FLT3', 'eSPN_DRD2', 'SST_NPY', 'VIP', 'CCK_EYA4', 'CHAT'))

pdf('PredictionPlot_JACCARD_ToCHIMP_FromHUMAN_FINAL.pdf', width = 10, height =10)
ggscatter(toplot, x = 'Chimpanzee', y = 'Human', color = 'Jaccard_Index', size = 10) +
scale_color_gradient2(low = "white", high = "blue", mid = "white", midpoint = 0, limit = c(0,1)) +
xlab('Chimpanzee') + ylab('Human') +
theme(axis.text.x = element_text(size=20),
	axis.text.y = element_text(size=20),
	axis.title = element_text(size=20)) +
	rotate_x_text(90) + grids(linetype = "dashed", color = 'grey') +
geom_text(aes(label = is_best), vjust = 0.8, colour = "red", size = 15 )
dev.off()



mapnames = setNames(toplot$Chimpanzee %>% as.character,
		      toplot$Annotation %>% as.character)

chimpSeur[["newannot"]] = mapnames[chimpSeur[["seurat_clusters"]][,1]]
Idents(chimpSeur) = chimpSeur$newannot

# PLOTS
pdf(paste0(pref, "_Clusters_Annotated.pdf"))
DimPlot(chimpSeur, group.by = 'newannot', label = T, raster = T, repel = T, label.size = 7) +
NoLegend() +
ggtitle('')
dev.off()

pdf(paste0(pref, "_Depth_Gene_Annotated.pdf"))
ggboxplot(chimpSeur[[]], x = 'newannot', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_Depth_UMI_Annotated.pdf"))
ggboxplot(chimpSeur[[]], x = 'newannot', y = 'nCount_RNA', ylim = c(0,80000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_NuclearFraction_Annotated.pdf"))
ggboxplot(chimpSeur[[]], x = 'newannot', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()



marksToPlot = c('FOSL2', 'NR4A1', 'NPAS4', 'TAC3', 'TH', 'TRH', 'TRHDE', 'SOX6', 'PTK2B', 'PRDM1', 'EYA4', 'FRMD7', 'SHD', 'CCK', 'PVALB', 'PTHLH', 'PDGFD', 'EYA2', 'NTNG1', 'KLHL1', 'NDNF', 'SST', 'NPY', 'VIP', 'CHAT', 'SULF1', 'ANGPT1', 'RBFOX3', 'SYT1', 'SNAP25', 'CASZ1', 'OTOF', 'CACNG5', 'PCDH8', 'FOXP2', 'DRD1', 'DRD2', 'TAC1', 'PENK', 'MOG', 'MBP', 'PCDH15', 'NRGN', 'GAD1', 'SLC1A2', 'AQP4', 'APBB1IP', 'FLT1', 'SLC17A7', 'SATB2', 'MALAT1')

DefaultAssay(chimpSeur) = 'RNA'
chimpSeur = NormalizeData(chimpSeur)
pdf(paste0(pref, "_MarkersDotPlot.pdf"), width = 25, height = 10)
DotPlot(chimpSeur, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

stackedbarplot(chimpSeur[[]], groupfill = 'newannot', groupx = 'orig.ident', fn = paste0(pref, '_Stacked_Annotated'))

stackedbarplot(chimpSeur[[]], groupfill = 'orig.ident', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_2'))

stackedbarplot(chimpSeur[[]], groupfill = 'Tissue', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Tissue'))

write_rds(chimpSeur, '~/workdir/01_CHIMP/03_CLUSTER_ANNOTATE/Chimp_CaudatePutamen_NEURONS_Annotated.RDS')
