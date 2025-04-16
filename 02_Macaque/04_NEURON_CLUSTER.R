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
caudSeur <- readRDS("/home2/gkonop/workdir/01_MACAQUE/INTEGRATED/MACAQUE_OUR_HE_INTEGRATED_GB_Annotated.RDS")

putSeur = read_rds('/home2/gkonop/workdir/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/MACAQUE_OUR_HE_INTEGRATED_Annotated_Putamen.RDS')
putSeur$tissue = gsub('Putamen', 'putamen', putSeur$tissue)

# Subset neurons
caudSeur = subset(caudSeur, subset = newannot %in% c('SPN', 'Non_SPN'))
putSeur = subset(putSeur, subset = newannot %in% c('SPN', 'Non_SPN'))
seurM = merge(caudSeur, putSeur)
seurM$Tissue = seurM$tissue
seurM$tissue = NULL
rm(caudSeur, putSeur)
gc()

####
## SET VARIABLES FOR CLUSTERING
####

pref = 'Macaque_CaudatePutamen_NEURONS'

markerGenes = c('RBFOX3', 'SATB2', 'SLC17A7', 'FOSL2', 'NPAS4', 'NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CFTR', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA')

####
## CLUSTER-1
####

setwd("/home2/gkonop/workdir/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN_NEURONS/CLUSTERING_1")

seurM_corr = clusterBatchCorrect(seurM, npcs = 30, res =1, pref = pref, sampleName = 'orig.ident')
stackedbarplot(seurM_corr[[]], 'seurat_clusters', 'Tissue', paste0(pref, '_Tissue'), wd = 10, hg = 10)
stackedbarplot(seurM_corr[[]], 'seurat_clusters', 'type', paste0(pref, '_Dataset'), wd = 10, hg = 10)

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
write_rds(seurM_corr, '~/workdir/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/Macaque_CaudatePutamen_NEURONS_Clustering1.RDS')
seurM_corr = read_rds('~/workdir/01_MACAQUE/03_CLUSTER_ANNOTATE/Macaque_CaudatePutamen_NEURONS_Clustering1.RDS')

qw = FindMarkers(seurM_corr, ident.1 = 21, only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
pdf('Cluster_4_Markers.pdf', width = 12)
DotPlot(seurM_corr, features = rownames(qw)[1:30]) + rotate_x_text(45)
dev.off()



####
## SRIATUM NEURON-GLIA DOUBLET SUSPECTOR IN NEURONS
####
#cluster_number,cell_type
#0,iSPN
#1,dSPN
#2,iSPN
#3,dSPN
#4,iSPN
#5,dSPN
#6,ePSN
#7,dSPN
#8,iSPN
#9,dSPN
#10,iSPN
#11*,iSPN+MOL?
#12,iSPN
#13,iSPN
#14*,SPN+MOL
#15,iSPN
#16*,iSPN+Ast?
#17,Non_SPN_FOXP2?
#18,Non_SPN?
#19,iSPN
#20*,SPN+OPC
#21*,Neu+Microglia
#22,Non_SPN
#23,Non_SPN
#24,dSPN

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
annots[annots$Cluster %in% c(21), 'CellType'] = 'Microglia'

####
## CLUSTER-2
####
setwd("/home2/gkonop/workdir/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN_NEURONS/CLUSTERING_2/")
# Remove doublets and recluster
doublets = annots[annots$CellType != 'Not_Doublet', 'Cluster'] %>% as.character %>% as.numeric
empty_droplets = c()

seurM = subset(seurM_corr, subset = seurat_clusters %in% c(doublets, empty_droplets), invert = T)
seurM_corr = clusterBatchCorrect(seurM, npcs = 30, res = 1, pref = pref, sampleName = 'orig.ident')

stackedbarplot(seurM_corr[[]], 'seurat_clusters', 'orig.ident', paste0(pref, '_Samples'), wd = 10, hg = 10)
stackedbarplot(seurM_corr[[]], 'seurat_clusters', 'Tissue', paste0(pref, '_Tissue'), wd = 10, hg = 10)

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
write_rds(seurM_corr, '~/workdir/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/Macaque_CaudatePutamen_NEURONS_Clustering2.RDS')
seurM_corr = read_rds('~/workdir/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/Macaque_CaudatePutamen_NEURONS_Clustering2.RDS')


qw = FindMarkers(seurM_corr, ident.1 = 13, only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
pdf('Cluster_4_Markers.pdf', width = 12)
DotPlot(seurM_corr, features = rownames(qw)[1:30]) + rotate_x_text(45)
dev.off()

####
## CLUSTER-3
####
setwd("/home2/gkonop/workdir/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN_NEURONS/CLUSTERING_3/")

#cluster_number,cell_type
#0,dSPN
#1,iSPN
#2,iSPN
#3,dSPN
#4,iSPN
#5,ePSN
#6,dSPN
#7,iSPN
#8,iSPN
#9,dSPN
#10,iSPN
#11,iSPN
#12,dSPN
#13,Non_SPN_FOXP2
#14,Non_SPN?
#15,iSPN
#16,iSPN
#17*,iSPN+MOL?
#18*,Exc_neu+OPC?
#19,dSPN
#20*,SPN+OPC


# Remove doublets and recluster
doublets = c(17,18,20)
empty_droplets = c()

seurM = subset(seurM_corr, subset = seurat_clusters %in% c(doublets, empty_droplets), invert = T)
seurM_corr = clusterBatchCorrect(seurM, npcs = 30, res = 0.5, pref = pref, sampleName = 'orig.ident')

stackedbarplot(seurM_corr[[]], 'seurat_clusters', 'Tissue', paste0(pref, '_Tissue'), wd = 10, hg = 10)
stackedbarplot(seurM_corr[[]], 'seurat_clusters', 'type', paste0(pref, '_Dataset'), wd = 10, hg = 10)

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
saveRDS(seurM_corr, '~/workdir/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/GB_Macaque_CaudatePutamen_NEURONS_Clustering3.RDS')
seurM_corr = read_rds('~/workdir/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/GB_Macaque_CaudatePutamen_NEURONS_Clustering3.RDS')

qw = FindMarkers(seurM_corr, ident.1 = 15, only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
pdf('Cluster_4_Markers.pdf', width = 12)
DotPlot(seurM_corr, features = rownames(qw)[1:30]) + rotate_x_text(45)
dev.off()

# PLOT NEURONAL CELL TYPE MARKERS (KRIENEN 2020)
marksToPlot = c('TAC3', 'TH', 'TRH', 'TRHDE', 'SOX6', 'PTK2B', 'PRDM1', 'EYA4', 'FRMD7', 'SHD', 'CCK', 'PVALB', 'PTHLH', 'PDGFD', 'EYA2', 'NTNG1', 'KLHL1', 'NDNF', 'SST', 'NPY', 'VIP', 'CHAT', 'SULF1', 'ANGPT1', 'RBFOX3', 'SYT1', 'SNAP25', 'CASZ1', 'OTOF', 'CACNG5', 'PCDH8', 'FOXP2', 'DRD1', 'DRD2', 'TAC1', 'PENK', 'MOG', 'MBP', 'PCDH15', 'NRGN', 'GAD1', 'SLC1A2', 'AQP4', 'APBB1IP', 'FLT1', 'SLC17A7', 'SATB2', 'MALAT1')

pdf(paste0(pref, "_GB_NEURON_MarkersDotPlot.pdf"), width = 20, height = 10)
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
macaqueSeur = read_rds('~/workdir/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/Macaque_CaudatePutamen_NEURONS_Clustering3.RDS')

# Remove doublet clusters
macaqueSeur = subset(macaqueSeur, subset = seurat_clusters %in% c(16,17), invert = T)

# Remove samples with drastically uneven distribution of clusters
macaqueSeur = subset(macaqueSeur, subset = orig.ident %in% c('SRR13808461', 'SRR13808459'), invert = T)

# Merge clusters from the same cell type
macaqueSeur$seurat_clusters = macaqueSeur$seurat_clusters %>% as.character
macaqueSeur$seurat_clusters = gsub('^7$|^2$', '0', macaqueSeur$seurat_clusters)
macaqueSeur$seurat_clusters = gsub('^3$|^8$', '1', macaqueSeur$seurat_clusters)

# Update some annotations in human
humanSeur$newannot = gsub('Exc_.*', 'Excitatory', humanSeur$newannot)

# Check how many genes will be used for label transfer
DefaultAssay(humanSeur) = 'SCT'
DefaultAssay(macaqueSeur) = 'SCT'
sum(VariableFeatures(macaqueSeur) %in% rownames(humanSeur@assays$SCT))

# Find transfer anchors
transfer_anchors = FindTransferAnchors(reference = humanSeur, query = macaqueSeur,
					dims = 1:20,
					reference.assay = "SCT",
					query.assay = "SCT",
					reduction = "cca")

# Make cell type predictions based on anchors
pred = TransferData(anchorset = transfer_anchors,
			refdata = humanSeur$newannot, 
			weight.reduction = macaqueSeur[["harmony"]],
			dims = 1:20)

# Add predictions to meta data
macaqueSeur$pred = pred$predicted.id
macaqueSeur$pred_score = pred$prediction.score.max

pdf("MaxPredictionScores_ToMACAQUE_FromHUMAN.pdf")
hist(macaqueSeur$pred_score)
dev.off()

# Plot label transfer
# Make the colors match
#newrna_int$pred <- factor(newrna_int$pred, levels = levels(allen_mtg_seur$cell_type))
Idents(macaqueSeur) = macaqueSeur$pred
Idents(humanSeur) = humanSeur$newannot
transferPlot(macaqueSeur, humanSeur, "LABEL_TRANSFER_ToMACAQUE_FromHUMAN.pdf")

# Examine label transfer per cluster. Increase resolution if cluster has two distinct labels.
macaqueSeur$predicted.id = macaqueSeur$pred
ltpred(macaqueSeur, vars = c('predicted.id', 'seurat_clusters'), fn = "LABEL_TRANSFER_CLARITY_ToMACAQUE_FromHUMAN.pdf")



# Jaccard similarity between prediction and annotation
meta = macaqueSeur@meta.data
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
toplot$Annotation = factor(toplot$Annotation, levels = c(0,1,4:6,9,10:15))
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

pdf('PredictionPlot_JACCARD_ToMACAQUE_FromHUMAN.pdf', width = 10, height =10)
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

toplot$Macaque = toplot$Annotation
toplot$Human = toplot$Prediction
toplot[["Macaque"]] = mapnames[toplot[["Annotation"]] %>% as.character]
toplot[toplot$Annotation == 11, 'Macaque'] = 'FOXP2_EYA2'
toplot[toplot$Macaque == 'VIP', 'Macaque'] = 'VIP_CCK_EYA4'

toplot$Human = factor(toplot$Human, levels = c('dSPN', 'iSPN', 'Excitatory', 'eSPN_FOXP2', 'PVALB_PDGFD', 'TAC3', 'PVALB_FLT3', 'eSPN_DRD2', 'SST_NPY', 'VIP', 'CCK_EYA4', 'CHAT'))
toplot$Macaque = factor(toplot$Macaque, levels = c('dSPN', 'iSPN', 'Excitatory', 'eSPN_FOXP2', 'PVALB_PDGFD', 'TAC3', 'PVALB_FLT3', 'eSPN_DRD2', 'SST_NPY', 'FOXP2_EYA2', 'VIP_CCK_EYA4', 'CHAT'))

pdf('PredictionPlot_JACCARD_ToMACAQUE_FromHUMAN_FINAL.pdf', width = 10, height =10)
ggscatter(toplot, x = 'Macaque', y = 'Human', color = 'Jaccard_Index', size = 10) +
scale_color_gradient2(low = "white", high = "blue", mid = "white", midpoint = 0, limit = c(0,1)) +
xlab('Macaque') + ylab('Human') +
theme(axis.text.x = element_text(size=20),
	axis.text.y = element_text(size=20),
	axis.title = element_text(size=20)) +
	rotate_x_text(90) + grids(linetype = "dashed", color = 'grey') +
geom_text(aes(label = is_best), vjust = 0.8, colour = "red", size = 15 )
dev.off()



mapnames = setNames(toplot$Macaque %>% as.character,
		      toplot$Annotation %>% as.character)

macaqueSeur[["newannot"]] = mapnames[macaqueSeur[["seurat_clusters"]][,1]]
Idents(macaqueSeur) = macaqueSeur$newannot

# PLOTS
pdf(paste0(pref, "_Clusters_Annotated.pdf"))
DimPlot(macaqueSeur, group.by = 'newannot', label = T, raster = T, repel = T, label.size = 7) +
NoLegend() +
ggtitle('')
dev.off()

pdf(paste0(pref, "_Depth_Gene_Annotated.pdf"))
ggboxplot(macaqueSeur[[]], x = 'newannot', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_Depth_UMI_Annotated.pdf"))
ggboxplot(macaqueSeur[[]], x = 'newannot', y = 'nCount_RNA', ylim = c(0,80000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_NuclearFraction_Annotated.pdf"))
ggboxplot(macaqueSeur[[]], x = 'newannot', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()



marksToPlot = c('FOSL2', 'NR4A1', 'NPAS4', 'TAC3', 'TH', 'TRH', 'TRHDE', 'SOX6', 'PTK2B', 'PRDM1', 'EYA4', 'FRMD7', 'SHD', 'CCK', 'PVALB', 'PTHLH', 'PDGFD', 'EYA2', 'NTNG1', 'KLHL1', 'NDNF', 'SST', 'NPY', 'VIP', 'CHAT', 'SULF1', 'ANGPT1', 'RBFOX3', 'SYT1', 'SNAP25', 'CASZ1', 'OTOF', 'CACNG5', 'PCDH8', 'FOXP2', 'DRD1', 'DRD2', 'TAC1', 'PENK', 'MOG', 'MBP', 'PCDH15', 'NRGN', 'GAD1', 'SLC1A2', 'AQP4', 'APBB1IP', 'FLT1', 'SLC17A7', 'SATB2', 'MALAT1')

DefaultAssay(macaqueSeur) = 'RNA'
macaqueSeur = NormalizeData(macaqueSeur)
pdf(paste0(pref, "_MarkersDotPlot.pdf"), width = 25, height = 10)
DotPlot(macaqueSeur, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

stackedbarplot(macaqueSeur[[]], groupfill = 'newannot', groupx = 'orig.ident', fn = paste0(pref, '_Stacked_Annotated'))

stackedbarplot(macaqueSeur[[]], groupfill = 'orig.ident', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_2'))

stackedbarplot(macaqueSeur[[]], groupfill = 'Tissue', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Tissue'))

write_rds(macaqueSeur, '~/workdir/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/Macaque_CaudatePutamen_NEURONS_Annotated.RDS')
