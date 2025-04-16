library(dplyr)
library(plyr)
library(Seurat)
library(ggpubr)
library(reshape2)
library(tidyr)
library(tidyverse)
library(DropletQC)
source('~/SCRIPTS/utility_functions.R')
set.seed(1234)


#####
## LOAD COUNT MATRIX
#####

# Set the samples and read them with QC
samplePath = '/home2/gkonop/project/03_MATRIX_FROM_CELLBENDER/FERRET_KRIENEN'
(samps = dir(samplePath))
bam_file_dirs <- c("/home2/gkonop/project/00_BAM_DOWNLOADED/KRIENEN_2020/FERRET/SRR11921037/cellranger_count/repaired_krienen_ferret_rxn1", "/home2/gkonop/project/00_BAM_DOWNLOADED/KRIENEN_2020/FERRET/SRR11921038/cellranger_count/repaired_krienen_ferret_rxn2")

allSeurL = list()
for(i in 1:length(samps)){

	#  Read CellBender output
	mat = Read10X_h5(paste0(samplePath, '/', samps[i], '/', samps[i], '_cellbender_output_top80000_filtered.h5'))

	# Soft filtering. Remove lower than 200 UMI
	mat = mat[, colSums(mat) > 200]
	cbBarcode = colnames(mat)
	colnames(mat) = gsub('-1', '', colnames(mat)) %>% paste0(., '_', samps[i])

	# Create seurat object
	tmpSeur = CreateSeuratObject(counts = mat, min.cells = 10, min.features = 1)

	# Calculate intronic read ratio. Skip this since we did not keep the BAMs after assigning to liftoffed genes
	nf2 = nuclear_fraction_annotation(
		annotation_path = '/home2/gkonop/workdir/REFERENCES/GTF_ONLY/GCF_011764305.1_ASM1176430v1.1_genomic.gtf',
		bam = paste0(bam_file_dirs[i], '/outs/possorted_genome_bam.bam'),
		barcodes = cbBarcode,
		cores = 23, verbose = T)

	nf2$barc = colnames(tmpSeur)
	tmpSeur$intronRat = nf2$nuclear_fraction
	tmpSeur$orig.ident = samps[i]
	tmpSeur$Species = 'Ferret'
	allSeurL[[i]] = tmpSeur

	print(i)
}

allseur = Reduce(merge, allSeurL)
saveRDS(allseur, '/home2/gkonop/workdir/01_FERRET/03_CLUSTER_ANNOTATE/KRIENEN_FERRET_SEURAT_RAW.RDS')
allseur = readRDS('/home2/gkonop/workdir/01_FERRET/03_CLUSTER_ANNOTATE/KRIENEN_FERRET_SEURAT_RAW.RDS')

# Add meta data
# ferret striatum data are only from caudate: https://www.ncbi.nlm.nih.gov/bioproject?LinkName=biosample_bioproject&from_uid=15096894

####
## CLUSTERING 1
####
table(allseur$orig.ident)

#SRR11921037 SRR11921038 
#      14245       14405 

pref = 'Ferret_Caudate_Krienen'
seurM = allseur
seurM = subset(seurM, subset = nCount_RNA > 0)
seurM_corr = clusterBatchCorrect(seurM, npcs = 30, res = 0.5, pref = pref, sampleName = 'orig.ident')
allseur$Tissue <- rep("Caudate", times = nrow(allseur[[]]))

pref = 'Ferret_Caudate_Krienen_res_1'
seurM = allseur
seurM = subset(seurM, subset = nCount_RNA > 0)
seurM_corr = clusterBatchCorrect(seurM, npcs = 30, res = 1, pref = pref, sampleName = 'orig.ident')
allseur$Tissue <- rep("Caudate", times = nrow(allseur[[]]))


# QC PLOTS
setwd("/home2/gkonop/workdir/01_FERRET/03_CLUSTER_ANNOTATE/CLUSTERING_1")
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
# set variables
markerGenes = c('NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CFTR', 'CHAT', 'ADARB2', 'LHX6', 'MEIS2', 'TAC3', 'VIP', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA', 'SATB2', 'SLC17A7', 'EYA2', 'PDGFD', 'PTHLH', 'PVALB', 'SPARCL1', 'RSPO2', 'RMST')



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

# Save (optional)
saveRDS(seurM_corr, paste0(pref, '_clustering1.RDS'))


####
## CLUSTERING 2
####
#CLUSTER#,Annotation
#0,MOL
#1,Non_SPN
#2,SPN
#3,SPN
#4,Microglia
#5,MOL
#6,Astrocyte
#7,SPN
#8,Non_SPN
#9*,empty
#10*,empty
#11,OPC
#12,Non_SPN
#13*,Astrocyte-neu doublet?
#14,Non_SPN
#15*,MOL+neu doublet?
#16,SPN
#17*,empty
#18*,SPN+glia doublet
#19*,Endo
#20*,Endo
#21,MOL
#22,Non_SPN
#23,MOL
#24,Non_SPN
#25*,OPC+neu doublet
#26,SPN
#27,SPN
#28,Non_SPN
#29,OPC
#30,Microglia
#31,OPC
#32,SPN
#################
setwd("/home2/gkonop/workdir/01_FERRET/03_CLUSTER_ANNOTATE/CLUSTERING_2")
pref = 'Ferret_Caudate_Krienen'
vasculature = c(19,20)
empty_drop = c(9,10,17)
doublet = c(13,15,18,25)
seurM = subset(seurM_corr, subset = seurat_clusters %in% c(vasculature, empty_drop, doublet), invert = T)
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
markerGenes = c('NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA')

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
 table(seurM_corr$orig.ident)

#SRR11921037 SRR11921038 
#      11489       11698 
# Save (optional)
saveRDS(seurM_corr, paste0(pref, '_clustering2.RDS'))


####
## ANNOTATE
####
#Cluster#,annotation
#0,MOL
#1,dSPN
#2,iSPN
#3,Microglia
#4,dSPN
#5,Astrocyte
#6,iSPN
#7,dSPN
#8,OPC
#9,iSPN
#10,Non_SPN
#11,eSPN
#12,dSPN
#13,MOL
#14,Non_SPN
#15,MOL
#16,MOL
#17,Non_SPN
#18,eSPN
#19,OPC
#20,Microglia
#21,dSPN
#22,OPC

setwd("/home2/gkonop/workdir/01_FERRET/03_CLUSTER_ANNOTATE/ANNOTATED")
mapnames <- (setNames(c("MOL","dSPN","iSPN","Microglia","dSPN","Astrocyte","iSPN",
"dSPN","OPC","iSPN","Non_SPN","eSPN","dSPN","MOL","Non_SPN","MOL","MOL","Non_SPN",
"eSPN","OPC","Microglia","dSPN","OPC"), c(0:22)))
seurM_corr[["newannot"]] = mapnames[seurM_corr[["seurat_clusters"]][,1]]

table(seurM_corr$newannot)
#Astrocyte      dSPN      eSPN      iSPN Microglia       MOL   Non_SPN       OPC 
#     1916      6197       698      4828      2203      5235      1083      1027 
# broad_annotations
mapnames <- (setNames(c("MOL","SPN","SPN","Microglia","SPN","Astrocyte","SPN","SPN","OPC","SPN",
"Non_SPN","SPN","SPN","MOL","Non_SPN","MOL","MOL","Non_SPN","SPN","OPC","Microglia","SPN","OPC"), c(0:22)))
seurM_corr[["broad_annot"]] = mapnames[seurM_corr[["seurat_clusters"]][,1]]
table(seurM_corr$broad_annot)

#Astrocyte Microglia       MOL   Non_SPN       OPC       SPN 
#     1916      2203      5235      1083      1027     11723 


# PLOTS
pdf(paste0(pref, "_Clusters_Annotated.pdf"))
DimPlot(seurM_corr, group.by = 'broad_annot', label = T, raster = T, repel = T, label.size = 7) +
NoLegend() +
ggtitle('')
dev.off()

pdf(paste0(pref, "_Depth_Gene_Annotated.pdf"))
ggboxplot(seurM_corr[[]], x = 'broad_annot', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_Depth_UMI_Annotated.pdf"))
ggboxplot(seurM_corr[[]], x = 'broad_annot', y = 'nCount_RNA', ylim = c(0,80000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_NuclearFraction_Annotated.pdf"))
ggboxplot(seurM_corr[[]], x = 'broad_annot', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()


DefaultAssay(seurM_corr) = 'RNA'
seurM_corr = NormalizeData(seurM_corr)

pdf(paste0(pref, "_MarkersDotPlot.pdf"), width = 25, height = 10)
DotPlot(seurM_corr, features = markerGenes, group.by = "broad_annot") +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

stackedbarplot(seurM_corr[[]], groupfill = 'broad_annot', groupx = 'Species', fn = paste0(pref, '_Stacked_Annotated_Species'))

stackedbarplot(seurM_corr[[]], groupfill = 'broad_annot', groupx = 'orig.ident', fn = paste0(pref, '_Stacked_Annotated_Samples'))

stackedbarplot(seurM_corr[[]], groupfill = 'orig.ident', groupx = 'broad_annot', fn = paste0(pref, '_Stacked_AnnotatedSamples2'), wd = 20)

stackedbarplot(seurM_corr[[]], groupfill = 'Tissue', groupx = 'broad_annot', fn = paste0(pref, '_Stacked_Annotated_Tissue'))


# Save (optional)
saveRDS(seurM_corr, paste0(pref, '_ANNOTATED.RDS'))
