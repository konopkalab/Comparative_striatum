# source ~/load_modules.sh
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
library(ggrepel)
source('~/SCRIPTS/utility_functions.R')
set.seed(1234)
# Define fixed jitter position with seed
fixed_jitter <- position_jitter(width = 0.2, height = 0, seed = 123)

################################################################
####
## LOAD NEURON DATASETS
####

### integrate interneurons of: Human, chimp, macaque, marmoset, bat, mouse, and ferret
# human
human <- readRDS("/home2/gkonop/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN/ANNOTATED/human_integrated_cuadate_putamen_ANNOTATED.RDS")

# chimp
chimp <- readRDS("/home2/gkonop/workdir/01_CHIMP/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN/ANNOTATED/chimp_integrated_cuadate_putamen_ANNOTATED.RDS")
chimp$Species <- "Chimp"

# macaque
macaque <- readRDS("/home2/gkonop/workdir/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/Caudate_Putamen/ANNOTATED/macaque_integrated_cuadate_putamen_ANNOTATED.RDS")

# marmoset
marmoset <- readRDS("/home2/gkonop/workdir/01_MARMOSET/MARMOSET_INTEGRATION/Caudate_putamen/ANNOTATED/marmoset_integrated_cuadate_putamen_ANNOTATED.RDS")

# mouse
mouse <- readRDS("/home2/gkonop/workdir/01_MOUSE/03_CLUSTER_ANNOTATE/Mouse_Caudate_Annotated_FINAL.RDS")
mouse$Species <- "Mouse"
table(mouse$newannot)
mouse$Tissue <- "Caudoputamen"

# bat
bat <- readRDS("/home2/gkonop/workdir/01_BAT/03_CLUSTER_ANNOTATE/Caudate_Putamen/ANNOTATED/bat_integrated_cuadate_putamen_ANNOTATED.RDS")

# ferret
ferret <- readRDS("/home2/gkonop/workdir/01_FERRET/03_CLUSTER_ANNOTATE/ANNOTATED/Ferret_Caudate_Krienen_ANNOTATED.RDS")
# create "id" metadata
ferret$id <- ferret$orig.ident
ferret$newannot <- ferret$broad_annot
########################################################################
# Found ortholog genes from human pr coding genes and extracted the genes which are ortholog in all 7 species using ncbi datasets tool
ortho_genes <- read.table("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/orthologs_in_7_species_human_pr_codingOnly_2.csv",  header = TRUE)$symbol #15055 pr coding ortho genes

### subset the genes and the interneurons from each Species
## For Human
# subset cells and metadata
human_non_spns <- subset(human, subset = newannot == "Non_SPN")
meta_human <- human_non_spns[[]]
# subset count matrix
mat <- human_non_spns@assays$RNA@counts
new_mat <- mat[rownames(mat) %in% ortho_genes,]
# generate new Seurat object
human_non_spns_sub <- CreateSeuratObject(count = new_mat, meta.data = meta_human)

## For Chimp
# subset cells and metadata
chimp_non_spns <- subset(chimp, subset = newannot == "Non_SPN")
meta_chimp <- chimp_non_spns[[]]
# subset count matrix
mat_chimp <- chimp_non_spns@assays$RNA@counts
new_mat_chimp <- mat_chimp[rownames(mat_chimp) %in% ortho_genes,]
# generate new Seurat object
chimp_non_spns_sub <- CreateSeuratObject(count = new_mat_chimp, meta.data = meta_chimp)

## For Macaque
# subset cells and metadata
macaque_non_spns <- subset(macaque, subset = newannot == "Non_SPN")
meta_macaque <- macaque_non_spns[[]]
# subset count matrix
mat_macaque <- macaque_non_spns@assays$RNA@counts
new_mat_macaque <- mat_macaque[rownames(mat_macaque) %in% ortho_genes,]
# generate new Seurat object
macaque_non_spns_sub <- CreateSeuratObject(count = new_mat_macaque, meta.data = meta_macaque)

## For Marmoset
# subset cells and metadata
marmoset_non_spns <- subset(marmoset, subset = newannot == "Non_SPN")
meta_marmoset <- marmoset_non_spns[[]]
# subset count matrix
mat_marmoset <- marmoset_non_spns@assays$RNA@counts
new_mat_marmoset <- mat_marmoset[rownames(mat_marmoset) %in% ortho_genes,]
# generate new Seurat object
marmoset_non_spns_sub <- CreateSeuratObject(count = new_mat_marmoset, meta.data = meta_marmoset)

## For Mouse
# subset cells and metadata
mouse_non_spns <- subset(mouse, subset = newannot == "Non_SPN")
meta_mouse <- mouse_non_spns[[]]
# subset count matrix
mat_mouse <- mouse_non_spns@assays$RNA@counts
rownames(mat_mouse) <- toupper(rownames(mat_mouse))
new_mat_mouse <- mat_mouse[rownames(mat_mouse) %in% ortho_genes,]
# generate new Seurat object
mouse_non_spns_sub <- CreateSeuratObject(count = new_mat_mouse, meta.data = meta_mouse)

## For Bat
# subset cells and metadata
bat_non_spns <- subset(bat, subset = newannot == "Non_SPN")
meta_bat <- bat_non_spns[[]]
# subset count matrix
mat_bat <- bat_non_spns@assays$RNA@counts
new_mat_bat <- mat_bat[rownames(mat_bat) %in% ortho_genes,]
# generate new Seurat object
bat_non_spns_sub <- CreateSeuratObject(count = new_mat_bat, meta.data = meta_bat)

## For Ferret
# subset cells and metadata
ferret_non_spns <- subset(ferret, subset = newannot == "Non_SPN")
meta_ferret <- ferret_non_spns[[]]
# subset count matrix
mat_ferret <- ferret_non_spns@assays$RNA@counts
new_mat_ferret <- mat_ferret[rownames(mat_ferret) %in% ortho_genes,]
# generate new Seurat object
ferret_non_spns_sub <- CreateSeuratObject(count = new_mat_ferret, meta.data = meta_ferret)


# Merge
# Step 1: List all the Seurat objects you want to merge
seurat_list <- list(human_non_spns_sub, 
                    chimp_non_spns_sub, 
                    macaque_non_spns_sub, 
                    marmoset_non_spns_sub, 
                    mouse_non_spns_sub, 
                    bat_non_spns_sub, 
                    ferret_non_spns_sub)

# Step 2: Merge all Seurat objects into one
seurM <- Reduce(merge, seurat_list)

saveRDS(seurM, '/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/NonSPN_seurM_AllSpecies_ALLTissues_human_pr_coding_orthologs.RDS')

####
## SET VARIABLES
####
# PPP1R1B IS AN SPN MARKER!
# TSHZ1 IS A PATCH MARKER
marksToPlot <- c('NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CFTR', 'CHAT', 'ADARB2', 'LHX6', 'MEIS2', 'TAC3', 'VIP', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA', 'SATB2', 'SLC17A7', 'EYA2', 'PDGFD', 'PTHLH', 'PVALB', 'SPARCL1', 'RSPO2', 'RMST', 'LMO3', 'TSHZ2', "CALB1", "CALB2", "NOS1", "TSHZ1", "OPRM1", "CHST9","GRM8", "PPP1R1B", 'GRIK3', 'CXCL14')

pref = 'GB_Non_SPN_AllSpecies_AllTissues_human_pr_coding_orthologs'

####
## INTEGRATE ACROSS SPECIES
####

# Split data
seurM$id = paste0(seurM$orig.ident, '_', seurM$Tissue)
seurML = SplitObject(seurM, split.by = "id")

# normalize and identify variable features for each dataset independently
seurML = lapply(X = seurML, FUN = function(x) {
    x = NormalizeData(x)
    x = FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features = SelectIntegrationFeatures(object.list = seurML)
print(paste0('Number of genes to use for integration: ', length(features)))
seurML = lapply(X = seurML, FUN = function(x) {
    x = ScaleData(x, features = features, verbose = FALSE)
    x = RunPCA(x, features = features, verbose = FALSE)
})

anchors = FindIntegrationAnchors(object.list = seurML, anchor.features = features,  normalization.method = "LogNormalize", reduction = "rpca")
saveRDS(anchors, '/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/Anchors_NonSPN_AllSpecies_ALLTissues_human_pr_coding_orthologs.RDS')

# save memory prior integration
gc()
###################################################3
## troubleshooting
## Remove the HE et al macaque sample SRR13808463_Caudate which has < 100 interneurons
# Iterate over the Seurat list (seurML)
seurML <- lapply(1:length(seurML), function(i) {
  # Check the number of cells in the current Seurat object
  num_cells <- dim(seurML[[i]]@assays$RNA@data)[2]
  
  # If the number of cells is less than 100, return NULL to remove it from the list
  if (num_cells < 100) {
    return(NULL)
  } else {
    return(seurML[[i]])
  }
})

# Remove NULL entries (those with less than 100 cells)
seurML <- seurML[!sapply(seurML, is.null)]

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features = SelectIntegrationFeatures(object.list = seurML)
print(paste0('Number of genes to use for integration: ', length(features)))
seurML = lapply(X = seurML, FUN = function(x) {
    x = ScaleData(x, features = features, verbose = FALSE)
    x = RunPCA(x, features = features, verbose = FALSE)
})

anchors = FindIntegrationAnchors(object.list = seurML, anchor.features = features,  normalization.method = "LogNormalize", reduction = "rpca")
saveRDS(anchors, '/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/Anchors_NonSPN_AllSpecies_ALLTissues_human_pr_coding_orthologs.RDS')

# save memory prior integration
gc()
# Check the features (gene names) across all Seurat objects
feature_list <- lapply(seurML, function(x) rownames(x@assays$RNA@data))

# Find the intersection of features across all Seurat objects
common_features <- Reduce(intersect, feature_list)

# Make sure all datasets have the same features
lapply(seurML, function(x) {
  missing_features <- setdiff(common_features, rownames(x@assays$RNA@data))
  if(length(missing_features) > 0) {
    print(paste("Missing features in object: ", x$orig.ident))
    return(missing_features)
  }
})

# Check the number of cells (columns) in each Seurat object
cell_counts <- sapply(seurML, function(x) ncol(x@assays$RNA@data))
print(cell_counts)
# Check metadata for missing or inconsistent values
lapply(seurML, function(x) {
  if (any(is.na(x@meta.data))) {
    print(paste("NA values in metadata for object: ", x$orig.ident))
    return(which(is.na(x@meta.data)))
  }
})

# Check if each Seurat object has variable features
lapply(seurML, function(x) length(x@assays$RNA@var.features))

# Function to integrate Seurat objects and catch errors during merging
integrate_with_error_handling <- function(seurML, anchors) {
  for (i in 1:(length(seurML)-1)) {
    tryCatch({
      # Print the dataset merge step to debug
      print(paste("Merging dataset", i, "into", i+1))
      
      # Integrate data
      integrated_data <- IntegrateData(anchorset = anchors, normalization.method = "LogNormalize")
      
      # If successful, print success message
      print(paste("Integration succeeded for dataset pair", i, "and", i+1))
      
    }, error = function(e) {
      # Print which step caused the error
      print(paste("Error during integration for dataset pair", i, "and", i+1, ":", e$message))
    })
  }
}

# Call the function to handle integration and print errors for each dataset merge
integrate_with_error_handling(seurML, anchors)
################################################################
# this command creates an 'integrated' data assay
#(kweg = floor(min(table(seurM$orig.ident))/10)*10) #90
options(future.globals.maxSize = 16000 * 1024^2)
allseur_integrated = IntegrateData(anchorset = anchors, k.weight = 60, normalization.method = 'LogNormalize')
saveRDS(allseur_integrated, '/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/Integrated_rpca_NonSPN_AllSpecies_ALLTissues_human_pr_coding_orthologs.RDS')


# Clustering
allseur_integrated = ScaleData(allseur_integrated, verbose = FALSE)
allseur_integrated = RunPCA(allseur_integrated, verbose = FALSE)
allseur_integrated = RunUMAP(allseur_integrated, dims = 1:20, reduction = 'pca')

# Basic Plots
setwd("~/workdir/03_INTEGRATE_ALL/Non_SPN/CLUSTERING_1")
pdf(paste0(pref, "_INTEGRATED_UMAP.pdf"))
DimPlot(allseur_integrated, group.by = 'Species', raster = T)
dev.off()

DefaultAssay(allseur_integrated) <- "integrated"
allseur_integrated = FindNeighbors(allseur_integrated, dims = 1:20, reduction = 'pca')
allseur_integrated = FindClusters(allseur_integrated, resolution = 1)

pdf(paste0(pref, "_INTEGRATED_CLUSTERS.pdf"))
DimPlot(allseur_integrated, label = T, raster = T) + NoLegend()
dev.off()

pdf(paste0(pref, "_INTEGRATED_PRIMATE_ANNOT.pdf"))
DimPlot(allseur_integrated, label = T, raster = T, group.by = 'newannot') + NoLegend()
dev.off()

stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_INTEGRATED_STACK_SAMPLES'), wd = 20)
stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'Species', paste0(pref, '_INTEGRATED_STACK_SPECIES'))

# Plot previously identified markers
DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)

pdf(paste0(pref, "_INTEGRATED_DOTPLOT.pdf"), width = 20, height = 10)
DotPlot(allseur_integrated, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "_INTEGRATED_DEPTH.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'black') + NoLegend() +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

saveRDS(allseur_integrated, paste0('/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/CLUSTERING_1/', pref, '_integrated_CLUSTERING1.RDS'))
allseur_integrated <- readRDS(paste0('/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/CLUSTERING_1/', pref, '_integrated_CLUSTERING1.RDS'))
###########################################################################
####
## FILTER AND REPLOT
####
#Cluster_#,Annotation
#0,PDGFD_PTHLH_PVALB+
#1,SST_NPY
#2,PDGFD_PTHLH_PVALB+
#3,LMO3
#4,TAC3
#5*,empty
#6,TH
#7,TAC3
#8,CHAT
#9,LMO3
#10,CHAT
#11*,Exc_neu
#12,LMO3_BAT
#13*,MOL
#14,CCK_VIP-
#15*,Ast
#16,CCK_VIP+
#17,PDGFD_PVALB-
#18,SST_NPY
#19,TAC3+MOL
#20*,microglia+MOL
#21*,PDGFD+MOL
#22,TAC3
#23,FOXP2_TSHZ2
#24,LMO3?
#25*,Exc_neu
#26,RSPO2
#27,PDGFD_PTHLH_PVALB+
#28,FOXP2_EYA2
#29,SST_NPY
#30*,OPC
#31,FOXP2_EYA2
#32,CHAT
#33*,MOL doublet
#34*,empty

(toremove = c(5,11,13,15,20,21,25,30,33,34))
allseur_integrated = subset(allseur_integrated, subset = seurat_clusters %in% toremove, invert = T)

# Clustering
DefaultAssay(allseur_integrated) = 'integrated'
allseur_integrated = ScaleData(allseur_integrated, verbose = FALSE)
allseur_integrated = RunPCA(allseur_integrated, verbose = FALSE)
allseur_integrated = RunUMAP(allseur_integrated, dims = 1:20, reduction = 'pca')

setwd("/home2/gkonop/workdir/03_INTEGRATE_ALL/Non_SPN/CLUSTERING_2")

# Basic Plots
pdf(paste0(pref, "_INTEGRATED_UMAP.pdf"))
DimPlot(allseur_integrated, group.by = 'Species', raster = T, label.size = 7)
dev.off()

allseur_integrated = FindNeighbors(allseur_integrated, dims = 1:20, reduction = 'pca')
allseur_integrated = FindClusters(allseur_integrated, resolution = 2)

pdf(paste0(pref, "_INTEGRATED_CLUSTERS.pdf"))
DimPlot(allseur_integrated, label = T, raster = T, label.size = 7) + NoLegend()
dev.off()

stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_INTEGRATED_STACK_SAMPLES'), wd = 20)
stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'Species', paste0(pref, '_INTEGRATED_STACK_SPECIES'))

# Plot previously identified markers
DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)

pdf(paste0(pref, "_INTEGRATED_DOTPLOT.pdf"), width = 25, height = 10)
DotPlot(allseur_integrated, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "_INTEGRATED_DEPTH.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'black') + NoLegend() +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

saveRDS(allseur_integrated, paste0('/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/CLUSTERING_2/', pref, '_integrated_CLUSTERING2.RDS'))
allseur_integrated <- readRDS(paste0('/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/CLUSTERING_2/', pref, '_integrated_CLUSTERING2.RDS'))
###############################################################################################
####
## CLUSTER and FILTER
####
#Cluster_#,Annotation
#0,SST_NPY
#1,PDGFD_PTHLH_PVALB-
#2,TAC3
#3,PDGFD_PTHLH_PVALB+
#4,TAC3
#5,TAC3
#6,LMO3
#7,PDGFD_PTHLH_PVALB+
#8,LMO3
#9,PDGFD_PTHLH_PVALB+
#10,CHAT
#11,PDGFD
#12*,MOL_doublet
#13,CCK_VIP-
#14,TAC3
#15,SST_NPY
#16,SST_NPY
#17,CCK_VIP+
#18,CHAT
#19*,doublet
#20,CHAT
#21,LMO3
#22*,MOL
#23,PDGFD_PTHLH_PVALB+
#24,LMO3
#25*,doublet
#26*,empty
#27,FOXP2_TSHZ2_BAT
#28,LMO3
#29,SST_NPY
#30,CHAT
#31,FOXP2_TSHZ2_BAT
#32,FOXP2_EYA2
#33,RSPO2
#34,FOXP2_EYA2
#35,FOXP2_TSHZ2_BAT
#36*,doublet
#37,RSPO2
#38,CCK_VIP-
#39*,MOL_DOUBLET
#40,CCK_VIP-
#41,?
#42,TH
#43,LMO3

(toremove = c(12,19,22,25,26,36,39))
allseur_integrated = subset(allseur_integrated, subset = seurat_clusters %in% toremove, invert = T)

# Clustering
DefaultAssay(allseur_integrated) = 'integrated'
allseur_integrated = ScaleData(allseur_integrated, verbose = FALSE)
allseur_integrated = RunPCA(allseur_integrated, verbose = FALSE)
allseur_integrated = RunUMAP(allseur_integrated, dims = 1:20, reduction = 'pca')

setwd("/home2/gkonop/workdir/03_INTEGRATE_ALL/Non_SPN/CLUSTERING_3")

# Basic Plots
pdf(paste0(pref, "_INTEGRATED_UMAP.pdf"))
DimPlot(allseur_integrated, group.by = 'Species', raster = T, label.size = 7)
dev.off()

allseur_integrated = FindNeighbors(allseur_integrated, dims = 1:20, reduction = 'pca')
allseur_integrated = FindClusters(allseur_integrated, resolution = 2)

pdf(paste0(pref, "_INTEGRATED_CLUSTERS.pdf"))
DimPlot(allseur_integrated, label = T, raster = T, label.size = 7) + NoLegend()
dev.off()

stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_INTEGRATED_STACK_SAMPLES'), wd = 20)
stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'Species', paste0(pref, '_INTEGRATED_STACK_SPECIES'))

# Plot previously identified markers
DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)

pdf(paste0(pref, "_INTEGRATED_DOTPLOT.pdf"), width = 25, height = 10)
DotPlot(allseur_integrated, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "_INTEGRATED_DEPTH.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'black') + NoLegend() +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

saveRDS(allseur_integrated, paste0('/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/CLUSTERING_3/', pref, '_integrated_CLUSTERING3.RDS'))
allseur_integrated <- readRDS(paste0('/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/CLUSTERING_3/', pref, '_integrated_CLUSTERING3.RDS'))
######################################################################################
####
## CLUSTER & FILTER
####
#Cluster_#,Annotation
#0,SST_NPY
#1,PDGFD_PTHLH_PVALB-
#2,PDGFD_PTHLH_PVALB+
#3,TAC3
#4,PDGFD
#5,LMO3_CALB1-
#6,TAC3
#7,PDGFD_PTHLH_PVALB+
#8,SST_NPY
#9,LMO3_CALB1+
#10,CHAT
#11,TAC3
#12,CCK_VIP-
#13,TAC3
#14,PDGFD_PTHLH_PVALB+
#15,TAC3
#16,TAC3
#17,SST_NPY
#18,CCK_VIP+
#19,PDGFD
#20,FOXP2_TSHZ2_BAT
#21,CHAT
#22,PDGFD
#23*,EXC_NEU
#24,CHAT
#25,PDGFD
#26,LMO3_CALB1-
#27,LMO3_CALB1-
#28*,AST_DOUBLET
#29,FOXP2_EYA2
#30,SST_NPY
#31,LMO3_CALB1-
#32,FOXP2_TSHZ2_BAT
#33,RSPO2
#34,FOXP2_EYA2
#35,RSPO2
#36,CHAT
#37,CCK_VIP-
#38,CCK_VIP-
#39*,MOL
#40,TH
#41,TH
#42,LMO3_CALB1+
#43,TAC3

(toremove = c(23, 28,39))
allseur_integrated = subset(allseur_integrated, subset = seurat_clusters %in% toremove, invert = T)

# Clustering
DefaultAssay(allseur_integrated) = 'integrated'
allseur_integrated = ScaleData(allseur_integrated, verbose = FALSE)
allseur_integrated = RunPCA(allseur_integrated, verbose = FALSE)
allseur_integrated = RunUMAP(allseur_integrated, dims = 1:20, reduction = 'pca')

setwd("/home2/gkonop/workdir/03_INTEGRATE_ALL/Non_SPN/CLUSTERING_4")

# Basic Plots
pdf(paste0(pref, "_INTEGRATED_UMAP.pdf"))
DimPlot(allseur_integrated, group.by = 'Species', raster = T, label.size = 7)
dev.off()

allseur_integrated = FindNeighbors(allseur_integrated, dims = 1:20, reduction = 'pca')
allseur_integrated = FindClusters(allseur_integrated, resolution = 2)

pdf(paste0(pref, "_INTEGRATED_CLUSTERS.pdf"))
DimPlot(allseur_integrated, label = T, raster = T, label.size = 7) + NoLegend()
dev.off()

stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_INTEGRATED_STACK_SAMPLES'), wd = 20)
stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'Species', paste0(pref, '_INTEGRATED_STACK_SPECIES'))

# Plot previously identified markers
DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)

pdf(paste0(pref, "_INTEGRATED_DOTPLOT.pdf"), width = 25, height = 10)
DotPlot(allseur_integrated, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "_INTEGRATED_DEPTH.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'black') + NoLegend() +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

saveRDS(allseur_integrated, paste0('/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/CLUSTERING_4/', pref, '_integrated_CLUSTERING4.RDS'))
allseur_integrated <- readRDS(paste0('/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/CLUSTERING_4/', pref, '_integrated_CLUSTERING4.RDS'))
####################################################################################################
####
## CLUSTER & FILTER
####
#Cluster#,Annotation
#0,PDGFD_PTHLH_PVALB+
#1,TAC3
#2,SST_NPY
#3,TAC3
#4,LMO3_PPP1R1B
#5,PDGFD
#6,PDGFD_PTHLH_PVALB-
#7,TAC3
#8,CHAT
#9,PDGFD_PTHLH_PVALB+
#10*,SPN
#11,SST_NPY
#12,TAC3
#13,PDGFD_PTHLH_PVALB+
#14,SST_NPY
#15,TAC3
#16,CCK_VIP-
#17,CCK_VIP+
#18,PDGFD
#19,CHAT
#20,FOXP2_TSHZ2_BAT
#21,SST_NPY
#22,PDGFD
#23,PDGFD
#24,LMO3_BAT
#25,CHAT
#26,FOXP2_EYA2
#27,LMO3_BAT
#28,RSPO2
#29,FOXP2_TSHZ2_BAT
#30,LMO3_BAT
#31,FOXP2_EYA2
#32,RSPO2
#33,SST_NPY
#34,CCK_VIP-
#35,CCK_VIP-
#36*,MOL
#37,PDGFD_PTHLH_PVALB+
#38,CHAT
#39,TH

DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)

# find markers of some clusters
Idents(allseur_integrated) <- "seurat_clusters"
qw = FindMarkers(allseur_integrated, ident.1 = 20, only.pos = T, logfc.threshold = 0.5, min.pct = 0.25)
DotPlot(allseur_integrated, features = rownames(qw)[1:30]) + rotate_x_text(90)


(toremove = c(10,36))
allseur_integrated = subset(allseur_integrated, subset = seurat_clusters %in% toremove, invert = T)

# Clustering
DefaultAssay(allseur_integrated) = 'integrated'
allseur_integrated = ScaleData(allseur_integrated, verbose = FALSE)
allseur_integrated = RunPCA(allseur_integrated, verbose = FALSE)
allseur_integrated = RunUMAP(allseur_integrated, dims = 1:20, reduction = 'pca')

setwd("/home2/gkonop/workdir/03_INTEGRATE_ALL/Non_SPN/CLUSTERING_5")

# Basic Plots
pdf(paste0(pref, "_INTEGRATED_UMAP.pdf"))
DimPlot(allseur_integrated, group.by = 'Species', raster = T, label.size = 7)
dev.off()

allseur_integrated = FindNeighbors(allseur_integrated, dims = 1:20, reduction = 'pca')
allseur_integrated = FindClusters(allseur_integrated, resolution = 2)

pdf(paste0(pref, "_INTEGRATED_CLUSTERS.pdf"))
DimPlot(allseur_integrated, label = T, raster = T, label.size = 7) + NoLegend()
dev.off()

stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_INTEGRATED_STACK_SAMPLES'), wd = 20)
stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'Species', paste0(pref, '_INTEGRATED_STACK_SPECIES'))

# Plot previously identified markers
DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)

pdf(paste0(pref, "_INTEGRATED_DOTPLOT.pdf"), width = 25, height = 10)
DotPlot(allseur_integrated, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "_INTEGRATED_DEPTH.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'black') + NoLegend() +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

saveRDS(allseur_integrated, paste0('/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/CLUSTERING_5/', pref, '_integrated_CLUSTERING5.RDS'))
allseur_integrated <- readRDS(paste0('/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/CLUSTERING_5/', pref, '_integrated_CLUSTERING5.RDS'))
###################################################################
####
## CLUSTER & FILTER
####
#Cluster_#,Annotation
#0,TAC3_TH
#1,SST_NPY
#2,PDGFD_PTHLH_PVALB-
#3,TAC3_TH
#4*,SPN
#5,PDGFD_PTHLH_PVALB+
#6,PDGFD
#7,PDGFD_PTHLH_PVALB+
#8,TAC3_TH
#9,CHAT
#10,SST_NPY
#11,PDGFD_PTHLH_PVALB+
#12,TAC3_TH
#13,SST_NPY
#14,TAC3_TH
#15,CCK_VIP-
#16,CCK_VIP+
#17,PDGFD
#18,CHAT
#19,LMO3_BAT
#20,PDGFD_PTHLH_PVALB-
#21,SST_NPY
#22,PDGFD
#23,CHAT
#24,LMO3_BAT
#25,LMO3_BAT
#26,FOXP2_EYA2
#27,LMO3_BAT
#28,FOXP2_TSHZ2_BAT
#29*,RSPO2
#30,FOXP2_EYA2
#31*,RSPO2
#32,SST_NPY
#33,CHAT
#34,CCK_VIP+
#35,CCK_VIP-
#36,PDGFD_PTHLH_PVALB+
#37,TH_macaque

(toremove = c(4,29,31))
allseur_integrated = subset(allseur_integrated, subset = seurat_clusters %in% toremove, invert = T)

# Clustering
DefaultAssay(allseur_integrated) = 'integrated'
allseur_integrated = ScaleData(allseur_integrated, verbose = FALSE)
allseur_integrated = RunPCA(allseur_integrated, verbose = FALSE)
allseur_integrated = RunUMAP(allseur_integrated, dims = 1:20, reduction = 'pca')

setwd("/home2/gkonop/workdir/03_INTEGRATE_ALL/Non_SPN/CLUSTERING_6")

# Basic Plots
pdf(paste0(pref, "_INTEGRATED_UMAP.pdf"))
DimPlot(allseur_integrated, group.by = 'Species', raster = T, label.size = 7)
dev.off()

allseur_integrated = FindNeighbors(allseur_integrated, dims = 1:20, reduction = 'pca')
allseur_integrated = FindClusters(allseur_integrated, resolution = 2)

pdf(paste0(pref, "_INTEGRATED_CLUSTERS.pdf"))
DimPlot(allseur_integrated, label = T, raster = T, label.size = 7) + NoLegend()
dev.off()

stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_INTEGRATED_STACK_SAMPLES'), wd = 20)
stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'Species', paste0(pref, '_INTEGRATED_STACK_SPECIES'))

# Plot previously identified markers
DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)

pdf(paste0(pref, "_INTEGRATED_DOTPLOT.pdf"), width = 25, height = 10)
DotPlot(allseur_integrated, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "_INTEGRATED_DEPTH.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'black') + NoLegend() +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

saveRDS(allseur_integrated, paste0('/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/CLUSTERING_6/', pref, '_integrated_CLUSTERING6.RDS'))
allseur_integrated <- readRDS(paste0('/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/CLUSTERING_6/', pref, '_integrated_CLUSTERING6.RDS'))
###############################################################################
####
## CLUSTER & FILTER
####
#Cluster_#,Annotation
#0,PDGFD_PTHLH_PVALB-
#1,SST_NPY
#2,PDGFD_PTHLH_PVALB+
#3,TAC3
#4,PDGFD
#5,PDGFD_PTHLH_PVALB-
#6,TAC3
#7,TAC3_TH
#8,SST_NPY
#9,TAC3
#10,TAC3
#11,PDGFD_PTHLH_PVALB-
#12,SST_NPY
#13,TAC3
#14,CHAT
#15,CCK_VIP-
#16,CCK_VIP+
#17,SST_NPY
#18,PDGFD
#19,CHAT
#20,CHAT
#21,LMO3_BAT
#22,PDGFD
#23,LMO3_BAT
#24,PDGFD
#25,LMO3_BAT
#26,FOXP2_EYA2
#27,FOXP2_TSHZ2_BAT
#28,LMO3_BAT
#29,FOXP2_EYA2
#30,CHAT
#31,CCK_VIP-
#32,SST_NPY
#33,CCK_VIP-
#34,CHAT
#35,PDGFD_PTHLH_PVALB+
#36,TH
#37*,DOUBLET


(toremove = c(37))
allseur_integrated = subset(allseur_integrated, subset = seurat_clusters %in% toremove, invert = T)

# Clustering
DefaultAssay(allseur_integrated) = 'integrated'
allseur_integrated = ScaleData(allseur_integrated, verbose = FALSE)
allseur_integrated = RunPCA(allseur_integrated, verbose = FALSE)
allseur_integrated = RunUMAP(allseur_integrated, dims = 1:20, reduction = 'pca')

setwd("/home2/gkonop/workdir/03_INTEGRATE_ALL/Non_SPN/CLUSTERING_7")

# Basic Plots
pdf(paste0(pref, "_INTEGRATED_UMAP.pdf"))
DimPlot(allseur_integrated, group.by = 'Species', raster = T, label.size = 7)
dev.off()

allseur_integrated = FindNeighbors(allseur_integrated, dims = 1:20, reduction = 'pca')
allseur_integrated = FindClusters(allseur_integrated, resolution = 2)

pdf(paste0(pref, "_INTEGRATED_CLUSTERS.pdf"))
DimPlot(allseur_integrated, label = T, raster = T, label.size = 7) + NoLegend()
dev.off()

stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_INTEGRATED_STACK_SAMPLES'), wd = 20)
stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'Species', paste0(pref, '_INTEGRATED_STACK_SPECIES'))

# Plot previously identified markers
DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)

pdf(paste0(pref, "_INTEGRATED_DOTPLOT.pdf"), width = 25, height = 10)
DotPlot(allseur_integrated, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "_INTEGRATED_DEPTH.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'black') + NoLegend() +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

saveRDS(allseur_integrated, paste0('/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/CLUSTERING_7/', pref, '_integrated_CLUSTERING7.RDS'))
allseur_integrated <- readRDS(paste0('/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/CLUSTERING_7/', pref, '_integrated_CLUSTERING7.RDS'))
#############################################################################
####
## CLUSTER & ANNOTATE
####
#0,PDGFD_PTHLH_PVALB-
#1,SST_NPY
#2,PDGFD_PTHLH_PVALB+
#3,TAC3_TH
#4,PDGFD_PTHLH_PVALB+
#5,TAC3_TH
#6,SST_NPY
#7,PDGFD
#8,TAC3_TH
#9,PDGFD_PTHLH_PVALB-
#10,TAC3_TH
#11,TAC3_TH
#12,SST_NPY
#13,TAC3_TH
#14,CCK_VIP-
#15,CCK_VIP+
#16,PDGFD
#17,CHAT
#18,CHAT
#19,CHAT
#20,SST_NPY
#21,FOXP2_EYA2
#22,PDGFD
#23,LMO3_BAT
#24,PDGFD_PTHLH_PVALB-
#25,LMO3_BAT
#26,CHAT
#27,FOXP2_EYA2
#28,FOXP2_TSHZ2_BAT
#29,LMO3_BAT
#30,FOXP2_EYA2
#31,SST_NPY
#32,CCK_VIP-
#33,CCK_VIP-
#34,CHAT
#35,TAC3_TH
#36,TAC3_TH

meta <- allseur_integrated[[]]
mapnames = setNames(c("PDGFD_PTHLH_PVALB-","SST_NPY","PDGFD_PTHLH_PVALB+",
"TAC3_TH","PDGFD_PTHLH_PVALB+","TAC3_TH","SST_NPY","PDGFD","TAC3_TH",
"PDGFD_PTHLH_PVALB-","TAC3_TH","TAC3_TH","SST_NPY","TAC3_TH","CCK_VIP-",
"CCK_VIP+","PDGFD","CHAT","CHAT","CHAT","SST_NPY","FOXP2_EYA2","PDGFD",
"LMO3_BAT","PDGFD_PTHLH_PVALB-","LMO3_BAT","CHAT","FOXP2_EYA2",
"FOXP2_TSHZ2_BAT","LMO3_BAT","FOXP2_EYA2","SST_NPY","CCK_VIP-",
"CCK_VIP-","CHAT","TAC3_TH","TAC3_TH"),c(0:36))

meta$newannot <- mapnames[meta$seurat_clusters]
meta$Species <- factor(meta$Species, levels = c('Human', 'Chimp', 'Macaque', 'Marmoset', 'Mouse', 'Bat', 'Ferret'))

setwd("/endosome/work/Neuroinformatics_Core/gkonop/03_INTEGRATE_ALL/Non_SPN/ANNOTATION/")

# Create unique id (cb_sample)
meta$cb = gsub('_.*|-.*', '', rownames(meta))
meta$id = paste0(meta$cb, '_', meta$orig.ident)

# Define the unique list of annotations without duplicates
annotation_levels <- c("PDGFD","PDGFD_PTHLH_PVALB-","PDGFD_PTHLH_PVALB+",
"CHAT","SST_NPY","CCK_VIP-","CCK_VIP+","TAC3_TH",
"FOXP2_EYA2","FOXP2_TSHZ2_BAT","LMO3_BAT")

# Assign these unique annotations as factor levels to `newannot`
meta$newannot <- factor(meta$newannot, levels = annotation_levels)

# Reassign meta data
allseur_integrated@meta.data = meta
################################
# convert the other TH cells to TAC3 as in the dotplot they show a lot more TAC3
ferret_sub <- subset(allseur_integrated, Species == "Ferret")
DefaultAssay(ferret_sub) <- "RNA"
ferret_sub <- NormalizeData(ferret_sub)
DotPlot(ferret_sub, features = c("TH", "TAC3", "PDGFD", "CHAT"), group.by = "seurat_clusters")

marm_sub <- subset(allseur_integrated, Species == "Marmoset")
DefaultAssay(marm_sub) <- "RNA"
marm_sub <- NormalizeData(marm_sub)
DotPlot(marm_sub, features = c("TH", "TAC3", "PDGFD", "CHAT"), group.by = "seurat_clusters")

DefaultAssay(allseur_integrated) <- "RNA"
allseur_integrated <- NormalizeData(allseur_integrated)
DotPlot(allseur_integrated, features = c("TH", "TAC3", "PDGFD", "CHAT"), group.by = "seurat_clusters")
####################################
# PLOTS
pdf(paste0(pref, "_Clusters_Annotated.pdf"))
DimPlot(allseur_integrated, group.by = 'newannot', label = T, raster = T, repel = T, label.size = 7) +
NoLegend() +
ggtitle('')
dev.off()

pdf(paste0(pref, "_Depth_Gene_Annotated.pdf"))
ggboxplot(allseur_integrated@meta.data, x = 'newannot', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_Depth_UMI_Annotated.pdf"))
ggboxplot(allseur_integrated@meta.data, x = 'newannot', y = 'nCount_RNA', ylim = c(0,80000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()


options(warn = 0) # prevent warnings to become errors
pdf(paste0(pref, "_NuclearFraction_Annotated.pdf"))
ggboxplot(allseur_integrated@meta.data, x = 'newannot', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()


DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)
pdf(paste0(pref, "_MarkersDotPlot.pdf"), width = 25, height = 10)
DotPlot(allseur_integrated, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

stackedbarplot(allseur_integrated@meta.data, groupfill = 'newannot', groupx = 'Species', fn = paste0(pref, '_Stacked_Annotated_Species'))

stackedbarplot(allseur_integrated@meta.data, groupfill = 'Species', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Species2'))

stackedbarplot(allseur_integrated@meta.data, groupfill = 'newannot', groupx = 'orig.ident', fn = paste0(pref, '_Stacked_Annotated_Samples'))

stackedbarplot(allseur_integrated@meta.data, groupfill = 'orig.ident', groupx = 'newannot', fn = paste0(pref, '_Stacked_AnnotatedSamples2'), wd = 20)

stackedbarplot(allseur_integrated@meta.data, groupfill = 'Tissue', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Tissue'))


# Primate vs bat mouse vs ferret
meta = allseur_integrated@meta.data
meta$Species_Broad = 'Primate'
meta[meta$Species == 'Bat', 'Species_Broad'] = 'Bat'
meta[meta$Species == 'Mouse', 'Species_Broad'] = 'Mouse'
meta[meta$Species == 'Ferret', 'Species_Broad'] = 'Ferret'

stackedbarplot(meta, groupfill = 'Species_Broad', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Species_Broad'))

# stacked bar plot of cell types for only primates
sub = subset(meta, subset = Species_Broad == c('Primate'))
table(sub$Species)
stackedbarplot(sub, groupfill = 'Species', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Species_Primates'))

# save
saveRDS(allseur_integrated, file = paste0("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/ANNOTATED/", pref, "_ANNOTATED.RDS"))
allseur_integrated <- readRDS(file = "/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/ANNOTATED/GB_Non_SPN_AllSpecies_AllTissues_human_pr_coding_orthologs_ANNOTATED.RDS")

DotPlot(allseur_integrated, features = c("TH", "TAC3"), split.by = "Species", group.by = "newannot", cols = c('blue', 'violetred', 'deepskyblue2', 'darkorchid', 'chocolate2', 'orange', 'purple'))


# Extract data from Seurat DotPlot
dotplot_data <- DotPlot(
  object = allseur_integrated,
  features = c("TH", "TAC3"),
  group.by = "newannot",
split.by = "Species",
  cols = c('blue', 'violetred', 'deepskyblue2', 'darkorchid', 'chocolate2', 'orange', 'purple')
)$data

library(stringr)

# Extract the last part of the `id` string as `Species`
dotplot_data$Species <- sapply(strsplit(x=as.character(dotplot_data$id), split= "_"), function(x) tail(x, 1))

# Extract cell types by splitting 'id' (ensure it's treated as character) and removing the last part (species)
dotplot_data$CellType <- sapply(strsplit(as.character(dotplot_data$id), "_"), function(x) paste(head(x, -1), collapse = "_"))

# Create a plot with facet_wrap
ggplot(dotplot_data, aes(x = CellType, y = features.plot, size = avg.exp, color = pct.exp)) +
  geom_point() +
  scale_color_gradientn(colors = c('blue', 'violetred', 'deepskyblue2', 'darkorchid', 'chocolate2', 'orange', 'purple')) +
  facet_wrap(~ Species, nrow = 1, scales = "free_x") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    strip.text = element_text(size = 12)
  )


dim(allseur_integrated) #[1] 14893 20927
table(allseur_integrated@meta.data$newannot, allseur_integrated@meta.data$Species)
                    
                     Human Chimp Macaque Marmoset Mouse  Bat Ferret
  SST_NPY              254   413     668      468   508 1026    127
  PDGFD_PTHLH_PVALB-   339   318    1221      279   243  163    137
  PDGFD                313   354     460      465    58  264     81
  PDGFD_PTHLH_PVALB+   102    97     288      197   458  556    159
  TAC3                 574   711    1496      794    69  237     52
  TH                    41    58     186      306   216  227      4
  CCK_VIP-             221    89     205      111    37  300      7
  CCK_VIP+             208   115     136        3     4   94     17
  CHAT                 152   304     410      184   296  536     59
  LMO3_BAT             146    48     181       12    49  754      6
  FOXP2_TSHZ2_BAT       27     5       7        9     0  229      0
  FOXP2_EYA2            46   144      67        1     7   97    144
  RSPO2                245    32     137       10     2   46      1

##############################################################################
## Troubleshooting
all.equal(rownames(allseur_integrated@meta.data), colnames(allseur_integrated))
# after manipulating metadata, the number of cells in the count matrix and the number of cells in the metadata do not match! Here is the fix:
# check dims
dim(allseur_integrated@meta.data) #[1] 21076    40
dim(allseur_integrated@assays$RNA@counts) #[1] 14893 21216
# find the missing cells
setdiff(rownames(allseur_integrated@meta.data), colnames(allseur_integrated))
setdiff(colnames(allseur_integrated), rownames(allseur_integrated@meta.data))
# remove those cells from the meta.data
allseur_integrated@meta.data <- allseur_integrated@meta.data[colnames(allseur_integrated), , drop = FALSE]
# subset the seurat object as well
allseur_integrated <- subset(allseur_integrated, cells = rownames(allseur_integrated@meta.data))
# check
all.equal(rownames(allseur_integrated@meta.data), colnames(allseur_integrated))
###########################################################################
## TAC3-TH sub-clustering
tac3_th <- subset(allseur_integrated, subset = newannot == "TAC3_TH")


# Clustering
DefaultAssay(tac3_th) = 'integrated'
tac3_th = ScaleData(tac3_th, verbose = FALSE)
tac3_th = RunPCA(tac3_th, verbose = FALSE)
tac3_th = RunUMAP(tac3_th, dims = 1:20, reduction = 'pca')

setwd("/home2/gkonop/workdir/03_INTEGRATE_ALL/Non_SPN/CLUSTERING_7")
pref <- "TAC3_TH_subsclustering_Non_SPN_AllSpecies_AllTissues_human_pr_coding_orthologs"

# Basic Plots
tac3_th = FindNeighbors(tac3_th, dims = 1:10, reduction = 'pca')
tac3_th = FindClusters(tac3_th, resolution = 1)

pdf(paste0(pref, "_INTEGRATED_CLUSTERS.pdf"))
DimPlot(tac3_th, label = T, raster = T, label.size = 7) + NoLegend()
dev.off()

stackedbarplot(tac3_th[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_INTEGRATED_STACK_SAMPLES'), wd = 20)
stackedbarplot(tac3_th[[]], groupx = 'seurat_clusters', groupfill = 'Species', paste0(pref, '_INTEGRATED_STACK_SPECIES'))

# plot cluster 15 separately
tac3_th_cl15 <- subset(tac3_th, seurat_clusters == 15)
stackedbarplot(tac3_th_cl15[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_cluster15_INTEGRATED_STACK_SAMPLES'), wd = 20)
	
# Plot previously identified markers
DefaultAssay(tac3_th) = 'RNA'
tac3_th = NormalizeData(tac3_th)

pdf(paste0(pref, "_INTEGRATED_DOTPLOT.pdf"), width = 25, height = 10)
DotPlot(tac3_th, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "_INTEGRATED_DEPTH.pdf"))
ggboxplot(tac3_th[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'black') + NoLegend() +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(tac3_th[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# save
saveRDS(tac3_th, paste0('/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/CLUSTERING_7/', pref, '_integrated_ANNOTATED.RDS'))
tac3_th <- readRDS("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/CLUSTERING_7/TAC3_TH_subsclustering_Non_SPN_AllSpecies_AllTissues_human_pr_coding_orthologs_integrated_CLUSTERING7.RDS")

allseur_integrated <- readRDS("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/ANNOTATED/GB_Non_SPN_AllSpecies_AllTissues_human_pr_coding_orthologs_ANNOTATED.RDS")
###############################################################
####
## CLUSTER & FILTER
####
#Cluster_#,Annotation
#0,TAC3
#1,TAC3
#2,TH
#3,TAC3
#4,TAC3
#5,TAC3
#6,TAC3
#7,TAC3
#8,TAC3
#9,TAC3
#10,TAC3
#11,TAC3
#12,TAC3
#13,TAC3
#14,TH
#15*,TH

tac3_th <- subset(tac3_th, subset = seurat_clusters == 15, invert = T)
setwd("/home2/gkonop/workdir/03_INTEGRATE_ALL/Non_SPN/CLUSTERING_8")
pref <- "TAC3_TH_subsclustering_Non_SPN_AllSpecies_AllTissues_human_pr_coding_orthologs"

# Clustering
DefaultAssay(tac3_th) = 'integrated'
tac3_th = ScaleData(tac3_th, verbose = FALSE)
tac3_th = RunPCA(tac3_th, verbose = FALSE)
tac3_th = RunUMAP(tac3_th, dims = 1:10, reduction = 'pca')


# Basic Plots
tac3_th = FindNeighbors(tac3_th, dims = 1:10, reduction = 'pca')
tac3_th = FindClusters(tac3_th, resolution = 1)

pdf(paste0(pref, "_INTEGRATED_CLUSTERS.pdf"))
DimPlot(tac3_th, label = T, raster = T, label.size = 7) + NoLegend()
dev.off()

stackedbarplot(tac3_th[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_INTEGRATED_STACK_SAMPLES'), wd = 20)
stackedbarplot(tac3_th[[]], groupx = 'seurat_clusters', groupfill = 'Species', paste0(pref, '_INTEGRATED_STACK_SPECIES'))

# Plot previously identified markers
DefaultAssay(tac3_th) = 'RNA'
tac3_th = NormalizeData(tac3_th)

pdf(paste0(pref, "_INTEGRATED_DOTPLOT.pdf"), width = 25, height = 10)
DotPlot(tac3_th, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "_INTEGRATED_DEPTH.pdf"))
ggboxplot(tac3_th[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'black') + NoLegend() +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(tac3_th[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()
############################################################
####
## CLUSTER & ANNOTATE
####
#Cluster_#,Annotation
#0,TAC3
#1,TAC3
#2,TH
#3,TAC3
#4,TAC3
#5,TAC3
#6,TAC3
#7,TAC3
#8,TAC3
#9,TAC3
#10,TAC3
#11,TAC3
#12,TH
#13,TAC3
#14,TAC3

(mapnames = setNames(c("TAC3","TAC3","TH","TAC3","TAC3",
"TAC3","TAC3","TAC3","TAC3","TAC3","TAC3","TAC3",
"TH","TAC3","TAC3"),c(0:14)))

tac3_th[["newannot"]] = mapnames[tac3_th[["seurat_clusters"]][,1]]
Idents(tac3_th) = tac3_th$newannot

setwd("/home2/gkonop/workdir/03_INTEGRATE_ALL/Non_SPN/ANNOTATION")
pref <- "TAC3_TH_subsclustering_Annotated_Non_SPN_AllSpecies_AllTissues_human_pr_coding_orthologs"

pdf(paste0(pref, "_INTEGRATED_CLUSTERS.pdf"))
DimPlot(tac3_th, label = T, raster = T, label.size = 7) + NoLegend()
dev.off()

stackedbarplot(tac3_th[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_INTEGRATED_STACK_SAMPLES'), wd = 20)
stackedbarplot(tac3_th[[]], groupx = 'seurat_clusters', groupfill = 'Species', paste0(pref, '_INTEGRATED_STACK_SPECIES'))

# Plot previously identified markers
DefaultAssay(tac3_th) = 'RNA'
tac3_th = NormalizeData(tac3_th)

# prevent warnings to become errors
options(warn = -1)
pdf(paste0(pref, "_INTEGRATED_DOTPLOT.pdf"), width = 25, height = 10)
DotPlot(tac3_th, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "_INTEGRATED_DEPTH.pdf"))
ggboxplot(tac3_th[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'black') + NoLegend() +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(tac3_th[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()
######################################################
## add the tac3 th labels to the interneuron seurat obj
allseur_integrated$newannot <- ifelse(names(allseur_integrated$id) %in% names(tac3_th$id), as.character(tac3_th$newannot), as.character(allseur_integrated$newannot))

# remove the 49 cells which were also removed while subclustering of tac3_th
allseur_integrated <- subset(allseur_integrated, newannot == "TAC3_TH", invert = T)

table(allseur_integrated$newannot)

          CCK_VIP-           CCK_VIP+               CHAT         FOXP2_EYA2 
               976                584               1947                961 
   FOZP2_TSHZ2_BAT           LMO3_BAT              PDGFD PDGFD_PTHLH_PVALB- 
               274                958               1813               2595 
PDGFD_PTHLH_PVALB+            SST_NPY               TAC3                 TH 
              2191               3460               4513                511

Idents(allseur_integrated) <- allseur_integrated$newannot
DimPlot(allseur_integrated, label = T, raster = T, label.size = 7) + NoLegend()

pref <- "Non_SPN_AllSpecies_AllTissues_human_pr_coding_orthologs_TAC3_TH_clustered"
# PLOTS
pdf(paste0(pref, "_Clusters_Annotated.pdf"))
DimPlot(allseur_integrated, group.by = 'newannot', label = T, raster = T, repel = T, label.size = 7) +
NoLegend() +
ggtitle('')
dev.off()

pdf(paste0(pref, "_Depth_Gene_Annotated.pdf"))
ggboxplot(allseur_integrated@meta.data, x = 'newannot', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_Depth_UMI_Annotated.pdf"))
ggboxplot(allseur_integrated@meta.data, x = 'newannot', y = 'nCount_RNA', ylim = c(0,80000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

options(warn = 0) # prevent warnings to become errors
pdf(paste0(pref, "_NuclearFraction_Annotated.pdf"))
ggboxplot(allseur_integrated@meta.data, x = 'newannot', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)
pdf(paste0(pref, "_MarkersDotPlot.pdf"), width = 25, height = 10)
DotPlot(allseur_integrated, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

stackedbarplot(allseur_integrated@meta.data, groupfill = 'newannot', groupx = 'Species', fn = paste0(pref, '_Stacked_Annotated_Species'))

# re-order species
allseur_integrated$Species <- factor(allseur_integrated$Species, levels = rev(c("Human", "Chimp", "Macaque", 
"Marmoset", "Mouse", "Bat", "Ferret")))
stackedbarplot(allseur_integrated@meta.data, groupfill = 'Species', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Species2'))

stackedbarplot(allseur_integrated@meta.data, groupfill = 'newannot', groupx = 'orig.ident', fn = paste0(pref, '_Stacked_Annotated_Samples'))

stackedbarplot(allseur_integrated@meta.data, groupfill = 'orig.ident', groupx = 'newannot', fn = paste0(pref, '_Stacked_AnnotatedSamples2'), wd = 20)

stackedbarplot(allseur_integrated@meta.data, groupfill = 'Tissue', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Tissue'))


# Primate vs bat mouse vs ferret
meta = allseur_integrated@meta.data
meta$Species_Broad = 'Primate'
meta[meta$Species == 'Bat', 'Species_Broad'] = 'Bat'
meta[meta$Species == 'Mouse', 'Species_Broad'] = 'Mouse'
meta[meta$Species == 'Ferret', 'Species_Broad'] = 'Ferret'

# re-order
meta$Species_Broad <- factor(meta$Species_Broad, levels = rev(c("Primate", "Mouse", "Bat", "Ferret")))

stackedbarplot(meta, groupfill = 'Species_Broad', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Species_Broad'))

# stacked bar plot of cell types for only primates
sub = subset(meta, subset = Species_Broad == c('Primate'))
table(sub$Species)
# re-order
sub$Species <- factor(sub$Species, levels = rev(c("Human", "Chimp", "Macaque", "Marmoset")))
stackedbarplot(sub, groupfill = 'Species', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Species_Primates'))

#save
saveRDS(tac3_th, paste0('/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/ANNOTATED/tac3_th_integrated_ANNOTATED.RDS'))


# generate unique sample id
allseur_integrated$id2 <- paste0(allseur_integrated$orig.ident, "_", allseur_integrated$Tissue)
saveRDS(allseur_integrated, "/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/ANNOTATED/GB_Non_SPN_AllSpecies_AllTissues_tac3_th_subclustered_human_pr_coding_orthologs_ANNOTATED.RDS")
allseur_integrated <- readRDS("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/ANNOTATED/GB_Non_SPN_AllSpecies_AllTissues_tac3_th_subclustered_human_pr_coding_orthologs_ANNOTATED.RDS")
###########################################################################
## check for TH cells in old seu obj (with non-coding genes)
# old seu obj
old_interneurons <- readRDS("/home2/gkonop/workdir/03_INTEGRATE_ALL/Non_SPN/ANNOTATION/GB_Non_SPN_AllSpecies_integrated_annotated_resolution_2.RDS")

# check the species distribution
stackedbarplot(old_interneurons@meta.data, groupfill = 'Species', groupx = 'newannot', fn = paste0("test", '_Stacked_Annotated_Species2'))

## create unique id for each cell
old_interneurons$cellbarc = gsub('_.*|-.*', '', rownames(old_interneurons[[]]))

# add id
old_interneurons$id2 <- paste0(old_interneurons$cellbarc, "_", old_interneurons$orig.ident)

# subset TH cells
old_th_cells <- subset(old_interneurons, subset = newannot == "TH")

# read the first clustering seu obj
first_non_spn <- readRDS("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/Integrated_rpca_NonSPN_AllSpecies_ALLTissues_human_pr_coding_orthologs.RDS")

# Create unique id (cb_sample)
first_non_spn$cellbarc = gsub('_.*|-.*', '', rownames(first_non_spn[[]]))
first_non_spn$id2 = paste0(first_non_spn$cellbarc, '_', first_non_spn$orig.ident)
#########################################################
### TROUBLESHOOTING
## The newannot labels are not correctly transferred: There are 471 LMO3_BAT macaque cells where there should be only 24!!
# Find matching indices
matched_indices <- match(first_non_spn$id2, allseur_integrated$id)

# Assign new annotations: match returns NA for non-matching values, so assign "Unknown" for them
first_non_spn$newannot <- ifelse(!is.na(matched_indices), allseur_integrated$newannot[matched_indices], "Unknown")

# check
table(first_non_spn$newannot, first_non_spn$Species)

stackedbarplot(first_non_spn[[]], groupx = 'newannot', groupfill = 'Species', paste0(pref, '_First_Non_SPN_INTEGRATED_STACK_SPECIES2'))
##########################################################################
# find old TH cells in new seu obj
# Find matching indices
matched_indices <- match(first_non_spn$id2, old_th_cells$id2)

# Assign new annotation where matches are found, otherwise keep the original annotation
first_non_spn$newannot <- ifelse(!is.na(matched_indices), "TH", first_non_spn$newannot)

# check
stackedbarplot(first_non_spn[[]], groupx = 'newannot', groupfill = 'Species', paste0(pref, '_First_Non_SPN_INTEGRATED_STACK_SPECIES2'))


# remove unknown and re-cluster
second_non_spn <- subset(first_non_spn, newannot == "Unknown", invert = T)
##########################################
# check 
table(second_non_spn$newannot, second_non_spn$Species)


                      Bat Chimp Ferret Human Macaque Marmoset Mouse
  CCK_VIP-            302    90      8   223     205      111    37
  CCK_VIP+             95   116     17   209     140        3     4
  CHAT                537   307     59   152     412      184   296
  FOXP2_EYA2          306   165    151    98     151       12    46
  FOZP2_TSHZ2_BAT     229     4      0    26       6        9     0
  LMO3_BAT            832    21      0    72      24        8     1
  PDGFD               168   314     47   263     615      373    26
  PDGFD_PTHLH_PVALB-  293   350    107   348     999      406    90
  PDGFD_PTHLH_PVALB+  526   105    227   142     353      163   644
  SST_NPY            1025   412    127   253     667      468   508
  TAC3                355   622     58   600    1430      974    97
  TH                  122   162      2    34     284      134   328

### PLOTS with Seurat Clusters
# Clustering
DefaultAssay(second_non_spn) = 'integrated'
second_non_spn = ScaleData(second_non_spn, verbose = FALSE)
second_non_spn = RunPCA(second_non_spn, verbose = FALSE)
second_non_spn = RunUMAP(second_non_spn, dims = 1:20, reduction = 'pca')

# set variables
pref <- "Non_SPN_AllSpecies_AllTissues_human_pr_coding_orthologs_TAC3_TH_clustered_Seurat_Clusters"
setwd("/endosome/work/Neuroinformatics_Core/gkonop/03_INTEGRATE_ALL/Non_SPN/ANNOTATION")

# Basic Plots
pdf(paste0(pref, "_INTEGRATED_UMAP.pdf"))
DimPlot(second_non_spn, group.by = 'Species', raster = T, label.size = 7)
dev.off()

second_non_spn = FindNeighbors(second_non_spn, dims = 1:20, reduction = 'pca')
second_non_spn = FindClusters(second_non_spn, resolution = 2)

pdf(paste0(pref, "_INTEGRATED_CLUSTERS.pdf"))
DimPlot(second_non_spn, label = T, raster = T, label.size = 7) + NoLegend()
dev.off()

stackedbarplot(second_non_spn[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_INTEGRATED_STACK_SAMPLES'), wd = 20)
stackedbarplot(second_non_spn[[]], groupx = 'seurat_clusters', groupfill = 'Species', paste0(pref, '_INTEGRATED_STACK_SPECIES'))

# Plot previously identified markers
DefaultAssay(second_non_spn) = 'RNA'
second_non_spn = NormalizeData(second_non_spn)

pdf(paste0(pref, "_INTEGRATED_DOTPLOT.pdf"), width = 25, height = 10)
DotPlot(second_non_spn, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "_INTEGRATED_DEPTH.pdf"))
ggboxplot(second_non_spn[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'black') + NoLegend() +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(second_non_spn[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()
##########################################
#### ANNOTATE !IGNORE!
#0,PDGFD_PTHLH_PVALB+
#1,PDGFD_PTHLH_PVALB-
#2,PDGFD_PTHLH_PVALB-
#3,TAC3
#4,TAC3
#5,SST_NPY
#6,SST_NPY
#7,TAC3
#8,CHAT
#9,CCK_VIP+
#10,CHAT
#11,TAC3
#12,CCK_VIP-
#13,SST_NPY
#14,PDGFD
#15,LMO3_BAT
#16,PDGFD_PTHLH_PVALB-
#17,LMO3_BAT
#18,LMO3_BAT
#19,FOXP2_EYA2
#20,FOXP2_TSHZ2
#21,FOXP2_EYA2
#22,CCK_VIP-
#23,CHAT
#24,TH

mapnames <- setNames(c("PDGFD_PTHLH_PVALB+","PDGFD_PTHLH_PVALB-","PDGFD_PTHLH_PVALB-",
"TAC3","TAC3","SST_NPY","SST_NPY","TAC3","CHAT","CCK_VIP+","CHAT","TAC3",
"CCK_VIP-","SST_NPY","PDGFD","LMO3_BAT","PDGFD_PTHLH_PVALB-","LMO3_BAT","LMO3_BAT",
"FOXP2_EYA2","FOXP2_TSHZ2","FOXP2_EYA2","CCK_VIP-","CHAT","TH"), 0:24)

second_non_spn[["newannot3"]] = mapnames[second_non_spn[["seurat_clusters"]][,1]]
Idents(second_non_spn) = second_non_spn$newannot3

setwd("/home2/gkonop/workdir/03_INTEGRATE_ALL/Non_SPN/ANNOTATION")
pref <- "Non_SPN_AllSpecies_AllTissues_human_pr_coding_orthologs_TAC3_TH_clustered_Annotated_OnItsOwn"

pdf(paste0(pref, "_INTEGRATED_CLUSTERS.pdf"))
DimPlot(second_non_spn, label = T, raster = T, label.size = 7) + NoLegend()
dev.off()

stackedbarplot(second_non_spn[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_INTEGRATED_STACK_SAMPLES'), wd = 20)
stackedbarplot(second_non_spn[[]], groupx = 'seurat_clusters', groupfill = 'Species', paste0(pref, '_INTEGRATED_STACK_SPECIES'))

# Plot previously identified markers
DefaultAssay(second_non_spn) = 'RNA'
second_non_spn = NormalizeData(second_non_spn)

# prevent warnings to become errors
options(warn = -1)
pdf(paste0(pref, "_INTEGRATED_DOTPLOT.pdf"), width = 25, height = 10)
DotPlot(second_non_spn, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "_INTEGRATED_DEPTH.pdf"))
ggboxplot(second_non_spn[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'black') + NoLegend() +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(second_non_spn[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()
#################################################################
# PLOTS with newannot
# set identity
Idents(second_non_spn) <- second_non_spn$newannot

# fix typo
second_non_spn$newannot <- gsub("FOZP2_TSHZ2_BAT", "FOXP2_TSHZ2_BAT", second_non_spn$newannot)

# re-order species
second_non_spn$Species <- factor(second_non_spn$Species, levels = rev(c("Human", "Chimp", "Macaque", 
"Marmoset", "Mouse", "Bat", "Ferret")))

# re-order newannot
# Define the unique list of annotations without duplicates
annotation_levels <- c("PDGFD","PDGFD_PTHLH_PVALB-","PDGFD_PTHLH_PVALB+",
"CHAT","SST_NPY","CCK_VIP-","CCK_VIP+","TAC3", "TH",
"FOXP2_EYA2","FOXP2_TSHZ2_BAT","LMO3_BAT")

# Assign these unique annotations as factor levels to `newannot`
second_non_spn$newannot <- factor(second_non_spn$newannot, levels = annotation_levels)

pref <- "Non_SPN_AllSpecies_AllTissues_human_pr_coding_orthologs_TAC3_TH_clustered_AnnoatedByAllseurIntegrated_and_OldSeuratObj"

pdf(paste0(pref, "_Clusters_Annotated.pdf"))
DimPlot(second_non_spn, group.by = 'newannot', label = T, raster = T, repel = T, label.size = 7) +
NoLegend() +
ggtitle('')
dev.off()

pdf(paste0(pref, "_Depth_Gene_Annotated.pdf"))
ggboxplot(second_non_spn@meta.data, x = 'newannot', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_Depth_UMI_Annotated.pdf"))
ggboxplot(second_non_spn@meta.data, x = 'newannot', y = 'nCount_RNA', ylim = c(0,80000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

options(warn = 0) # prevent warnings to become errors
pdf(paste0(pref, "_NuclearFraction_Annotated.pdf"))
ggboxplot(second_non_spn@meta.data, x = 'newannot', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# normalize gene exp
DefaultAssay(second_non_spn) = 'RNA'
second_non_spn = NormalizeData(second_non_spn)

# plot gene exp
pdf(paste0(pref, "_MarkersDotPlot.pdf"), width = 25, height = 10)
DotPlot(second_non_spn, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

### plot developmental otigin markers
interesting_marks <- c("NKX2-1", "LHX6", "FOXP2", "PAX6","MEIS2", "SCGN","SP8", "NR2F1","NR2F2", "ADARB2", "LAMP5", "CXCL14", "MAF","CRABP1", "LHX8", "TAC2", "TAC3", "TSHZ1", "TSHZ2", "MKI67")




devo_marks <- c("SATB1", "PROX1","SHH","GBX2 ", "OTX","CNTNAP2","MAF","LHX8", "NKX2-1", "LHX6", "ZIC1", "ZIC4", "FOXP2", "PAX6","MEIS2","ISL1","EBF1", "ZNF503","NR2F1","NR2F2","ADARB2", "LAMP5", "CXCL14", "SIX3", "HTT", "CHD7", "DLX2", "ASCL1", "HES1", "NUSAP1","DRD1", "DRD2", "OTOF", "CASZ1", "HTR7")

SPN_marks <- c("DRD1", "DRD2", "OTOF", "CASZ1", "HTR7", "ADORA2 ", "ISL1", "TRANK1", "AGTR1A", 
"PCDH8", "PPP1R1B")
#### markers
## ISL1: DRD1 devo marker
## SIX3: DRD2 devo marker
## HTT: huntintin
## CHD7: olfactory bulb interneuron
## DLX2 : interneuron precursor
## ASCL1: interneuron precursor 
## HES1 : neuron progenitor marker
## NUSAP1 : neuron progenitor marker
## MAF: MGE
## CNTNAP2: MGE
## ZIC1 : ventral MGE
## ZIC4 : ventral MGE
#MGE-derived interneuron precursors (MGE1 and MGE2) highly expressed LHX6, LHX8 and NKX2-1; CGE-derived interneuron precursors abundantly expressed NR2F1 (also known as COUP-TFI) and NR2F2 (also known as COUP-TFII); and LGE-derived interneuron precursors (LGE1, LGE2 and LGE3) mainly expressed EBF1, ISL1 and ZNF503 (Extended Data Fig. 2d?f). https://www.nature.com/articles/s41593-021-00940-3?fromPaywallRec=false
# For instance, loss of Satb1 results in distinct deficits in SST+ versus PV+ interneurons [58], and Prox1 removal differentially affects distinct CGE-derived interneuron subtype https://pmc.ncbi.nlm.nih.gov/articles/PMC5699457/
#in brain sections, and we detected higher expression levels of MEIS2 and NR2F1 in the dorsal MGE and higher expression levels of LHX8 and ZIC4 in the ventral MGE (Extended Data Fig. 7d,e).https://www.nature.com/articles/s41593-021-00940-3#Fig2

second_non_spn <- NormalizeData(second_non_spn) 

# plot gene exp
pdf(paste0(pref, "_SPNMarkersDotPlot.pdf"), width = 25, height = 10)
DotPlot(second_non_spn, features = SPN_marks,cols = c("blue", "red")) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()


# plot gene exp
pdf(paste0(pref, "_InterestingMarkersDotPlot.pdf"), width = 25, height = 10)
DotPlot(second_non_spn, features = devo_marks,cols = c("blue", "red")) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()


# plot gene exp
pdf(paste0(pref, "_DevoMarkersDotPlot.pdf"), width = 25, height = 10)
DotPlot(second_non_spn, features = devo_marks,cols = c("blue", "red")) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "_DevoMarkersDotPlot.pdf"), width = 25, height = 10)
bat_filt <- NormalizeData(bat_filt) 
DotPlot(bat_filt, features = devo_marks) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)


stackedbarplot(second_non_spn@meta.data, groupfill = 'newannot', groupx = 'Species', fn = paste0(pref, '_Stacked_Annotated_Species'))

stackedbarplot(second_non_spn@meta.data, groupfill = 'Species', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Species2'))

stackedbarplot(second_non_spn@meta.data, groupfill = 'newannot', groupx = 'orig.ident', fn = paste0(pref, '_Stacked_Annotated_Samples'))

stackedbarplot(second_non_spn@meta.data, groupfill = 'orig.ident', groupx = 'newannot', fn = paste0(pref, '_Stacked_AnnotatedSamples2'), wd = 20)

stackedbarplot(second_non_spn@meta.data, groupfill = 'Tissue', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Tissue'))


# Primate vs bat mouse vs ferret
meta = second_non_spn@meta.data
meta$Species_Broad = 'Primate'
meta[meta$Species == 'Bat', 'Species_Broad'] = 'Bat'
meta[meta$Species == 'Mouse', 'Species_Broad'] = 'Mouse'
meta[meta$Species == 'Ferret', 'Species_Broad'] = 'Ferret'

# re-order
meta$Species_Broad <- factor(meta$Species_Broad, levels = rev(c("Primate", "Mouse", "Bat", "Ferret")))

stackedbarplot(meta, groupfill = 'Species_Broad', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Species_Broad'))

# stacked bar plot of cell types for only primates
sub = subset(meta, subset = Species_Broad == c('Primate'))
table(sub$Species)
# re-order
sub$Species <- factor(sub$Species, levels = rev(c("Human", "Chimp", "Macaque", "Marmoset")))
stackedbarplot(sub, groupfill = 'Species', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Species_Primates'))

# save
saveRDS(second_non_spn, "/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/ANNOTATED/GB_Paper_Non_SPN_AllSpecies_AllTissues_tac3_th_subclustered_th_recovered_human_pr_coding_orthologs_ANNOTATED.RDS")
second_non_spn <- readRDS("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/Non_SPN/ANNOTATED/GB_Paper_Non_SPN_AllSpecies_AllTissues_tac3_th_subclustered_th_recovered_human_pr_coding_orthologs_ANNOTATED.RDS")
###########################################################################
## dendogram
# change default assay to normalized assay
DefaultAssay(second_non_spn) = 'integrated'

# build the dendogram of hiearchical clustering
second_non_spn <- BuildClusterTree(object = second_non_spn) # Calculates the average distance of the PCs for average cell of each cell type (identity class).
# check
str(second_non_spn@tools)

# plot
pdf(paste0(pref,'NonSPN_AllSpecies_Cluster_Tree.pdf'))
PlotClusterTree(object = second_non_spn)
dev.off()
##################################################
#### Heatmap 
# normalize the count matrix
DefaultAssay(second_non_spn) <- "RNA"
second_non_spn <- NormalizeData(second_non_spn)
# centralize the normalized counts 
second_non_spn <- ScaleData(second_non_spn)

# select genes to be mapped
marker_genes = c('CHAT', 'CPNE4',
		 'SST', 'NPY', 'NOS1', 'DACH1', 'LHX6',
		 'TAC3', 'TMEM163', 'PTPRK', 'CALB2', 
		 'TRHDE', 'TH',
		 'LMO3', 'ELAVL4', 'TAFA1', 'SASH1', 'DCC',
		 'FOXP2', 'MEIS2', 'TSHZ2', 'EYA2',  
                 'PDGFD', 'PTHLH', 'PVALB', 
                 'CCK','ADARB2', 'VIP')

# re-order newannot
annotation_levels <- c("CHAT","SST_NPY","TAC3", "TH","LMO3_BAT",
"FOXP2_EYA2","FOXP2_TSHZ2_BAT","PDGFD","PDGFD_PTHLH_PVALB-","PDGFD_PTHLH_PVALB+",
"CCK_VIP-","CCK_VIP+")

# Assign these unique annotations as factor levels to `newannot`
second_non_spn$newannot <- factor(second_non_spn$newannot, levels = annotation_levels)

library(viridis)
# Create the heatmap using DoHeatmap
pdf("Interneuron_Heatmap.pdf")
DoHeatmap(second_non_spn, 
          features = marker_genes,            # Genes to display
          group.by = "newannot",              # Grouping by 'newannot' column                       
          size = 3,                           # Adjust size of text on the heatmap
          disp.max = 3,                       # Maximum value for color scaling
          disp.min = -3,                      # Minimum value for color scaling
          label = TRUE,                       # Show labels on the heatmap
)  + scale_fill_viridis()
dev.off()
######################################################
### test the cell type proportions across species
### perform stat test to compare the cell type proportions
meta <- second_non_spn[[]]

# Calculate the size of each 'newannot' cluster for each combination of 'orig.ident', 'Tissue', and 'Species'
df <- meta %>%
  group_by(orig.ident, Tissue, Species) %>%
  summarize(
    totsize = n(), 
    chat_size = sum(newannot == "CHAT"),  
    sst_npy_size = sum(newannot == "SST_NPY"),  
    tac3_size = sum(newannot == "TAC3"),  
    th_size = sum(newannot == "TH"), 
    lmo3_bat_size = sum(newannot == "LMO3_BAT"),  
    foxp2_eya2_size = sum(newannot == "FOXP2_EYA2"),  
    foxp2_tshz2_bat_size = sum(newannot == "FOXP2_TSHZ2_BAT"), 
    pdgfd_size = sum(newannot == "PDGFD"), 
    pdgfd_pthlh_pvalb_minus_size = sum(newannot == "PDGFD_PTHLH_PVALB-"), 
    pdgfd_pthlh_pvalb_plus_size = sum(newannot == "PDGFD_PTHLH_PVALB+"), 
    cck_vip_minus_size = sum(newannot == "CCK_VIP-"), 
    cck_vip_plus_size = sum(newannot == "CCK_VIP+"), 
    .groups = 'drop'  # To drop the grouping after summarization
  )

# Compute ratios
df <- df %>%
  mutate(
    chat_size_to_allInterneu = chat_size / totsize,
    sst_npy_size_to_allInterneu = sst_npy_size / totsize,
    tac3_size_to_allInterneu = tac3_size / totsize,
    th_size_to_allInterneu = th_size / totsize,
    lmo3_bat_size_to_allInterneu = lmo3_bat_size / totsize,
    foxp2_eya2_size_to_allInterneu = foxp2_eya2_size / totsize,
    foxp2_tshz2_bat_size_to_allInterneu = foxp2_tshz2_bat_size / totsize,
    pdgfd_size_to_allInterneu = pdgfd_size / totsize,
    pdgfd_pthlh_pvalb_minus_size_to_allInterneu = pdgfd_pthlh_pvalb_minus_size / totsize,
    pdgfd_pthlh_pvalb_plus_size_to_allInterneu = pdgfd_pthlh_pvalb_plus_size / totsize,
    cck_vip_minus_size_to_allInterneu = cck_vip_minus_size / totsize,
    cck_vip_plus_size_to_allInterneu = cck_vip_plus_size / totsize
  )

df$Species <- factor(df$Species, levels = c('Bat', 'Ferret', 'Mouse', 'Marmoset', 'Macaque', 'Chimp', 'Human'))
###################################################################3
### save in a table
write.csv(df, file = "/home2/gkonop/workdir/03_INTEGRATE_ALL/Non_SPN/ANNOTATION/SuppTable3_interneuron_props.csv", quote = F, row.names = F)

### generate the proportions primate vs non-primate and save in the same table
meta = second_non_spn@meta.data
meta$Species_Broad = 'Primate'
meta[meta$Species %in% c('Bat','Ferret','Mouse'), 'Species_Broad'] = 'Non_primate'

# Calculate the size of each 'newannot' cluster for each combination of 'orig.ident', 'Tissue', and 'Species'
df2 <- meta %>%
  group_by(Tissue, Species_Broad) %>%
  summarize(
    totsize = n(), 
    chat_size = sum(newannot == "CHAT"),  
    sst_npy_size = sum(newannot == "SST_NPY"),  
    tac3_size = sum(newannot == "TAC3"),  
    th_size = sum(newannot == "TH"), 
    lmo3_bat_size = sum(newannot == "LMO3_BAT"),  
    foxp2_eya2_size = sum(newannot == "FOXP2_EYA2"),  
    foxp2_tshz2_bat_size = sum(newannot == "FOXP2_TSHZ2_BAT"), 
    pdgfd_size = sum(newannot == "PDGFD"), 
    pdgfd_pthlh_pvalb_minus_size = sum(newannot == "PDGFD_PTHLH_PVALB-"), 
    pdgfd_pthlh_pvalb_plus_size = sum(newannot == "PDGFD_PTHLH_PVALB+"), 
    cck_vip_minus_size = sum(newannot == "CCK_VIP-"), 
    cck_vip_plus_size = sum(newannot == "CCK_VIP+"), 
    .groups = 'drop'  # To drop the grouping after summarization
  )

# Compute ratios
df2 <- df2 %>%
  mutate(
    chat_size_to_allInterneu = chat_size / totsize,
    sst_npy_size_to_allInterneu = sst_npy_size / totsize,
    tac3_size_to_allInterneu = tac3_size / totsize,
    th_size_to_allInterneu = th_size / totsize,
    lmo3_bat_size_to_allInterneu = lmo3_bat_size / totsize,
    foxp2_eya2_size_to_allInterneu = foxp2_eya2_size / totsize,
    foxp2_tshz2_bat_size_to_allInterneu = foxp2_tshz2_bat_size / totsize,
    pdgfd_size_to_allInterneu = pdgfd_size / totsize,
    pdgfd_pthlh_pvalb_minus_size_to_allInterneu = pdgfd_pthlh_pvalb_minus_size / totsize,
    pdgfd_pthlh_pvalb_plus_size_to_allInterneu = pdgfd_pthlh_pvalb_plus_size / totsize,
    cck_vip_minus_size_to_allInterneu = cck_vip_minus_size / totsize,
    cck_vip_plus_size_to_allInterneu = cck_vip_plus_size / totsize
  )

df2$Species_Broad <- factor(df2$Species_Broad, levels = rev(c('Primate', 'Non_primate')))

## save primate vs non-primate and all species ratio calculations in separate sheets
write.csv(df2, file = "/home2/gkonop/workdir/03_INTEGRATE_ALL/Non_SPN/ANNOTATION/SuppTable3_interneuron_props2.csv", quote = F, row.names = F)
#############################################################
### generate the proportions across Species and Tissue and save in the same table
meta = second_non_spn@meta.data

# Calculate the size of each 'newannot' cluster for each combination of 'orig.ident', 'Tissue', and 'Species'
df3 <- meta %>%
  group_by(Tissue, Species) %>%
  summarize(
    totsize = n(), 
    chat_size = sum(newannot == "CHAT"),  
    sst_npy_size = sum(newannot == "SST_NPY"),  
    tac3_size = sum(newannot == "TAC3"),  
    th_size = sum(newannot == "TH"), 
    lmo3_bat_size = sum(newannot == "LMO3_BAT"),  
    foxp2_eya2_size = sum(newannot == "FOXP2_EYA2"),  
    foxp2_tshz2_bat_size = sum(newannot == "FOXP2_TSHZ2_BAT"), 
    pdgfd_size = sum(newannot == "PDGFD"), 
    pdgfd_pthlh_pvalb_minus_size = sum(newannot == "PDGFD_PTHLH_PVALB-"), 
    pdgfd_pthlh_pvalb_plus_size = sum(newannot == "PDGFD_PTHLH_PVALB+"), 
    cck_vip_minus_size = sum(newannot == "CCK_VIP-"), 
    cck_vip_plus_size = sum(newannot == "CCK_VIP+"), 
    .groups = 'drop'  # To drop the grouping after summarization
  )

# Compute ratios
df3 <- df3 %>%
  mutate(
    chat_size_to_allInterneu = chat_size / totsize,
    sst_npy_size_to_allInterneu = sst_npy_size / totsize,
    tac3_size_to_allInterneu = tac3_size / totsize,
    th_size_to_allInterneu = th_size / totsize,
    lmo3_bat_size_to_allInterneu = lmo3_bat_size / totsize,
    foxp2_eya2_size_to_allInterneu = foxp2_eya2_size / totsize,
    foxp2_tshz2_bat_size_to_allInterneu = foxp2_tshz2_bat_size / totsize,
    pdgfd_size_to_allInterneu = pdgfd_size / totsize,
    pdgfd_pthlh_pvalb_minus_size_to_allInterneu = pdgfd_pthlh_pvalb_minus_size / totsize,
    pdgfd_pthlh_pvalb_plus_size_to_allInterneu = pdgfd_pthlh_pvalb_plus_size / totsize,
    cck_vip_minus_size_to_allInterneu = cck_vip_minus_size / totsize,
    cck_vip_plus_size_to_allInterneu = cck_vip_plus_size / totsize
  )

df3$Species <- factor(df3$Species, levels = c('Bat', 'Ferret', 'Mouse', 'Marmoset', 'Macaque', 'Chimp', 'Human'))
## save 
write.csv(df3, file = "/home2/gkonop/workdir/03_INTEGRATE_ALL/Non_SPN/ANNOTATION/SuppTable3_interneuron_props3.csv", quote = F, row.names = F)
#########################################################################


# ratio for Caudate for all species and for mouse CaudoPutamen 
# subset the dataset
df_no_put <- df[df$Tissue != "Putamen",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_put %>%
  group_by(Species) %>%
  summarize(
    mean_chat_to_allInterneu = mean(chat_size_to_allInterneu, na.rm = TRUE),
    sd_chat_to_allInterneu = sd(chat_size_to_allInterneu, na.rm = TRUE),
    n_chat_to_allInterneu = sum(!is.na(chat_size_to_allInterneu)),
    sem_chat_to_allInterneu = ifelse(sd_chat_to_allInterneu == 0 | n_chat_to_allInterneu <= 1, 0, sd_chat_to_allInterneu / sqrt(n_chat_to_allInterneu)),
    
    mean_sst_npy_to_allInterneu = mean(sst_npy_size_to_allInterneu, na.rm = TRUE),
    sd_sst_npy_to_allInterneu = sd(sst_npy_size_to_allInterneu, na.rm = TRUE),
    n_sst_npy_to_allInterneu = sum(!is.na(sst_npy_size_to_allInterneu)),
    sem_sst_npy_to_allInterneu = ifelse(sd_sst_npy_to_allInterneu == 0 | n_sst_npy_to_allInterneu <= 1, 0, sd_sst_npy_to_allInterneu / sqrt(n_sst_npy_to_allInterneu)),
    
    mean_tac3_to_allInterneu = mean(tac3_size_to_allInterneu, na.rm = TRUE),
    sd_tac3_to_allInterneu = sd(tac3_size_to_allInterneu, na.rm = TRUE),
    n_tac3_to_allInterneu = sum(!is.na(tac3_size_to_allInterneu)),
    sem_tac3_to_allInterneu = ifelse(sd_tac3_to_allInterneu == 0 | n_tac3_to_allInterneu <= 1, 0, sd_tac3_to_allInterneu / sqrt(n_tac3_to_allInterneu)),
    
    mean_th_to_allInterneu = mean(th_size_to_allInterneu, na.rm = TRUE),
    sd_th_to_allInterneu = sd(th_size_to_allInterneu, na.rm = TRUE),
    n_th_to_allInterneu = sum(!is.na(th_size_to_allInterneu)),
    sem_th_to_allInterneu = ifelse(sd_th_to_allInterneu == 0 | n_th_to_allInterneu <= 1, 0, sd_th_to_allInterneu / sqrt(n_th_to_allInterneu)),
    
    mean_lmo3_bat_to_allInterneu = mean(lmo3_bat_size_to_allInterneu, na.rm = TRUE),
    sd_lmo3_bat_to_allInterneu = sd(lmo3_bat_size_to_allInterneu, na.rm = TRUE),
    n_lmo3_bat_to_allInterneu = sum(!is.na(lmo3_bat_size_to_allInterneu)),
    sem_lmo3_bat_to_allInterneu = ifelse(sd_lmo3_bat_to_allInterneu == 0 | n_lmo3_bat_to_allInterneu <= 1, 0, sd_lmo3_bat_to_allInterneu / sqrt(n_lmo3_bat_to_allInterneu)),
    
    mean_foxp2_eya2_to_allInterneu = mean(foxp2_eya2_size_to_allInterneu, na.rm = TRUE),
    sd_foxp2_eya2_to_allInterneu = sd(foxp2_eya2_size_to_allInterneu, na.rm = TRUE),
    n_foxp2_eya2_to_allInterneu = sum(!is.na(foxp2_eya2_size_to_allInterneu)),
    sem_foxp2_eya2_to_allInterneu = ifelse(sd_foxp2_eya2_to_allInterneu == 0 | n_foxp2_eya2_to_allInterneu <= 1, 0, sd_foxp2_eya2_to_allInterneu / sqrt(n_foxp2_eya2_to_allInterneu)),
    
    mean_foxp2_tshz2_bat_to_allInterneu = mean(foxp2_tshz2_bat_size_to_allInterneu, na.rm = TRUE),
    sd_foxp2_tshz2_bat_to_allInterneu = sd(foxp2_tshz2_bat_size_to_allInterneu, na.rm = TRUE),
    n_foxp2_tshz2_bat_to_allInterneu = sum(!is.na(foxp2_tshz2_bat_size_to_allInterneu)),
    sem_foxp2_tshz2_bat_to_allInterneu = ifelse(sd_foxp2_tshz2_bat_to_allInterneu == 0 | n_foxp2_tshz2_bat_to_allInterneu <= 1, 0, sd_foxp2_tshz2_bat_to_allInterneu / sqrt(n_foxp2_tshz2_bat_to_allInterneu)),
    
    mean_pdgfd_to_allInterneu = mean(pdgfd_size_to_allInterneu, na.rm = TRUE),
    sd_pdgfd_to_allInterneu = sd(pdgfd_size_to_allInterneu, na.rm = TRUE),
    n_pdgfd_to_allInterneu = sum(!is.na(pdgfd_size_to_allInterneu)),
    sem_pdgfd_to_allInterneu = ifelse(sd_pdgfd_to_allInterneu == 0 | n_pdgfd_to_allInterneu <= 1, 0, sd_pdgfd_to_allInterneu / sqrt(n_pdgfd_to_allInterneu)),
    
    mean_pdgfd_pthlh_pvalb_minus_to_allInterneu = mean(pdgfd_pthlh_pvalb_minus_size_to_allInterneu, na.rm = TRUE),
    sd_pdgfd_pthlh_pvalb_minus_to_allInterneu = sd(pdgfd_pthlh_pvalb_minus_size_to_allInterneu, na.rm = TRUE),
    n_pdgfd_pthlh_pvalb_minus_to_allInterneu = sum(!is.na(pdgfd_pthlh_pvalb_minus_size_to_allInterneu)),
    sem_pdgfd_pthlh_pvalb_minus_to_allInterneu = ifelse(sd_pdgfd_pthlh_pvalb_minus_to_allInterneu == 0 | n_pdgfd_pthlh_pvalb_minus_to_allInterneu <= 1, 0, sd_pdgfd_pthlh_pvalb_minus_to_allInterneu / sqrt(n_pdgfd_pthlh_pvalb_minus_to_allInterneu)),
    
    mean_pdgfd_pthlh_pvalb_plus_to_allInterneu = mean(pdgfd_pthlh_pvalb_plus_size_to_allInterneu, na.rm = TRUE),
    sd_pdgfd_pthlh_pvalb_plus_to_allInterneu = sd(pdgfd_pthlh_pvalb_plus_size_to_allInterneu, na.rm = TRUE),
    n_pdgfd_pthlh_pvalb_plus_to_allInterneu = sum(!is.na(pdgfd_pthlh_pvalb_plus_size_to_allInterneu)),
    sem_pdgfd_pthlh_pvalb_plus_to_allInterneu = ifelse(sd_pdgfd_pthlh_pvalb_plus_to_allInterneu == 0 | n_pdgfd_pthlh_pvalb_plus_to_allInterneu <= 1, 0, sd_pdgfd_pthlh_pvalb_plus_to_allInterneu / sqrt(n_pdgfd_pthlh_pvalb_plus_to_allInterneu)),
    
    mean_cck_vip_minus_to_allInterneu = mean(cck_vip_minus_size_to_allInterneu, na.rm = TRUE),
    sd_cck_vip_minus_to_allInterneu = sd(cck_vip_minus_size_to_allInterneu, na.rm = TRUE),
    n_cck_vip_minus_to_allInterneu = sum(!is.na(cck_vip_minus_size_to_allInterneu)),
    sem_cck_vip_minus_to_allInterneu = ifelse(sd_cck_vip_minus_to_allInterneu == 0 | n_cck_vip_minus_to_allInterneu <= 1, 0, sd(cck_vip_minus_size_to_allInterneu) / sqrt(n_cck_vip_minus_to_allInterneu)),
    
    mean_cck_vip_plus_to_allInterneu = mean(cck_vip_plus_size_to_allInterneu, na.rm = TRUE),
    sd_cck_vip_plus_to_allInterneu = sd(cck_vip_plus_size_to_allInterneu, na.rm = TRUE),
    n_cck_vip_plus_to_allInterneu = sum(!is.na(cck_vip_plus_size_to_allInterneu)),
    sem_cck_vip_plus_to_allInterneu = ifelse(sd(cck_vip_plus_size_to_allInterneu) == 0 | n_cck_vip_plus_to_allInterneu <= 1, 0, sd(cck_vip_plus_size_to_allInterneu) / sqrt(n_cck_vip_plus_to_allInterneu)),
    
    .groups = 'drop'
  )
### RUN source("~/SCRIPTS/03_INTEGRATION_ALL/interneu_prop_calc.R")

## Perform t-tests for each species combination
# get species
species_levels <- c("Bat", "Ferret", "Mouse", "Marmoset", "Macaque", "Chimp","Human")

# List of variables (statistics) to compare
stat_vars <- c("chat_size_to_allInterneu", 
               "sst_npy_size_to_allInterneu",
               "tac3_size_to_allInterneu",
               "th_size_to_allInterneu",
               "lmo3_bat_size_to_allInterneu",
               "foxp2_eya2_size_to_allInterneu", 
               "foxp2_tshz2_bat_size_to_allInterneu", 
               "pdgfd_size_to_allInterneu", 
	       "pdgfd_pthlh_pvalb_minus_size_to_allInterneu", 
	       "pdgfd_pthlh_pvalb_plus_size_to_allInterneu", 
	       "cck_vip_minus_size_to_allInterneu",
               "cck_vip_plus_size_to_allInterneu")
################################################################################
#### TROUBLESHOOTING
##t.test(species1_cells, species2_cells)and stat_compare_means(comparisons = comps, method = 't.test', paired = FALSE) are giving very different p vals
## Could it be equal/unequal var assumption?
## Could it be multi-test correction?
### First method
species_pairs <- combn(species_levels, 2, simplify = FALSE)  # This generates unique pairs

species_pair <- species_pairs[[6]]
species1 <- species_pair[1]
species2 <- species_pair[2]

stat <- "chat_size_to_allInterneu"

(cell_type <- toupper(sub("_size.*", "", stat)))  # This removes everything after "_size" and converts to uppercase

# Subset the cells by species and cell type for the current statistic
(species1_cells <- df_no_put[[stat]][df_no_put$Species == species1])
(species2_cells <- df_no_put[[stat]][df_no_put$Species == species2])

if (length(species1_cells) > 1 && length(species2_cells) > 1) {
      t_test <- t.test(species1_cells, species2_cells)
      
      # Store the results with the statistic, cell type, and species comparison
      comparison_name <- paste(stat, species1, species2, sep = "_vs_")
      t_test_results[[comparison_name]] <- list(
        Comparison = comparison_name,
        Cell_Type = cell_type,  # Add Cell_Type (extracted from statistic name)
        Species1 = species1,
        Species2 = species2,
        p_value = t_test$p.value
      )
    }

### Second method
# CHAT
meta <- second_non_spn[[]]
# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species) %>%
  summarize(
    totsize = n(),
    CHATsize = sum(newannot == 'CHAT'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    CHAT_to_allInterneu = CHATsize / totsize   
  )

df$Species <- factor(df$Species, levels = c('Bat', 'Ferret', 'Mouse', 'Marmoset', 'Macaque', 'Chimp', 'Human'))

# ratio for Caudate for all species and for mouse CaudoPutamen 
# subset the dataset
df_no_put <- df[df$Tissue != "Putamen",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_put %>%
  group_by(Species) %>%
  summarize(
    mean_CHAT_to_allInterneu = mean(chat_size_to_allInterneu, na.rm = TRUE),
    sd_CHAT_to_allInterneu = sd(chat_size_to_allInterneu, na.rm = TRUE),
    n = sum(!is.na(chat_size_to_allInterneu)),
    sem_CHAT_to_allInterneu = sd_CHAT_to_allInterneu / sqrt(n),
    .groups = 'drop'
  )
    
# Comparison with mouse (caudoputamen)
comps <- list(
  c('Human', 'Bat')
)

df_no_put <- as.data.frame(df_no_put)
# Subset the dataframe and create a new dataframe with the specified columns
b <- df_no_put[, c("orig.ident", "Tissue", "Species", "chat_size_to_allInterneu")]

(a <-
  ggplot(data = b, aes(x = Species, y = chat_size_to_allInterneu, color = Species)) +
  
  labs(
    x = NULL,
    y = paste('CHAT / Interneuron Ratio in Caudate'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black")+ 
  stat_compare_means(comparisons = comps, method = 't.test', paired = FALSE))
### t_test

#	Welch Two Sample t-test

#data:  species1_cells and species2_cells
#t = 9.9263, df = 7.2797, p-value = 1.726e-05
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
# 0.06991451 0.11319785
#sample estimates:
#mean of x mean of y 
#0.1252533 0.0336971 

### ggplot t-test p val = 1.7e-05
#Both p vals are equal. SUCCESS!
################################################################################
# Initialize an empty list to store t-test results
t_test_results <- list()

# Loop through each statistic and each combination of species
for (stat in stat_vars) {
  # Generate unique species pairs
  species_pairs <- combn(species_levels, 2, simplify = FALSE)  # This generates unique pairs
  
  for (species_pair in species_pairs) {
    species1 <- species_pair[1]
    species2 <- species_pair[2]

    # Extract the Cell_Type from the statistic name (before "_size" and in uppercase)
    cell_type <- toupper(sub("_size.*", "", stat))  # This removes everything after "_size" and converts to uppercase
    
    # Subset the cells by species and cell type for the current statistic
    species1_cells <- df_no_put[[stat]][df_no_put$Species == species1]
    species2_cells <- df_no_put[[stat]][df_no_put$Species == species2]
    
    # Perform t-test only if there is data for both species
    if (length(species1_cells) > 1 && length(species2_cells) > 1) {
      t_test <- t.test(species1_cells, species2_cells)
      
      # Store the results with the statistic, cell type, and species comparison
      comparison_name <- paste(stat, species1, species2, sep = "_vs_")
      t_test_results[[comparison_name]] <- list(
        Comparison = comparison_name,
        Cell_Type = cell_type,  # Add Cell_Type (extracted from statistic name)
        Species1 = species1,
        Species2 = species2,
        p_value = t_test$p.value
      )
    }
  }
}

# Convert the t_test_results list to a data frame
t_test_df <- do.call(rbind, lapply(t_test_results, data.frame, stringsAsFactors = FALSE))

# Add a new column 'is_signif' to indicate if p-value < 0.05
t_test_df$is_signif_pLessThan0_05 <- t_test_df$p_value < 0.05
t_test_df$is_signif_pLessThan0_01 <- t_test_df$p_value < 0.01

# Write the results to a CSV file
write.csv(t_test_df, "Caudate_t_test_results_with_celltype_comparison.csv", row.names = FALSE)

# Write the results to a CSV file
write.csv(df, "df.csv", row.names = FALSE)
###############################################################
##### PUTAMEN AND CAUDOPUTAMEN
# ratio for Caudate for all species and for mouse CaudoPutamen 
# subset the dataset
df_no_caud <- df[df$Tissue != "Caudate",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_caud %>%
  group_by(Species) %>%
  summarize(
    mean_chat_to_allInterneu = mean(chat_size_to_allInterneu, na.rm = TRUE),
    sd_chat_to_allInterneu = sd(chat_size_to_allInterneu, na.rm = TRUE),
    n_chat_to_allInterneu = sum(!is.na(chat_size_to_allInterneu)),
    sem_chat_to_allInterneu = ifelse(sd_chat_to_allInterneu == 0 | n_chat_to_allInterneu <= 1, 0, sd_chat_to_allInterneu / sqrt(n_chat_to_allInterneu)),
    
    mean_sst_npy_to_allInterneu = mean(sst_npy_size_to_allInterneu, na.rm = TRUE),
    sd_sst_npy_to_allInterneu = sd(sst_npy_size_to_allInterneu, na.rm = TRUE),
    n_sst_npy_to_allInterneu = sum(!is.na(sst_npy_size_to_allInterneu)),
    sem_sst_npy_to_allInterneu = ifelse(sd_sst_npy_to_allInterneu == 0 | n_sst_npy_to_allInterneu <= 1, 0, sd_sst_npy_to_allInterneu / sqrt(n_sst_npy_to_allInterneu)),
    
    mean_tac3_to_allInterneu = mean(tac3_size_to_allInterneu, na.rm = TRUE),
    sd_tac3_to_allInterneu = sd(tac3_size_to_allInterneu, na.rm = TRUE),
    n_tac3_to_allInterneu = sum(!is.na(tac3_size_to_allInterneu)),
    sem_tac3_to_allInterneu = ifelse(sd_tac3_to_allInterneu == 0 | n_tac3_to_allInterneu <= 1, 0, sd_tac3_to_allInterneu / sqrt(n_tac3_to_allInterneu)),
    
    mean_th_to_allInterneu = mean(th_size_to_allInterneu, na.rm = TRUE),
    sd_th_to_allInterneu = sd(th_size_to_allInterneu, na.rm = TRUE),
    n_th_to_allInterneu = sum(!is.na(th_size_to_allInterneu)),
    sem_th_to_allInterneu = ifelse(sd_th_to_allInterneu == 0 | n_th_to_allInterneu <= 1, 0, sd_th_to_allInterneu / sqrt(n_th_to_allInterneu)),
    
    mean_lmo3_bat_to_allInterneu = mean(lmo3_bat_size_to_allInterneu, na.rm = TRUE),
    sd_lmo3_bat_to_allInterneu = sd(lmo3_bat_size_to_allInterneu, na.rm = TRUE),
    n_lmo3_bat_to_allInterneu = sum(!is.na(lmo3_bat_size_to_allInterneu)),
    sem_lmo3_bat_to_allInterneu = ifelse(sd_lmo3_bat_to_allInterneu == 0 | n_lmo3_bat_to_allInterneu <= 1, 0, sd_lmo3_bat_to_allInterneu / sqrt(n_lmo3_bat_to_allInterneu)),
    
    mean_foxp2_eya2_to_allInterneu = mean(foxp2_eya2_size_to_allInterneu, na.rm = TRUE),
    sd_foxp2_eya2_to_allInterneu = sd(foxp2_eya2_size_to_allInterneu, na.rm = TRUE),
    n_foxp2_eya2_to_allInterneu = sum(!is.na(foxp2_eya2_size_to_allInterneu)),
    sem_foxp2_eya2_to_allInterneu = ifelse(sd_foxp2_eya2_to_allInterneu == 0 | n_foxp2_eya2_to_allInterneu <= 1, 0, sd_foxp2_eya2_to_allInterneu / sqrt(n_foxp2_eya2_to_allInterneu)),
    
    mean_foxp2_tshz2_bat_to_allInterneu = mean(foxp2_tshz2_bat_size_to_allInterneu, na.rm = TRUE),
    sd_foxp2_tshz2_bat_to_allInterneu = sd(foxp2_tshz2_bat_size_to_allInterneu, na.rm = TRUE),
    n_foxp2_tshz2_bat_to_allInterneu = sum(!is.na(foxp2_tshz2_bat_size_to_allInterneu)),
    sem_foxp2_tshz2_bat_to_allInterneu = ifelse(sd_foxp2_tshz2_bat_to_allInterneu == 0 | n_foxp2_tshz2_bat_to_allInterneu <= 1, 0, sd_foxp2_tshz2_bat_to_allInterneu / sqrt(n_foxp2_tshz2_bat_to_allInterneu)),
    
    mean_pdgfd_to_allInterneu = mean(pdgfd_size_to_allInterneu, na.rm = TRUE),
    sd_pdgfd_to_allInterneu = sd(pdgfd_size_to_allInterneu, na.rm = TRUE),
    n_pdgfd_to_allInterneu = sum(!is.na(pdgfd_size_to_allInterneu)),
    sem_pdgfd_to_allInterneu = ifelse(sd_pdgfd_to_allInterneu == 0 | n_pdgfd_to_allInterneu <= 1, 0, sd_pdgfd_to_allInterneu / sqrt(n_pdgfd_to_allInterneu)),
    
    mean_pdgfd_pthlh_pvalb_minus_to_allInterneu = mean(pdgfd_pthlh_pvalb_minus_size_to_allInterneu, na.rm = TRUE),
    sd_pdgfd_pthlh_pvalb_minus_to_allInterneu = sd(pdgfd_pthlh_pvalb_minus_size_to_allInterneu, na.rm = TRUE),
    n_pdgfd_pthlh_pvalb_minus_to_allInterneu = sum(!is.na(pdgfd_pthlh_pvalb_minus_size_to_allInterneu)),
    sem_pdgfd_pthlh_pvalb_minus_to_allInterneu = ifelse(sd_pdgfd_pthlh_pvalb_minus_to_allInterneu == 0 | n_pdgfd_pthlh_pvalb_minus_to_allInterneu <= 1, 0, sd_pdgfd_pthlh_pvalb_minus_to_allInterneu / sqrt(n_pdgfd_pthlh_pvalb_minus_to_allInterneu)),
    
    mean_pdgfd_pthlh_pvalb_plus_to_allInterneu = mean(pdgfd_pthlh_pvalb_plus_size_to_allInterneu, na.rm = TRUE),
    sd_pdgfd_pthlh_pvalb_plus_to_allInterneu = sd(pdgfd_pthlh_pvalb_plus_size_to_allInterneu, na.rm = TRUE),
    n_pdgfd_pthlh_pvalb_plus_to_allInterneu = sum(!is.na(pdgfd_pthlh_pvalb_plus_size_to_allInterneu)),
    sem_pdgfd_pthlh_pvalb_plus_to_allInterneu = ifelse(sd_pdgfd_pthlh_pvalb_plus_to_allInterneu == 0 | n_pdgfd_pthlh_pvalb_plus_to_allInterneu <= 1, 0, sd_pdgfd_pthlh_pvalb_plus_to_allInterneu / sqrt(n_pdgfd_pthlh_pvalb_plus_to_allInterneu)),
    
    mean_cck_vip_minus_to_allInterneu = mean(cck_vip_minus_size_to_allInterneu, na.rm = TRUE),
    sd_cck_vip_minus_to_allInterneu = sd(cck_vip_minus_size_to_allInterneu, na.rm = TRUE),
    n_cck_vip_minus_to_allInterneu = sum(!is.na(cck_vip_minus_size_to_allInterneu)),
    sem_cck_vip_minus_to_allInterneu = ifelse(sd_cck_vip_minus_to_allInterneu == 0 | n_cck_vip_minus_to_allInterneu <= 1, 0, sd(cck_vip_minus_size_to_allInterneu) / sqrt(n_cck_vip_minus_to_allInterneu)),
    
    mean_cck_vip_plus_to_allInterneu = mean(cck_vip_plus_size_to_allInterneu, na.rm = TRUE),
    sd_cck_vip_plus_to_allInterneu = sd(cck_vip_plus_size_to_allInterneu, na.rm = TRUE),
    n_cck_vip_plus_to_allInterneu = sum(!is.na(cck_vip_plus_size_to_allInterneu)),
    sem_cck_vip_plus_to_allInterneu = ifelse(sd(cck_vip_plus_size_to_allInterneu) == 0 | n_cck_vip_plus_to_allInterneu <= 1, 0, sd(cck_vip_plus_size_to_allInterneu) / sqrt(n_cck_vip_plus_to_allInterneu)), .groups = 'drop')
### RUN source("~/SCRIPTS/03_INTEGRATION_ALL/interneu_prop_calcPutamen.R")

## Perform t-tests for each species combination
# get species
species_levels <- c("Bat", "Ferret", "Mouse", "Marmoset", "Macaque", "Chimp","Human")

# List of variables (statistics) to compare
stat_vars <- c("chat_size_to_allInterneu", 
               "sst_npy_size_to_allInterneu",
               "tac3_size_to_allInterneu",
               "th_size_to_allInterneu",
               "lmo3_bat_size_to_allInterneu",
               "foxp2_eya2_size_to_allInterneu", 
               "foxp2_tshz2_bat_size_to_allInterneu", 
               "pdgfd_size_to_allInterneu", 
	       "pdgfd_pthlh_pvalb_minus_size_to_allInterneu", 
	       "pdgfd_pthlh_pvalb_plus_size_to_allInterneu", 
	       "cck_vip_minus_size_to_allInterneu",
               "cck_vip_plus_size_to_allInterneu")


# Initialize an empty list to store t-test results
t_test_results <- list()

# Loop through each statistic and each combination of species
for (stat in stat_vars) {
  # Generate unique species pairs
  species_pairs <- combn(species_levels, 2, simplify = FALSE)  # This generates unique pairs
  
  for (species_pair in species_pairs) {
    species1 <- species_pair[1]
    species2 <- species_pair[2]
    
    # Extract the Cell_Type from the statistic name (before "_size" and in uppercase)
    cell_type <- toupper(sub("_size.*", "", stat))  # This removes everything after "_size" and converts to uppercase
    
    # Subset the cells by species and cell type for the current statistic
    species1_cells <- df_no_caud[[stat]][df_no_caud$Species == species1]
    species2_cells <- df_no_caud[[stat]][df_no_caud$Species == species2]
    
    # Perform t-test only if there is data for both species
    if (length(species1_cells) > 1 && length(species2_cells) > 1) {
      t_test <- t.test(species1_cells, species2_cells)
      
      # Store the results with the statistic, cell type, and species comparison
      comparison_name <- paste(stat, species1, species2, sep = "_vs_")
      t_test_results[[comparison_name]] <- list(
        Comparison = comparison_name,
        Cell_Type = cell_type,  # Add Cell_Type (extracted from statistic name)
        Species1 = species1,
        Species2 = species2,
        p_value = t_test$p.value
      )
    }
  }
}

# Convert the t_test_results list to a data frame
t_test_df <- do.call(rbind, lapply(t_test_results, data.frame, stringsAsFactors = FALSE))

# Add a new column 'is_signif' to indicate if p-value < 0.05
t_test_df$is_signif_pLessThan0_05 <- t_test_df$p_value < 0.05
t_test_df$is_signif_pLessThan0_01 <- t_test_df$p_value < 0.01

# Write the results to a CSV file
write.csv(t_test_df, "Putamen_t_test_results_with_celltype_comparison.csv", row.names = FALSE)
################################################################
#### COMPARISONS
# change directory
setwd("/endosome/work/Neuroinformatics_Core/gkonop/03_INTEGRATE_ALL/Non_SPN/ANNOTATION")

### Human specific comparisons
## Caudate
# CHAT
meta <- second_non_spn[[]]
# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species) %>%
  summarize(
    totsize = n(),
    CHATsize = sum(newannot == 'CHAT'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    CHAT_to_allInterneu = CHATsize / totsize   
  )

df$Species <- factor(df$Species, levels = c('Bat', 'Ferret', 'Mouse', 'Marmoset', 'Macaque', 'Chimp', 'Human'))

# ratio for Caudate for all species and for mouse CaudoPutamen 
# subset the dataset
df_no_put <- df[df$Tissue != "Putamen",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_put %>%
  group_by(Species) %>%
  summarize(
    mean_CHAT_to_allInterneu = mean(CHAT_to_allInterneu, na.rm = TRUE),
    sd_CHAT_to_allInterneu = sd(CHAT_to_allInterneu, na.rm = TRUE),
    n = sum(!is.na(CHAT_to_allInterneu)),
    sem_CHAT_to_allInterneu = sd_CHAT_to_allInterneu / sqrt(n),
    .groups = 'drop'
  )
    
# Comparison with mouse (caudoputamen)
comps <- list(
  c('Ferret', 'Human'),
  c('Mouse', 'Human'), 
  c('Marmoset', 'Human'), 
  c('Macaque', 'Human'), 
  c('Chimp', 'Human'), 
  c('Bat', 'Human')
)


# Plot with p-values
pdf(file = "GB_barplot_CHAT_to_all_interneurons_Ratio_Caudate_with_mouse_Caudoputamen_t-test_pvals_Across_Species_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put, aes(x = Species, y = CHAT_to_allInterneu, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_CHAT_to_allInterneu, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_CHAT_to_allInterneu - sem_CHAT_to_allInterneu, ymax = mean_CHAT_to_allInterneu + sem_CHAT_to_allInterneu),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('CHAT / Interneuron Ratio in Caudate'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black")+ 
  stat_compare_means(comparisons = comps, method = 't.test', paired = FALSE) + coord_flip()
)
dev.off()
  
# Plot without p-values
pdf(file = "GB_barplot_CHAT_to_all_interneurons_Ratio_Caudate_with_mouse_Caudoputamen_Across_Species_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put, aes(x = Species, y = CHAT_to_allInterneu, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_CHAT_to_allInterneu, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_CHAT_to_allInterneu - sem_CHAT_to_allInterneu, ymax = mean_CHAT_to_allInterneu + sem_CHAT_to_allInterneu),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('CHAT / Interneuron Ratio in Caudate'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black")+ coord_flip()
)
dev.off()
############################################################
### Human specific comparisons
## Caudate
# CCK_VIP+
meta <- second_non_spn[[]]
# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species) %>%
  summarize(
    totsize = n(),
    CCK_VIP_pos_size = sum(newannot == 'CCK_VIP+'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    CCK_VIP_pos_to_allInterneu = CCK_VIP_pos_size / totsize   
  )

df$Species <- factor(df$Species, levels = c('Bat', 'Ferret', 'Mouse', 'Marmoset', 'Macaque', 'Chimp', 'Human'))

# ratio for Caudate for all species and for mouse CaudoPutamen 
# subset the dataset
df_no_put <- df[df$Tissue != "Putamen",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_put %>%
  group_by(Species) %>%
  summarize(
    mean_CCK_VIP_pos_to_allInterneu = mean(CCK_VIP_pos_to_allInterneu, na.rm = TRUE),
    sd_CCK_VIP_pos_to_allInterneu = sd(CCK_VIP_pos_to_allInterneu, na.rm = TRUE),
    n = sum(!is.na(CCK_VIP_pos_to_allInterneu)),
    sem_CCK_VIP_pos_to_allInterneu = sd_CCK_VIP_pos_to_allInterneu / sqrt(n),
    .groups = 'drop'
  )
    
# Comparison with mouse (caudoputamen)
comps <- list(
  c('Ferret', 'Human'),
  c('Mouse', 'Human'), 
  c('Marmoset', 'Human'), 
  c('Macaque', 'Human'), 
  c('Chimp', 'Human'), 
  c('Bat', 'Human')
)

# Plot with p-values
pdf(file = "GB_barplot_CCK_VIP_pos_to_all_interneurons_Ratio_Caudate_with_mouse_Caudoputamen_t-test_pvals_Across_Species_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put, aes(x = Species, y = CCK_VIP_pos_to_allInterneu, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_CCK_VIP_pos_to_allInterneu, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_CCK_VIP_pos_to_allInterneu - sem_CCK_VIP_pos_to_allInterneu, ymax = mean_CCK_VIP_pos_to_allInterneu + sem_CCK_VIP_pos_to_allInterneu),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('CCK_VIP_pos / Interneuron Ratio in Caudate'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black")+ 
  stat_compare_means(comparisons = comps, method = 't.test', paired = FALSE) + coord_flip()
)
dev.off()
  
# Plot without p-values
pdf(file = "GB_barplot_CCK_VIP_pos_to_all_interneurons_Ratio_Caudate_with_mouse_Caudoputamen_Across_Species_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put, aes(x = Species, y = CCK_VIP_pos_to_allInterneu, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_CCK_VIP_pos_to_allInterneu, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_CCK_VIP_pos_to_allInterneu - sem_CCK_VIP_pos_to_allInterneu, ymax = mean_CCK_VIP_pos_to_allInterneu + sem_CCK_VIP_pos_to_allInterneu),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('CCK_VIP_pos / Interneuron Ratio in Caudate'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black")+ coord_flip()
)
dev.off()
############################################################
### Human specific comparisons
## Putamen
# CCK_VIP+
meta <- second_non_spn[[]]
# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species) %>%
  summarize(
    totsize = n(),
    CCK_VIP_pos_size = sum(newannot == 'CCK_VIP+'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    CCK_VIP_pos_to_allInterneu = CCK_VIP_pos_size / totsize   
  )

df$Species <- factor(df$Species, levels = c('Bat', 'Mouse', 'Marmoset', 'Macaque', 'Chimp', 'Human'))



# ratio for Caudate for all species and for mouse CaudoPutamen 
# subset the dataset
df_no_caud <- df[df$Tissue != "Caudate",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_caud %>%
  group_by(Species) %>%
  summarize(
    mean_CCK_VIP_pos_to_allInterneu = mean(CCK_VIP_pos_to_allInterneu, na.rm = TRUE),
    sd_CCK_VIP_pos_to_allInterneu = sd(CCK_VIP_pos_to_allInterneu, na.rm = TRUE),
    n = sum(!is.na(CCK_VIP_pos_to_allInterneu)),
    sem_CCK_VIP_pos_to_allInterneu = sd_CCK_VIP_pos_to_allInterneu / sqrt(n),
    .groups = 'drop'
  )
    
# Comparison with mouse (caudoputamen)
comps <- list(
  c('Mouse', 'Human'), 
  c('Marmoset', 'Human'), 
  c('Macaque', 'Human'), 
  c('Chimp', 'Human'), 
  c('Bat', 'Human')
)


# Plot with p-values
pdf(file = "GB_barplot_CCK_VIP_pos_to_all_interneurons_Ratio_Putamen_with_mouse_Caudoputamen_t-test_pvals_Across_Species_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_caud, aes(x = Species, y = CCK_VIP_pos_to_allInterneu, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_CCK_VIP_pos_to_allInterneu, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_CCK_VIP_pos_to_allInterneu - sem_CCK_VIP_pos_to_allInterneu, ymax = mean_CCK_VIP_pos_to_allInterneu + sem_CCK_VIP_pos_to_allInterneu),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('CCK_VIP_pos / Interneuron Ratio in Putamen'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") + 
  stat_compare_means(comparisons = comps, method = 't.test', paired = FALSE) + coord_flip()
)
dev.off()
  
# Plot without p-values
pdf(file = "GB_barplot_CCK_VIP_pos_to_all_interneurons_Ratio_Putamen_with_mouse_Caudoputamen_Across_Species_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_caud, aes(x = Species, y = CCK_VIP_pos_to_allInterneu, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_CCK_VIP_pos_to_allInterneu, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_CCK_VIP_pos_to_allInterneu - sem_CCK_VIP_pos_to_allInterneu, ymax = mean_CCK_VIP_pos_to_allInterneu + sem_CCK_VIP_pos_to_allInterneu),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('CCK_VIP_pos / Interneuron Ratio in Putamen'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black")+ coord_flip()
)
dev.off()
############################################################
### Mouse specific comparisons
## AllTissue
# TH
meta <- second_non_spn[[]]
# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species) %>%
  summarize(
    totsize = n(),
    TH_size = sum(newannot == 'TH'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    TH_to_allInterneu = TH_size / totsize   
  )

df$Species <- factor(df$Species, levels = c('Bat', 'Ferret', 'Mouse', 'Marmoset', 'Macaque', 'Chimp', 'Human'))


# Calculate summary statistics for the bar plots
df_summary <- df %>%
  group_by(Species) %>%
  summarize(
    mean_TH_to_allInterneu = mean(TH_to_allInterneu, na.rm = TRUE),
    sd_TH_to_allInterneu = sd(TH_to_allInterneu, na.rm = TRUE),
    n = sum(!is.na(TH_to_allInterneu)),
    sem_TH_to_allInterneu = sd_TH_to_allInterneu / sqrt(n),
    .groups = 'drop'
  )
    
# Comparison with mouse (caudoputamen)
comps <- list(
  c('Ferret', 'Mouse'),
  c('Mouse', 'Human'), 
  c('Marmoset', 'Mouse'), 
  c('Macaque', 'Mouse'), 
  c('Chimp', 'Mouse'), 
  c('Bat', 'Mouse')
)

# Plot with p-values
pdf(file = "GB_barplot_TH_to_all_interneurons_Ratio_AllTissue_t-test_pvals_Across_Species_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df, aes(x = Species, y = TH_to_allInterneu, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_TH_to_allInterneu, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_TH_to_allInterneu - sem_TH_to_allInterneu, ymax = mean_TH_to_allInterneu + sem_TH_to_allInterneu),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('TH / Interneuron Ratio in All Tissue'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black")+ 
  stat_compare_means(comparisons = comps, method = 't.test', paired = FALSE) + coord_flip()
)
dev.off()
  
# Plot with p-values
pdf(file = "GB_barplot_TH_to_all_interneurons_Ratio_AllTissue_Across_Species_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df, aes(x = Species, y = TH_to_allInterneu, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_TH_to_allInterneu, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_TH_to_allInterneu - sem_TH_to_allInterneu, ymax = mean_TH_to_allInterneu + sem_TH_to_allInterneu),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('TH / Interneuron Ratio in All Tissue'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black")+ coord_flip()
)
dev.off()
############################################################
### Mouse specific comparisons
## AllTissue
# PDGFD_PTHLH_PVALB+
meta <- second_non_spn[[]] 
# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species) %>%
  summarize(
    totsize = n(),
    PDGFD_PTHLH_PVALB_plus_size = sum(newannot == 'PDGFD_PTHLH_PVALB+'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    PDGFD_PTHLH_PVALB_plus_to_allInterneu = PDGFD_PTHLH_PVALB_plus_size / totsize   
  )

df$Species <- factor(df$Species, levels = c('Bat', 'Ferret', 'Mouse', 'Marmoset', 'Macaque', 'Chimp', 'Human'))

# Calculate summary statistics for the bar plots
df_summary <- df %>%
  group_by(Species) %>%
  summarize(
    mean_PDGFD_PTHLH_PVALB_plus_to_allInterneu = mean(PDGFD_PTHLH_PVALB_plus_to_allInterneu, na.rm = TRUE),
    sd_PDGFD_PTHLH_PVALB_plus_to_allInterneu = sd(PDGFD_PTHLH_PVALB_plus_to_allInterneu, na.rm = TRUE),
    n = sum(!is.na(PDGFD_PTHLH_PVALB_plus_to_allInterneu)),
    sem_PDGFD_PTHLH_PVALB_plus_to_allInterneu = sd_PDGFD_PTHLH_PVALB_plus_to_allInterneu / sqrt(n),
    .groups = 'drop'
  )
    
# Comparison with mouse (caudoputamen)
comps <- list(
  c('Ferret', 'Mouse'),
  c('Mouse', 'Human'), 
  c('Marmoset', 'Mouse'), 
  c('Macaque', 'Mouse'), 
  c('Chimp', 'Mouse'), 
  c('Bat', 'Mouse')
)

# Plot with p-values
pdf(file = "GB_barplot_PDGFD_PTHLH_PVALB_plus_to_all_interneurons_Ratio_AllTissue_t-test_pvals_Across_Species_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df, aes(x = Species, y = PDGFD_PTHLH_PVALB_plus_to_allInterneu, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_PDGFD_PTHLH_PVALB_plus_to_allInterneu, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_PDGFD_PTHLH_PVALB_plus_to_allInterneu - sem_PDGFD_PTHLH_PVALB_plus_to_allInterneu, ymax = mean_PDGFD_PTHLH_PVALB_plus_to_allInterneu + sem_PDGFD_PTHLH_PVALB_plus_to_allInterneu),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('PDGFD_PTHLH_PVALB+ / Interneuron Ratio in All Tissue'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black")+ 
  stat_compare_means(comparisons = comps, method = 't.test', paired = FALSE) + coord_flip()
)
dev.off()
  
# Plot without p-values
pdf(file = "GB_barplot_PDGFD_PTHLH_PVALB_plus_to_all_interneurons_Ratio_AllTissue_Across_Species_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df, aes(x = Species, y = PDGFD_PTHLH_PVALB_plus_to_allInterneu, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_PDGFD_PTHLH_PVALB_plus_to_allInterneu, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_PDGFD_PTHLH_PVALB_plus_to_allInterneu - sem_PDGFD_PTHLH_PVALB_plus_to_allInterneu, ymax = mean_PDGFD_PTHLH_PVALB_plus_to_allInterneu + sem_PDGFD_PTHLH_PVALB_plus_to_allInterneu),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('PDGFD_PTHLH_PVALB+ / Interneuron Ratio in All Tissue'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black")+ coord_flip()
)
dev.off()
############################################################
### Bat specific comparisons
## Putamen
# LMO3_BAT
meta <- second_non_spn[[]]

# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species) %>%
  summarize(
    totsize = n(),
    LMO3_BAT_size = sum(newannot == 'LMO3_BAT'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    LMO3_BAT_to_allInterneu = LMO3_BAT_size / totsize   
  )

df$Species <- factor(df$Species, levels = c('Bat', 'Mouse', 'Marmoset', 'Macaque', 'Chimp', 'Human'))

# ratio for Caudate for all species and for mouse CaudoPutamen 
# subset the dataset
df_no_caud <- df[df$Tissue != "Caudate",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_caud %>%
  group_by(Species) %>%
  summarize(
    mean_LMO3_BAT_to_allInterneu = mean(LMO3_BAT_to_allInterneu, na.rm = TRUE),
    sd_LMO3_BAT_to_allInterneu = sd(LMO3_BAT_to_allInterneu, na.rm = TRUE),
    n = sum(!is.na(LMO3_BAT_to_allInterneu)),
    sem_LMO3_BAT_to_allInterneu = sd_LMO3_BAT_to_allInterneu / sqrt(n),
    .groups = 'drop'
  )
    
# Comparison with mouse (caudoputamen)
comps <- list(
  c('Mouse', 'Bat'), 
  c('Marmoset', 'Bat'), 
  c('Macaque', 'Bat'), 
  c('Chimp', 'Bat'), 
  c('Bat', 'Human')
)


# Plot with p-values
pdf(file = "GB_barplot_LMO3_BAT_to_all_interneurons_Ratio_Putamen_with_mouse_Caudoputamen_t-test_pvals_Across_Species_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_caud, aes(x = Species, y = LMO3_BAT_to_allInterneu, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_LMO3_BAT_to_allInterneu, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_LMO3_BAT_to_allInterneu - sem_LMO3_BAT_to_allInterneu, ymax = mean_LMO3_BAT_to_allInterneu + sem_LMO3_BAT_to_allInterneu),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('LMO3_BAT / Interneuron Ratio in Putamen'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") + 
  stat_compare_means(comparisons = comps, method = 't.test', paired = FALSE) + coord_flip()
)
dev.off()
  
# Plot without p-values
pdf(file = "GB_barplot_LMO3_BAT_to_all_interneurons_Ratio_Putamen_with_mouse_Caudoputamen_Across_Species_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_caud, aes(x = Species, y = LMO3_BAT_to_allInterneu, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_LMO3_BAT_to_allInterneu, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_LMO3_BAT_to_allInterneu - sem_LMO3_BAT_to_allInterneu, ymax = mean_LMO3_BAT_to_allInterneu + sem_LMO3_BAT_to_allInterneu),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('LMO3_BAT / Interneuron Ratio in Putamen'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black")  + coord_flip()
)
dev.off()
############################################################
### Bat specific comparisons
## Putamen
# FOXP2_TSHZ2_BAT
meta <- second_non_spn[[]]

# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species) %>%
  summarize(
    totsize = n(),
    FOXP2_TSHZ2_BAT_size = sum(newannot == 'FOXP2_TSHZ2_BAT'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    FOXP2_TSHZ2_BAT_to_allInterneu = FOXP2_TSHZ2_BAT_size / totsize   
  )

df$Species <- factor(df$Species, levels = c('Bat', 'Mouse', 'Marmoset', 'Macaque', 'Chimp', 'Human'))

# Ratio for Caudate for all species and for mouse CaudoPutamen 
# Subset the dataset
df_no_caud <- df[df$Tissue != "Caudate",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_caud %>%
  group_by(Species) %>%
  summarize(
    mean_FOXP2_TSHZ2_BAT_to_allInterneu = mean(FOXP2_TSHZ2_BAT_to_allInterneu, na.rm = TRUE),
    sd_FOXP2_TSHZ2_BAT_to_allInterneu = sd(FOXP2_TSHZ2_BAT_to_allInterneu, na.rm = TRUE),
    n = sum(!is.na(FOXP2_TSHZ2_BAT_to_allInterneu)),
    sem_FOXP2_TSHZ2_BAT_to_allInterneu = sd_FOXP2_TSHZ2_BAT_to_allInterneu / sqrt(n),
    .groups = 'drop'
  )

# Comparison with mouse (caudoputamen)
comps <- list(
  c('Mouse', 'Bat'), 
  c('Marmoset', 'Bat'), 
  c('Macaque', 'Bat'), 
  c('Chimp', 'Bat'), 
  c('Bat', 'Human')
)


# Plot with p-values
pdf(file = "GB_barplot_FOXP2_TSHZ2_BAT_to_all_interneurons_Ratio_Putamen_with_mouse_Caudoputamen_t-test_pvals_Across_Species_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_caud, aes(x = Species, y = FOXP2_TSHZ2_BAT_to_allInterneu, color = Species)) +
    geom_col(data = df_summary, aes(x = Species, y = mean_FOXP2_TSHZ2_BAT_to_allInterneu, fill = Species), position = position_dodge()) +
    geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_FOXP2_TSHZ2_BAT_to_allInterneu - sem_FOXP2_TSHZ2_BAT_to_allInterneu, ymax = mean_FOXP2_TSHZ2_BAT_to_allInterneu + sem_FOXP2_TSHZ2_BAT_to_allInterneu),
                  position = position_dodge(width = 0.9), width = 0.25) +
    theme_minimal() +
    labs(
      x = NULL,
      y = paste('FOXP2_TSHZ2_BAT / Interneuron Ratio in Putamen'),
      fill = 'Species',
      color = 'Species'
    ) +
    theme(
      legend.position = "none",
      text = element_text(size = 20, face = 'bold'),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
    ) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") + 
    stat_compare_means(comparisons = comps, method = 't.test', paired = FALSE) + coord_flip()
)
dev.off()

# Plot without p-values
pdf(file = "GB_barplot_FOXP2_TSHZ2_BAT_to_all_interneurons_Ratio_Putamen_with_mouse_Caudoputamen_Across_Species_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_caud, aes(x = Species, y = FOXP2_TSHZ2_BAT_to_allInterneu, color = Species)) +
    geom_col(data = df_summary, aes(x = Species, y = mean_FOXP2_TSHZ2_BAT_to_allInterneu, fill = Species), position = position_dodge()) +
    geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_FOXP2_TSHZ2_BAT_to_allInterneu - sem_FOXP2_TSHZ2_BAT_to_allInterneu, ymax = mean_FOXP2_TSHZ2_BAT_to_allInterneu + sem_FOXP2_TSHZ2_BAT_to_allInterneu),
                  position = position_dodge(width = 0.9), width = 0.25) +
    theme_minimal() +
    labs(
      x = NULL,
      y = paste('FOXP2_TSHZ2_BAT / Interneuron Ratio in Putamen'),
      fill = 'Species',
      color = 'Species'
    ) +
    theme(
      legend.position = "none",
      text = element_text(size = 20, face = 'bold'),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
    ) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black")  + coord_flip()
)
dev.off()
############################################################
### Primate specific
## Caudate
# PDGFD_PTHLH_PVALB+
meta = second_non_spn@meta.data
meta$Species_Broad = 'Primate'
meta[meta$Species %in% c('Bat','Ferret','Mouse'), 'Species_Broad'] = 'Non_primate'

# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species_Broad) %>%
  summarize(
    totsize = n(),
    PDGFD_PTHLH_PVALB_possize = sum(newannot == 'PDGFD_PTHLH_PVALB+'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    PDGFD_PTHLH_PVALB_pos_to_allInterneu = PDGFD_PTHLH_PVALB_possize / totsize
  )

df$Species_Broad <- factor(df$Species_Broad, levels = rev(c('Primate', 'Non_primate')))

# Ratio for Caudate for all species and for mouse CaudoPutamen 
# Subset the dataset
df_no_put <- df[df$Tissue != "Putamen",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_put %>%
  group_by(Species_Broad) %>%
  summarize(
    mean_PDGFD_PTHLH_PVALB_pos_to_allInterneu = mean(PDGFD_PTHLH_PVALB_pos_to_allInterneu, na.rm = TRUE),
    sd_PDGFD_PTHLH_PVALB_pos_to_allInterneu = sd(PDGFD_PTHLH_PVALB_pos_to_allInterneu, na.rm = TRUE),
    n = sum(!is.na(PDGFD_PTHLH_PVALB_pos_to_allInterneu)),
    sem_PDGFD_PTHLH_PVALB_pos_to_allInterneu = sd_PDGFD_PTHLH_PVALB_pos_to_allInterneu / sqrt(n),
    .groups = 'drop'
  )
    
# Comparison with mouse (caudoputamen)
comps <- list(
  c('Non_primate', 'Primate')
)

# Plot with p-values
pdf(file = "GB_barplot_PDGFD_PTHLH_PVALB_pos_AllInterneuron_Ratio_Caudate_with_mouse_Caudoputamen_Across_BroadSpecies_with_t-testpvals_withFerret_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put, aes(x = Species_Broad, y = PDGFD_PTHLH_PVALB_pos_to_allInterneu, color = Species_Broad)) +
  geom_col(data = df_summary, aes(x = Species_Broad, y = mean_PDGFD_PTHLH_PVALB_pos_to_allInterneu, fill = Species_Broad), position = position_dodge(), width = 0.4) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species_Broad, ymin = mean_PDGFD_PTHLH_PVALB_pos_to_allInterneu - sem_PDGFD_PTHLH_PVALB_pos_to_allInterneu, ymax = mean_PDGFD_PTHLH_PVALB_pos_to_allInterneu + sem_PDGFD_PTHLH_PVALB_pos_to_allInterneu),
                position = position_dodge(width = 0.9), width = 0.1) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('PDGFD_PTHLH_PVALB+ / Interneuron Ratio in Caudate'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(position = fixed_jitter, size = 2, alpha = 0.7, color = "black") + 
  stat_compare_means(comparisons = comps, method = 't.test', paired = FALSE)
)
dev.off()

# Plot without p-values
pdf(file = "GB_barplot_PDGFD_PTHLH_PVALB_pos_AllInterneuron_Ratio_Caudate_with_mouse_Caudoputamen_Across_BroadSpecies_withFerret_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put, aes(x = Species_Broad, y = PDGFD_PTHLH_PVALB_pos_to_allInterneu, color = Species_Broad)) +
  geom_col(data = df_summary, aes(x = Species_Broad, y = mean_PDGFD_PTHLH_PVALB_pos_to_allInterneu, fill = Species_Broad), position = position_dodge(0), width = 0.4) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species_Broad, ymin = mean_PDGFD_PTHLH_PVALB_pos_to_allInterneu - sem_PDGFD_PTHLH_PVALB_pos_to_allInterneu, ymax = mean_PDGFD_PTHLH_PVALB_pos_to_allInterneu + sem_PDGFD_PTHLH_PVALB_pos_to_allInterneu),
                position = position_dodge(width = 0.9), width = 0.1) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('PDGFD_PTHLH_PVALB+ / Interneuron Ratio in Caudate'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(position = fixed_jitter, size = 2, alpha = 0.7, color = "black")
)
dev.off()
##########################################################
## Putamen
# PDGFD_PTHLH_PVALB+
meta = second_non_spn@meta.data
meta$Species_Broad = 'Primate'
meta[meta$Species %in% c('Bat','Ferret','Mouse'), 'Species_Broad'] = 'Non_primate'

# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species_Broad) %>%
  summarize(
    totsize = n(),
    PDGFD_PTHLH_PVALB_pos_size = sum(newannot == 'PDGFD_PTHLH_PVALB+'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    PDGFD_PTHLH_PVALB_pos_to_allInterneu = PDGFD_PTHLH_PVALB_pos_size / totsize
  )

df$Species_Broad <- factor(df$Species_Broad, levels = rev(c('Primate', 'Non_primate')))

# Ratio for Putamen for all species and for mouse CaudoPutamen 
# Subset the dataset
df_no_caud <- df[df$Tissue != "Caudate",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_caud %>%
  group_by(Species_Broad) %>%
  summarize(
    mean_PDGFD_PTHLH_PVALB_pos_to_allInterneu = mean(PDGFD_PTHLH_PVALB_pos_to_allInterneu, na.rm = TRUE),
    sd_PDGFD_PTHLH_PVALB_pos_to_allInterneu = sd(PDGFD_PTHLH_PVALB_pos_to_allInterneu, na.rm = TRUE),
    n = sum(!is.na(PDGFD_PTHLH_PVALB_pos_to_allInterneu)),
    sem_PDGFD_PTHLH_PVALB_pos_to_allInterneu = sd_PDGFD_PTHLH_PVALB_pos_to_allInterneu / sqrt(n),
    .groups = 'drop'
  )
    
# Comparison with mouse (caudoputamen)
comps <- list(
  c('Non_primate', 'Primate')
)


# Plot with p-values
pdf(file = "GB_barplot_PDGFD_PTHLH_PVALB_pos_AllInterneuron_Ratio_Putamen_with_mouse_Caudoputamen_Across_BroadSpecies_with_t-testpvals_withFerret_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_caud, aes(x = Species_Broad, y = PDGFD_PTHLH_PVALB_pos_to_allInterneu, color = Species_Broad)) +
  geom_col(data = df_summary, aes(x = Species_Broad, y = mean_PDGFD_PTHLH_PVALB_pos_to_allInterneu, fill = Species_Broad), position = position_dodge(), width = 0.4) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species_Broad, ymin = mean_PDGFD_PTHLH_PVALB_pos_to_allInterneu - sem_PDGFD_PTHLH_PVALB_pos_to_allInterneu, ymax = mean_PDGFD_PTHLH_PVALB_pos_to_allInterneu + sem_PDGFD_PTHLH_PVALB_pos_to_allInterneu),
                position = position_dodge(width = 0.9), width = 0.1) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('PDGFD_PTHLH_PVALB+ / Interneuron Ratio in Putamen'), 
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(position = fixed_jitter, size = 2, alpha = 0.7, color = "black") + 
  stat_compare_means(comparisons = comps, method = 't.test', paired = FALSE) +
 # Use ggrepel for automatic label placement with more robust settings
  geom_text_repel(data = df_no_caud %>% filter(Species_Broad == "Non_primate" & PDGFD_PTHLH_PVALB_pos_to_allInterneu < 0.1),
                  aes(label = orig.ident),
                  size = 5, 
                  color = "black", 
                  box.padding = 0.2,  # Space around the labels
                  max.overlaps = 10, # Allow a certain number of overlaps if needed
                  direction = "both",  # Allow labels to be repelled in both x and y directions
                  nudge_y = 0.01, # Nudge labels a bit vertically
		  nudge_x = 0.01,  # Nudge labels a bit horizontally
                  seed = 42) # Set seed for consistent labeling
)
dev.off()

# Plot without p-values
pdf(file = "GB_barplot_PDGFD_PTHLH_PVALB_pos_AllInterneuron_Ratio_Putamen_with_mouse_Caudoputamen_Across_BroadSpecies_withFerret_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_caud, aes(x = Species_Broad, y = PDGFD_PTHLH_PVALB_pos_to_allInterneu, color = Species_Broad)) +
  geom_col(data = df_summary, aes(x = Species_Broad, y = mean_PDGFD_PTHLH_PVALB_pos_to_allInterneu, fill = Species_Broad), position = position_dodge(), width = 0.4) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species_Broad, ymin = mean_PDGFD_PTHLH_PVALB_pos_to_allInterneu - sem_PDGFD_PTHLH_PVALB_pos_to_allInterneu, ymax = mean_PDGFD_PTHLH_PVALB_pos_to_allInterneu + sem_PDGFD_PTHLH_PVALB_pos_to_allInterneu),
                position = position_dodge(width = 0.9), width = 0.1) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('PDGFD_PTHLH_PVALB+ / Interneuron Ratio in Putamen'), 
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black")
)
dev.off()
##############################################################
######
## Caudate
# TAC3
meta = second_non_spn@meta.data
meta$Species_Broad = 'Primate'
meta[meta$Species %in% c('Bat','Ferret','Mouse'), 'Species_Broad'] = 'Non_primate'

# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species_Broad) %>%
  summarize(
    totsize = n(),
    TAC3size = sum(newannot == 'TAC3'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    TAC3_to_allInterneu = TAC3size / totsize
  )

df$Species_Broad <- factor(df$Species_Broad, levels = rev(c('Primate', 'Non_primate')))

# Ratio for Caudate for all species and for mouse CaudoPutamen 
# Subset the dataset
df_no_put <- df[df$Tissue != "Putamen",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_put %>%
  group_by(Species_Broad) %>%
  summarize(
    mean_TAC3_to_allInterneu = mean(TAC3_to_allInterneu, na.rm = TRUE),
    sd_TAC3_to_allInterneu = sd(TAC3_to_allInterneu, na.rm = TRUE),
    n = sum(!is.na(TAC3_to_allInterneu)),
    sem_TAC3_to_allInterneu = sd_TAC3_to_allInterneu / sqrt(n),
    .groups = 'drop'
  )
    
# Comparison with mouse (caudoputamen)
comps <- list(
  c('Non_primate', 'Primate')
)

# Plot with p-values
pdf(file = "GB_barplot_TAC3_AllInterneuron_Ratio_Caudate_with_mouse_Caudoputamen_Across_BroadSpecies_with_t-testpvals_withFerret_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put, aes(x = Species_Broad, y = TAC3_to_allInterneu, color = Species_Broad)) +
  geom_col(data = df_summary, aes(x = Species_Broad, y = mean_TAC3_to_allInterneu, fill = Species_Broad), position = position_dodge(), width = 0.4) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species_Broad, ymin = mean_TAC3_to_allInterneu - sem_TAC3_to_allInterneu, ymax = mean_TAC3_to_allInterneu + sem_TAC3_to_allInterneu),
                position = position_dodge(width = 0.9), width = 0.1) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('TAC3 / Interneuron Ratio in Caudate'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") + 
  stat_compare_means(comparisons = comps, method = 't.test', paired = FALSE)
)
dev.off()

# Plot without p-values
pdf(file = "GB_barplot_TAC3_AllInterneuron_Ratio_Caudate_with_mouse_Caudoputamen_Across_BroadSpecies_withFerret_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put, aes(x = Species_Broad, y = TAC3_to_allInterneu, color = Species_Broad)) +
  geom_col(data = df_summary, aes(x = Species_Broad, y = mean_TAC3_to_allInterneu, fill = Species_Broad), position = position_dodge(), width = 0.4) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species_Broad, ymin = mean_TAC3_to_allInterneu - sem_TAC3_to_allInterneu, ymax = mean_TAC3_to_allInterneu + sem_TAC3_to_allInterneu),
                position = position_dodge(width = 0.9), width = 0.1) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('TAC3 / Interneuron Ratio in Caudate'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black")
)
dev.off()
##########################################################
## Putamen
# TAC3
meta = second_non_spn@meta.data
meta$Species_Broad = 'Primate'
meta[meta$Species %in% c('Bat','Ferret','Mouse'), 'Species_Broad'] = 'Non_primate'

# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species_Broad) %>%
  summarize(
    totsize = n(),
    TAC3size = sum(newannot == 'TAC3'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    TAC3_to_allInterneu = TAC3size / totsize
  )

df$Species_Broad <- factor(df$Species_Broad, levels = rev(c('Primate', 'Non_primate')))

# Ratio for Putamen for all species and for mouse CaudoPutamen 
# Subset the dataset
df_no_caud <- df[df$Tissue != "Caudate",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_caud %>%
  group_by(Species_Broad) %>%
  summarize(
    mean_TAC3_to_allInterneu = mean(TAC3_to_allInterneu, na.rm = TRUE),
    sd_TAC3_to_allInterneu = sd(TAC3_to_allInterneu, na.rm = TRUE),
    n = sum(!is.na(TAC3_to_allInterneu)),
    sem_TAC3_to_allInterneu = sd_TAC3_to_allInterneu / sqrt(n),
    .groups = 'drop'
  )
    
# Comparison with mouse (caudoputamen)
comps <- list(
  c('Non_primate', 'Primate')
)

# Plot with p-values
pdf(file = "GB_barplot_TAC3_AllInterneuron_Ratio_Putamen_with_mouse_Caudoputamen_Across_BroadSpecies_with_t-testpvals_withFerret_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_caud, aes(x = Species_Broad, y = TAC3_to_allInterneu, color = Species_Broad)) +
  geom_col(data = df_summary, aes(x = Species_Broad, y = mean_TAC3_to_allInterneu, fill = Species_Broad), position = position_dodge(), width = 0.4) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species_Broad, ymin = mean_TAC3_to_allInterneu - sem_TAC3_to_allInterneu, ymax = mean_TAC3_to_allInterneu + sem_TAC3_to_allInterneu),
                position = position_dodge(width = 0.9), width = 0.1) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('TAC3 / Interneuron Ratio in Putamen'), 
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") + 
  stat_compare_means(comparisons = comps, method = 't.test', paired = FALSE)
)
dev.off()

# Plot without p-values
pdf(file = "GB_barplot_TAC3_AllInterneuron_Ratio_Putamen_with_mouse_Caudoputamen_Across_BroadSpecies_withFerret_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_caud, aes(x = Species_Broad, y = TAC3_to_allInterneu, color = Species_Broad)) +
  geom_col(data = df_summary, aes(x = Species_Broad, y = mean_TAC3_to_allInterneu, fill = Species_Broad), position = position_dodge(), width = 0.4) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species_Broad, ymin = mean_TAC3_to_allInterneu - sem_TAC3_to_allInterneu, ymax = mean_TAC3_to_allInterneu + sem_TAC3_to_allInterneu),
                position = position_dodge(width = 0.9), width = 0.1) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('TAC3 / Interneuron Ratio in Putamen'), 
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black")
)
dev.off()
##############################################################
## Caudate
# PDGFD
meta = second_non_spn@meta.data
meta$Species_Broad = 'Primate'
meta[meta$Species %in% c('Bat','Ferret','Mouse'), 'Species_Broad'] = 'Non_primate'

# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species_Broad) %>%
  summarize(
    totsize = n(),
    PDGFDsize = sum(newannot == 'PDGFD'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    PDGFD_to_allInterneu = PDGFDsize / totsize
  )

df$Species_Broad <- factor(df$Species_Broad, levels = rev(c('Primate', 'Non_primate')))

# Ratio for Caudate for all species and for mouse CaudoPutamen 
# Subset the dataset
df_no_put <- df[df$Tissue != "Putamen",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_put %>%
  group_by(Species_Broad) %>%
  summarize(
    mean_PDGFD_to_allInterneu = mean(PDGFD_to_allInterneu, na.rm = TRUE),
    sd_PDGFD_to_allInterneu = sd(PDGFD_to_allInterneu, na.rm = TRUE),
    n = sum(!is.na(PDGFD_to_allInterneu)),
    sem_PDGFD_to_allInterneu = sd_PDGFD_to_allInterneu / sqrt(n),
    .groups = 'drop'
  )
    
# Comparison with mouse (caudoputamen)
comps <- list(
  c('Non_primate', 'Primate')
)

# Plot with p-values
pdf(file = "GB_barplot_PDGFD_AllInterneuron_Ratio_Caudate_with_mouse_Caudoputamen_Across_BroadSpecies_with_t-testpvals_withFerret_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put, aes(x = Species_Broad, y = PDGFD_to_allInterneu, color = Species_Broad)) + 
  geom_col(data = df_summary, aes(x = Species_Broad, y = mean_PDGFD_to_allInterneu, fill = Species_Broad), position = position_dodge(), width = 0.4) + 
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species_Broad, ymin = mean_PDGFD_to_allInterneu - sem_PDGFD_to_allInterneu, ymax = mean_PDGFD_to_allInterneu + sem_PDGFD_to_allInterneu),  
                position = position_dodge(width = 0.9), width = 0.1) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('PDGFD / Interneuron Ratio in Caudate'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") + 
  stat_compare_means(comparisons = comps, method = 't.test', paired = FALSE)
)
dev.off()

# Plot without p-values
pdf(file = "GB_barplot_PDGFD_AllInterneuron_Ratio_Caudate_with_mouse_Caudoputamen_Across_BroadSpecies_withFerret_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put, aes(x = Species_Broad, y = PDGFD_to_allInterneu, color = Species_Broad)) + 
  geom_col(data = df_summary, aes(x = Species_Broad, y = mean_PDGFD_to_allInterneu, fill = Species_Broad), position = position_dodge(), width = 0.4) + 
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species_Broad, ymin = mean_PDGFD_to_allInterneu - sem_PDGFD_to_allInterneu, ymax = mean_PDGFD_to_allInterneu + sem_PDGFD_to_allInterneu),  
                position = position_dodge(width = 0.9), width = 0.1) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('PDGFD / Interneuron Ratio in Caudate'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black")
)
dev.off()
##########################################################
## Putamen
# PDGFD
meta = second_non_spn@meta.data
meta$Species_Broad = 'Primate'
meta[meta$Species %in% c('Bat','Ferret','Mouse'), 'Species_Broad'] = 'Non_primate'

# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species_Broad) %>%
  summarize(
    totsize = n(),
    PDGFDsize = sum(newannot == 'PDGFD'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    PDGFD_to_allInterneu = PDGFDsize / totsize
  )

df$Species_Broad <- factor(df$Species_Broad, levels = rev(c('Primate', 'Non_primate')))

# Ratio for Putamen for all species and for mouse CaudoPutamen 
# Subset the dataset
df_no_caud <- df[df$Tissue != "Caudate",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_caud %>%
  group_by(Species_Broad) %>%
  summarize(
    mean_PDGFD_to_allInterneu = mean(PDGFD_to_allInterneu, na.rm = TRUE),
    sd_PDGFD_to_allInterneu = sd(PDGFD_to_allInterneu, na.rm = TRUE),
    n = sum(!is.na(PDGFD_to_allInterneu)),
    sem_PDGFD_to_allInterneu = sd_PDGFD_to_allInterneu / sqrt(n),
    .groups = 'drop'
  )
    
# Comparison with mouse (caudoputamen)
comps <- list(
  c('Non_primate', 'Primate')
)

# Plot with p-values
pdf(file = "GB_barplot_PDGFD_AllInterneuron_Ratio_Putamen_with_mouse_Caudoputamen_Across_BroadSpecies_with_t-testpvals_withFerret_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_caud, aes(x = Species_Broad, y = PDGFD_to_allInterneu, color = Species_Broad)) + 
  geom_col(data = df_summary, aes(x = Species_Broad, y = mean_PDGFD_to_allInterneu, fill = Species_Broad), position = position_dodge(), width = 0.4) + 
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species_Broad, ymin = mean_PDGFD_to_allInterneu - sem_PDGFD_to_allInterneu, ymax = mean_PDGFD_to_allInterneu + sem_PDGFD_to_allInterneu),  
                position = position_dodge(width = 0.9), width = 0.1) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('PDGFD / Interneuron Ratio in Putamen'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") + 
  stat_compare_means(comparisons = comps, method = 't.test', paired = FALSE)
)
dev.off()

# Plot without p-values
pdf(file = "GB_barplot_PDGFD_AllInterneuron_Ratio_Putamen_with_mouse_Caudoputamen_Across_BroadSpecies_withFerret_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_caud, aes(x = Species_Broad, y = PDGFD_to_allInterneu, color = Species_Broad)) + 
  geom_col(data = df_summary, aes(x = Species_Broad, y = mean_PDGFD_to_allInterneu, fill = Species_Broad), position = position_dodge(), width = 0.4) + 
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species_Broad, ymin = mean_PDGFD_to_allInterneu - sem_PDGFD_to_allInterneu, ymax = mean_PDGFD_to_allInterneu + sem_PDGFD_to_allInterneu),  
                position = position_dodge(width = 0.9), width = 0.1) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('PDGFD / Interneuron Ratio in Putamen'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black")
)
dev.off()
############################################################
## Caudate
# PDGFD_PTHLH_PVALBminus
meta = second_non_spn@meta.data
meta$Species_Broad = 'Primate'
meta[meta$Species %in% c('Bat','Ferret','Mouse'), 'Species_Broad'] = 'Non_primate'

# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species_Broad) %>%
  summarize(
    totsize = n(),
    PDGFD_PTHLH_PVALBminussize = sum(newannot == 'PDGFD_PTHLH_PVALB-'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    PDGFD_PTHLH_PVALBminus_to_allInterneu = PDGFD_PTHLH_PVALBminussize / totsize
  )

df$Species_Broad <- factor(df$Species_Broad, levels = rev(c('Primate', 'Non_primate')))

# Ratio for Caudate for all species and for mouse CaudoPutamen 
# Subset the dataset
df_no_put <- df[df$Tissue != "Putamen",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_put %>%
  group_by(Species_Broad) %>%
  summarize(
    mean_PDGFD_PTHLH_PVALBminus_to_allInterneu = mean(PDGFD_PTHLH_PVALBminus_to_allInterneu, na.rm = TRUE),
    sd_PDGFD_PTHLH_PVALBminus_to_allInterneu = sd(PDGFD_PTHLH_PVALBminus_to_allInterneu, na.rm = TRUE),
    n = sum(!is.na(PDGFD_PTHLH_PVALBminus_to_allInterneu)),
    sem_PDGFD_PTHLH_PVALBminus_to_allInterneu = sd_PDGFD_PTHLH_PVALBminus_to_allInterneu / sqrt(n),
    .groups = 'drop'
  )
    
# Comparison with mouse (caudoputamen)
comps <- list(
  c('Non_primate', 'Primate')
)

# Plot with p-values
pdf(file = "GB_barplot_PDGFD_PTHLH_PVALBminus_AllInterneuron_Ratio_Caudate_with_mouse_Caudoputamen_Across_BroadSpecies_with_t-testpvals_withFerret_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put, aes(x = Species_Broad, y = PDGFD_PTHLH_PVALBminus_to_allInterneu, color = Species_Broad)) + 
  geom_col(data = df_summary, aes(x = Species_Broad, y = mean_PDGFD_PTHLH_PVALBminus_to_allInterneu, fill = Species_Broad), position = position_dodge(), width = 0.4) + 
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species_Broad, ymin = mean_PDGFD_PTHLH_PVALBminus_to_allInterneu - sem_PDGFD_PTHLH_PVALBminus_to_allInterneu, ymax = mean_PDGFD_PTHLH_PVALBminus_to_allInterneu + sem_PDGFD_PTHLH_PVALBminus_to_allInterneu),  
                position = position_dodge(width = 0.9), width = 0.1) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('PDGFD_PTHLH_PVALBminus / Interneuron Ratio in Caudate'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") + 
  stat_compare_means(comparisons = comps, method = 't.test', paired = FALSE)
)
dev.off()

# Plot without p-values
pdf(file = "GB_barplot_PDGFD_PTHLH_PVALBminus_AllInterneuron_Ratio_Caudate_with_mouse_Caudoputamen_Across_BroadSpecies_withFerret_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put, aes(x = Species_Broad, y = PDGFD_PTHLH_PVALBminus_to_allInterneu, color = Species_Broad)) + 
  geom_col(data = df_summary, aes(x = Species_Broad, y = mean_PDGFD_PTHLH_PVALBminus_to_allInterneu, fill = Species_Broad), position = position_dodge(), width = 0.4) + 
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species_Broad, ymin = mean_PDGFD_PTHLH_PVALBminus_to_allInterneu - sem_PDGFD_PTHLH_PVALBminus_to_allInterneu, ymax = mean_PDGFD_PTHLH_PVALBminus_to_allInterneu + sem_PDGFD_PTHLH_PVALBminus_to_allInterneu),  
                position = position_dodge(width = 0.9), width = 0.1) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('PDGFD_PTHLH_PVALBminus / Interneuron Ratio in Caudate'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black")
)
dev.off()
##########################################################
## Putamen
# PDGFD_PTHLH_PVALBminus
meta = second_non_spn@meta.data
meta$Species_Broad = 'Primate'
meta[meta$Species %in% c('Bat','Ferret','Mouse'), 'Species_Broad'] = 'Non_primate'

# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species_Broad) %>%
  summarize(
    totsize = n(),
    PDGFD_PTHLH_PVALBminussize = sum(newannot == 'PDGFD_PTHLH_PVALB-'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    PDGFD_PTHLH_PVALBminus_to_allInterneu = PDGFD_PTHLH_PVALBminussize / totsize
  )

df$Species_Broad <- factor(df$Species_Broad, levels = rev(c('Primate', 'Non_primate')))

# Ratio for Putamen for all species and for mouse CaudoPutamen 
# Subset the dataset
df_no_caud <- df[df$Tissue != "Caudate",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_caud %>%
  group_by(Species_Broad) %>%
  summarize(
    mean_PDGFD_PTHLH_PVALBminus_to_allInterneu = mean(PDGFD_PTHLH_PVALBminus_to_allInterneu, na.rm = TRUE),
    sd_PDGFD_PTHLH_PVALBminus_to_allInterneu = sd(PDGFD_PTHLH_PVALBminus_to_allInterneu, na.rm = TRUE),
    n = sum(!is.na(PDGFD_PTHLH_PVALBminus_to_allInterneu)),
    sem_PDGFD_PTHLH_PVALBminus_to_allInterneu = sd_PDGFD_PTHLH_PVALBminus_to_allInterneu / sqrt(n),
    .groups = 'drop'
  )
    
# Comparison with mouse (caudoputamen)
comps <- list(
  c('Non_primate', 'Primate')
)

# Plot with p-values
pdf(file = "GB_barplot_PDGFD_PTHLH_PVALBminus_AllInterneuron_Ratio_Putamen_with_mouse_Caudoputamen_Across_BroadSpecies_with_t-testpvals_withFerret_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_caud, aes(x = Species_Broad, y = PDGFD_PTHLH_PVALBminus_to_allInterneu, color = Species_Broad)) + 
  geom_col(data = df_summary, aes(x = Species_Broad, y = mean_PDGFD_PTHLH_PVALBminus_to_allInterneu, fill = Species_Broad), position = position_dodge(), width = 0.4) + 
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species_Broad, ymin = mean_PDGFD_PTHLH_PVALBminus_to_allInterneu - sem_PDGFD_PTHLH_PVALBminus_to_allInterneu, ymax = mean_PDGFD_PTHLH_PVALBminus_to_allInterneu + sem_PDGFD_PTHLH_PVALBminus_to_allInterneu),  
                position = position_dodge(width = 0.9), width = 0.1) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('PDGFD_PTHLH_PVALBminus / Interneuron Ratio in Putamen'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") + 
  stat_compare_means(comparisons = comps, method = 't.test', paired = FALSE)
)
dev.off()

# Plot without p-values
pdf(file = "GB_barplot_PDGFD_PTHLH_PVALBminus_AllInterneuron_Ratio_Putamen_with_mouse_Caudoputamen_Across_BroadSpecies_withFerret_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_caud, aes(x = Species_Broad, y = PDGFD_PTHLH_PVALBminus_to_allInterneu, color = Species_Broad)) + 
  geom_col(data = df_summary, aes(x = Species_Broad, y = mean_PDGFD_PTHLH_PVALBminus_to_allInterneu, fill = Species_Broad), position = position_dodge(), width = 0.4) + 
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species_Broad, ymin = mean_PDGFD_PTHLH_PVALBminus_to_allInterneu - sem_PDGFD_PTHLH_PVALBminus_to_allInterneu, ymax = mean_PDGFD_PTHLH_PVALBminus_to_allInterneu + sem_PDGFD_PTHLH_PVALBminus_to_allInterneu),  
                position = position_dodge(width = 0.9), width = 0.1) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('PDGFD_PTHLH_PVALBminus / Interneuron Ratio in Putamen'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black")
)
dev.off()
