# load packages
library(patchwork)
library(Seurat)
library(rhdf5)
library(DropletUtils)
library(DropletQC)
library(dplyr)
library(Matrix)
library(ggplot2)
library(plyr)
library(tidyverse)
library(tidyr)
library(ggpubr)
library(reshape2)
library(rio)
library(data.table)
library(harmony)
library(ggrepel)
library(edgeR)
library(variancePartition)
library(Matrix.utils)
library(SingleCellExperiment)
library(RColorBrewer)
set.seed(1234)
source("/project/Neuroinformatics_Core/Konopka_lab/s422071/SCRIPTS_pr/SCRIPTS/utility_functions.R")
library(curl)
#conda install bioconda::r-wgcna
#BiocManager::install('WGCNA')
library(WGCNA)
library(flashClust)
#remotes::install_url("https://cran.r-project.org/src/contrib/Archive/matrixStats/matrixStats_1.1.0.tar.gz")
library(matrixStats)
library(RColorBrewer)

####
## PREPARE DATASETS
####
# read the data
human_str <- readRDS(paste0(file_dir,"human_integrated_caudate_putamen_ANNOTATED.RDS"))
chimp_str = readRDS(paste0(file_dir,"chimp_integrated_caudate_putamen_ANNOTATED.RDS"))
macaque_str = readRDS(paste0(file_dir,"macaque_integrated_caudate_putamen_ANNOTATED.RDS"))
marmoset_str = readRDS(paste0(file_dir,"marmoset_integrated_caudate_putamen_ANNOTATED.RDS"))
mouse_str = readRDS(paste0(file_dir,"Mouse_Caudate_Annotated_FINAL.RDS"))
mouse_str$Tissue = rep("Caudoputamen", nrow(mouse_str[[]]))
mouse_str$Species = rep("Mouse", nrow(mouse_str[[]]))
mouse_str$id = paste0(mouse_str$orig.ident, mouse_str$Tissue)
bat_str = readRDS(paste0(file_dir,"bat_integrated_caudate_putamen_ANNOTATED.RDS"))
ferret_str = readRDS(paste0(file_dir,"Ferret_Caudate_Krienen_ANNOTATED.RDS"))
ferret_str$newannot_2 = ferret_str$newannot
ferret_str$newannot = ferret_str$broad_annot


# Reshape and combine metadata
human_meta = human_str@meta.data
human_meta$Species = 'Human'

chimp_meta = chimp_str@meta.data
chimp_meta$Species = 'Chimp'

macaque_meta = macaque_str@meta.data
macaque_meta$Species = 'Macaque'

marmoset_meta = marmoset_str@meta.data
marmoset_meta$Species = 'Marmoset'

bat_meta = bat_str@meta.data
bat_meta$Species = 'Bat'

sub_obj_new_list <- list(
  Human    = human_str,
  Chimp    = chimp_str,
  Macaque  = macaque_str,
  Marmoset = marmoset_str,
  Mouse    = mouse_str,
  Bat      = bat_str,
  Ferret   = ferret_str
)

saveRDS(sub_obj_new_list, file = "/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/seu_objs/sub_obj_new_list_11_12_2025.RDS")

# make RNA as active assay
a <- lapply(sub_obj_new_list, function(x) {
  DefaultAssay(x) <- "RNA"
  x@assays$SCT = NULL
  x[["RNA"]] <- as(object = x[["RNA"]], Class = "Assay")
  x
})

b_new = Reduce(merge, a)

### microglia cleanup
seurM = subset(b_new, subset = newannot_2 == "Microglia")
               
# Found ortholog genes from human pr coding genes and extracted the genes which are ortholog in all 7 species using ncbi datasets tool
ortho_genes <- read.table("/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/metadata/orthologs_in_7_species_human_pr_codingOnly_2.csv",  header = TRUE)$symbol #15055 pr coding ortho genes

### subset the genes and the Microglia from each Species
# subset cells and metadata
# Extract the count matrix from the Seurat object (default is "RNA" assay)
mat = seurM@assays$RNA$counts
new_mat <- mat[rownames(mat) %in% ortho_genes,]

# generate new Seurat object
glia_merged <- CreateSeuratObject(count = new_mat, meta.data = seurM[[]])

####
## SET VARIABLES
####
marksToPlot <- c('NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'APOE', 'GPC5', 'GLI3', 'AQP4', 'CSF1R', 'ARHGAP15', 'PLDL1', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CFTR', 'CHAT', 'ADARB2', 'LHX6', 'MEIS2', 'TAC3', 'VIP', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA', 'SATB2', 'SLC17A7', 'EYA2', 'PDGFD', 'PTHLH', 'PVALB', 'SPARCL1', 'RSPO2', 'RMST', 'LMO3', 'TSHZ2', "CALB1", "CALB2", "NOS1", "TSHZ1", "OPRM1", "CHST9","GRM8", "PPP1R1B", 'GRIK3', 'CXCL14')

pref = 'Microglia_AllSpecies_AllTissues_human_pr_coding_orthologs'

####
## INTEGRATE ACROSS SPECIES
####
# Split data
seurM = glia_merged
seurM[["RNA"]] = as(object = seurM[["RNA"]], Class = "Assay")
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
    x = RunPCA(x, features = features, verbose = FALSE, npcs = 30)
})

anchors = FindIntegrationAnchors(object.list = seurML, anchor.features = features,  normalization.method = "LogNormalize", reduction = "rpca")

# save memory prior integration
gc()

# this command creates an 'integrated' data assay
options(future.globals.maxSize = 16000 * 1024^2)
allseur_integrated = IntegrateData(anchorset = anchors, k.weight = 30, normalization.method = 'LogNormalize')

saveRDS(allseur_integrated, '/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/seu_objs/Integrated_rpca_Microglia_AllSpecies_ALLTissues_human_pr_coding_orthologs.RDS')

# Clustering
allseur_integrated = ScaleData(allseur_integrated, verbose = FALSE)
allseur_integrated = RunPCA(allseur_integrated, verbose = FALSE)
allseur_integrated = RunUMAP(allseur_integrated, dims = 1:20, reduction = 'pca')

# Basic Plots
setwd("/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/01_MICROGLIA/CLUSTERING_1")
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

saveRDS(allseur_integrated, paste0('/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/01_MICROGLIA/CLUSTERING_1/', pref, '_integrated_CLUSTERING1.RDS'))
allseur_integrated = readRDS( paste0('/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/01_MICROGLIA/CLUSTERING_1/', pref, '_integrated_CLUSTERING1.RDS'))
################# CLUSTERING_2 ##########################
##Cluster#,CellType
#12*,MOL+micro
#19*,Neu+Micro
#21*,Endo+micro

setwd("/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/01_MICROGLIA/CLUSTERING_2/")
toremove = c(12,19,21)
allseur_integrated = subset(allseur_integrated, subset = seurat_clusters %in% toremove, invert = T)

# Clustering
DefaultAssay(allseur_integrated) = 'integrated'
allseur_integrated = ScaleData(allseur_integrated, verbose = FALSE)
allseur_integrated = RunPCA(allseur_integrated, verbose = FALSE)
allseur_integrated = RunUMAP(allseur_integrated, dims = 1:20, reduction = 'pca')

# Basic Plots
pdf(paste0(pref, "_UMAP.pdf"))
DimPlot(allseur_integrated, group.by = 'Species', raster = T)
dev.off()

allseur_integrated = FindNeighbors(allseur_integrated, dims = 1:20, reduction = 'pca')
allseur_integrated = FindClusters(allseur_integrated, resolution = 1)

pdf(paste0(pref, "_CLUSTERS.pdf"))
DimPlot(allseur_integrated, label = T, raster = T) + NoLegend()
dev.off()

stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_STACK_SAMPLES'), wd = 20)
stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'Species', paste0(pref, '_STACK_SPECIES'))

# Plot previously identified markers
DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)

pdf(paste0(pref, "INTEGRATED_DOTPLOT.pdf"), width = 20, height = 10)
DotPlot(allseur_integrated, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "INTEGRATED_DEPTH.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'black') + NoLegend() +
rotate_x_text(90)
dev.off()


pdf(paste0(pref, "INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# save
saveRDS(allseur_integrated, paste0(pref,"_integrated_CLUSTERING_2.RDS"))
allseur_integrated = readRDS(paste0(pref,"_integrated_CLUSTERING_2.RDS"))
##################### CLUSTERING_3 ################# 
#14*, MOL+micro
#19*, MOL+micro

setwd("/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/01_MICROGLIA/CLUSTERING_3/")
toremove = c(14,19)
allseur_integrated = subset(allseur_integrated, subset = seurat_clusters %in% toremove, invert = T)

# Clustering
DefaultAssay(allseur_integrated) = 'integrated'
allseur_integrated = ScaleData(allseur_integrated, verbose = FALSE)
allseur_integrated = RunPCA(allseur_integrated, verbose = FALSE)
allseur_integrated = RunUMAP(allseur_integrated, dims = 1:20, reduction = 'pca')

# Basic Plots
pdf(paste0(pref, "_UMAP.pdf"))
DimPlot(allseur_integrated, group.by = 'Species', raster = T)
dev.off()


allseur_integrated = FindNeighbors(allseur_integrated, dims = 1:20, reduction = 'pca')
allseur_integrated = FindClusters(allseur_integrated, resolution = 1)

pdf(paste0(pref, "_CLUSTERS.pdf"))
DimPlot(allseur_integrated, label = T, raster = T) + NoLegend()
dev.off()

stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_STACK_SAMPLES'), wd = 20)
stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'Species', paste0(pref, '_STACK_SPECIES'))

# Plot previously identified markers
DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)

pdf(paste0(pref, "INTEGRATED_DOTPLOT.pdf"), width = 20, height = 10)
DotPlot(allseur_integrated, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "INTEGRATED_DEPTH.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'black') + NoLegend() +
rotate_x_text(90)
dev.off()


pdf(paste0(pref, "INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# save
saveRDS(allseur_integrated, paste0(pref,"_integrated_CLUSTERING_3.RDS"))
allseur_integrated = readRDS(paste0(pref,"_integrated_CLUSTERING_3.RDS"))
################### CLUSTERING_4 ###################
#15*,Endo

setwd("/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/01_MICROGLIA/CLUSTERING_4/")
toremove = c(15)
allseur_integrated = subset(allseur_integrated, subset = seurat_clusters %in% toremove, invert = T)

# Clustering
DefaultAssay(allseur_integrated) = 'integrated'
allseur_integrated = ScaleData(allseur_integrated, verbose = FALSE)
allseur_integrated = RunPCA(allseur_integrated, verbose = FALSE)
allseur_integrated = RunUMAP(allseur_integrated, dims = 1:20, reduction = 'pca')

# Basic Plots
pdf(paste0(pref, "_UMAP.pdf"))
DimPlot(allseur_integrated, group.by = 'Species', raster = T)
dev.off()


allseur_integrated = FindNeighbors(allseur_integrated, dims = 1:20, reduction = 'pca')
allseur_integrated = FindClusters(allseur_integrated, resolution = 1)

pdf(paste0(pref, "_CLUSTERS.pdf"))
DimPlot(allseur_integrated, label = T, raster = T) + NoLegend()
dev.off()

stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_STACK_SAMPLES'), wd = 20)
stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'Species', paste0(pref, '_STACK_SPECIES'))

# Plot previously identified markers
DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)

pdf(paste0(pref, "INTEGRATED_DOTPLOT.pdf"), width = 20, height = 10)
DotPlot(allseur_integrated, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "INTEGRATED_DEPTH.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'black') + NoLegend() +
rotate_x_text(90)
dev.off()


pdf(paste0(pref, "INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# save
saveRDS(allseur_integrated, paste0(pref,"_integrated_CLUSTERING_4.RDS"))

## annotate
(mapnames <-setNames(rep("Microglia", 27),c(0:26)))

# Create the broadannot
allseur_integrated$Micro_annot = unname(mapnames[allseur_integrated[["seurat_clusters"]][,1]]) # to extract the annotations and not the cluster #s

table(allseur_integrated$Micro_annot)
#Microglia 
#    32961 

setwd("/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/01_MICROGLIA/ANNOTATED")

# QC PLOTS
pdf(paste0(pref, "_Clusters.pdf"))
DimPlot(allseur_integrated, group.by = 'Micro_annot', label = T, raster = T, label.size = 7) +
NoLegend() +
ggtitle('')
dev.off()

pdf(paste0(pref, "_Depth_Gene.pdf"))
ggboxplot(allseur_integrated[[]], x = 'Micro_annot', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_Depth_UMI.pdf"))
ggboxplot(allseur_integrated[[]], x = 'Micro_annot', y = 'nCount_RNA', ylim = c(0,20000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_NuclearFraction.pdf"))
ggboxplot(allseur_integrated[[]], x = 'Micro_annot', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# MARKER GENE PLOTS
DefaultAssay(allseur_integrated) = 'RNA'
seurM_corr = NormalizeData(allseur_integrated)

pdf(paste0(pref, "_MarkersDotPlot.pdf"), width = 20, height = 10)
DotPlot(allseur_integrated, group.by = "Micro_annot", features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

stackedbarplot(allseur_integrated[[]], groupx = 'Micro_annot', groupfill = 'orig.ident', paste0(pref, '_STACK_SAMPLES'), wd = 20)
stackedbarplot(allseur_integrated[[]], groupx = 'Micro_annot', groupfill = 'Species', paste0(pref, '_STACK_SPECIES'))

# save
saveRDS(allseur_integrated, paste0(pref,"_integrated_ANNOTATED.RDS"))
## add these annotations to the original seu obj
# extract the newest annotations
b$detailed_annot <- as.character(b$newannot_2)

# Match cell names between objects
common_cells <- intersect(rownames(allseur_integrated@meta.data), rownames(b@meta.data))

# Update annotations for those cells
b$detailed_annot[common_cells] <- allseur_integrated$Micro_annot[common_cells]
b$is_micro = ifelse(rownames(b@meta.data) %in% rownames(allseur_integrated@meta.data), "yes", "no")

# Identify microglia cells in b that are not in the integrated object
b_keep = subset(b, subset = newannot_2 != "Microglia" | is_micro == "yes")
############################################
####### MOL + OPC
################################################
### subset the genes and the MOL+OPC from each Species
# subset cells and metadata
subset_b = subset(b_keep, subset = newannot_2 %in% c("MOL", "OPC"))  # already orthologs only
  
# Extract the count matrix from the Seurat object (default is "RNA" assay)
mat = subset_b@assays$RNA$counts
new_mat <- mat[rownames(mat) %in% ortho_genes,]

# generate new Seurat object
glia_merged <- CreateSeuratObject(count = new_mat, meta.data = subset_b[[]])

####
## SET VARIABLES
####
# PPP1R1B IS AN SPN MARKER!
# TSHZ1 IS A PATCH MARKER
marksToPlot <- c('NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'ARHGAP15', 'PLD1','CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CFTR', 'CHAT', 'ADARB2', 'LHX6', 'MEIS2', 'TAC3', 'VIP', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA', 'SATB2', 'SLC17A7', 'EYA2', 'PDGFD', 'PTHLH', 'PVALB', 'SPARCL1', 'RSPO2', 'RMST', 'LMO3', 'TSHZ2', "CALB1", "CALB2", "NOS1", "TSHZ1", "OPRM1", "CHST9","GRM8", "PPP1R1B", 'GRIK3', 'CXCL14')

pref = 'MOL_OPC_AllSpecies_AllTissues_human_pr_coding_orthologs'

####
## INTEGRATE ACROSS SPECIES
####
#### use old seurat version
# Split data
seurM = glia_merged
seurM[["RNA"]] = as(object = seurM[["RNA"]], Class = "Assay")
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
    x = RunPCA(x, features = features, verbose = FALSE, npcs = 30)
})

anchors = FindIntegrationAnchors(object.list = seurML, anchor.features = features,  normalization.method = "LogNormalize", reduction = "rpca")

# save memory prior integration
gc()

# this command creates an 'integrated' data assay
#(kweg = floor(min(table(seurM$orig.ident))/10)*10) #90
options(future.globals.maxSize = 16000 * 1024^2)
allseur_integrated = IntegrateData(anchorset = anchors, k.weight = 30, normalization.method = 'LogNormalize')

saveRDS(allseur_integrated, '/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/seu_objs/Integrated_rpca_MOL_OPC_AllSpecies_ALLTissues_human_pr_coding_orthologs.RDS')

# Clustering
allseur_integrated = ScaleData(allseur_integrated, verbose = FALSE)
allseur_integrated = RunPCA(allseur_integrated, verbose = FALSE)
allseur_integrated = RunUMAP(allseur_integrated, dims = 1:20, reduction = 'pca')

# Basic Plots
setwd("/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/01_MOL_OPC/CLUSTERING_1")
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

saveRDS(allseur_integrated, paste0('/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/01_MOL_OPC/CLUSTERING_1/', pref, '_integrated_CLUSTERING1.RDS'))
################# CLUSTERING_2 ##########################
##Cluster#,CellType
#17*,Neu+MOL
#19*,MOL+Ast
#25*,Neu+MOL

setwd("/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/01_MOL_OPC/CLUSTERING_2/")
toremove = c(17,19,25)
allseur_integrated = subset(allseur_integrated, subset = seurat_clusters %in% toremove, invert = T)

# Clustering
DefaultAssay(allseur_integrated) = 'integrated'
allseur_integrated = ScaleData(allseur_integrated, verbose = FALSE)
allseur_integrated = RunPCA(allseur_integrated, verbose = FALSE)
allseur_integrated = RunUMAP(allseur_integrated, dims = 1:20, reduction = 'pca')

# Basic Plots
pdf(paste0(pref, "_UMAP.pdf"))
DimPlot(allseur_integrated, group.by = 'Species', raster = T)
dev.off()

allseur_integrated = FindNeighbors(allseur_integrated, dims = 1:20, reduction = 'pca')
allseur_integrated = FindClusters(allseur_integrated, resolution = 1)

pdf(paste0(pref, "_CLUSTERS.pdf"))
DimPlot(allseur_integrated, label = T, raster = T) + NoLegend()
dev.off()

stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_STACK_SAMPLES'), wd = 20)
stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'Species', paste0(pref, '_STACK_SPECIES'))

# Plot previously identified markers
DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)

pdf(paste0(pref, "INTEGRATED_DOTPLOT.pdf"), width = 20, height = 10)
DotPlot(allseur_integrated, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "INTEGRATED_DEPTH.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'black') + NoLegend() +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# save
saveRDS(allseur_integrated, paste0(pref,"_integrated_CLUSTERING_2.RDS"))
allseur_integrated= readRDS(paste0(pref,"_integrated_CLUSTERING_2.RDS"))
##################### CLUSTERING_3 ################# 
#21*, MOL+micro
#24*, neu+MOL

setwd("/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/01_MOL_OPC/CLUSTERING_3/")
toremove = c(21,24)
allseur_integrated = subset(allseur_integrated, subset = seurat_clusters %in% toremove, invert = T)

# Clustering
DefaultAssay(allseur_integrated) = 'integrated'
allseur_integrated = ScaleData(allseur_integrated, verbose = FALSE)
allseur_integrated = RunPCA(allseur_integrated, verbose = FALSE)
allseur_integrated = RunUMAP(allseur_integrated, dims = 1:20, reduction = 'pca')

# Basic Plots
pdf(paste0(pref, "_UMAP.pdf"))
DimPlot(allseur_integrated, group.by = 'Species', raster = T)
dev.off()

allseur_integrated = FindNeighbors(allseur_integrated, dims = 1:20, reduction = 'pca')
allseur_integrated = FindClusters(allseur_integrated, resolution = 1)

pdf(paste0(pref, "_CLUSTERS.pdf"))
DimPlot(allseur_integrated, label = T, raster = T) + NoLegend()
dev.off()

stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_STACK_SAMPLES'), wd = 20)
stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'Species', paste0(pref, '_STACK_SPECIES'))

# Plot previously identified markers
DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)

pdf(paste0(pref, "INTEGRATED_DOTPLOT.pdf"), width = 20, height = 10)
DotPlot(allseur_integrated, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "INTEGRATED_DEPTH.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'black') + NoLegend() +
rotate_x_text(90)
dev.off()


pdf(paste0(pref, "INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# save
saveRDS(allseur_integrated, paste0(pref,"_integrated_CLUSTERING_3.RDS"))
############### CLUSTERING_4 #################
#21*,MOL+OPC

setwd("/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/01_MOL_OPC/CLUSTERING_4/")
toremove = c(21)
allseur_integrated = subset(allseur_integrated, subset = seurat_clusters %in% toremove, invert = T)

# Clustering
DefaultAssay(allseur_integrated) = 'integrated'
allseur_integrated = ScaleData(allseur_integrated, verbose = FALSE)
allseur_integrated = RunPCA(allseur_integrated, verbose = FALSE)
allseur_integrated = RunUMAP(allseur_integrated, dims = 1:20, reduction = 'pca')

# Basic Plots
pdf(paste0(pref, "_UMAP.pdf"))
DimPlot(allseur_integrated, group.by = 'Species', raster = T)
dev.off()


allseur_integrated = FindNeighbors(allseur_integrated, dims = 1:20, reduction = 'pca')
allseur_integrated = FindClusters(allseur_integrated, resolution = 1)

pdf(paste0(pref, "_CLUSTERS.pdf"))
DimPlot(allseur_integrated, label = T, raster = T) + NoLegend()
dev.off()

stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_STACK_SAMPLES'), wd = 20)
stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'Species', paste0(pref, '_STACK_SPECIES'))

# Plot previously identified markers
DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)

pdf(paste0(pref, "INTEGRATED_DOTPLOT.pdf"), width = 20, height = 10)
DotPlot(allseur_integrated, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "INTEGRATED_DEPTH.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'black') + NoLegend() +
rotate_x_text(90)
dev.off()


pdf(paste0(pref, "INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# save
saveRDS(allseur_integrated, paste0(pref,"_integrated_CLUSTERING_4.RDS"))
######################### ANNOTATE #################
###Custer#,cell_type_annotation
#0,MOL
#1,MOL
#2,MOL
#3,MOL
#4,MOL
#5,MOL
#6,MOL
#7,MOL
#8,MOL
#9,MOL
#10,OPC
#11,MOL
#12,MOL
#13,MOL
#14,OPC
#15,MOL
#16,MOL
#17,OPC
#18,OPC
#19,MOL
#20,MOL
#21,COP
#22,OPC
#23,OPC

## annotate
(mapnames <-setNames(c("MOL","MOL","MOL","MOL","MOL","MOL","MOL","MOL","MOL",
"MOL","OPC","MOL","MOL","MOL","OPC","MOL","MOL","OPC","OPC","MOL",
"MOL","COP","OPC","OPC"),c(0:23)))

# Create the broadannot
allseur_integrated$MOL_OPC_annot = unname(mapnames[allseur_integrated[["seurat_clusters"]][,1]]) # to extract the annotations and not the cluster #s

table(allseur_integrated$MOL_OPC_annot)

#   COP    MOL    OPC 
#   447 263193  28394 

setwd("/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/01_MOL_OPC/ANNOTATED")

# QC PLOTS
pdf(paste0(pref, "_Clusters.pdf"))
DimPlot(allseur_integrated, group.by = 'MOL_OPC_annot', label = T, raster = T, label.size = 7) +
NoLegend() +
ggtitle('')
dev.off()

pdf(paste0(pref, "_Depth_Gene.pdf"))
ggboxplot(allseur_integrated[[]], x = 'MOL_OPC_annot', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_Depth_UMI.pdf"))
ggboxplot(allseur_integrated[[]], x = 'MOL_OPC_annot', y = 'nCount_RNA', ylim = c(0,20000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_NuclearFraction.pdf"))
ggboxplot(allseur_integrated[[]], x = 'MOL_OPC_annot', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# MARKER GENE PLOTS
DefaultAssay(allseur_integrated) = 'RNA'
seurM_corr = NormalizeData(allseur_integrated)

pdf(paste0(pref, "_MarkersDotPlot.pdf"), width = 20, height = 10)
DotPlot(allseur_integrated, group.by = "MOL_OPC_annot", features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

stackedbarplot(allseur_integrated[[]], groupx = 'MOL_OPC_annot', groupfill = 'orig.ident', paste0(pref, '_STACK_SAMPLES'), wd = 20)
stackedbarplot(allseur_integrated[[]], groupx = 'MOL_OPC_annot', groupfill = 'Species', paste0(pref, '_STACK_SPECIES'))

# save
saveRDS(allseur_integrated, paste0(pref,"_integrated_ANNOTATED.RDS"))
### put these annots to the original seurat obj
b = readRDS("/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/seu_objs/merged_clean_AllCellTypes_AllTissues_seu.RDS")

# Match cell names between objects
common_cells <- intersect(rownames(allseur_integrated@meta.data), rownames(b@meta.data))

# Update annotations for those cells
b$detailed_annot = b$newannot_2
b$detailed_annot[common_cells] <- allseur_integrated$MOL_OPC_annot[common_cells]

b$is_mol = ifelse(rownames(b@meta.data) %in% common_cells, "yes", "no")

# Identify cells in b that are not in the integrated object
b_keep <- subset(
  b,
  subset = !(newannot_2 %in% c("MOL", "OPC", "COP")) | is_mol == "yes")
)

# check
Idents(b_keep) = b_keep$detailed_annot
DotPlot(b_keep, features = marksToPlot) + rotate_x_text(45)
# save
saveRDS(b_keep,"/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/seu_objs/merged_clean_AllCellTypes_AllTissues_seu.RDS")
############################################
#######  Astrocyte
################################################
### subset the genes and the MOL+OPC from each Species
# subset cells and metadata
# Subset Seurat object to only include 'Microglia' cells and exclude the specified 'id' values bcs they have <10 microglial cells
seurM = subset(b_keep, subset = newannot_2 %in% c("Astrocyte"))  # already orthologs only
  
# Extract the count matrix from the Seurat object (default is "RNA" assay)
mat = seurM@assays$RNA$counts
new_mat <- mat[rownames(mat) %in% ortho_genes,]

# generate new Seurat object
glia_merged <- CreateSeuratObject(count = new_mat, meta.data = subset_b[[]])

####
## SET VARIABLES
####
marksToPlot <- c('NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'ARHGAP15', 'PLD1','CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CFTR', 'CHAT', 'ADARB2', 'LHX6', 'MEIS2', 'TAC3', 'VIP', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA', 'SATB2', 'SLC17A7', 'EYA2', 'PDGFD', 'PTHLH', 'PVALB', 'SPARCL1', 'RSPO2', 'RMST', 'LMO3', 'TSHZ2', "CALB1", "CALB2", "NOS1", "TSHZ1", "OPRM1", "CHST9","GRM8", "PPP1R1B", 'GRIK3', 'CXCL14')

pref = 'Astrocyte_AllSpecies_AllTissues_human_pr_coding_orthologs'

####
## INTEGRATE ACROSS SPECIES
####
#### use old seurat version
# Split data
seurM = glia_merged
seurM[["RNA"]] = as(object = seurM[["RNA"]], Class = "Assay")
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
    x = RunPCA(x, features = features, verbose = FALSE, npcs = 30)
})

anchors = FindIntegrationAnchors(object.list = seurML, anchor.features = features,  normalization.method = "LogNormalize", reduction = "rpca")

# save memory prior integration
gc()

# this command creates an 'integrated' data assay
#(kweg = floor(min(table(seurM$orig.ident))/10)*10) #90
options(future.globals.maxSize = 16000 * 1024^2)
allseur_integrated = IntegrateData(anchorset = anchors, k.weight = 30, normalization.method = 'LogNormalize')

saveRDS(allseur_integrated, '/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/seu_objs/Integrated_rpca_Astrocyte_AllSpecies_ALLTissues_human_pr_coding_orthologs.RDS')

# Clustering
allseur_integrated = ScaleData(allseur_integrated, verbose = FALSE)
allseur_integrated = RunPCA(allseur_integrated, verbose = FALSE)
allseur_integrated = RunUMAP(allseur_integrated, dims = 1:20, reduction = 'pca')

# Basic Plots
setwd("/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/01_ASTROCYTE/CLUSTERING_1")
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

saveRDS(allseur_integrated, paste0('/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/01_MOL_OPC/CLUSTERING_1/', pref, '_integrated_CLUSTERING1.RDS'))
################# CLUSTERING_2 ##########################
##Cluster#,CellType
#10*,MOL+Ast
#12*,MOL+Ast
#13*,Neu+Ast
#15*,MOL+Ast
#21*,MOL+Ast

setwd("/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/01_ASTROCYTE/CLUSTERING_2/")
toremove = c(10,12,13,15,21)
allseur_integrated = subset(allseur_integrated, subset = seurat_clusters %in% toremove, invert = T)

# Clustering
DefaultAssay(allseur_integrated) = 'integrated'
allseur_integrated = ScaleData(allseur_integrated, verbose = FALSE)
allseur_integrated = RunPCA(allseur_integrated, verbose = FALSE)
allseur_integrated = RunUMAP(allseur_integrated, dims = 1:20, reduction = 'pca')

# Basic Plots
pdf(paste0(pref, "_UMAP.pdf"))
DimPlot(allseur_integrated, group.by = 'Species', raster = T)
dev.off()

allseur_integrated = FindNeighbors(allseur_integrated, dims = 1:20, reduction = 'pca')
allseur_integrated = FindClusters(allseur_integrated, resolution = 1)

pdf(paste0(pref, "_CLUSTERS.pdf"))
DimPlot(allseur_integrated, label = T, raster = T) + NoLegend()
dev.off()

stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_STACK_SAMPLES'), wd = 20)
stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'Species', paste0(pref, '_STACK_SPECIES'))

# Plot previously identified markers
DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)

pdf(paste0(pref, "INTEGRATED_DOTPLOT.pdf"), width = 20, height = 10)
DotPlot(allseur_integrated, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "INTEGRATED_DEPTH.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'black') + NoLegend() +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# save
saveRDS(allseur_integrated, paste0(pref,"_integrated_CLUSTERING_2.RDS"))
################# CLUSTERING_2 ##########################
##Cluster#,CellType
#14*,Neu+Ast
#15*,MOL+Ast
#17*,Endo+Ast

setwd("/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/01_ASTROCYTE/CLUSTERING_3/")
toremove = c(14,15,17)
allseur_integrated = subset(allseur_integrated, subset = seurat_clusters %in% toremove, invert = T)

# Clustering
DefaultAssay(allseur_integrated) = 'integrated'
allseur_integrated = ScaleData(allseur_integrated, verbose = FALSE)
allseur_integrated = RunPCA(allseur_integrated, verbose = FALSE)
allseur_integrated = RunUMAP(allseur_integrated, dims = 1:20, reduction = 'pca')

# Basic Plots
pdf(paste0(pref, "_UMAP.pdf"))
DimPlot(allseur_integrated, group.by = 'Species', raster = T)
dev.off()

allseur_integrated = FindNeighbors(allseur_integrated, dims = 1:20, reduction = 'pca')
allseur_integrated = FindClusters(allseur_integrated, resolution = 1)

pdf(paste0(pref, "_CLUSTERS.pdf"))
DimPlot(allseur_integrated, label = T, raster = T) + NoLegend()
dev.off()

stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_STACK_SAMPLES'), wd = 20)
stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'Species', paste0(pref, '_STACK_SPECIES'))

# Plot previously identified markers
DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)

pdf(paste0(pref, "INTEGRATED_DOTPLOT.pdf"), width = 20, height = 10)
DotPlot(allseur_integrated, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

pdf(paste0(pref, "INTEGRATED_DEPTH.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'black') + NoLegend() +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# save
saveRDS(allseur_integrated, paste0(pref,"_integrated_CLUSTERING_3.RDS"))
############ ANNOTATE #############
(mapnames <-setNames(rep("Astrocyte", 22),c(0:21)))

# Create the broadannot
allseur_integrated$Ast_annot = unname(mapnames[allseur_integrated[["seurat_clusters"]][,1]]) # to extract the annotations and not the cluster #s

table(allseur_integrated$Ast_annot)
#Astrocyte 
#    64660 

setwd("/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/01_ASTROCYTE/ANNOTATED")

# QC PLOTS
pdf(paste0(pref, "_Clusters.pdf"))
DimPlot(allseur_integrated, group.by = 'Ast_annot', label = T, raster = T, label.size = 7) +
NoLegend() +
ggtitle('')
dev.off()

pdf(paste0(pref, "_Depth_Gene.pdf"))
ggboxplot(allseur_integrated[[]], x = 'Ast_annot', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_Depth_UMI.pdf"))
ggboxplot(allseur_integrated[[]], x = 'Ast_annot', y = 'nCount_RNA', ylim = c(0,20000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_NuclearFraction.pdf"))
ggboxplot(allseur_integrated[[]], x = 'Ast_annot', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# MARKER GENE PLOTS
DefaultAssay(allseur_integrated) = 'RNA'
seurM_corr = NormalizeData(allseur_integrated)

pdf(paste0(pref, "_MarkersDotPlot.pdf"), width = 20, height = 10)
DotPlot(allseur_integrated, group.by = "Ast_annot", features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

stackedbarplot(allseur_integrated[[]], groupx = 'Ast_annot', groupfill = 'orig.ident', paste0(pref, '_STACK_SAMPLES'), wd = 20)
stackedbarplot(allseur_integrated[[]], groupx = 'Ast_annot', groupfill = 'Species', paste0(pref, '_STACK_SPECIES'))

# save
saveRDS(allseur_integrated, paste0(pref,"_integrated_ANNOTATED.RDS"))

## add these annotations to the original seu obj
# extract the newest annotations
# Match cell names between objects
b_keep$is_ast = ifelse(rownames(b_keep@meta.data) %in% rownames(allseur_integrated@meta.data), "yes", "no")

# Identify cells in b that are not in the integrated object
b_keep_2 = subset(b_keep, subset = newannot_2 != "Astrocyte" | is_ast == "yes")

# check
Idents(b_keep_2) = b_keep_2$detailed_annot
DotPlot(b_keep_2, features = marksToPlot) + rotate_x_text(45)
# save
saveRDS(b_keep_2,"/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/seu_objs/merged_clean_AllCellTypes_AllTissues_seu.RDS")
####################### POST GLIA CLEANUP ###########################
file_outdir = "/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/seu_objs/"

species_list <- c("human", "chimp", "macaque", "marmoset", "mouse", "bat", "ferret")
primates_list <- c("human", "chimp", "macaque", "marmoset")

#### new obj AFTER GLIA CLEANUP
b = readRDS("/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/seu_objs/merged_clean_AllCellTypes_AllTissues_seu.RDS")
b$CellType = b$detailed_annot
b$CellType = gsub("COP", "OPC", b$CellType)

# put the interneuron subtype annotations to this seu obj
interneurons_seu = readRDS("/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/seu_objs/Interneurons_AllSpecies_AllTissues_ANNOTATED.RDS")

## add these annotations to the original seu obj
# extract the newest annotations
# Match cell names between objects
b$is_interneu = ifelse(rownames(b@meta.data) %in% rownames(interneurons_seu@meta.data), "yes", "no")

# Identify cells in b that are not in the integrated object
b_keep = subset(b, subset = detailed_annot != "Non_SPN" | is_interneu == "yes")

# order the rownames of smaller seu obj according to the bigger one
ordered_meta_interneurons_seu = interneurons_seu@meta.data[match(rownames(b_keep@meta.data), rownames(interneurons_seu@meta.data)),]

b_keep$is_interneu = ordered_meta_interneurons_seu$newannot
b_keep$a = ifelse(is.na(b_keep$is_interneu), as.character(b_keep$detailed_annot), as.character(b_keep$is_interneu))

# check
table(b_keep$detailed_annot, b_keep$a)
# re-assign
b_keep$broadannot_2 = b_keep$detailed_annot
b_keep$detailed_annot = b_keep$a

# orgnize metadata
b_keep$is_micro = NULL
b_keep$is_mol = NULL
b_keep$keep_cell = NULL
b_keep$is_ast = NULL
b_keep$is_interneu = NULL
b_keep$a = NULL
b_keep$broad_annot = NULL
b_keep$tissue = NULL
b_keep$id = paste0(b_keep$orig.ident, "_", b_keep$Tissue) 

### remove human caudate sample 
b_keep_2 = subset(b_keep, subset = id %in% "Sample_242999_Caudate", invert = T)

#### find the sex and age for all 
orig_meta <- b_keep_2@meta.data           # safe copy

### add info to the metadata
a = orig_meta
a[grep("v324", a$orig.ident),]$Sex = "Male"
a[grep("v321", a$orig.ident),]$Sex = "Female"
a$Sex = sub("female", "Female", a$Sex)
a$Sex = sub("male", "Male", a$Sex)
a[grep("Bat", a$Species),]$Sex = "Male"
a[grep("Bat", a$Species),]$Age = "3y"
a[grep("marm027", a$orig.ident),]$Sex = "Male"
a[grep("marm028", a$orig.ident),]$Sex = "Female"
a[grep("marm029", a$orig.ident),]$Sex = "Male"
a[grep("marm027", a$orig.ident),]$Age = "2y4m"
a[grep("marm028", a$orig.ident),]$Age = "3y2m"
a[grep("marm029", a$orig.ident),]$Age = "2y6m"
a[grep("v324", a$orig.ident),]$Age = "2035d"
a[grep("v321", a$orig.ident),]$Age = "2051d"
#a[grep("SRR13808459", a$orig.ident),]$Sex = "Female"
a[grep("SRR13808462", a$orig.ident),]$Sex = "Female"
a[grep("SRR13808463", a$orig.ident),]$Sex = "Female"
#a[grep("SRR13808466", a$orig.ident),]$Sex = "Female"
#a[grep("SRR13808467", a$orig.ident),]$Sex = "Female"
a[grep("SRR11921037", a$orig.ident),]$Sex = "Female"
a[grep("SRR11921038", a$orig.ident),]$Sex = "Female"
a[grep("SRR11921037", a$orig.ident),]$Age = "42d"
a[grep("SRR11921038", a$orig.ident),]$Age = "42d"
patterns <- c(paste0("SRR1192100", 5:9), paste0("SRR119210", 10:12))
pattern_all <- paste(patterns, collapse = "|")
idx <- grep(pattern_all, a$orig.ident)
a[idx, ]$Age <- "70d"
a[idx, ]$Sex = "Male"
a$sex = NULL
a$age = NULL

a$cellbarc = str_split_i(rownames(a), "-", 1)
a$cellID = paste0(a$cellbarc, "_", a$id)
a$Sex = sub("FeMale", "Female", a$Sex)

### add humanized ages
# read species traits taken from https://genomics.senescence.info/species/entry.php?species=Mus_musculus, https://genomics.senescence.info/species/entry.php?species=Mustela_nigripes, https://genomics.senescence.info/species/entry.php?species=Phyllostomus_hastatus, https://animaldiversity.org/accounts/Phyllostomus_hastatus/, https://genomics.senescence.info/species/entry.php?species=Phyllostomus_discolor
outdir = "/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/"
species_traits = read.csv(paste0(outdir, "primate_lifeTraits.csv"))
sample_metadata = read.csv(paste0(outdir, "comparative_striatum_sample_metadata.csv"), sep = "\t")
sample_metadata$id = as.character(paste0(sample_metadata$Barcode.ID, "_", sample_metadata$Brain.Area))

a$id = as.character(a$id)

# add the ages in years to the single cell-level metadata
meta_df <- merge(a, sample_metadata, by = "id", all.x = TRUE)
meta_df$Age <- meta_df$Age.x
meta_df$Sex <- meta_df$Sex.x
meta_df$Species <- meta_df$Species.x

meta_df <- meta_df[, !grepl("\\.x$|\\.y$", names(meta_df))]
meta_df$Brain.Area = NULL

# Add orig.ident only if rownames are not already "cellbarc-orig.ident"
rownames(meta_df) <- ifelse(
  grepl("[-_]", meta_df$cellbarc),  # check for presence of "-" or "_"
  meta_df$cellbarc,                  # if present, keep as is
  paste0(meta_df$cellbarc, "-", meta_df$orig.ident)  # else, append orig.ident
)

# restore original cell order (this prevents scrambling)
meta_df <- meta_df[rownames(a), ]

# Check if same set but different order
setequal(rownames(meta_df), rownames(a))

humanize_all_ages <- function(meta_df, species_traits) {
  # Initialize vector for humanized ages (same length and order)
  humanized_ages <- numeric(nrow(meta_df))
  
  # Get unique species in the metadata
  species_list <- unique(meta_df$Species)
  
  for (species in species_list) {
    # Convert ages to numerical
    species_ages <- meta_df$Age_in_years[meta_df$Species == species]
    species_ages_num <- as.numeric(species_ages)
    if (species %in% c("Human")) {
      # For humans, humanized age = original age
      humanized_ages[meta_df$Species == species] <- species_ages_num
    } else {
      # Check if species exists in lifetraits
      if (species %in% colnames(species_traits)) {
        # Build linear model: species ~ Human (lifetraits)
        model <- lm(species_traits[[species]] ~ species_traits$Human)
        
        # Calculate humanized ages using the model coefficients
        humanized <- (species_ages_num - model$coefficients[1]) / model$coefficients[2]
        
        # Assign back to correct rows in humanized_ages vector
        humanized_ages[meta_df$Species == species] <- humanized
      } else {
        warning(paste("Species", species, "not found in lifetraits. Assigning NA"))
        humanized_ages[meta_df$Species == species] <- NA
      }
    }
  }
  
  return(humanized_ages)
}

meta_df$Humanized_age = humanize_all_ages(meta_df, species_traits)

# update the metadata of the seu obj
# Join by cell ID (ensure cell names match exactly)
b_keep_2@meta.data = meta_df

# check
marksToPlot <- c('NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'APOE', 'GPC5', 'GLI3', 'AQP4', 'CSF1R', 'ARHGAP15', 'PLDL1', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CFTR', 'CHAT', 'ADARB2', 'LHX6', 'MEIS2', 'TAC3', 'VIP', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA', 'SATB2', 'SLC17A7', 'EYA2', 'PDGFD', 'PTHLH', 'PVALB', 'SPARCL1', 'RSPO2', 'RMST', 'LMO3', 'TSHZ2', "CALB1", "CALB2", "NOS1", "TSHZ1", "OPRM1", "CHST9","GRM8", "PPP1R1B", 'GRIK3', 'CXCL14')

# Plot previously identified markers
DefaultAssay(b_keep_2) = 'RNA'
b_keep_2 = NormalizeData(b_keep_2)

pdf(paste0("merged_clean_AllCellTypes_AllTissues_seu_detailedAnnot_DOTPLOT.pdf"), width = 20, height = 10)
DotPlot(b_keep_2, features = marksToPlot, group.by = "detailed_annot") +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

# save
saveRDS(b_keep_2,"/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/seu_objs/merged_clean_AllCellTypes_AllTissues_seu.RDS")

##### generate pseudobulk counts
file_outdir = "/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/seu_objs/"

species_list <- c("human", "chimp", "macaque", "marmoset", "mouse", "bat", "ferret")
primates_list <- c("human", "chimp", "macaque", "marmoset")

#### new obj AFTER GLIA CLEANUP
b = readRDS("/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/seu_objs/merged_clean_AllCellTypes_AllTissues_seu.RDS")

sub_obj_new_list = SplitObject(b, split.by = "Species")
names(sub_obj_new_list) = tolower(names(sub_obj_new_list))
meta_df = b@meta.data

### run function to generate pseudobulk count matrices for each species and celltype. Keep cellTypes which are present in at least 3 of the samples for a given Species-Tissue as well as genes with at least 1 copy of at least 3 of the samples for a given Species-Tissue. Generate the pseudobulk count matrix for all species-cellTypes except ferret which has only 2 samples using utilities function: pseudobulk_species = function(seurObj, ctype, features = c())

pb_count_list = list()
for (sp in species_list) {
  sub_obj_new = sub_obj_new_list[[sp]]
  celltypes = unique(sub_obj_new$CellType)

  for (ctype in celltypes) {
    if (sp != "ferret") {
      pb_count_list[[sp]][[ctype]] <- pseudobulk_species(
        seurObj = sub_obj_new,
        ctype = ctype,
        n_copy = 1,
        sample_number = 3
      )
    } else {
      pb_count_list[["ferret"]][[ctype]] <- NULL
    }
  }
}

#save
saveRDS(pb_count_list, file = "/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/seu_objs/pseudobulk_count_matrices_AllSpecies_AllCellTypes_broadannot_PostgliaCleanup.RDS")

### Genes were filtered to prep the pseudobulk counts (> 1 copy in >= 3 samples of that Species-Tissue) To continue with downstream analysis, get only the common genes across primates
common_genes_per_celltype = list()

primate_pb_count_list = pb_count_list[primates_list]
for (ctype in names(primate_pb_count_list$human)) {
  # Extract the gene vectors for this cell type across all species
  genes_across_species <- lapply(primate_pb_count_list, function(species_count) {
    rownames(species_count[[ctype]])
  })
  
  # Remove any NULLs (if some species lack this cell type)
  genes_across_species <- genes_across_species[!sapply(genes_across_species, is.null)]
  
  # Take intersection across species
  common_genes_per_celltype[[ctype]] <- Reduce(intersect, genes_across_species)
}

# Inspect results
str(common_genes_per_celltype)

# save
saveRDS(common_genes_per_celltype, file = "/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/seu_objs/pseudobulk_common_genes_withinPrimates_per_celltype_broadannot_postgliaCleanup.RDS")		
common_genes_per_celltype = readRDS("/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/seu_objs/pseudobulk_common_genes_withinPrimates_per_celltype_broadannot_postgliaCleanup.RDS")

##### subset each count matrix with these genes
sub_seu_pb_count_list = list()
pb_meta_list = list()
meta_df_list = list()
## Generate singleCellExperiment object
for (sp in primates_list){
  for (t in unique(sub_obj_new_list[[sp]]$Tissue)){
      # change COP to OPC
      sub_obj_new_list[[sp]]$CellType <- gsub("COP", "OPC", sub_obj_new_list[[sp]]$CellType)
      for (ctype in unique(sub_obj_new_list[[sp]]$CellType)){
      # get the seu obj
      seu = sub_obj_new_list[[sp]]
      
      # subset for the given cell type and tissue
      sub_seu = subset(seu, subset = CellType %in% ctype & Tissue %in% t)
      
      # subset pseudobulk count mat with common genes
      sub_seu_pb_count = primate_pb_count_list[[sp]][[ctype]][, grepl(t, colnames(primate_pb_count_list[[sp]][[ctype]]))]
      # subset the pseudocount matrix with common genes for each celltype
      sub_seu_pb_count_list[[t]][[ctype]][[sp]] = sub_seu_pb_count[common_genes_per_celltype[[ctype]],]
      
      # generate pseudobulk metadata
      sub_seu_meta <- sub_seu@meta.data %>%
            as.data.frame() %>%
            dplyr::select(Species, Sample, Tissue, Sex, Age, Age_in_years, Humanized_age, CellType, id, nCount_RNA)
            
      ## calculate the sum of the total UMI at the log scale for each sample
      tibble_meta <- as_tibble(sub_seu_meta) %>%
            dplyr::group_by(Sample) %>%
            dplyr::summarise(
                  log10_sum_ncount_RNA = log10(sum(nCount_RNA, na.rm = TRUE)),
                  n = dplyr::n(),
                  .groups = "drop"
            )
      tibble_meta = as.data.frame(tibble_meta)
      
      # aggregate based on sample
      a <- sub_seu_meta %>% distinct(Sample, .keep_all = TRUE)
      rownames(a) = a$Sample

      # merge two dataframes
      meta_df_sample = merge(a,tibble_meta ,by = "Sample", all.x = TRUE)
      rownames(meta_df_sample) = meta_df_sample$Sample
       
      # Check matching of matrix columns and metadata rows
      all(colnames(sub_seu_pb_count) == rownames(meta_df_sample))
      
      # if it says FALSE due to order not being the same. Order: 
      # order the metadata rownames based on the col names of the count matrix
      #rownames(meta_df_sample) <- rownames(meta_df_sample)[order(match(rownames(meta_df_sample), colnames(sub_seu_pb_count)))]

      # Check matching of matrix columns and metadata rows
      #all(colnames(sub_seu_pb_count) == rownames(meta_df_sample))
      # save in a list
      pb_meta_list[[t]][[ctype]][[sp]] = meta_df_sample
      meta_df_list[[t]][[ctype]][[sp]] = sub_seu@meta.data
      #rownames(meta_df_list[[t]][[ctype]][[sp]]) <- meta_df_list[[t]][[ctype]][[sp]]$Sample
      # get sub_seu_pb_count_list for the next step
}}}
saveRDS(pb_meta_list, file = paste0(file_outdir,"primates_pb_meta_list_PostgliaCleanup.RDS"))
saveRDS(meta_df_list, file = paste0(file_outdir,"primates_meta_df_list_PostgliaCleanup.RDS"))
saveRDS(sub_seu_pb_count_list, file = paste0(file_outdir,"primates_sub_seu_pb_count_list_PostgliaCleanup.RDS"))
  

### generate a list of metadata and pseudobulk counts for the same tissue and celltype 
merged_pb_meta <- list()
merged_sub_seu_pb_count <- list()
merged_meta_df <- list()

for (t in names(pb_meta_list)){
    for (ctype in names(pb_meta_list[[t]])){
      # Combine all species' metadata rows
      merged_pb_meta[[t]][[ctype]] <- do.call(rbind, pb_meta_list[[t]][[ctype]])
      merged_sub_seu_pb_count[[t]][[ctype]] <- do.call(cbind, sub_seu_pb_count_list[[t]][[ctype]])
      merged_meta_df[[t]][[ctype]] <- do.call(rbind, meta_df_list[[t]][[ctype]])
}}

# save
saveRDS(merged_pb_meta, file = paste0(file_outdir,"primates_merged_pb_meta_PostgliaCleanup.RDS"))
saveRDS(merged_sub_seu_pb_count, file = paste0(file_outdir,"primates_merged_sub_seu_pb_count_PostgliaCleanup.RDS"))
saveRDS(merged_meta_df, file = paste0(file_outdir,"primates_merged_meta_df_PostgliaCleanup.RDS"))

############################################################
sessionInfo()
R version 4.2.3 (2023-03-15)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux Server 7.9 (Maipo)

Matrix products: default
BLAS/LAPACK: /cm/shared/apps/rstudio-desktop/2022.12.0/lib/libopenblasp-r0.3.21.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] fastDummies_1.7.3           flashClust_1.01-2          
 [3] WGCNA_1.73                  fastcluster_1.3.0          
 [5] dynamicTreeCut_1.63-1       curl_5.2.2                 
 [7] Matrix.utils_0.9.8          variancePartition_1.28.9   
 [9] BiocParallel_1.32.5         edgeR_3.40.2               
[11] limma_3.54.0                ggrepel_0.9.3              
[13] DESeq2_1.38.3               harmony_1.1.0              
[15] Rcpp_1.0.10                 data.table_1.14.8          
[17] rio_1.0.1                   reshape2_1.4.4             
[19] ggpubr_0.6.0                lubridate_1.9.3            
[21] forcats_1.0.0               stringr_1.5.0              
[23] purrr_1.0.1                 readr_2.1.4                
[25] tidyr_1.3.0                 tibble_3.2.1               
[27] tidyverse_2.0.0             plyr_1.8.9                 
[29] ggplot2_3.4.4               Matrix_1.6-5               
[31] dplyr_1.1.4                 DropletQC_0.0.0.9000       
[33] DropletUtils_1.18.1         SingleCellExperiment_1.20.0
[35] SummarizedExperiment_1.28.0 Biobase_2.58.0             
[37] GenomicRanges_1.50.0        GenomeInfoDb_1.34.9        
[39] IRanges_2.32.0              S4Vectors_0.36.0           
[41] BiocGenerics_0.44.0         MatrixGenerics_1.10.0      
[43] matrixStats_0.63.0          rhdf5_2.42.1               
[45] BPCells_0.3.1               Seurat_5.3.0               
[47] SeuratObject_5.1.0          sp_1.6-0                   
[49] patchwork_1.3.0.9000       

loaded via a namespace (and not attached):
  [1] scattermore_1.2           R.methodsS3_1.8.2        
  [3] knitr_1.48                bit64_4.0.5              
  [5] irlba_2.3.5.1             DelayedArray_0.24.0      
  [7] R.utils_2.12.2            rpart_4.1.19             
  [9] KEGGREST_1.38.0           RCurl_1.98-1.12          
 [11] doParallel_1.0.17         generics_0.1.3           
 [13] preprocessCore_1.60.2     RhpcBLASctl_0.23-42      
 [15] cowplot_1.1.1             RSQLite_2.3.1            
 [17] RANN_2.6.1                future_1.33.1            
 [19] bit_4.0.5                 tzdb_0.3.0               
 [21] spatstat.data_3.0-0       httpuv_1.6.9             
 [23] xfun_0.47                 hms_1.1.2                
 [25] evaluate_0.20             promises_1.2.0.1         
 [27] fansi_1.0.6               progress_1.2.2           
 [29] caTools_1.18.2            igraph_1.4.2             
 [31] DBI_1.1.3                 geneplotter_1.76.0       
 [33] htmlwidgets_1.6.2         spatstat.geom_3.0-6      
 [35] ellipsis_0.3.2            RSpectra_0.16-1          
 [37] backports_1.4.1           annotate_1.76.0          
 [39] aod_1.3.3                 deldir_1.0-6             
 [41] sparseMatrixStats_1.10.0  vctrs_0.6.5              
 [43] ROCR_1.0-11               abind_1.4-5              
 [45] cachem_1.0.7              withr_2.5.2              
 [47] grr_0.9.5                 progressr_0.13.0         
 [49] checkmate_2.3.0           sctransform_0.4.1        
 [51] prettyunits_1.1.1         mclust_6.0.0             
 [53] goftest_1.2-3             cluster_2.1.4            
 [55] dotCall64_1.1-0           lazyeval_0.2.2           
 [57] crayon_1.5.2              spatstat.explore_3.0-6   
 [59] pkgconfig_2.0.3           nlme_3.1-162             
 [61] nnet_7.3-18               rlang_1.1.4              
 [63] globals_0.16.2            lifecycle_1.0.3          
 [65] miniUI_0.1.1.1            polyclip_1.10-4          
 [67] RcppHNSW_0.4.1            lmtest_0.9-40            
 [69] carData_3.0-5             Rhdf5lib_1.20.0          
 [71] boot_1.3-28.1             zoo_1.8-11               
 [73] base64enc_0.1-3           ggridges_0.5.4           
 [75] png_0.1-8                 viridisLite_0.4.1        
 [77] bitops_1.0-7              R.oo_1.25.0              
 [79] KernSmooth_2.23-20        spam_2.10-0              
 [81] rhdf5filters_1.10.1       Biostrings_2.66.0        
 [83] blob_1.2.3                DelayedMatrixStats_1.20.0
 [85] parallelly_1.35.0         spatstat.random_3.1-3    
 [87] remaCor_0.0.16            rstatix_0.7.2            
 [89] ggsignif_0.6.4            beachmat_2.14.0          
 [91] scales_1.3.0              memoise_2.0.1            
 [93] magrittr_2.0.3            ica_1.0-3                
 [95] gplots_3.1.3              zlibbioc_1.44.0          
 [97] compiler_4.2.3            dqrng_0.3.0              
 [99] RColorBrewer_1.1-3        lme4_1.1-35.1            
[101] fitdistrplus_1.1-8        cli_3.6.2                
[103] XVector_0.38.0            listenv_0.9.0            
[105] pbapply_1.7-0             htmlTable_2.4.2          
[107] Formula_1.2-5             MASS_7.3-58.3            
[109] tidyselect_1.2.0          stringi_1.7.12           
[111] locfit_1.5-9.8            grid_4.2.3               
[113] tools_4.2.3               timechange_0.2.0         
[115] future.apply_1.10.0       parallel_4.2.3           
[117] rstudioapi_0.14           foreign_0.8-85           
[119] foreach_1.5.2             gridExtra_2.3            
[121] EnvStats_2.8.1            farver_2.1.1             
[123] Rtsne_0.16                digest_0.6.31            
[125] shiny_1.7.4               car_3.1-2                
[127] broom_1.0.3               scuttle_1.8.0            
[129] later_1.3.0               RcppAnnoy_0.0.20         
[131] httr_1.4.5                AnnotationDbi_1.60.2     
[133] Rdpack_2.6                colorspace_2.1-0         
[135] XML_3.99-0.14             tensor_1.5               
[137] reticulate_1.42.0         splines_4.2.3            
[139] uwot_0.1.14               spatstat.utils_3.1-4     
[141] plotly_4.10.1             xtable_1.8-4             
[143] jsonlite_1.8.4            nloptr_2.0.3             
[145] R6_2.5.1                  Hmisc_5.1-1              
[147] pillar_1.9.0              htmltools_0.5.8.1        
[149] mime_0.12                 glue_1.6.2               
[151] fastmap_1.1.1             minqa_1.2.6              
[153] codetools_0.2-19          mvtnorm_1.2-3            
[155] utf8_1.2.4                lattice_0.21-8           
[157] spatstat.sparse_3.0-0     pbkrtest_0.5.2           
[159] gtools_3.9.4              GO.db_3.16.0             
[161] survival_3.5-3            rmarkdown_2.28           
[163] munsell_0.5.0             GenomeInfoDbData_1.2.9   
[165] iterators_1.0.14          impute_1.72.3            
[167] HDF5Array_1.26.0          gtable_0.3.3             
[169] rbibutils_2.2.16        
