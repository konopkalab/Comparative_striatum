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
library(DropletQC)
library(xlsx)
library(biomaRt)
library(ggbreak)
source('~/SCRIPTS/utility_functions.R')
set.seed(1234)


####
## LOAD NEURON DATASETS
####

### integrate SPNs of: Human, chimp, macaque, marmoset, bat, mouse, and ferret
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
mouse_caud_put <- readRDS("/home2/gkonop/workdir/01_MOUSE/03_CLUSTER_ANNOTATE/Mouse_Caudate_Annotated_FINAL.RDS")
mouse_caud_put$Species <- "Mouse"
table(mouse_caud_put$newannot)
mouse_caud_put$Tissue <- "Caudoputamen"

# bat
bat <- readRDS("/home2/gkonop/workdir/01_BAT/03_CLUSTER_ANNOTATE/Caudate_Putamen/ANNOTATED/bat_integrated_cuadate_putamen_ANNOTATED.RDS")

# ferret
ferret <- readRDS("/home2/gkonop/workdir/01_FERRET/03_CLUSTER_ANNOTATE/ANNOTATED/Ferret_Caudate_Krienen_ANNOTATED.RDS")
# create "id" metadata
ferret$id <- ferret$orig.ident
ferret$newannot <- ferret$broad_annot
########################################################
# put all the genes in each of the count matrices in a vector
all_genes <- c(rownames(human@assays$RNA@counts),
		rownames(chimp@assays$RNA@counts),
		rownames(macaque@assays$RNA@counts),
		rownames(marmoset@assays$RNA@counts),
		toupper(rownames(mouse_caud_put@assays$RNA@counts)),
		rownames(bat@assays$RNA@counts),
		rownames(ferret@assays$RNA@counts))

sorted_all_genes <- sort(all_genes)
# get the unique ones
all_unique_genes <- unique(sorted_all_genes)

## get the protein coding human genes only
# get all human genes
gene_names <- rownames(human@assays$RNA@counts)

library(biomaRt)

# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get the list of protein-coding genes
protein_coding_genes <- getBM(attributes = c("external_gene_name", "gene_biotype"),
                              filters = "biotype",
                              values = "protein_coding",
                              mart = ensembl)

# Extract the gene names of protein-coding genes
protein_coding_gene_names <- protein_coding_genes$external_gene_name

# write all the genes in a txt file
write.table(protein_coding_gene_names, "/endosome/work/Neuroinformatics_Core/gkonop/03_INTEGRATE_ALL/SPN/human_pr_coding_genes.txt", quote = F, row.names = F, col.names = F)
write.table(all_unique_genes, "/endosome/work/Neuroinformatics_Core/gkonop/03_INTEGRATE_ALL/SPN/all_genes.txt", quote = F, row.names = F, col.names = F)
write.table(rownames(human@assays$RNA@counts), "/endosome/work/Neuroinformatics_Core/gkonop/03_INTEGRATE_ALL/SPN/human_genes.txt", quote = F, row.names = F, col.names = F)

### find the orthologue genes between any other species and human
## chimp, macaque, and marmoset genomes were "lifted off" to the human genome
################################################
# Found ortholog genes from human pr coding genes and extracted the genes which are ortholog in all 7 species using ncbi datasets tool
ortho_genes <- read.table("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/orthologs_in_7_species_human_pr_codingOnly_2.csv",  header = TRUE)$symbol #15055 pr coding ortho genes

## some markers are not present in the common orthologs so add these genes to the list
# markers
marksToPlot = c('SLC17A7', 'SATB2', 'RBFOX3', 'NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA')

# Genes in marksToPlot that are not in ortho_genes$symbol
missing_genes <- marksToPlot[!marksToPlot %in% ortho_genes]
print(missing_genes) #APBB1IP gene is missing macaque ortho info

### Create new seurat objects for each species with the common orthologue genes
## ferret
# extract the count matrix
ferret_old_count_mat <- ferret@assays$RNA@counts
# Subset the expression matrix to include only common orthologous genes
ferret_humanized_common_mat <- ferret_old_count_mat[rownames(ferret_old_count_mat) %in% c(ortho_genes,missing_genes), ]
# create new seurat object
ferret_humanized = CreateSeuratObject(ferret_humanized_common_mat, meta.data = ferret@meta.data)
# subset SPNs
ferret_humanized_SPN <- subset(ferret_humanized, subset = newannot == "SPN")

## marmoset
# extract the count matrix
marmoset_old_count_mat <- marmoset@assays$RNA@counts
# Subset the expression matrix to include only common orthologous genes
marmoset_humanized_common_mat <- marmoset_old_count_mat[rownames(marmoset_old_count_mat) %in% c(ortho_genes,missing_genes),]
# create new seurat object
marmoset_humanized = CreateSeuratObject(marmoset_humanized_common_mat, meta.data = marmoset@meta.data)
# subset SPNs and caudate
marmoset_humanized_SPN <- subset(marmoset_humanized, subset = newannot == "SPN")

## macaque
# extract the count matrix
macaque_old_count_mat <- macaque@assays$RNA@counts
# Subset the expression matrix to include only common orthologous genes
macaque_humanized_common_mat <- macaque_old_count_mat[rownames(macaque_old_count_mat) %in% c(ortho_genes,missing_genes),]
# create new seurat object
macaque_humanized = CreateSeuratObject(macaque_humanized_common_mat, meta.data = macaque@meta.data)
# subset SPNs
macaque_humanized_SPN <- subset(macaque_humanized, subset = newannot == "SPN")

## mouse
## change the gene names to human annotaiton
# extract the count matrix
mouse_old_count_mat <- mouse_caud_put@assays$RNA@counts
# change
rownames(mouse_old_count_mat) <- toupper(rownames(mouse_old_count_mat))
# Subset the expression matrix to include only common orthologous genes
mouse_humanized_common_mat <- mouse_old_count_mat[rownames(mouse_old_count_mat) %in% c(ortho_genes,missing_genes),]
# create new seurat object
mouse_humanized = CreateSeuratObject(mouse_humanized_common_mat, meta.data = mouse_caud_put@meta.data)
# subset the spns
mouse_humanized_SPN <- subset(mouse_humanized, subset = newannot =="SPN")

## chimp
# extract the count matrix
chimp_old_count_mat <- chimp@assays$RNA@counts
# Subset the expression matrix to include only common orthologous genes
chimp_humanized_common_mat <- chimp_old_count_mat[rownames(chimp_old_count_mat) %in% c(ortho_genes,missing_genes),  ]
# create new seurat object
chimp_humanized = CreateSeuratObject(chimp_humanized_common_mat, meta.data = chimp@meta.data)
# subset the spns
chimp_humanized_SPN <- subset(chimp_humanized, subset = newannot =="SPN")

## bat
bat_count <- bat@assays$RNA@counts
# Subset the expression matrix to include only common orthologous genes
bat_humanized_common_mat <- bat_count[rownames(bat_count) %in% c(ortho_genes,missing_genes),  ]
# create new seurat object
bat_humanized = CreateSeuratObject(bat_humanized_common_mat, meta.data = bat@meta.data)
# subset the spns
bat_humanized_SPN <- subset(bat_humanized, subset = newannot =="SPN")

# human
human_count <- human@assays$RNA@counts
# Subset the expression matrix to include only common orthologous genes
human_common_mat <- human_count[rownames(human_count) %in% c(ortho_genes,missing_genes),]
new_human <- CreateSeuratObject(human_common_mat, meta.data = human@meta.data)
# subset the spns
human_SPN <- subset(new_human, subset = newannot =="SPN")

# merge SPN seurat objects
seurM = Reduce(merge, list(human_SPN, chimp_humanized_SPN, macaque_humanized_SPN, marmoset_humanized_SPN, bat_humanized_SPN, mouse_humanized_SPN, ferret_humanized_SPN))
dim(seurM) #[1]  14894 293822

# save
saveRDS(seurM, file = "/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/merged_seu_obj_SPN_primates_bat_mouse_ferret_common_ncbi_database_human_pr_coding_orthologs.RDS")
seurM <- readRDS("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/merged_seu_obj_SPN_primates_bat_mouse_ferret_common_ncbi_database_human_pr_coding_orthologs.RDS")
######################################################
# Define the sample to remove
# Set assay to RNA
DefaultAssay(seurM) <- "RNA"

#seurM_no_mouse <- subset(seurM, subset = Species == "Mouse", invert = T) 

# check
table(seurM$Species)
#     Bat    Chimp   Ferret    Human  Macaque Marmoset    Mouse 
#   63799    28334    11723    35690    55198    43087    55991 

# Split data by sample and tissue
seurM$id <- paste0(seurM$orig.ident, '_', seurM$Tissue)
seurM_split <- SplitObject(seurM, split.by = "id")

# Normalize and identify variable features for each dataset independently
seurM_split <- lapply(X = seurM_split, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    return(x) # Ensure each processed object is returned
})

# Check the result
lapply(seurM_split, function(x) dim(x@assays$RNA@data)) # Check dimensions after processing

# find the samples with < 500 SPNs
lapply(seurM_split, function(x) {
	spn_count <- dim(x@assays$RNA@data)[2]
	if(spn_count < 600)
		print(paste0(x$Species,", ", x$Tissue,", ", spn_count," is < 600"))
	else
		print("All good")
})
#$SRR13808461_Putamen [1] "Macaque, Putamen, 146 is < 600"
#$SRR13808459_Caudate [1] "Macaque, Caudate, 561 is < 600"
#Sample_242939_Putamen [1] "Human, Putamen, 437 is < 600"


# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features = SelectIntegrationFeatures(object.list = seurM_split)
print(paste0('Number of genes to use for integration: ', length(features)))
seurML = lapply(X = seurM_split, FUN = function(x) {
    x = ScaleData(x, features = features, verbose = FALSE)
    x = RunPCA(x, features = features, verbose = FALSE)
})

## cca gives a mmerory problem try integrating with rpca
anchors = FindIntegrationAnchors(object.list = seurML, anchor.features = features, reduction = "rpca")

saveRDS(anchors, "/project/Neuroinformatics_Core/Konopka_lab/shared/For_Gozde/03_INTEGRATED_ALL/SPN/GB_SPN_anchors_rpca_primates_bat_mouse_ferret_common_ncbi_database_human_pr_coding_orthologs.RDS")

# integrate
# set this option when analyzing large datasets
options(future.globals.maxSize = 3e+09)
# List all objects in the global environment
all_objects <- ls()

# Specify the objects you want to keep
keep_objects <- c("anchors")

# Delete all objects except the specified one
rm(list = setdiff(all_objects, keep_objects))
gc() # make space

allseur_integrated = IntegrateData(anchorset = anchors, k.weight  = 90) #!!!! Put the kweight less than your the least cell containg sample!!!!
# in projectdir
saveRDS(allseur_integrated, "/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/GB_Integrated_rpca_SPN_primates_bat_mouse_ferret_common_ncbi_database_human_pr_coding_orthologs.RDS")
allseur_integrated <- readRDS("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/GB_Integrated_rpca_SPN_primates_bat_mouse_ferret_common_ncbi_database_human_pr_coding_orthologs.RDS")


# Clustering
marksToPlot = c('SLC17A7', 'SATB2', 'RBFOX3', 'NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA')

DefaultAssay(allseur_integrated) = 'integrated'
allseur_integrated = ScaleData(allseur_integrated, verbose = FALSE)
allseur_integrated = RunPCA(allseur_integrated, verbose = FALSE)
allseur_integrated = RunUMAP(allseur_integrated, dims = 1:20, reduction = 'pca')

# Basic Plots
pref = 'GB_AllTissues_Primates_bat_mouse_ferret_Integrated_SPN_ncbi_human_pr_coding_orthologs'
setwd("/home2/gkonop/workdir/03_INTEGRATE_ALL/SPN/CLUSTERING_1")

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

pdf(paste0(pref, "_INTEGRATED_DEPTH.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'black') + NoLegend() +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_INTEGRATED_FEATUREPLOT.pdf"), width = 15, height = 35)
FeaturePlot(allseur_integrated, features = marksToPlot, sort = T, raster = T, pt.size = 2)
dev.off()

pdf(paste0(pref, "_INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()
# save
saveRDS(allseur_integrated,"/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/CLUSTERING_1/GB_Caudate_Primates_bat_ferret_Integrated_SPN_ncbi_human_pr_coding_orthologs_CLUSTERING_1.RDS")
allseur_integrated <- readRDS("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/CLUSTERING_1/GB_Caudate_Primates_bat_ferret_Integrated_SPN_ncbi_human_pr_coding_orthologs_CLUSTERING_1.RDS")
##########################################
###
## FILTER AND REPLOT
####
#Cluster#,Cell_type
#0,iSPN
#1,dSPN
#2,dSPN
#3,dSPN
#4,dSPN
#5,dSPN
#6,iSPN
#7,iSPN
#8,iSPN
#9,iSPN
#10,dSPN
#11,iSPN
#12*,dSPN+iSPN
#13*,MOL+SPN
#14,eSPN
#15*,SPN+MOL+OPC+Ast
#16,iSPN
#17,iSPN
#18,?
#19,iSPN
#20*,SPN+MOL+Ast
#21,iSPN
#22*,SPN+MOL
#23*,SPN+MOL+microglia
#24,eSPN
#25*,doublet
#26*,doublet

setwd("/home2/gkonop/workdir/03_INTEGRATE_ALL/SPN/CLUSTERING_2")
toremove = c(12,13,15,20,22,23,25,26)
allseur_integrated = subset(allseur_integrated, subset = seurat_clusters %in% toremove, invert = T)

# Clustering
DefaultAssay(allseur_integrated) = 'integrated'
allseur_integrated = ScaleData(allseur_integrated, verbose = FALSE)
allseur_integrated = RunPCA(allseur_integrated, verbose = FALSE)
allseur_integrated = RunUMAP(allseur_integrated, dims = 1:20, reduction = 'pca')

# Basic Plots
setwd("/home2/gkonop/workdir/03_INTEGRATE_ALL/SPN/CLUSTERING_2")

pdf(paste0(pref, "_UMAP.pdf"))
DimPlot(allseur_integrated, group.by = 'Species', raster = T)
dev.off()


allseur_integrated = FindNeighbors(allseur_integrated, dims = 1:20, reduction = 'pca')
allseur_integrated = FindClusters(allseur_integrated, resolution = 1.5)

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

pdf(paste0(pref, "INTEGRATED_FEATUREPLOT.pdf"), width = 15, height = 35)
FeaturePlot(allseur_integrated, features = marksToPlot, sort = T, raster = T, pt.size = 2)
dev.off()

pdf(paste0(pref, "INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# save
saveRDS(allseur_integrated, paste0("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/CLUSTERING_2/",pref, "_CLUSTERING_2.RDS"))
allseur_integrated <- readRDS(paste0("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/CLUSTERING_2/GB_AllTissues_Primates_bat_ferret_Integrated_SPN_ncbi_human_pr_coding_orthologs_CLUSTERING_2.RDS"))
#######################################################################################################
###
## FILTER AND REPLOT
####
#Cluster#,Cell_type
#0,iSPN
#1,dSPN
#2,iSPN
#3,dSPN
#4,dSPN
#5,dSPN
#6,iSPN
#7,iSPN
#8,dSPN
#9,iSPN
#10,iSPN
#11,dSPN
#12,dSPN
#13,eSPN
#14,iSPN
#15,dSPN
#16,dSPN
#17,iSPN
#18*,Neu+OPC
#19,eSPN
#20,iSPN
#21,eSPN
#22*,MOL+eSPN

setwd("/home2/gkonop/workdir/03_INTEGRATE_ALL/SPN/CLUSTERING_3")
toremove = c(18,22)
allseur_integrated = subset(allseur_integrated, subset = seurat_clusters %in% toremove, invert = T)

# Clustering
DefaultAssay(allseur_integrated) = 'integrated'
allseur_integrated = ScaleData(allseur_integrated, verbose = FALSE)
allseur_integrated = RunPCA(allseur_integrated, verbose = FALSE)
allseur_integrated = RunUMAP(allseur_integrated, dims = 1:10, reduction = 'pca')

# Basic Plots
pdf(paste0(pref, "_UMAP.pdf"))
DimPlot(allseur_integrated, group.by = 'Species', raster = T)
dev.off()


allseur_integrated = FindNeighbors(allseur_integrated, dims = 1:10, reduction = 'pca')
allseur_integrated = FindClusters(allseur_integrated, resolution = 0.5)

pdf(paste0(pref, "_CLUSTERS.pdf"))
DimPlot(allseur_integrated, label = T, raster = T) + NoLegend()
dev.off()

stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_STACK_SAMPLES'), wd = 20)
stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'Species', paste0(pref, '_STACK_SPECIES'))

# Plot previously identified markers
DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)

## label for striosome and matrix
# markers
marksToPlot = c('SLC17A7', 'SATB2', 'RBFOX3', 'NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1','COL11A1', 'PAPC', 'DNER', 'DRD2', 'SPON1', 'BACH2','KCNIP1', 'EPHAA5', 'HTR2A', 'RXRG', 'SEMA5B', 'MFGE8', 'KREMEN1','PDE1C','OPRM1', 'CALB1', 'CRYM', 'ID4', 'ZFHX3', 'SGK1','SV2B','FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA')


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

pdf(paste0(pref, "INTEGRATED_FEATUREPLOT.pdf"), width = 15, height = 35)
FeaturePlot(allseur_integrated, features = marksToPlot, sort = T, raster = T, pt.size = 2)
dev.off()

pdf(paste0(pref, "INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# save
saveRDS(allseur_integrated, paste0("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/CLUSTERING_3/",pref, "_CLUSTERING_3.RDS"))
allseur_integrated <- readRDS(paste0("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/CLUSTERING_3/",pref, "_CLUSTERING_3.RDS"))
#######################################################################################################
###
## CLUSTERING_4 AND REPLOT
####
#Cluster#,Cell_type (res 1.5)
#0,dSPN_matrix
#1,dSPN_matrix
#2,iSPN_matrix
#3,iSPN_matrix
#4,dSPN_matrix
#5,dSPN_matrix
#6,iSPN_matrix
#7,iSPN_matrix
#8,dSPN_matrix
#9,iSPN_matrix
#10,iSPN_matrix
#11,dSPN_striosome
#12,dSPN_striosome
#13,dSPN_striosome
#14,iSPN_striosome
#15,dSPN_matrix
#16,iSPN_matrix
#17,iSPN_matrix
#18,dSPN_matrix
#19,iSPN_striosome
#20,iSPN_striosome
#21,eSPN
#22,iSPN_matrix
#23,iSPN_matrix
#24,dSPN_matrix
#25,eSPN
#26,iSPN_striosome
#27,dSPN_matrix
#28,eSPN
#29,iSPN_striosome
#30,iSPN_matrix
#31,eSPN
#32*,SPN+MOL
setwd("/home2/gkonop/workdir/03_INTEGRATE_ALL/SPN/CLUSTERING_4")
toremove = c(32)
allseur_integrated = subset(allseur_integrated, subset = seurat_clusters %in% toremove, invert = T)

#### (res 1)
#Cluster#,Cell_type (res 1)
#0,iSPN_matrix
#1,iSPN_matrix
#2,dSPN_matrix
#3,dSPN_matrix
#4,dSPN_matrix
#5,iSPN_matrix
#6,iSPN_striosome
#7,dSPN_matrix
#8,iSPN_matrix
#9,dSPN_matrix
#10,iSPN_matrix
#11,dSPN_striosome
#12,dSPN_matrix
#13,dSPN_matrix
#14,eSPN
#15,dSPN_matrix
#16,iSPN_matrix
#17,iSPN_matrix
#18,iSPN_striosome
#19,eSPN
#20,iSPN_matrix
#21,iSPN_striosome
#22,eSPN
#23,iSPN_matrix

##### CREATE res 0.2 SPN seu obj

#### (res 0.2)
#Cluster#,Cell_type
#0,iSPN_matrix
#1,dSPN_matrix
#2,dSPN_matrix
#3,dSPN_striosome
#4,iSPN_striosome
#5,eSPN
#6,iSPN_matrix
#7,iSPN_striosome

(mapnames <- setNames(c("iSPN_matrix","dSPN_matrix","dSPN_matrix",
"dSPN_striosome","iSPN_striosome","eSPN","iSPN_matrix","iSPN_striosome"), 0:7))

allseur_integrated[["newannot"]] = mapnames[allseur_integrated[["seurat_clusters"]][,1]]
Idents(allseur_integrated) = allseur_integrated$newannot
allseur_integrated[["newannot_2"]] = str_split_i(allseur_integrated$newannot, "_", 1)

table(allseur_integrated$newannot_2)

  dSPN   eSPN   iSPN 
134486  10033 111186
 
dim(allseur_integrated) #[1]  14894 255705
saveRDS(allseur_integrated, paste0("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/ANNOTATED/GB_AllTissues_Primates_bat_mouse_ferret_Integrated_SPN_ncbi_human_pr_coding_orthologs_res1_res0_2_ANNOTATED.RDS"))

#############################################################
###
# Clustering
DefaultAssay(allseur_integrated) = 'integrated'
allseur_integrated = ScaleData(allseur_integrated, verbose = FALSE)
allseur_integrated = RunPCA(allseur_integrated, verbose = FALSE)
allseur_integrated = RunUMAP(allseur_integrated, dims = 1:10, reduction = 'pca')

# Basic Plots
pdf(paste0(pref, "_UMAP.pdf"))
DimPlot(allseur_integrated, group.by = 'Species', raster = T)
dev.off()


allseur_integrated = FindNeighbors(allseur_integrated, dims = 1:10, reduction = 'pca')
allseur_integrated = FindClusters(allseur_integrated, resolution = 0.2)

pdf(paste0(pref, "_CLUSTERS.pdf"))
DimPlot(allseur_integrated, label = T, raster = T) + NoLegend()
dev.off()

stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', paste0(pref, '_STACK_SAMPLES'), wd = 20)
stackedbarplot(allseur_integrated[[]], groupx = 'seurat_clusters', groupfill = 'Species', paste0(pref, '_STACK_SPECIES'))

# Plot previously identified markers
DefaultAssay(allseur_integrated) = 'RNA'
allseur_integrated = NormalizeData(allseur_integrated)

## label for striosome and matrix
# markers
marksToPlot = c('SLC17A7', 'SATB2', 'RBFOX3', 'NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1','COL11A1', 'PAPC', 'DNER', 'DRD2', 'SPON1', 'BACH2','KCNIP1', 'EPHAA5', 'HTR2A', 'RXRG', 'SEMA5B', 'MFGE8', 'KREMEN1','PDE1C','OPRM1', 'CALB1', 'CRYM', 'ID4', 'ZFHX3', 'SGK1','SV2B','FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA')


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

pdf(paste0(pref, "INTEGRATED_FEATUREPLOT.pdf"), width = 15, height = 35)
FeaturePlot(allseur_integrated, features = marksToPlot, sort = T, raster = T, pt.size = 2)
dev.off()

pdf(paste0(pref, "INTEGRATED_NUCLEAR_FRACTION.pdf"))
ggboxplot(allseur_integrated[[]], x = 'seurat_clusters', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# save
saveRDS(allseur_integrated, paste0("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/CLUSTERING_4/",pref, "CLUSTERING_4.RDS"))
allseur_integrated <- readRDS(paste0("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/CLUSTERING_4/",pref, "CLUSTERING_4.RDS"))

###########################################################################################
###
## ANNOTATE AND REPLOT
####
#Cluster#,Cell_type
#0,iSPN_matrix
#1,dSPN_matrix
#2,dSPN_matrix
#3,dSPN_striosome
#4,iSPN_matrix
#5,iSPN_striosome
#6,dSPN_matrix
#7,iSPN_matrix
#8,eSPN
#9,dSPN_matrix
#10,iSPN_matrix
#11,iSPN_striosome
#12,iSPN_matrix


setwd("/home2/gkonop/workdir/03_INTEGRATE_ALL/SPN/ANNOTATED")
(mapnames <- setNames(c("iSPN_matrix","dSPN_matrix","dSPN_matrix","dSPN_striosome","iSPN_matrix",
"iSPN_striosome","dSPN_matrix","iSPN_matrix","eSPN","dSPN_matrix","iSPN_matrix",
"iSPN_striosome","iSPN_matrix"), 0:12))

allseur_integrated[["newannot"]] = mapnames[allseur_integrated[["seurat_clusters"]][,1]]
Idents(allseur_integrated) = allseur_integrated$newannot

# PLOTS for striosome and matrix
pdf(paste0(pref, "_Clusters_Annotated.pdf"))
DimPlot(allseur_integrated, group.by = 'newannot', cols = c('iSPN_matrix' = 'springgreen3', 'iSPN_striosome' = 'violetred', 'dSPN_matrix' = 'deepskyblue2', 'dSPN_striosome' = 'darkorchid', 'eSPN' = 'chocolate2'), label = T, raster = T, repel = T, label.size = 7) +
NoLegend()
ggtitle('')
dev.off()

pdf(paste0(pref, "_Depth_Gene_Annotated.pdf"))
ggboxplot(allseur_integrated[[]], x = 'newannot', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_Depth_UMI_Annotated.pdf"))
ggboxplot(allseur_integrated[[]], x = 'newannot', y = 'nCount_RNA', ylim = c(0,80000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_NuclearFraction_Annotated.pdf"))
ggboxplot(allseur_integrated[[]], x = 'newannot', y = 'intronRat') +
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

allseur_integrated$newannot = factor(allseur_integrated$newannot, levels = c("dSPN_matrix", "dSPN_striosome","iSPN_matrix", "iSPN_striosome","eSPN"))

allseur_integrated$Species = factor(allseur_integrated$Species, levels = c('Ferret', 'Bat', 'Mouse', 'Marmoset', 'Macaque', 'Chimp', 'Human'))

stackedbarplot(allseur_integrated[[]], groupfill = 'newannot', groupx = 'Species', fn = paste0(pref, '_Stacked_Annotated_Species'))

stackedbarplot(allseur_integrated[[]], groupfill = 'Species', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Species2'))

stackedbarplot(allseur_integrated[[]], groupfill = 'newannot', groupx = 'orig.ident', fn = paste0(pref, '_Stacked_Annotated_Samples'), hg = 20)

stackedbarplot(allseur_integrated[[]], groupfill = 'orig.ident', groupx = 'newannot', fn = paste0(pref, '_Stacked_AnnotatedSamples2'), wd = 20)

stackedbarplot(allseur_integrated[[]], groupfill = 'Tissue', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Tissue'))

# Primate vs bat vs ferret
meta = allseur_integrated@meta.data
meta$Species_Broad = 'Primate'
meta[meta$Species == 'Bat', 'Species_Broad'] = 'Bat'
meta[meta$Species == 'Ferret', 'Species_Broad'] = 'Ferret'

stackedbarplot(meta, groupfill = 'Species_Broad', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Species_Primate_vsRest'))

# stacked bar plot of cell types for only primates
sub = subset(meta, subset = Species_Broad == c('Primate'))
table(sub$Species)
stackedbarplot(sub, groupfill = 'Species', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Species_Primates'))

######################################################################################
# PLOTS for dSPN, iSPN, and eSPNs
allseur_integrated[["newannot_2"]]  = str_split_i(allseur_integrated$newannot, "_",1)
Idents(allseur_integrated) = allseur_integrated$newannot_2

## remove any sample which has < 500 SPNs
# extract sample names
samples <- table(allseur_integrated$orig.ident)
bad_samples <- names(samples[samples < 500])
filtered_spn <- subset(allseur_integrated, orig.ident %in% bad_samples, invert = T)
filtered_spn$newannot_2 <- factor(filtered_spn$newannot_2, levels = c("dSPN", "iSPN", "eSPN"))
filtered_spn$newannot <- factor(filtered_spn$newannot, levels = c("dSPN_striosome", "dSPN_matrix", "iSPN_striosome", "iSPN_matrix","eSPN"))

pref = 'GB_AllTissues_Primates_bat_mouse_ferret_filtered_Integrated_SPN_ncbi_human_pr_coding_orthologs'
setwd("/endosome/work/Neuroinformatics_Core/gkonop/03_INTEGRATE_ALL/SPN/ANNOTATED")

# PLOTS for striosome and matrix
pdf(paste0(pref, "_Clusters_Annotated_BroadAnnot.pdf"))
DimPlot(filtered_spn, group.by = 'newannot_2', label = T, raster = T, repel = T, label.size = 7) +
NoLegend() +
ggtitle('')
dev.off()

pdf(paste0(pref, "_Depth_Gene_Annotated_BroadAnnot.pdf"))
ggboxplot(filtered_spn[[]], x = 'newannot_2', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_Depth_UMI_Annotated_BroadAnnot.pdf"))
ggboxplot(filtered_spn[[]], x = 'newannot_2', y = 'nCount_RNA', ylim = c(0,80000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# prevent warnings to become errors
options(warn = -1)

pdf(paste0(pref, "_NuclearFraction_Annotated_BroadAnnot.pdf"))
ggboxplot(filtered_spn[[]], x = 'newannot_2', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()


DefaultAssay(filtered_spn) = 'RNA'
filtered_spn = NormalizeData(filtered_spn)

pdf(paste0(pref, "_MarkersDotPlot_BroadAnnot.pdf"), width = 25, height = 10)
DotPlot(filtered_spn, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()


filtered_spn$Species = factor(filtered_spn$Species, levels = c('Ferret', 'Bat', 'Mouse', 'Marmoset', 'Macaque', 'Chimp', 'Human'))

stackedbarplot(filtered_spn[[]], groupfill = 'newannot_2', groupx = 'Species', fn = paste0(pref, '_Stacked_Annotated_Species_BroadAnnot'))

stackedbarplot(filtered_spn[[]], groupfill = 'Species', groupx = 'newannot_2', fn = paste0(pref, '_Stacked_Annotated_Species2_BroadAnnot'))

stackedbarplot(filtered_spn[[]], groupfill = 'newannot_2', groupx = 'orig.ident', fn = paste0(pref, '_Stacked_Annotated_Samples_BroadAnnot'))

stackedbarplot(filtered_spn[[]], groupfill = 'newannot_2', groupx = 'orig.ident', fn = paste0(pref, '_Stacked_Annotated_Samples_BroadAnnot'), hg = 20)

stackedbarplot(filtered_spn[[]], groupfill = 'orig.ident', groupx = 'newannot_2', fn = paste0(pref, '_Stacked_AnnotatedSamples2_BroadAnnot'), wd = 20)


stackedbarplot(filtered_spn[[]], groupfill = 'Tissue', groupx = 'newannot_2', fn = paste0(pref, '_Stacked_Annotated_Tissue_BroadAnnot'))

# Primate vs bat vs ferret
meta = filtered_spn@meta.data
meta$Species_Broad = 'Primate'
meta[meta$Species == 'Bat', 'Species_Broad'] = 'Bat'
meta[meta$Species == 'Ferret', 'Species_Broad'] = 'Ferret'

stackedbarplot(meta, groupfill = 'Species_Broad', groupx = 'newannot_2', fn = paste0(pref, '_Stacked_Annotated_Species_Broad_BroadAnnot'))


# stacked bar plot of cell types for only primates
sub = subset(meta, subset = Species_Broad == c('Primate'))
table(sub$Species)

stackedbarplot(sub, groupfill = 'Species', groupx = 'newannot_2', fn = paste0(pref, '_Stacked_Annotated_Species_Primates_BroadAnnot'))

table(filtered_spn$Species) # dims are [1]  14894 254583

#  Ferret      Bat    Mouse Marmoset  Macaque    Chimp    Human 
#   11049    59004    50605    41894    46442    20461    25128 


# save
saveRDS(filtered_spn, paste0("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/ANNOTATED/",pref, "_ANNOTATED.RDS"))
filtered_spn <- readRDS("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/ANNOTATED/GB_AllTissues_Primates_bat_mouse_ferret_filtered_Integrated_SPN_ncbi_human_pr_coding_orthologs_ANNOTATED.RDS")

# check if eSPNs show similar markers in all the species
# Plot only necessary markers across species
DefaultAssay(filtered_spn) = 'RNA'
split_colors = c('blue', 'orange', 'darkgreen', 'pink', 'deepskyblue2', 'tan3', 'violetred')

dotout = DotPlot(filtered_spn, features = c("DRD1", "DRD2", "CASZ1"), group.by = 'newannot_2', split.by = 'Species', cols = split_colors)
dotout = dotout$data
dotout$Species = gsub('.*_', '', dotout$id)
dotout$Species = factor(dotout$Species, levels = c('Human', 'Chimp', 'Macaque', 'Marmoset', 'Mouse',  'Ferret', 'Bat'))
dotout$CellType = gsub('_Human|_Chimp|_Macaque|_Marmoset|_Mouse|_Ferret|_Bat', '', dotout$id)


pdf(paste0("SPN_AllSpecies_SelectMarkersDotPlot_human_pr_coding_orthologs_filtered.pdf"), width = 30, height = 10)
ggscatter(dotout, y = 'features.plot', x = 'Species', color = 'Species',
			palette = split_colors, size = 'pct.exp') +
scale_size(range = c(-0.5,6)) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
theme(text=element_text(size=20, face = 'bold')) +
facet_wrap(~CellType, nrow = 1) +
rotate_x_text(45)
dev.off()
######################################################################################################
##!!! Remove ferret data!!!
filtered_spn_noFerret  <- subset(filtered_spn, subset = Species == "Ferret", invert = T)
filtered_spn_noFerret$Species = factor(filtered_spn_noFerret$Species, levels = c('Bat', 'Mouse', 'Marmoset', 'Macaque', 'Chimp', 'Human'))
# dims are [1]  14894 244258

# set variables
pref = 'GB_AllTissues_Primates_bat_mouse_filtered_Integrated_SPN_ncbi_human_pr_coding_orthologs'
setwd("/home2/gkonop/workdir/03_INTEGRATE_ALL/SPN/ANNOTATED")

# Basic Plots
## label for striosome and matrix
# markers
marksToPlot = c('SLC17A7', 'SATB2', 'RBFOX3', 'NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1','COL11A1', 'PAPC', 'DNER', 'DRD2', 'SPON1', 'BACH2','KCNIP1', 'EPHAA5', 'HTR2A', 'RXRG', 'SEMA5B', 'MFGE8', 'KREMEN1','PDE1C','OPRM1', 'CALB1', 'CRYM', 'ID4', 'ZFHX3', 'SGK1','SV2B','FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA')

Idents(filtered_spn_noFerret) = filtered_spn_noFerret$newannot

# PLOTS for striosome and matrix
pdf(paste0(pref, "_Clusters_Annotated.pdf"))
DimPlot(filtered_spn_noFerret, group.by = 'newannot', cols = c('iSPN_matrix' = 'springgreen3', 'iSPN_striosome' = 'violetred', 'dSPN_matrix' = 'deepskyblue2', 'dSPN_striosome' = 'darkorchid', 'eSPN' = 'chocolate2'), label = T, raster = T, repel = T, label.size = 7) +
NoLegend()
ggtitle('')
dev.off()

pdf(paste0(pref, "_Depth_Gene_Annotated.pdf"))
ggboxplot(filtered_spn_noFerret[[]], x = 'newannot', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_Depth_UMI_Annotated.pdf"))
ggboxplot(filtered_spn_noFerret[[]], x = 'newannot', y = 'nCount_RNA', ylim = c(0,80000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_NuclearFraction_Annotated.pdf"))
ggboxplot(filtered_spn_noFerret[[]], x = 'newannot', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()


DefaultAssay(filtered_spn_noFerret) = 'RNA'
filtered_spn_noFerret = NormalizeData(filtered_spn_noFerret)

pdf(paste0(pref, "_MarkersDotPlot.pdf"), width = 25, height = 10)
DotPlot(filtered_spn_noFerret, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

filtered_spn_noFerret$newannot = factor(filtered_spn_noFerret$newannot, levels = c("dSPN_matrix", "dSPN_striosome","iSPN_matrix", "iSPN_striosome","eSPN"))

filtered_spn_noFerret$newannot_2 = factor(filtered_spn_noFerret$newannot_2, levels = c("dSPN", "iSPN_matrix", "iSPN","eSPN"))

filtered_spn_noFerret$Species = factor(filtered_spn_noFerret$Species, levels = c('Bat', 'Mouse', 'Marmoset', 'Macaque', 'Chimp', 'Human'))

stackedbarplot(filtered_spn_noFerret[[]], groupfill = 'newannot', groupx = 'Species', fn = paste0(pref, '_Stacked_Annotated_Species'))

stackedbarplot(filtered_spn_noFerret[[]], groupfill = 'Species', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Species2'))

stackedbarplot(filtered_spn_noFerret[[]], groupfill = 'Species', groupx = 'newannot_2', fn = paste0(pref, '_Stacked_Annotated_Species2'))

stackedbarplot(filtered_spn_noFerret[[]], groupfill = 'newannot', groupx = 'orig.ident', fn = paste0(pref, '_Stacked_Annotated_Samples'), hg = 20)

stackedbarplot(filtered_spn_noFerret[[]], groupfill = 'orig.ident', groupx = 'newannot', fn = paste0(pref, '_Stacked_AnnotatedSamples2'), wd = 20)

stackedbarplot(filtered_spn_noFerret[[]], groupfill = 'Tissue', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Tissue'))

# Primate vs bat vs mouse
meta = filtered_spn_noFerret@meta.data
meta$Species_Broad = 'Primate'
meta[meta$Species == 'Bat', 'Species_Broad'] = 'Bat'
meta[meta$Species == 'Mouse', 'Species_Broad'] = 'Mouse'

stackedbarplot(meta, groupfill = 'Species_Broad', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Species_Primate_vsRest'))

# stacked bar plot of cell types for only primates
sub = subset(meta, subset = Species_Broad == c('Primate'))
table(sub$Species)
stackedbarplot(sub, groupfill = 'Species', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Species_Primates'))
######################################################################################
####@@@@@@@
# PLOTS for dSPN, iSPN, and eSPNs
filtered_spn_noFerret[["newannot_2"]]  = str_split_i(allseur_integrated_noFerret$newannot, "_",1)
Idents(allseur_integrated_noFerret) = allseur_integrated_noFerret$newannot_2

# PLOTS for striosome and matrix
pdf(paste0(pref, "_Clusters_Annotated_BroadAnnot.pdf"))
DimPlot(allseur_integrated_noFerret, group.by = 'newannot_2', label = T, raster = T, repel = T, label.size = 7) +
NoLegend() +
ggtitle('')
dev.off()

pdf(paste0(pref, "_Depth_Gene_Annotated_BroadAnnot.pdf"))
ggboxplot(allseur_integrated_noFerret[[]], x = 'newannot_2', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_Depth_UMI_Annotated_BroadAnnot.pdf"))
ggboxplot(allseur_integrated_noFerret[[]], x = 'newannot_2', y = 'nCount_RNA', ylim = c(0,80000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_NuclearFraction_Annotated_BroadAnnot.pdf"))
ggboxplot(allseur_integrated_noFerret[[]], x = 'newannot_2', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()


DefaultAssay(allseur_integrated_noFerret) = 'RNA'
allseur_integrated_noFerret = NormalizeData(allseur_integrated_noFerret)

pdf(paste0(pref, "_MarkersDotPlot_BroadAnnot.pdf"), width = 25, height = 10)
DotPlot(allseur_integrated_noFerret, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

allseur_integrated_noFerret$newannot = factor(allseur_integrated_noFerret$newannot_2, levels = c("dSPN", "iSPN", "eSPN"))

allseur_integrated_noFerret$Species = factor(allseur_integrated_noFerret$Species, levels = c('Bat', 'Mouse', 'Marmoset', 'Macaque', 'Chimp', 'Human'))

stackedbarplot(allseur_integrated_noFerret[[]], groupfill = 'newannot_2', groupx = 'Species', fn = paste0(pref, '_Stacked_Annotated_Species_BroadAnnot'))

stackedbarplot(allseur_integrated_noFerret[[]], groupfill = 'Species', groupx = 'newannot_2', fn = paste0(pref, '_Stacked_Annotated_Species2_BroadAnnot'))

stackedbarplot(allseur_integrated_noFerret[[]], groupfill = 'newannot_2', groupx = 'orig.ident', fn = paste0(pref, '_Stacked_Annotated_Samples_BroadAnnot'))

stackedbarplot(allseur_integrated_noFerret[[]], groupfill = 'newannot_2', groupx = 'orig.ident', fn = paste0(pref, '_Stacked_Annotated_Samples_BroadAnnot'), hg = 20)

stackedbarplot(allseur_integrated_noFerret[[]], groupfill = 'orig.ident', groupx = 'newannot_2', fn = paste0(pref, '_Stacked_AnnotatedSamples2_BroadAnnot'), wd = 20)


stackedbarplot(allseur_integrated_noFerret[[]], groupfill = 'Tissue', groupx = 'newannot_2', fn = paste0(pref, '_Stacked_Annotated_Tissue_BroadAnnot'))

# Primate vs bat vs mouse
meta = allseur_integrated_noFerret@meta.data
meta$Species_Broad = 'Primate'
meta[meta$Species == 'Bat', 'Species_Broad'] = 'Bat'
meta[meta$Species == 'Mouse', 'Species_Broad'] = 'Mouse'

stackedbarplot(meta, groupfill = 'Species_Broad', groupx = 'newannot_2', fn = paste0(pref, '_Stacked_Annotated_Species_Broad_BroadAnnot'))


# stacked bar plot of cell types for only primates
sub = subset(meta, subset = Species_Broad == c('Primate'))
table(sub$Species)

stackedbarplot(sub, groupfill = 'Species', groupx = 'newannot_2', fn = paste0(pref, '_Stacked_Annotated_Species_Primates_BroadAnnot'))

table(allseur_integrated_noFerret$Species)


# save
saveRDS(allseur_integrated_noFerret, paste0("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/ANNOTATED/",pref, "_ANNOTATED.RDS"))
##################################################################################################
##!!! Remove putamen data!!!
allseur_integrated_noPutamen  <- subset(allseur_integrated, subset = Tissue == "Putamen", invert = T)
allseur_integrated_noPutamen$Species = factor(allseur_integrated_noPutamen$Species, levels = c('Ferret','Bat', 'Mouse', 'Marmoset', 'Macaque', 'Chimp', 'Human'))
# dims are [1]  14894 162148

# Basic Plots
pref = 'GB_Caudate_Primates_bat_mouse_ferret_Integrated_SPN_ncbi_human_pr_coding_orthologs'

setwd("/home2/gkonop/workdir/03_INTEGRATE_ALL/SPN/ANNOTATED")

## label for striosome and matrix
# markers
marksToPlot = c('SLC17A7', 'SATB2', 'RBFOX3', 'NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1','COL11A1', 'PAPC', 'DNER', 'DRD2', 'SPON1', 'BACH2','KCNIP1', 'EPHAA5', 'HTR2A', 'RXRG', 'SEMA5B', 'MFGE8', 'KREMEN1','PDE1C','OPRM1', 'CALB1', 'CRYM', 'ID4', 'ZFHX3', 'SGK1','SV2B','FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA')

# change idents prior plotting
Idents(allseur_integrated_noPutamen) = allseur_integrated_noPutamen$newannot

# PLOTS for striosome and matrix
pdf(paste0(pref, "_Clusters_Annotated.pdf"))
DimPlot(allseur_integrated_noPutamen, group.by = 'newannot', cols = c('iSPN_matrix' = 'springgreen3', 'iSPN_striosome' = 'violetred', 'dSPN_matrix' = 'deepskyblue2', 'dSPN_striosome' = 'darkorchid', 'eSPN' = 'chocolate2'), label = T, raster = T, repel = T, label.size = 7) +
NoLegend()
ggtitle('')
dev.off()

pdf(paste0(pref, "_Depth_Gene_Annotated.pdf"))
ggboxplot(allseur_integrated_noPutamen[[]], x = 'newannot', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_Depth_UMI_Annotated.pdf"))
ggboxplot(allseur_integrated_noPutamen[[]], x = 'newannot', y = 'nCount_RNA', ylim = c(0,80000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_NuclearFraction_Annotated.pdf"))
ggboxplot(allseur_integrated_noPutamen[[]], x = 'newannot', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()


DefaultAssay(allseur_integrated_noPutamen) = 'RNA'
allseur_integrated_noPutamen = NormalizeData(allseur_integrated_noPutamen)

pdf(paste0(pref, "_MarkersDotPlot.pdf"), width = 25, height = 10)
DotPlot(allseur_integrated_noPutamen, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

allseur_integrated_noPutamen$newannot = factor(allseur_integrated_noPutamen$newannot, levels = c("dSPN_matrix", "dSPN_striosome","iSPN_matrix", "iSPN_striosome","eSPN"))

stackedbarplot(allseur_integrated_noPutamen[[]], groupfill = 'newannot', groupx = 'Species', fn = paste0(pref, '_Stacked_Annotated_Species'))

stackedbarplot(allseur_integrated_noPutamen[[]], groupfill = 'Species', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Species2'))

stackedbarplot(allseur_integrated_noPutamen[[]], groupfill = 'newannot', groupx = 'orig.ident', fn = paste0(pref, '_Stacked_Annotated_Samples'), hg = 20)

stackedbarplot(allseur_integrated_noPutamen[[]], groupfill = 'orig.ident', groupx = 'newannot', fn = paste0(pref, '_Stacked_AnnotatedSamples2'), wd = 20)

stackedbarplot(allseur_integrated_noPutamen[[]], groupfill = 'Tissue', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Tissue'))

# Primate vs bat vs mouse vs ferret
meta = allseur_integrated_noPutamen@meta.data
meta$Species_Broad = 'Primate'
meta[meta$Species == 'Bat', 'Species_Broad'] = 'Bat'
meta[meta$Species == 'Mouse', 'Species_Broad'] = 'Mouse'
meta[meta$Species == 'Ferret', 'Species_Broad'] = 'Ferret'

stackedbarplot(meta, groupfill = 'Species_Broad', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Species_Primate_vsRest'))

# stacked bar plot of cell types for only primates
sub = subset(meta, subset = Species_Broad == c('Primate'))
table(sub$Species)
stackedbarplot(sub, groupfill = 'Species', groupx = 'newannot', fn = paste0(pref, '_Stacked_Annotated_Species_Primates'))
######################################################################################
# PLOTS for dSPN, iSPN, and eSPNs
allseur_integrated_noPutamen[["newannot_2"]]  = str_split_i(allseur_integrated_noPutamen$newannot, "_",1)
Idents(allseur_integrated_noPutamen) = allseur_integrated_noPutamen$newannot_2

# PLOTS for striosome and matrix
pdf(paste0(pref, "_Clusters_Annotated_BroadAnnot.pdf"))
DimPlot(allseur_integrated_noPutamen, group.by = 'newannot_2', label = T, raster = T, repel = T, label.size = 7) +
NoLegend() +
ggtitle('')
dev.off()

pdf(paste0(pref, "_Depth_Gene_Annotated_BroadAnnot.pdf"))
ggboxplot(allseur_integrated_noPutamen[[]], x = 'newannot_2', y = 'nFeature_RNA') +
rotate_x_text(90)
dev.off()

pdf(paste0(pref, "_Depth_UMI_Annotated_BroadAnnot.pdf"))
ggboxplot(allseur_integrated_noPutamen[[]], x = 'newannot_2', y = 'nCount_RNA', ylim = c(0,80000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

pdf(paste0(pref, "_NuclearFraction_Annotated_BroadAnnot.pdf"))
ggboxplot(allseur_integrated_noPutamen[[]], x = 'newannot_2', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()


DefaultAssay(allseur_integrated_noPutamen) = 'RNA'
allseur_integrated_noPutamen = NormalizeData(allseur_integrated_noPutamen)

pdf(paste0(pref, "_MarkersDotPlot_BroadAnnot.pdf"), width = 25, height = 10)
DotPlot(allseur_integrated_noPutamen, features = marksToPlot) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

allseur_integrated_noPutamen$newannot = factor(allseur_integrated_noPutamen$newannot_2, levels = c("dSPN", "iSPN", "eSPN"))

allseur_integrated_noPutamen$Species = factor(allseur_integrated_noPutamen$Species, levels = c('Ferret','Bat', 'Mouse', 'Marmoset', 'Macaque', 'Chimp', 'Human'))

stackedbarplot(allseur_integrated_noPutamen[[]], groupfill = 'newannot_2', groupx = 'Species', fn = paste0(pref, '_Stacked_Annotated_Species_BroadAnnot'))

stackedbarplot(allseur_integrated_noPutamen[[]], groupfill = 'Species', groupx = 'newannot_2', fn = paste0(pref, '_Stacked_Annotated_Species2_BroadAnnot'))

stackedbarplot(allseur_integrated_noPutamen[[]], groupfill = 'newannot_2', groupx = 'orig.ident', fn = paste0(pref, '_Stacked_Annotated_Samples_BroadAnnot'))

stackedbarplot(allseur_integrated_noPutamen[[]], groupfill = 'newannot_2', groupx = 'orig.ident', fn = paste0(pref, '_Stacked_Annotated_Samples_BroadAnnot'), hg = 20)

stackedbarplot(allseur_integrated_noPutamen[[]], groupfill = 'orig.ident', groupx = 'newannot_2', fn = paste0(pref, '_Stacked_AnnotatedSamples2_BroadAnnot'), wd = 20)


stackedbarplot(allseur_integrated_noPutamen[[]], groupfill = 'Tissue', groupx = 'newannot_2', fn = paste0(pref, '_Stacked_Annotated_Tissue_BroadAnnot'))

# Primate vs bat vs mouse
meta = allseur_integrated_noPutamen@meta.data
meta$Species_Broad = 'Primate'
meta[meta$Species == 'Bat', 'Species_Broad'] = 'Bat'
meta[meta$Species == 'Mouse', 'Species_Broad'] = 'Mouse'
meta[meta$Species == 'Ferret', 'Species_Broad'] = 'Ferret'

stackedbarplot(meta, groupfill = 'Species_Broad', groupx = 'newannot_2', fn = paste0(pref, '_Stacked_Annotated_Species_Broad_BroadAnnot'))


# stacked bar plot of cell types for only primates
sub = subset(meta, subset = Species_Broad == c('Primate'))
table(sub$Species)

stackedbarplot(sub, groupfill = 'Species', groupx = 'newannot_2', fn = paste0(pref, '_Stacked_Annotated_Species_Primates_BroadAnnot'))

table(allseur_integrated_noPutamen$Species)
#  Ferret      Bat    Mouse Marmoset  Macaque    Chimp    Human 
#   11049    39751    50605    25923    17711     7840     9269 

# save
saveRDS(allseur_integrated_noPutamen, paste0("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/ANNOTATED/",pref, "_ANNOTATED.RDS"))

#############################################################################
## Compare and plot 
# change plot dir
setwd("/endosome/work/Neuroinformatics_Core/gkonop/03_INTEGRATE_ALL/SPN/ANNOTATED")

### perform stat test to compare the cell type proportions
## remove any sample which has < 100 SPNs
# extract sample names
samples <- table(allseur_integrated$orig.ident)
bad_samples <- names(samples[samples < 100])
filtered_spn <- subset(spn, orig.ident %in% bad_samples, invert = T)

meta <- filtered_spn[[]]
# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species) %>%
  summarize(
    totsize = n(),
    eSPNsize = sum(newannot_2 == 'eSPN'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    eSPN_to_allSPN = eSPNsize / totsize   
  )

df$Species <- factor(df$Species, levels = rev(c('Bat', 'Ferret', 'Mouse', 'Marmoset', 'Macaque', 'Chimp', 'Human')))

# ratio for Caudate for all species and for mouse CaudoPutamen 
# subset the dataset
df_no_put <- df[df$Tissue != "Putamen",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_put %>%
  group_by(Species) %>%
  summarize(
    mean_eSPN_to_allSPN = mean(eSPN_to_allSPN, na.rm = TRUE),
    sd_eSPN_to_allSPN = sd(eSPN_to_allSPN, na.rm = TRUE),
    n = sum(!is.na(eSPN_to_allSPN)),
    sem_eSPN_to_allSPN = sd_eSPN_to_allSPN / sqrt(n),
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

# change directory
setwd("/endosome/work/Neuroinformatics_Core/gkonop/03_INTEGRATE_ALL/SPN/ANNOTATED")
# Plot with p-values
pdf(file = "GB_barplot_eSPN_to_AllSPN_Ratio_Caudate_with_mouse_Caudoputamen_Across_Species_with_Wilcoxpvals_withFerret_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put, aes(x = Species, y = eSPN_to_allSPN, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_eSPN_to_allSPN, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_eSPN_to_allSPN - sem_eSPN_to_allSPN, ymax = mean_eSPN_to_allSPN + sem_eSPN_to_allSPN),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('eSPN / SPN Ratio in Caudate'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") + 
  stat_compare_means(comparisons = comps, method = 'wilcox.test', paired = FALSE)
)
dev.off()
  
# Plot with p-values
pdf(file = "GB_barplot_eSPN_to_AllSPN_Ratio_Caudate_with_mouse_Caudoputamen_Across_Species_with_t-testpvals_withFerret_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put, aes(x = Species, y = eSPN_to_allSPN, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_eSPN_to_allSPN, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_eSPN_to_allSPN - sem_eSPN_to_allSPN, ymax = mean_eSPN_to_allSPN + sem_eSPN_to_allSPN),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('eSPN / SPN Ratio in Caudate'),
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
pdf(file = "GB_barplot_eSPN_to_AllSPN_Ratio_Caudate_with_mouse_Caudoputamen_Across_Species_withFerret_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put, aes(x = Species, y = eSPN_to_allSPN, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_eSPN_to_allSPN, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_eSPN_to_allSPN - sem_eSPN_to_allSPN, ymax = mean_eSPN_to_allSPN + sem_eSPN_to_allSPN),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('eSPN / SPN Ratio in Caudate'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") + coord_flip()
)
dev.off()

####
# Primates vs non-primates
meta = allseur_integrated@meta.data
meta$Species_Broad = 'Primate'
meta[meta$Species %in% c('Bat','Ferret','Mouse'), 'Species_Broad'] = 'Non_primate'

# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species_Broad) %>%
  summarize(
    totsize = n(),
    eSPNsize = sum(newannot_2 == 'eSPN'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    eSPN_to_allSPN = eSPNsize / totsize   
  )

df$Species_Broad <- factor(df$Species_Broad, levels = rev(c('Primate', 'Non_primate')))

# ratio for Caudate for all species and for mouse CaudoPutamen 
# subset the dataset
df_no_put <- df[df$Tissue != "Putamen",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_put %>%
  group_by(Species_Broad) %>%
  summarize(
    mean_eSPN_to_allSPN = mean(eSPN_to_allSPN, na.rm = TRUE),
    sd_eSPN_to_allSPN = sd(eSPN_to_allSPN, na.rm = TRUE),
    n = sum(!is.na(eSPN_to_allSPN)),
    sem_eSPN_to_allSPN = sd_eSPN_to_allSPN / sqrt(n),
    .groups = 'drop'
  )
    
# Comparison with mouse (caudoputamen)
comps <- list(
  c('Non_primate', 'Primate')
)

# change directory
setwd("/endosome/work/Neuroinformatics_Core/gkonop/03_INTEGRATE_ALL/SPN/ANNOTATED")
# Plot with p-values
pdf(file = "GB_barplot_eSPN_to_AllSPN_Ratio_Caudate_with_mouse_Caudoputamen_Across_BroadSpecies_with_Wilcoxpvals_withFerret_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put, aes(x = Species_Broad, y = eSPN_to_allSPN, color = Species_Broad)) +
  geom_col(data = df_summary, aes(x = Species_Broad, y = mean_eSPN_to_allSPN, fill = Species_Broad), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species_Broad, ymin = mean_eSPN_to_allSPN - sem_eSPN_to_allSPN, ymax = mean_eSPN_to_allSPN + sem_eSPN_to_allSPN),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('eSPN / SPN Ratio in Caudate'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") + 
  stat_compare_means(comparisons = comps, method = 'wilcox.test', paired = FALSE)
)
dev.off()
  
# Plot with p-values
pdf(file = "GB_barplot_eSPN_to_AllSPN_Ratio_Caudate_with_mouse_Caudoputamen_Across_BroadSpecies_with_t-testpvals_withFerret_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put, aes(x = Species_Broad, y = eSPN_to_allSPN, color = Species_Broad)) +
  geom_col(data = df_summary, aes(x = Species_Broad, y = mean_eSPN_to_allSPN, fill = Species_Broad), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species_Broad, ymin = mean_eSPN_to_allSPN - sem_eSPN_to_allSPN, ymax = mean_eSPN_to_allSPN + sem_eSPN_to_allSPN),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('eSPN / SPN Ratio in Caudate'),
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
###############################################################
# ratio for Caudate for all species except Ferret and for mouse CaudoPutamen 
# extract the metadata
meta <- filtered_spn[[]]
# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species) %>%
  summarize(
    totsize = n(),
    eSPNsize = sum(newannot_2 == 'eSPN'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    eSPN_to_allSPN = eSPNsize / totsize   
  )

df$Species <- factor(df$Species, levels = c('Bat', 'Ferret', 'Mouse', 'Marmoset', 'Macaque', 'Chimp', 'Human'))
# subset the dataset
df_no_put <- df[df$Tissue != "Putamen",] 
df_no_put_noFerret <- df_no_put[df_no_put$Species != "Ferret",]

# Calculate summary statistics for the bar plots
df_summary <- df_no_put_noFerret %>%
  group_by(Species) %>%
  summarize(
    mean_eSPN_to_allSPN = mean(eSPN_to_allSPN, na.rm = TRUE),
    sd_eSPN_to_allSPN = sd(eSPN_to_allSPN, na.rm = TRUE),
    n = sum(!is.na(eSPN_to_allSPN)),
    sem_eSPN_to_allSPN = sd_eSPN_to_allSPN / sqrt(n),
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

# change directory
setwd("/endosome/work/Neuroinformatics_Core/gkonop/03_INTEGRATE_ALL/SPN/ANNOTATED")

# Plot with p-values
pdf(file = "GB_barplot_eSPN_to_AllSPN_Ratio_Caudate_with_mouse_Caudoputamen_Across_Species_with_Wilcoxpvals_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put_noFerret, aes(x = Species, y = eSPN_to_allSPN, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_eSPN_to_allSPN, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_eSPN_to_allSPN - sem_eSPN_to_allSPN, ymax = mean_eSPN_to_allSPN + sem_eSPN_to_allSPN),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('eSPN / SPN Ratio in Caudate'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") + 
  stat_compare_means(comparisons = comps, method = 'wilcox.test', paired = FALSE)
)
dev.off()
  
# Plot with p-values
pdf(file = "GB_barplot_eSPN_to_AllSPN_Ratio_Caudate_with_mouse_Caudoputamen_Across_Species_with_t-testpvals_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put_noFerret, aes(x = Species, y = eSPN_to_allSPN, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_eSPN_to_allSPN, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_eSPN_to_allSPN - sem_eSPN_to_allSPN, ymax = mean_eSPN_to_allSPN + sem_eSPN_to_allSPN),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('eSPN / SPN Ratio in Caudate'),
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
pdf(file = "GB_barplot_eSPN_to_AllSPN_Ratio_Caudate_with_mouse_Caudoputamen_Across_Species_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put_noFerret, aes(x = Species, y = eSPN_to_allSPN, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_eSPN_to_allSPN, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_eSPN_to_allSPN - sem_eSPN_to_allSPN, ymax = mean_eSPN_to_allSPN + sem_eSPN_to_allSPN),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('eSPN / SPN Ratio in Caudate'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") + coord_flip()
)
dev.off()

## primate vs non-primate (Caudate w/o ferret)
meta = allseur_integrated@meta.data
meta$Species_Broad = 'Primate'
meta[meta$Species %in% c('Mouse','Bat'), 'Species_Broad'] = 'Non_primate'


# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species_Broad) %>%
  summarize(
    totsize = n(),
    eSPNsize = sum(newannot_2 == 'eSPN'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    eSPN_to_allSPN = eSPNsize / totsize   
  )

df$Species_Broad <- factor(df$Species_Broad, levels = rev(c('Primate', 'Non_primate')))

# ratio for Caudate for all species and for mouse CaudoPutamen 
# subset the dataset
df_no_put <- df[df$Tissue != "Putamen",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_put %>%
  group_by(Species_Broad) %>%
  summarize(
    mean_eSPN_to_allSPN = mean(eSPN_to_allSPN, na.rm = TRUE),
    sd_eSPN_to_allSPN = sd(eSPN_to_allSPN, na.rm = TRUE),
    n = sum(!is.na(eSPN_to_allSPN)),
    sem_eSPN_to_allSPN = sd_eSPN_to_allSPN / sqrt(n),
    .groups = 'drop'
  )
    
# Comparison
comps <- list(
  c('Non_primate', 'Primate')
)

# change directory
setwd("/endosome/work/Neuroinformatics_Core/gkonop/03_INTEGRATE_ALL/SPN/ANNOTATED")
# Plot with p-values
pdf(file = "GB_barplot_eSPN_to_AllSPN_Ratio_Caudate_with_mouse_Caudoputamen_Across_BroadSpecies_withoutFerret_with_Wilcoxpvals_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put, aes(x = Species_Broad, y = eSPN_to_allSPN, color = Species_Broad)) +
  geom_col(data = df_summary, aes(x = Species_Broad, y = mean_eSPN_to_allSPN, fill = Species_Broad), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species_Broad, ymin = mean_eSPN_to_allSPN - sem_eSPN_to_allSPN, ymax = mean_eSPN_to_allSPN + sem_eSPN_to_allSPN),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('eSPN / SPN Ratio in Caudate'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") + 
  stat_compare_means(comparisons = comps, method = 'wilcox.test', paired = FALSE)
)
dev.off()
  
# Plot with p-values
pdf(file = "GB_barplot_eSPN_to_AllSPN_Ratio_Caudate_with_mouse_Caudoputamen_Across_BroadSpecies_withoutFerret_with_t-testpvals_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put, aes(x = Species_Broad, y = eSPN_to_allSPN, color = Species_Broad)) +
  geom_col(data = df_summary, aes(x = Species_Broad, y = mean_eSPN_to_allSPN, fill = Species_Broad), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species_Broad, ymin = mean_eSPN_to_allSPN - sem_eSPN_to_allSPN, ymax = mean_eSPN_to_allSPN + sem_eSPN_to_allSPN),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('eSPN / SPN Ratio in Caudate'),
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


## primate vs non-primate
meta = allseur_integrated@meta.data
meta$Species_Broad = 'Primate'
meta[meta$Species %in% c('Mouse','Bat', 'Ferret'), 'Species_Broad'] = 'Non_primate'


# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species_Broad) %>%
  summarize(
    totsize = n(),
    eSPNsize = sum(newannot_2 == 'eSPN'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    eSPN_to_allSPN = eSPNsize / totsize   
  )

df$Species_Broad <- factor(df$Species_Broad, levels = rev(c('Primate', 'Non_primate')))

# ratio for Caudate for all species and for mouse CaudoPutamen 
# subset the dataset
df_no_put <- df[df$Tissue != "Putamen",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_put %>%
  group_by(Species_Broad) %>%
  summarize(
    mean_eSPN_to_allSPN = mean(eSPN_to_allSPN, na.rm = TRUE),
    sd_eSPN_to_allSPN = sd(eSPN_to_allSPN, na.rm = TRUE),
    n = sum(!is.na(eSPN_to_allSPN)),
    sem_eSPN_to_allSPN = sd_eSPN_to_allSPN / sqrt(n),
    .groups = 'drop'
  )
    
# Comparison
comps <- list(
  c('Non_primate', 'Primate')
)

# change directory
setwd("/endosome/work/Neuroinformatics_Core/gkonop/03_INTEGRATE_ALL/SPN/ANNOTATED")
# Plot with p-values
pdf(file = "GB_barplot_eSPN_to_AllSPN_Ratio_Caudate_with_mouse_Caudoputamen_Across_BroadSpecies_withoutFerret_with_Wilcoxpvals_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put, aes(x = Species_Broad, y = eSPN_to_allSPN, color = Species_Broad)) +
  geom_col(data = df_summary, aes(x = Species_Broad, y = mean_eSPN_to_allSPN, fill = Species_Broad), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species_Broad, ymin = mean_eSPN_to_allSPN - sem_eSPN_to_allSPN, ymax = mean_eSPN_to_allSPN + sem_eSPN_to_allSPN),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('eSPN / SPN Ratio in Caudate'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") + 
  stat_compare_means(comparisons = comps, method = 'wilcox.test', paired = FALSE)
)
dev.off()
  
# Plot with p-values
pdf(file = "GB_barplot_eSPN_to_AllSPN_Ratio_Caudate_with_mouse_Caudoputamen_Across_BroadSpecies_withFerret_with_t-testpvals_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put, aes(x = Species_Broad, y = eSPN_to_allSPN, color = Species_Broad)) +
  geom_col(data = df_summary, aes(x = Species_Broad, y = mean_eSPN_to_allSPN, fill = Species_Broad), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species_Broad, ymin = mean_eSPN_to_allSPN - sem_eSPN_to_allSPN, ymax = mean_eSPN_to_allSPN + sem_eSPN_to_allSPN),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('eSPN / SPN Ratio in Caudate'),
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
#################################################################################
#### odds ratio
meta = filtered_spn@meta.data
meta$Species_Broad = 'Primate'
meta[meta$Species %in% c('Mouse','Bat','Ferret'), 'Species_Broad'] = 'Non_primate'

# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species_Broad) %>%
  summarize(
    totsize = n(),
    eSPNsize = sum(newannot_2 == 'eSPN'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    eSPN_to_allSPN = eSPNsize / totsize   
  )

df$Species_Broad <- factor(df$Species_Broad, levels = rev(c('Primate', 'Non_primate')))

### caudate with ferret
# subset the dataset
df_no_put <- df[df$Tissue != "Putamen",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_put %>%
  group_by(Species_Broad) %>%
  summarize(
    mean_eSPN_to_allSPN = mean(eSPN_to_allSPN, na.rm = TRUE),
    sd_eSPN_to_allSPN = sd(eSPN_to_allSPN, na.rm = TRUE),
    n = sum(!is.na(eSPN_to_allSPN)),
    sem_eSPN_to_allSPN = sd_eSPN_to_allSPN / sqrt(n),
    .groups = 'drop'
  )
# odds ratio
df_summary <- as.data.frame(df_summary)
df_summary[df_summary$Species_Broad == "Primate",]$mean_eSPN_to_allSPN / df_summary[df_summary$Species_Broad == "Non_primate",]$mean_eSPN_to_allSPN #[1] 2.357699

### caudate w/o ferret
no_ferret <- subset(filtered_spn, Species == "Ferret", invert = T)
meta = no_ferret@meta.data
meta$Species_Broad = 'Primate'
meta[meta$Species %in% c('Mouse','Bat'), 'Species_Broad'] = 'Non_primate'

# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species_Broad) %>%
  summarize(
    totsize = n(),
    eSPNsize = sum(newannot_2 == 'eSPN'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    eSPN_to_allSPN = eSPNsize / totsize   
  )

df$Species_Broad <- factor(df$Species_Broad, levels = rev(c('Primate', 'Non_primate')))


df_no_put <- df[df$Tissue != "Putamen",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_put %>%
  group_by(Species_Broad) %>%
  summarize(
    mean_eSPN_to_allSPN = mean(eSPN_to_allSPN, na.rm = TRUE),
    sd_eSPN_to_allSPN = sd(eSPN_to_allSPN, na.rm = TRUE),
    n = sum(!is.na(eSPN_to_allSPN)),
    sem_eSPN_to_allSPN = sd_eSPN_to_allSPN / sqrt(n),
    .groups = 'drop'
  )

# convert tibble to dataframe
df_summary <- as.data.frame(df_summary)
df_summary[df_summary$Species_Broad == "Primate",]$mean_eSPN_to_allSPN / df_summary[df_summary$Species_Broad == "Non_primate",]$mean_eSPN_to_allSPN #[1] 2.576812

### putamen
# ratio for Putamen for all species and for mouse CaudoPutamen 
# subset the dataset
df_no_caud <- df[df$Tissue != "Caudate",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_caud %>%
  group_by(Species_Broad) %>%
  summarize(
    mean_eSPN_to_allSPN = mean(eSPN_to_allSPN, na.rm = TRUE),
    sd_eSPN_to_allSPN = sd(eSPN_to_allSPN, na.rm = TRUE),
    n = sum(!is.na(eSPN_to_allSPN)),
    sem_eSPN_to_allSPN = sd_eSPN_to_allSPN / sqrt(n),
    .groups = 'drop'
  )
df_summary <- as.data.frame(df_summary)
df_summary[df_summary$Species_Broad == "Primate",]$mean_eSPN_to_allSPN / df_summary[df_summary$Species_Broad == "Non_primate",]$mean_eSPN_to_allSPN #[1] 2.227278
    
###############################################################
# ratio for Putamen for all species except Ferret 
# subset the dataset
df_no_caud_noFerret <- df[df$Tissue != "Caudate",] 

meta = df_no_caud_noFerret@meta.data
meta$Species_Broad = 'Primate'
meta[meta$Species %in% c('Mouse','Bat'), 'Species_Broad'] = 'Non_primate'

# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species) %>%
  summarize(
    totsize = n(),
    eSPNsize = sum(newannot_2 == 'eSPN'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    eSPN_to_allSPN = eSPNsize / totsize   
  )

# Calculate summary statistics for the bar plots
df_summary <- df %>%
  group_by(Species) %>%
  summarize(
    mean_eSPN_to_allSPN = mean(eSPN_to_allSPN, na.rm = TRUE),
    sd_eSPN_to_allSPN = sd(eSPN_to_allSPN, na.rm = TRUE),
    n = sum(!is.na(eSPN_to_allSPN)),
    sem_eSPN_to_allSPN = sd_eSPN_to_allSPN / sqrt(n),
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

# change directory
setwd("/endosome/work/Neuroinformatics_Core/gkonop/03_INTEGRATE_ALL/SPN/ANNOTATED")
# Plot with p-values
pdf(file = "GB_barplot_eSPN_to_AllSPN_Ratio_Putamen_with_mouse_Caudoputamen_Across_Species_with_Wilcoxpvals_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df, aes(x = Species, y = eSPN_to_allSPN, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_eSPN_to_allSPN, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_eSPN_to_allSPN - sem_eSPN_to_allSPN, ymax = mean_eSPN_to_allSPN + sem_eSPN_to_allSPN),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('eSPN / SPN Ratio in Putamen'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") +
    geom_text(aes(label = orig.ident), 
              position = position_jitter(width = 0.2, height = 0), 
              size = 3, 
              hjust = 0.5, 
              vjust = -0.5) +
  stat_compare_means(comparisons = comps, method = 'wilcox.test', paired = FALSE)
)
dev.off()
 
# t-test
pdf(file = "GB_barplot_eSPN_to_AllSPN_Ratio_Putamen_with_mouse_Caudoputamen_Across_Species_with_t-testpvals_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df, aes(x = Species, y = eSPN_to_allSPN, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_eSPN_to_allSPN, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_eSPN_to_allSPN - sem_eSPN_to_allSPN, ymax = mean_eSPN_to_allSPN + sem_eSPN_to_allSPN),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('eSPN / SPN Ratio in Putamen'),
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
###@@@@@@

# Plot without p-values
pdf(file = "GB_barplot_eSPN_to_AllSPN_Ratio_Putamen_with_mouse_Caudoputamen_Across_Species_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df, aes(x = Species, y = eSPN_to_allSPN, color = Species)) + ylim(0,0.77)+
  geom_col(data = df_summary, aes(x = Species, y = mean_eSPN_to_allSPN, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_eSPN_to_allSPN - sem_eSPN_to_allSPN, ymax = mean_eSPN_to_allSPN + sem_eSPN_to_allSPN),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('eSPN / SPN Ratio in Putamen'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  scale_y_cut(breaks=c(0.2, 0.76), which=c(2,0,2), scales=c(0.25,0,0.25))  + # introduce breaks to y axis
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black")
)
dev.off()

####
#Primates vs non-primates
meta = allseur_integrated@meta.data
meta$Species_Broad = 'Primate'
meta[meta$Species %in% c('Mouse','Bat'), 'Species_Broad'] = 'Non_primate'


# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species_Broad) %>%
  summarize(
    totsize = n(),
    eSPNsize = sum(newannot_2 == 'eSPN'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    eSPN_to_allSPN = eSPNsize / totsize   
  )

df$Species_Broad <- factor(df$Species_Broad, levels = rev(c('Primate', 'Non_primate')))

# ratio for Putamen for all species and for mouse CaudoPutamen 
# subset the dataset
df_no_caud <- df[df$Tissue != "Caudate",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_caud %>%
  group_by(Species_Broad) %>%
  summarize(
    mean_eSPN_to_allSPN = mean(eSPN_to_allSPN, na.rm = TRUE),
    sd_eSPN_to_allSPN = sd(eSPN_to_allSPN, na.rm = TRUE),
    n = sum(!is.na(eSPN_to_allSPN)),
    sem_eSPN_to_allSPN = sd_eSPN_to_allSPN / sqrt(n),
    .groups = 'drop'
  )
    
# Comparison with mouse (caudoputamen)
comps <- list(
  c('Non_primate', 'Primate')
)

# change directory
setwd("/endosome/work/Neuroinformatics_Core/gkonop/03_INTEGRATE_ALL/SPN/ANNOTATED")
# Plot with p-values
pdf(file = "GB_barplot_eSPN_to_AllSPN_Ratio_Putamen_with_mouse_Caudoputamen_Across_BroadSpecies_with_Wilcoxpvals_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_caud, aes(x = Species_Broad, y = eSPN_to_allSPN, color = Species_Broad)) +
  geom_col(data = df_summary, aes(x = Species_Broad, y = mean_eSPN_to_allSPN, fill = Species_Broad), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species_Broad, ymin = mean_eSPN_to_allSPN - sem_eSPN_to_allSPN, ymax = mean_eSPN_to_allSPN + sem_eSPN_to_allSPN),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('eSPN / SPN Ratio in Putamen'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") + 
  stat_compare_means(comparisons = comps, method = 'wilcox.test', paired = FALSE)
)
dev.off()
  
# Plot with p-values
pdf(file = "GB_barplot_eSPN_to_AllSPN_Ratio_Putamen_with_mouse_Caudoputamen_Across_BroadSpecies_with_t-testpvals_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_caud, aes(x = Species_Broad, y = eSPN_to_allSPN, color = Species_Broad)) +
  geom_col(data = df_summary, aes(x = Species_Broad, y = mean_eSPN_to_allSPN, fill = Species_Broad), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species_Broad, ymin = mean_eSPN_to_allSPN - sem_eSPN_to_allSPN, ymax = mean_eSPN_to_allSPN + sem_eSPN_to_allSPN),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('eSPN / SPN Ratio in Putamen'),
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
###############################################################
## dSPN/iSPN
setwd("/endosome/work/Neuroinformatics_Core/gkonop/03_INTEGRATE_ALL/SPN/ANNOTATED")
## with ferret ! NO PUTAMEN !
meta = allseur_integrated[[]]

# Calculate the ratio for tissues
df <- meta %>%
  group_by(orig.ident, Tissue, Species) %>%
  summarize(
    totsize = n(),
    dSPNsize = sum(newannot_2 == 'dSPN'),
    iSPNsize = sum(newannot_2 == 'iSPN'),
    .groups = 'drop'
  )

# Compute ratios
df <- df %>%
  mutate(
    dSPN_to_iSPN = dSPNsize / iSPNsize   
  )

df$Species <- factor(df$Species, levels = c('Bat','Ferret', 'Mouse','Marmoset', 'Macaque', 'Chimp','Human'))

## ratio for Caudate for all species and for mouse CaudoPutamen 
# subset the dataset
df_no_put <- df[df$Tissue != "Putamen",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_put %>%
  group_by(Species) %>%
  summarize(
    mean_dSPN_to_iSPN = mean(dSPN_to_iSPN, na.rm = TRUE),
    sd_dSPN_to_iSPN = sd(dSPN_to_iSPN, na.rm = TRUE),
    n = sum(!is.na(dSPN_to_iSPN)),
    sem_dSPN_to_iSPN = sd_dSPN_to_iSPN / sqrt(n),
    .groups = 'drop'
  )
    
# Comparison with mouse (caudate+putamen)
comps <- comps <- list(
  c('Ferret', 'Human'),
  c('Mouse', 'Human'), 
  c('Marmoset', 'Human'), 
  c('Macaque', 'Human'), 
  c('Chimp', 'Human'), 
  c('Bat', 'Human')
)


# Define fixed jitter position with seed
fixed_jitter <- position_jitter(width = 0.2, height = 0, seed = 123)

# Plot with p-values two tailed t-test
pdf(file = "GB_barplot_dSPN_to_iSPN_Ratio_Caudate_with_mouse_Caudoputamen_Across_Species_with_TwoTailedt-test_pvals_withFerret_comparingEverthingToHuman_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put, aes(x = Species, y = dSPN_to_iSPN, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_dSPN_to_iSPN, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_dSPN_to_iSPN - sem_dSPN_to_iSPN, ymax = mean_dSPN_to_iSPN + sem_dSPN_to_iSPN),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('dSPN / iSPN Ratio in Caudate'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(position = fixed_jitter, size = 2, alpha = 0.7, color = "black") +
  stat_compare_means(comparisons = comps, method = 't.test', paired = FALSE) + coord_flip()
)
dev.off()

# Plot with Wilcoxon rank-sum p-values
pdf(file = "GB_barplot_dSPN_to_iSPN_Ratio_Caudate_with_mouse_Caudoputamen_Across_Species_with_Wilcoxon_pvals_withFerret_comparingEverthingToHuman_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put, aes(x = Species, y = dSPN_to_iSPN, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_dSPN_to_iSPN, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_dSPN_to_iSPN - sem_dSPN_to_iSPN, ymax = mean_dSPN_to_iSPN + sem_dSPN_to_iSPN),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('dSPN / iSPN Ratio in Caudate'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(position = fixed_jitter, size = 2, alpha = 0.7, color = "black") +
  stat_compare_means(comparisons = comps, method = 'wilcox.test')  + coord_flip()
)
dev.off()

  
# Plot without p-values
pdf(file = "GB_barplot_dSPN_to_iSPN_Ratio_Caudate_with_mouse_Caudoputamen_Across_Species_withFerret.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_put, aes(x = Species, y = dSPN_to_iSPN, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_dSPN_to_iSPN, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_dSPN_to_iSPN - sem_dSPN_to_iSPN, ymax = mean_dSPN_to_iSPN + sem_dSPN_to_iSPN),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('dSPN / iSPN Ratio in Caudate'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(position = fixed_jitter, size = 2, alpha = 0.7, color = "black")  + coord_flip()
)
dev.off()

## ratio for Putamen for all species and for mouse CaudoPutamen 
# subset the dataset
df_no_caud <- df[df$Tissue != "Caudate",] 

# Calculate summary statistics for the bar plots
df_summary <- df_no_caud %>%
  group_by(Species) %>%
  summarize(
    mean_dSPN_to_iSPN = mean(dSPN_to_iSPN, na.rm = TRUE),
    sd_dSPN_to_iSPN = sd(dSPN_to_iSPN, na.rm = TRUE),
    n = sum(!is.na(dSPN_to_iSPN)),
    sem_dSPN_to_iSPN = sd_dSPN_to_iSPN / sqrt(n),
    .groups = 'drop'
  )
    
# Comparison with mouse (caudate+putamen)
comps <- comps <- list(
  c('Mouse', 'Human'), 
  c('Marmoset', 'Human'), 
  c('Macaque', 'Human'), 
  c('Chimp', 'Human'), 
  c('Bat', 'Human')
)


# Define fixed jitter position with seed
fixed_jitter <- position_jitter(width = 0.2, height = 0, seed = 123)

# Plot with p-values two tailed t-test
pdf(file = "GB_barplot_dSPN_to_iSPN_Ratio_Putamen_with_mouse_Caudoputamen_Across_Species_with_TwoTailedt-test_pvals_withFerret_comparingEverthingToHuman_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_caud, aes(x = Species, y = dSPN_to_iSPN, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_dSPN_to_iSPN, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_dSPN_to_iSPN - sem_dSPN_to_iSPN, ymax = mean_dSPN_to_iSPN + sem_dSPN_to_iSPN),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('dSPN / iSPN Ratio in Putamen'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(position = fixed_jitter, size = 2, alpha = 0.7, color = "black") +
  stat_compare_means(comparisons = comps, method = 't.test', paired = FALSE) + coord_flip()
)
dev.off()

# Plot with Wilcoxon rank-sum p-values
pdf(file = "GB_barplot_dSPN_to_iSPN_Ratio_Putamen_with_mouse_Caudoputamen_Across_Species_with_Wilcoxon_pvals_withFerret_comparingEverthingToHuman_NCBI_human_pr_coding_orthologs.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_caud, aes(x = Species, y = dSPN_to_iSPN, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_dSPN_to_iSPN, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_dSPN_to_iSPN - sem_dSPN_to_iSPN, ymax = mean_dSPN_to_iSPN + sem_dSPN_to_iSPN),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('dSPN / iSPN Ratio in Putamen'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(position = fixed_jitter, size = 2, alpha = 0.7, color = "black") +
  stat_compare_means(comparisons = comps, method = 'wilcox.test')  + coord_flip()
)
dev.off()

  
# Plot without p-values
pdf(file = "GB_barplot_dSPN_to_iSPN_Ratio_Putamen_with_mouse_Caudoputamen_Across_Species_withFerret.pdf", width = 6, height = 6)
print(
  ggplot(data = df_no_caud, aes(x = Species, y = dSPN_to_iSPN, color = Species)) +
  geom_col(data = df_summary, aes(x = Species, y = mean_dSPN_to_iSPN, fill = Species), position = position_dodge()) +
  geom_errorbar(inherit.aes = FALSE, data = df_summary, aes(x = Species, ymin = mean_dSPN_to_iSPN - sem_dSPN_to_iSPN, ymax = mean_dSPN_to_iSPN + sem_dSPN_to_iSPN),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = paste('dSPN / iSPN Ratio in Putamen'),
    fill = 'Species',
    color = 'Species'
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  geom_jitter(position = fixed_jitter, size = 2, alpha = 0.7, color = "black")  + coord_flip()
)
dev.off()
#############################################################################


