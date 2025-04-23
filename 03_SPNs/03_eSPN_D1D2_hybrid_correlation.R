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
library(stringr)
source('~/SCRIPTS/utility_functions.R')
set.seed(1234)
# prevent warnings become errors
old_ops <- options(warn = 1)

warning("this is a warning")
#> Error in eval(expr, envir, enclos): (converted from warning) this is a warning

x <- "a"
as.numeric(x)
#> Error in eval(expr, envir, enclos): (converted from warning) NAs introduced by coercion

options(old_ops)
###########################################################
## data
# our data
our_spn_seu <- readRDS("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/ANNOTATED/GB_SPN_primates_bat_mouse_and_ferret_striosome_matrix_ANNOTATED.RDS")

# He at al macaque SPN processed data
#wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167920/suppl/GSE167920%5FResults%5FMSNs%5Fprocessed%5Ffinal.rds.gz
#gunzip GSE167920_Results_MSNs_processed_final.rds.gz 

others_macaque_spn_seu <- readRDS("/home2/gkonop/workdir/03_INTEGRATE_ALL/SPN/He_etal_macaque_data/GSE167920_Results_MSNs_processed_final.rds")

# markers
marksToPlot = c('SLC17A7', 'SATB2', 'RBFOX3', 'NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1','COL11A1', 'PAPC', 'DNER', 'DRD2', 'SPON1', 'BACH2','KCNIP1', 'EPHAA5', 'HTR2A', 'RXRG', 'SEMA5B', 'MFGE8', 'KREMEN1','PDE1C','OPRM1', 'CALB1', 'CRYM', 'ID4', 'ZFHX3', 'SGK1','SV2B','FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA')
##############################################################################################
## perform correlation analysis between the markers of eSPNs and the D1/D2 hybrid cells

# subset their gene exp matrix with the same genes as ours 
our_genes <- rownames(our_spn_seu@assays$RNA@counts)
their_genes <- rownames(others_macaque_spn_seu@assays$RNA@counts)

# Find common genes
common_genes <- intersect(our_genes, their_genes)

# Subset expression matrices by common genes
our_counts_common <- our_spn_seu@assays$RNA@counts[common_genes, ]
their_counts_common <- others_macaque_spn_seu@assays$RNA@counts[common_genes, ]

# Find genes with non-zero expression in both datasets
common_genes_1 <- rownames(our_counts_common[rowSums(our_counts_common > 0) > 0, ])
common_genes_2 <- rownames(their_counts_common[rowSums(their_counts_common > 0) > 0, ])

# Check the lengths
length(common_genes_1) # 11966
length(common_genes_2) # 12030
# continue with common_genes_1
############################################
# ours
# remove lowly abundant genes to improve comparison
dim(our_spn_seu@assays$RNA@counts[common_genes_1,])#[1] 12039 244824
sum(rowSums(our_spn_seu@assays$RNA@counts[common_genes_1,]) > 0) #[1] 11966

# get the count matrix
counts_data_our <- our_spn_seu@assays$RNA@counts[common_genes_1,]
# remove genes with zero counts
counts_data_filtered_our <- counts_data_our[rowSums(counts_data_our) > 0,]

# remove the genes with zero count in our count matrix
genes_to_remove <- (setdiff(rownames(counts_data_filtered_our), rownames(counts_data_filtered_their)))

# Remove those genes from `counts_data_filtered_their`
counts_data_filtered_our_2 <- counts_data_filtered_our[!rownames(counts_data_filtered_our) %in% genes_to_remove, ]

# Create a new Seurat object with the filtered counts data
our_spn_seu_filtered <- CreateSeuratObject(counts = counts_data_filtered_our_2, meta.data = our_spn_seu@meta.data)

# Copy over other metadata 
our_spn_seu_filtered$newannot <- our_spn_seu$newannot
our_spn_seu_filtered$Tissue <- our_spn_seu$Tissue

## get normalized counts:
# normalize data
our_spn_seu_filtered <- NormalizeData(our_spn_seu_filtered)
# get normalized data
our_normalized_data <- GetAssayData(our_spn_seu_filtered, slot = "data")

# generate new column with only eSPN, dSPN, and iSPN
our_spn_seu_filtered$broadannot <- str_split_i(our_spn_seu_filtered$newannot, "_", 1)

# Group by clusters
Idents(our_spn_seu_filtered) <- our_spn_seu_filtered$broadannot
our_cluster_ids <- Idents(our_spn_seu_filtered)

# Calculate the mean expression for each gene across each cluster
# We'll calculate the row means per cluster for each gene
cls_2 = unique(our_cluster_ids)
our_cluster_means <- sapply(cls_2, function(cluster) {
  cluster_cells <- which(our_cluster_ids == cluster)
  rowMeans(our_normalized_data[, cluster_cells], na.rm = TRUE)
})

# Convert the result to a data frame for easier viewing
our_cluster_means_df <- as.data.frame(our_cluster_means)

# Add cluster names as rownames
colnames(our_cluster_means_df) <- cls_2
############################################
## theirs
# get the count matrix
counts_data_their <- others_macaque_spn_seu@assays$RNA@counts[common_genes_1,]
# remove genes with zero counts
counts_data_filtered_their <- counts_data_their[rowSums(counts_data_their) > 0,]

# Create a new Seurat object with the filtered counts data
their_spn_seu_filtered <- CreateSeuratObject(counts = counts_data_filtered_their, meta.data = others_macaque_spn_seu@meta.data)

# Copy over other metadata 
their_spn_seu_filtered$newannot <- others_macaque_spn_seu$MSN_type
their_spn_seu_filtered$Tissue <- others_macaque_spn_seu$region_name

## get normalized counts:
# normalize data
their_spn_seu_filtered <- NormalizeData(their_spn_seu_filtered)
# get normalized data
their_normalized_data <- GetAssayData(their_spn_seu_filtered, slot = "data")

# Group by clusters
Idents(their_spn_seu_filtered) <- their_spn_seu_filtered$newannot
their_cluster_ids <- Idents(their_spn_seu_filtered)

# Calculate the mean expression for each gene across each cluster
# We'll calculate the row means per cluster for each gene
their_cls = unique(their_cluster_ids)
their_cluster_means <- sapply(their_cls, function(cluster) {
  cluster_cells <- which(their_cluster_ids == cluster)
  rowMeans(their_normalized_data[, cluster_cells], na.rm = TRUE)
})

# Convert the result to a data frame for easier viewing
their_cluster_means_df <- as.data.frame(their_cluster_means)

# Add cluster names as rownames
colnames(their_cluster_means_df) <- their_cls
##################################################################
# Match the gene names (rowname) orders
cmns = intersect(rownames(our_cluster_means_df), rownames(their_cluster_means_df))
our_cluster_means_df = our_cluster_means_df[cmns,]
their_cluster_means_df = their_cluster_means_df[cmns,]

# Convert to a matrix
cluster_means_matrix <- cbind(our_cluster_means_df, "D1_D2_hybrid" = their_cluster_means_df[,"D1/D2-Hybrid"])

# Calculate the correlation between clusters
correlation_matrix <- cor(cluster_means_matrix)

# adjust the color for the map:
# Check the minimum and maximum values of the correlation matrix
(min_value <- min(correlation_matrix, na.rm = TRUE)) # 0.8
(max_value <- max(correlation_matrix, na.rm = TRUE)) # 1

# Set custom color range from 0.7 to 1
color_range <- seq(0.7, 1, length.out = 101)

# Visualize the correlation matrix with a heatmap
library(pheatmap)
pdf("/home2/gkonop/workdir/03_INTEGRATE_ALL/SPN/He_etal_macaque_data/correlation_map_eSPN_D1D2_hybrid.pdf")
pheatmap(correlation_matrix, cluster_rows = TRUE, cluster_cols = TRUE, 
         show_rownames = TRUE, show_colnames = TRUE, 
         main = "Cluster-wise Gene Expression Correlation", breaks = color_range)
dev.off()
