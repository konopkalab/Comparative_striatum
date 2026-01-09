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
library(DESeq2)
library(SingleCellExperiment)
library(apeglm)
library(ashr)
library(edgeR)
library(glmGamPoi)
library(ggrepel)

source('~/SCRIPTS/utility_functions.R')
set.seed(1234)
################################################################
bat_filt <- readRDS("/home2/gkonop/projectShared/For_Gozde/01_BAT/bat_filt.RDS")
###########################################################################
#### find correlation of gene expression profiles across interneurons
# subset Putamen
sub <- subset(bat_filt, Tissue == "Putamen")
# subset the interneurons
sub = subset(sub, subset = newannot3 %in% c("MOL", "OPC", "COP", "Astrocyte", "Microglia", "SPN"), invert = T)

# set identity
Idents(sub) <- sub$newannot3

# remove lowly abundant genes to improve comparison
dim(sub@assays$RNA@counts)#[1] 27343  2616
sum(rowSums(sub@assays$RNA@counts) > 0) #[1] 19606

# get counts
counts_data <- sub@assays$RNA@counts
# remove genes with zero counts
counts_data_filtered <- counts_data[rowSums(counts_data) > 0,]

# Create a new Seurat object with the filtered counts data
sub_filtered <- CreateSeuratObject(counts = counts_data_filtered, meta.data = sub@meta.data)

# Copy over other metadata 
sub_filtered$newannot3 <- sub$newannot3
sub_filtered$Tissue <- sub$Tissue

## get normalized counts:
# normalize data
sub_filtered <- NormalizeData(sub_filtered)
# get normalized data
normalized_data <- GetAssayData(sub_filtered, slot = "data")

# Group by clusters
Idents(sub_filtered) <- sub_filtered$newannot3
cluster_ids <- Idents(sub_filtered)

# Calculate the mean expression for each gene across each cluster
# We'll calculate the row means per cluster for each gene
cluster_means <- sapply(unique(cluster_ids), function(cluster) {
  cluster_cells <- which(cluster_ids == cluster)
  rowMeans(normalized_data[, cluster_cells], na.rm = TRUE)
})

# Convert the result to a data frame for easier viewing
cluster_means_df <- as.data.frame(cluster_means)

# Add cluster names as rownames
colnames(cluster_means_df) <- unique(cluster_ids)

# Convert to a matrix
cluster_means_matrix <- do.call(cbind, cluster_means_df)

# Calculate the correlation between clusters
correlation_matrix <- cor(cluster_means_matrix)

# Visualize the correlation matrix with a heatmap
library(pheatmap)
pdf("correlation_map_InterneuronswithinBat.pdf")
pheatmap(correlation_matrix, cluster_rows = TRUE, cluster_cols = TRUE, 
         show_rownames = TRUE, show_colnames = TRUE, 
         main = "Cluster-wise Gene Expression Correlation")
dev.off()

correlation_matrix["LMO3_BAT", "FOXP2_EYA2"]
correlation_matrix["LMO3_BAT", "FOXP2_TSHZ2_BAT"]
correlation_matrix["FOXP2_EYA2", "TH"]
correlation_matrix["FOXP2_TSHZ2_BAT", "TAC3"]
[1] 0.9578292
[1] 0.9276468
[1] 0.94205
[1] 0.8972202
correlation_matrix[c("LMO3_BAT","FOXP2_EYA2","FOXP2_TSHZ2_BAT"),
        c("TH","TAC3")]
                       TH      TAC3
LMO3_BAT        0.9258643 0.9142284
FOXP2_EYA2      0.9420500 0.9248304
FOXP2_TSHZ2_BAT 0.9089461 0.8972202
