# load packages
library(patchwork)
library(Seurat)
library(BPCells)
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
library(DESeq2)
library(ggrepel)
library(edgeR)
library(ashr)
library(variancePartition)
library(Matrix.utils)
set.seed(1234)
source("/project/Neuroinformatics_Core/Konopka_lab/s422071/SCRIPTS_pr/SCRIPTS/utility_functions.R")
library(curl)
#conda install bioconda::r-wgcna
#BiocManager::install('WGCNA')
library(WGCNA)
library(flashClust)
#BiocManager::install("UpSetR")
library(UpSetR)
#BiocManager::install("gprofiler2")
library(gprofiler2)
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(DT)
#devtools::install_github("krassowski/complex-upset")
library(ComplexUpset)
library(ComplexHeatmap)
library(purrr)


## load the data
species_list = c("human", "chimp", "macaque", "marmoset")

# load data
file_outdir = "/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/seu_objs/"
merged_pb_meta = readRDS(file = paste0(file_outdir,"primates_merged_pb_meta_PostgliaCleanup.RDS"))
merged_sub_seu_pb_count = readRDS(file = paste0(file_outdir,"primates_merged_sub_seu_pb_count_PostgliaCleanup.RDS"))
merged_meta_df = readRDS(file = paste0(file_outdir,"primates_merged_meta_df_PostgliaCleanup.RDS"))

# path to save files
file_dir = "/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/DEG_analysis/DESeq2/"

####
## Pseudobulk DESeq2 DEG proportions for celltypes
####
# Initialize output lists
all_species_results <- list()
top_10_all_species_results <- list()
sp_up <- list()
sp_down <- list()

#Tissue = "Putamen"
#cellType = "iSPN"

for(Tissue in names(merged_pb_meta)) {
	cellTypes <- setdiff(names(merged_pb_meta[[Tissue]]), "Non_SPN")
	for (cellType in cellTypes) {
	  species <- unique(merged_pb_meta[[Tissue]][[cellType]]$Species)

	  # List all relevant files for this cell type
	  files <- list.files(
	    file_dir,
	    pattern = paste0("^log10_sum_ncount_RNA_Humanized_age_covar_", Tissue, "_", cellType, ".*_all_genes\\.csv$"),
	    full.names = TRUE
	  )
	  
	  for (file in files) {
	    # Extract just the file name
	    base_file <- basename(file)

	    # Remove prefix and suffix to get species pair
	    species_pair <- sub(paste0("^log10_sum_ncount_RNA_Humanized_age_covar_", Tissue, "_", cellType, "_"), "", base_file)
	    species_pair <- sub("_all_genes\\.csv$", "", species_pair)

	    # Extract species names
	    pair <- strsplit(species_pair, "_vs_")[[1]]
	    if (length(pair) != 2) next

	    sp1 <- pair[1]
	    sp2 <- pair[2]

	    # Filter: Only process if the current species is part of the comparison
	    if (!(sp1 %in% species && sp2 %in% species)) next

	    # Read in result table
	    result_table <- read.table(file, header = TRUE, sep = ",")
	    
	    # Add to nested results list
	    if (is.null(all_species_results[[Tissue]][[cellType]])) {
	      all_species_results[[Tissue]][[cellType]] <- list()
	    }

	    pair_name <- paste(sp1, sp2, sep = "_vs_")
	    all_species_results[[Tissue]][[cellType]][[pair_name]] <- result_table
	    top_10_all_species_results[[Tissue]][[cellType]][[pair_name]] 

	    # Add up/downregulated genes
	    up_genes <- result_table[result_table$log2FoldChange > 0.6, "gene"]
	    down_genes <- result_table[result_table$log2FoldChange < -0.6, "gene"]

	    if (is.null(sp_up[[Tissue]][[cellType]])) sp_up[[Tissue]][[cellType]] <- list()
	    if (is.null(sp_down[[Tissue]][[cellType]])) sp_down[[Tissue]][[cellType]] <- list()

	    sp_up[[Tissue]][[cellType]][[pair_name]] <- up_genes
	    sp_down[[Tissue]][[cellType]][[pair_name]] <- down_genes
	  }
	}
}
# Save combined results
saveRDS(all_species_results, file = paste0(file_dir, "Primates_results_log10_sum_ncount_RNA_Humanized_age_covar_pseudobulk_DESeq2.RDS"))
all_species_results = readRDS(paste0(file_dir, "Primates_results_log10_sum_ncount_RNA_Humanized_age_covar_pseudobulk_DESeq2.RDS"))

### Fig 1: Plot of number of signif DEGs in each comparison
# Step 1: Extract the gene names for all combinations

# Create a named list where each key is a combination (like 'Caudate_Human_vs_Chimp')
gene_list <- list()

# Loop through the list to extract gene names
for (tissue in names(sp_up)) {
  for (cellType in names(sp_up[[tissue]])) {
    for (comparison in names(sp_up[[tissue]][[cellType]])) {
	
		# Extract gene names
		gene_names <- sp_up[[tissue]][[cellType]][[comparison]]
		combination_name <- paste(tissue, cellType,comparison, sep = "_")
		gene_list[[combination_name]] <- gene_names
		}
	}
}

# Check the gene list
head(gene_list)

# Step 2: Create the binary matrix
# Create a data frame of all genes
all_genes <- unique(unlist(gene_list))

# Create the binary matrix (1 if gene is in the set, 0 if not)
binary_matrix <- data.frame(
  gene = all_genes,
  matrix(unlist(lapply(gene_list, function(x) ifelse(all_genes %in% x, 1, 0))),
         ncol = length(gene_list),
         byrow = FALSE)
)

# Set column names for the binary matrix
colnames(binary_matrix)[-1] <- names(gene_list)
rownames(binary_matrix) <- binary_matrix$gene

# Check the binary matrix
head(binary_matrix)

### Caudate
# Convert to long format for ComplexUpset
caud_binary_matrix_long <- binary_matrix[, c("gene", grep("Caudate", colnames(binary_matrix), value = TRUE))]
rownames(caud_binary_matrix_long) <- caud_binary_matrix_long$gene
caud_binary_matrix_long$gene <- NULL


# Make sure all are 0/1 integers
caud_binary_matrix_long[] <- lapply(caud_binary_matrix_long, as.integer)

# Keep only rows that have at least one TRUE
caud_binary_matrix_long <- caud_binary_matrix_long[rowSums(caud_binary_matrix_long) > 0, ]

# Auto scale width
#pdf_width <- max(15, ncol(binary_matrix_long) * 0.5)

top_intersections <- colnames(caud_binary_matrix_long)


pdf(paste0(file_dir, "Primates_Caudate_upsetPlot_top_30_sets_complexUpset.pdf"), width = 30, height = 20)
upset(
  caud_binary_matrix_long,
  intersect = colnames(caud_binary_matrix_long),
  n_intersections = 18,
  min_size = 10,
  width_ratio = 0.08,
  sort_intersections = "descending",
  base_annotations = list(
    'Intersection size' = intersection_size(
      counts = TRUE,
      bar_number_threshold = 1
    )
  ),
  themes = upset_default_themes(
    text = element_text(size = 14)
  )
) +
theme(
  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, margin = margin(t = 5)),
  axis.text.y = element_text(size = 12, margin = margin(r = 5))
)

dev.off()

### Putamen
# Convert to long format for ComplexUpset
put_binary_matrix_long <- binary_matrix[, c("gene", grep("Putamen", colnames(binary_matrix), value = TRUE))]
rownames(put_binary_matrix_long) <- put_binary_matrix_long$gene
put_binary_matrix_long$gene <- NULL


# Make sure all are 0/1 integers
put_binary_matrix_long[] <- lapply(put_binary_matrix_long, as.integer)

# Keep only rows that have at least one TRUE
put_binary_matrix_long <- put_binary_matrix_long[rowSums(put_binary_matrix_long) > 0, ]

# Auto scale width
#pdf_width <- max(15, ncol(binary_matrix_long) * 0.5)

top_intersections <- colnames(put_binary_matrix_long)


pdf(paste0(file_dir, "Primates_Putamen_upsetPlot_top_30_sets_complexUpset.pdf"), width = 30, height = 20)
upset(
  put_binary_matrix_long,
  intersect = colnames(put_binary_matrix_long),
  n_intersections = 18,
  min_size = 10,
  width_ratio = 0.08,
  sort_intersections = "descending",
  base_annotations = list(
    'Intersection size' = intersection_size(
      counts = TRUE,
      bar_number_threshold = 1
    )
  ),
  themes = upset_default_themes(
    text = element_text(size = 14)
  )
) +
theme(
  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, margin = margin(t = 5)),
  axis.text.y = element_text(size = 12, margin = margin(r = 5))
)

dev.off()
##########################################
### Find the Human-specific DEGs
HCM_list <- list()
HC_list = list()
HM_list = list()
CM_list = list()
HMarm_list = list()
CMarm_list = list()
all_genes_tested = list()

#tissue = "Putamen"
#cellType = "iSPN"
# Loop through the list to extract gene names
for (tissue in names(all_species_results)) {
  for (cellType in names(all_species_results[[tissue]])) {
    for (comparison in names(all_species_results[[tissue]][[cellType]])) {
	if (comparison == "chimp_vs_Human") {
        prefix <- "HC_"
	HC_list[[tissue]][[cellType]] = all_species_results[[tissue]][[cellType]][[comparison]]
	colnames(HC_list[[tissue]][[cellType]]) <- 
        paste0(prefix, colnames(HC_list[[tissue]][[cellType]]))
      } else if (comparison == "Human_vs_Macaque") {
        prefix <- "HM_"
	HM_list[[tissue]][[cellType]] = all_species_results[[tissue]][[cellType]][[comparison]]
	colnames(HM_list[[tissue]][[cellType]]) <- 
        paste0(prefix, colnames(HM_list[[tissue]][[cellType]]))
        } else if (comparison == "Human_vs_Marmoset") {
        prefix <- "HMarm_"
	HMarm_list[[tissue]][[cellType]] = all_species_results[[tissue]][[cellType]][[comparison]]
	colnames(HMarm_list[[tissue]][[cellType]]) <- 
        paste0(prefix, colnames(HMarm_list[[tissue]][[cellType]]))
      } else if (comparison == "chimp_vs_Macaque") {
        prefix <- "CM_"
	CM_list[[tissue]][[cellType]] = all_species_results[[tissue]][[cellType]][[comparison]]
	colnames(CM_list[[tissue]][[cellType]]) <- 
        paste0(prefix, colnames(CM_list[[tissue]][[cellType]]))
      } else if (comparison == "chimp_vs_Marmoset") {
        prefix <- "CMarm_"
	CMarm_list[[tissue]][[cellType]] = all_species_results[[tissue]][[cellType]][[comparison]]
	colnames(CMarm_list[[tissue]][[cellType]]) <- 
        paste0(prefix, colnames(CMarm_list[[tissue]][[cellType]]))
      } else {
        prefix <- ""  # do nothing
      }
    #HCM_list[[tissue]][[cellType]] = cbind(HC_list[[tissue]][[cellType]], HM_list[[tissue]][[cellType]], CM_list[[tissue]][[cellType]])
    }
  }
}

### calculate the props of human changed DEGs over the human changed DEGs and chimp changed DEGs
props = list()
n_DEG_changes = list()
chi_list_1 = list()
human_up_gene_list = list()
human_down_gene_list = list()
all_genes_tested = list()

for (i in names(HC_list)) {
  for (j in names(HC_list[[i]])) {
	#i = "Putamen"
	#j = "Non_SPN"

	### Find Human and Chimpanzee specific up / down regulation
	## Human
	Human_down_1 = subset(HC_list[[i]][[j]], HC_padj < 0.05 & HC_log2FoldChange > 0.5)$HC_gene

	Human_down_2 = subset(HM_list[[i]][[j]], HM_padj < 0.05 & HM_log2FoldChange < -0.5)$HM_gene

	Human_down_3 = subset(HMarm_list[[i]][[j]], HMarm_padj < 0.05 & HMarm_log2FoldChange < -0.5)$HMarm_gene

	Human_down_4 = subset(CM_list[[i]][[j]],CM_padj > 0.1)$CM_gene

	Human_down = Reduce(intersect, list(Human_down_1,Human_down_2,Human_down_3,Human_down_4))
	
	stats = HC_list[[i]][[j]][HC_list[[i]][[j]]$HC_gene %in% Human_down,] # be careful the positive logFC means it's downregulated in Human because the DESeq2 comparison is chimp vs Human
	human_down_gene_list[[i]][[j]] = stats

	Human_up_1 = subset(HC_list[[i]][[j]], HC_padj < 0.05 & HC_log2FoldChange < -0.5)$HC_gene

	Human_up_2 = subset(HM_list[[i]][[j]], HM_padj < 0.05 & HM_log2FoldChange > 0.5)$HM_gene

	Human_up_3 = subset(HMarm_list[[i]][[j]], HMarm_padj < 0.05 & HMarm_log2FoldChange > 0.5)$HMarm_gene
	
	Human_up_4 = subset(CM_list[[i]][[j]],CM_padj > 0.1)$CM_gene

	Human_up = Reduce(intersect,list(Human_up_1,Human_up_2,Human_up_3,Human_up_4))

	stats = HC_list[[i]][[j]][HC_list[[i]][[j]]$HC_gene %in% Human_up,] # be careful the negative logFC means it's upregulated in Human because the DESeq2 comparison is chimp vs Human
	human_up_gene_list[[i]][[j]] = stats
# check human specifically upregulated genes in Putamen Non-SPNs 
#write.csv(Human_up, quote = F, file =paste0(file_dir, "Putamen_Non_SPN_Human_up.csv"))

	# combine the up and down regulated genes
	Human_change = union(Human_up, Human_down)
	
	all_genes_tested[[i]][[j]]  <- Reduce(union, list(HC_list[[i]][[j]]$HC_gene, 
				HM_list[[i]][[j]]$HM_gene, HMarm_list[[i]][[j]]$HMarm_gene))

	## proportions
	props[[i]][[j]] = length(Human_change) / length(all_genes_tested[[i]][[j]])

	n_DEG_changes[[j]][[i]]$Human_up <- length(Human_up)
	n_DEG_changes[[j]][[i]]$Human_down <- length(Human_down)
	n_DEG_changes[[j]][[i]]$Human_both_changes <- length(Human_change)
  }
}

# Create a helper function to flatten each nested list
flatten_gene_list <- function(gene_list, species, direction) {
  imap_dfr(gene_list, function(tissue_list, tissue_name) {
    imap_dfr(tissue_list, function(genes, celltype_name) {
      tibble(
        Species = species,
        Direction = direction,
        Tissue = tissue_name,
        CellType = celltype_name,
        Gene = genes
      )
    })
  })
}

# Flatten all 4 lists
df_human_up   <- flatten_gene_list(human_up_gene_list, "Human", "Up")
df_human_down <- flatten_gene_list(human_down_gene_list, "Human", "Down")

# Combine into one data.frame
all_genes_df <- bind_rows(df_human_up, df_human_down)

# Write to CSV
write.csv(all_genes_df, paste0(file_dir, "Human_specific_genes_logFC_0_5.csv"), row.names = FALSE, quote = T)
write.csv(df_human_up, paste0(file_dir, "Human_specific_Upgenes_logFC_0_5.csv"), row.names = FALSE, quote = T)
write.csv(df_human_down, paste0(file_dir, "Human_specific_Downgenes_logFC_0_5.csv"), row.names = FALSE, quote = T)

# save all genes tested per each CellType-Tissue
library(tidyverse)

all_genes_df <- map_dfr(
  names(all_genes_tested),
  function(region) {
    map_dfr(
      names(all_genes_tested[[region]]),
      function(cell_type) {
        tibble(
          region = region,
          cell_type = cell_type,
          gene = all_genes_tested[[region]][[cell_type]]
        )
      }
    )
  }
)

write.csv(
  all_genes_df,
  file = paste0(file_dir,"all_genes_tested_long.csv"),
  row.names = FALSE,
  quote = TRUE)


## test the # of HS DEG difference for SPN vs Glia in each tissue
# sum changes for all SPN types
(caud_SPN_human_changes_deg <- sum(
    n_DEG_changes[["iSPN"]][["Caudate"]]$Human_both_changes,
    n_DEG_changes[["eSPN"]][["Caudate"]]$Human_both_changes,
    n_DEG_changes[["dSPN"]][["Caudate"]]$Human_both_changes
))

(caud_SPN_total_genes_tested <- sum(
  sapply(
    list(all_genes_tested$Caudate[["dSPN"]],
         all_genes_tested$Caudate[["iSPN"]],
         all_genes_tested$Caudate[["eSPN"]]),
    length
  )
))

(put_SPN_human_changes_deg <- sum(
    n_DEG_changes[["iSPN"]][["Putamen"]]$Human_both_changes,
    n_DEG_changes[["eSPN"]][["Putamen"]]$Human_both_changes,
    n_DEG_changes[["dSPN"]][["Putamen"]]$Human_both_changes
))

(put_SPN_total_genes_tested <- sum(
  sapply(
    list(all_genes_tested$Putamen[["dSPN"]],
         all_genes_tested$Putamen[["iSPN"]],
         all_genes_tested$Putamen[["eSPN"]]),
    length
  )
))

chi_sq_list_2 = list()
# Now you can store the total SPN changes in the chi_sq_list_2
chi_sq_list_2[["SPN"]] <- list()
chi_sq_list_2[["SPN"]]$caud_SPN_human_changes_deg <- caud_SPN_human_changes_deg
chi_sq_list_2[["SPN"]]$put_SPN_human_changes_deg <- put_SPN_human_changes_deg
chi_sq_list_2[["SPN"]]$caud_SPN_total_genes_tested <- caud_SPN_total_genes_tested
chi_sq_list_2[["SPN"]]$put_SPN_total_genes_tested <- put_SPN_total_genes_tested
  
# sum changes for all Glia types
(caud_Glia_human_changes_deg <- sum(
    n_DEG_changes[["MOL"]][["Caudate"]]$Human_both_changes,
    n_DEG_changes[["OPC"]][["Caudate"]]$Human_both_changes,
    n_DEG_changes[["Astrocyte"]][["Caudate"]]$Human_both_changes,
    n_DEG_changes[["Microglia"]][["Caudate"]]$Human_both_changes
))

(caud_Glia_total_genes_tested <- sum(
  sapply(
    list(all_genes_tested$Caudate[["MOL"]],
         all_genes_tested$Caudate[["OPC"]],
         all_genes_tested$Caudate[["Astrocyte"]],
         all_genes_tested$Caudate[["Microglia"]]),
    length
  )
))

(put_Glia_human_changes_deg <- sum(
    n_DEG_changes[["MOL"]][["Putamen"]]$Human_both_changes,
    n_DEG_changes[["OPC"]][["Putamen"]]$Human_both_changes,
    n_DEG_changes[["Astrocyte"]][["Putamen"]]$Human_both_changes,
    n_DEG_changes[["Microglia"]][["Putamen"]]$Human_both_changes
))

(put_Glia_total_genes_tested <- sum(
  sapply(
    list(all_genes_tested$Putamen[["MOL"]],
         all_genes_tested$Putamen[["OPC"]],
         all_genes_tested$Putamen[["Astrocyte"]],
         all_genes_tested$Putamen[["Microglia"]]),
    length
  )
))

# Now you can store the total Glia changes in the chi_sq_list_2
chi_sq_list_2[["Glia"]] <- list()
chi_sq_list_2[["Glia"]]$caud_Glia_human_changes_deg <- caud_Glia_human_changes_deg
chi_sq_list_2[["Glia"]]$put_Glia_human_changes_deg <- put_Glia_human_changes_deg
chi_sq_list_2[["Glia"]]$caud_Glia_total_genes_tested <- caud_Glia_total_genes_tested
chi_sq_list_2[["Glia"]]$put_Glia_total_genes_tested <- put_Glia_total_genes_tested

# compare the change in ratios SPN vs Glia
caud_chisq_result <- prop.test(c(caud_SPN_human_changes_deg,caud_Glia_human_changes_deg), c(caud_SPN_total_genes_tested, caud_Glia_total_genes_tested))
    
put_chisq_result <- prop.test(c(put_SPN_human_changes_deg,put_Glia_human_changes_deg),
c(put_SPN_total_genes_tested, put_Glia_total_genes_tested))

# Store results
chi_sq_list_2 <- list(
    caud_SPN_Human_change = caud_SPN_human_changes_deg,
    put_SPN_Human_change = put_SPN_human_changes_deg,
    caud_SPN_Genes_tested = caud_SPN_total_genes_tested,
    put_SPN_Genes_tested = put_SPN_total_genes_tested,
    caud_glia_Human_change = caud_Glia_human_changes_deg,
    put_glia_Human_change  = put_Glia_human_changes_deg,
    caud_glia_Genes_tested = caud_Glia_total_genes_tested,
    put_glia_Genes_tested  = put_Glia_total_genes_tested,
    caud_chisq_result_change_pval = caud_chisq_result$p.value,
    put_chisq_result_change_pval = put_chisq_result$p.value
)

# save
saveRDS(chi_sq_list_2, file = paste0(file_dir, "SPN_vs_Glia_chi_sq_list_2_logFC_0_5.RDS"))

### plot comp for each tissue between SPNs & glia
# flatten the dataframe
(chi_sq_list_2_flat = cbind(
    data.frame(
      Tissue = factor(c("Caudate","Putamen")),
      SPN_Human_change = c(chi_sq_list_2$caud_SPN_Human_change,chi_sq_list_2$put_SPN_Human_change),
      Glia_Human_change = c(chi_sq_list_2$caud_glia_Human_change,chi_sq_list_2$put_glia_Human_change),
      SPN_Genes_tested = c(chi_sq_list_2$caud_SPN_Genes_tested,chi_sq_list_2$put_SPN_Genes_tested),
      Glia_Genes_tested = c(chi_sq_list_2$caud_glia_Genes_tested,chi_sq_list_2$put_glia_Genes_tested),
      SPN_props_H_over_AllGenesTested = c(c(chi_sq_list_2$caud_SPN_Human_change / chi_sq_list_2$caud_SPN_Genes_tested), c(chi_sq_list_2$put_SPN_Human_change / chi_sq_list_2$put_SPN_Genes_tested)),
      Glia_props_H_over_AllGenesTested = c(
  c(chi_sq_list_2$caud_glia_Human_change / chi_sq_list_2$caud_glia_Genes_tested),
  c(chi_sq_list_2$put_glia_Human_change / chi_sq_list_2$put_glia_Genes_tested)),
  chi_sq_pval = c(chi_sq_list_2$caud_chisq_result_change_pval, chi_sq_list_2$put_chisq_result_change_pval),
  FDR = p.adjust(c(chi_sq_list_2$caud_chisq_result_change_pval, chi_sq_list_2$put_chisq_result_change_pval), method =  "fdr")
)))


# Reshape the data into a long format
chi_sq_list_2_flat_long <- chi_sq_list_2_flat %>%
  pivot_longer(
    cols = matches("^(SPN|Glia)_"),
    names_to = c("CellType", ".value"),
    names_pattern = "(SPN|Glia)_(.*)"
  ) %>%
  arrange(Tissue, CellType)
                                       
# Plotting: Comparing Caudate and Putamen proportions side-by-side for each CellType
pdf(paste0(file_dir, "Str_HS_up_and_down_genes_AllGenesTested_ratios_comparing_SPN_vs_Glia.pdf"), height = 5, width = 14)
ggplot(chi_sq_list_2_flat_long, aes(x = Tissue, y = props_H_over_AllGenesTested, fill = CellType)) +
    # Custom fill colors for Caudate and Putamen
  geom_bar(stat = "identity", position = position_dodge(width = 1)) +
    # Display the p-value above the bars (no line, just text on top)
  geom_text(aes(label = sprintf("p=%.2e", FDR)),
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3.5, color = "black") +
  theme_classic() +
  labs(y = "HS / AllGenesTested DEG Ratio") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for clarity
  ) +   # Facet by CellType (side-by-side bars for each)
  facet_wrap(~Tissue, nrow = 1, scales = "free_x")
dev.off()


#######
## test the # of HS DEG difference for each celltype-tissue
chi_sq_list_2 = list()

# Determine total number of genes (assuming all DEG tables have same genes)
# Loop over cell types and tissues
for (ctype in names(n_DEG_changes)) {
  chi_sq_list_2[[ctype]] <- list()  # initialize list
  # Extract DEG counts for Human and total genes tested for DEG
  caud_human_changes_deg <- n_DEG_changes[[ctype]][["Caudate"]]$Human_both_changes
  caud_all_genes_tested_deg <- all_genes_tested[["Caudate"]][[ctype]]
  put_human_changes_deg <- n_DEG_changes[[ctype]][["Putamen"]]$Human_both_changes
  put_all_genes_tested_deg <- all_genes_tested[["Putamen"]][[ctype]]
  # compare the change in ratios across tissues
  chisq_result <- prop.test(c(caud_human_changes_deg,put_human_changes_deg), c(length(caud_all_genes_tested_deg), length(put_all_genes_tested_deg)))
  # Store results
  chi_sq_list_2[[ctype]] <- list(
    caud_Human_change = caud_human_changes_deg,
    put_Human_change = put_human_changes_deg,
    caud_Genes_tested = length(caud_all_genes_tested_deg),
    put_Genes_tested = length(put_all_genes_tested_deg),
    chisq_result_change_pval = chisq_result$p.value)
}

### plot comp for each celltype between caud & putamen
# flatten the dataframe
(chi_sq_list_2_flat = rbindlist(lapply(names(chi_sq_list_2), function(ctype) {
data.frame(
      CellType = ctype,
      caud_Human_change = chi_sq_list_2[[ctype]]$caud_Human_change,
      put_Human_change = chi_sq_list_2[[ctype]]$put_Human_change,
      caud_Genes_tested = chi_sq_list_2[[ctype]]$caud_Genes_tested,
      put_Genes_tested = chi_sq_list_2[[ctype]]$put_Genes_tested,
      chi_sq_change_pval = chi_sq_list_2[[ctype]]$chisq_result_change_pval,
      caud_props_H_over_AllGenesTested = chi_sq_list_2[[ctype]]$caud_Human_change / chi_sq_list_2[[ctype]]$caud_Genes_tested,
      put_props_H_over_AllGenesTested = chi_sq_list_2[[ctype]]$put_Human_change / chi_sq_list_2[[ctype]]$put_Genes_tested)}), fill = TRUE))

# Reshape the data into a long format
chi_sq_list_2_flat_long <- chi_sq_list_2_flat %>%
  pivot_longer(cols = c(caud_props_H_over_AllGenesTested, put_props_H_over_AllGenesTested),
               names_to = "Tissue", 
               values_to = "props_H_over_AllGenesTested") %>%
  mutate(Tissue = recode(Tissue, 
                        caud_props_H_over_AllGenesTested = "Caudate", 
                        put_props_H_over_AllGenesTested = "Putamen"))

# plot all celltypes wihtin a tissue next to each other
pdf(paste0(file_dir, "Str_HS_up_and_down_genes_AllGenesTested_ratios_comparing_CellTypes.pdf"), height = 5, width = 14)
ggplot(chi_sq_list_2_flat_long, aes(x = CellType, y = props_H_over_AllGenesTested, fill = CellType)) +
  geom_bar(stat = "identity", position = position_dodge(width = 1)) +
  theme_classic() +
  labs(y = "HS / AllGenesTested DEG Ratio") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for clarity
  ) +   # Facet by CellType (side-by-side bars for each)
  facet_wrap(~Tissue, nrow = 1, scales = "free_x")
dev.off()
######################
#### plot HS up/down
# prep the data
# Count the number of genes per category
deg_summary <- all_genes_df %>%
  group_by(Species, Direction, Tissue, CellType) %>%
  summarise(DEG_count = n(), .groups = "drop") %>%
  # Make "Down" counts negative so they plot to the left
  mutate(DEG_count_signed = ifelse(Direction == "Down", -DEG_count, DEG_count),
         Species_Direction = paste(Species, Direction, sep = "_"))

# re-order celltypes
 deg_summary$CellType = factor(deg_summary$CellType, levels = rev(c("dSPN", "iSPN", "eSPN", "Non_SPN", "Microglia", "Astrocyte", "MOL", "OPC")))
  
  
# Order Species_Direction for consistent grouping
deg_summary$Species_Direction <- factor(
  deg_summary$Species_Direction,
  levels = c("Human_Down", "Human_Up")
)

# Define fill colors: darker for ups, lighter for downs
species_fill_map <- c(
  "Human_Up"   = "#1f77b4",   # deep blue
  "Human_Down" = "#aec7e8"   # light blue
)

# --- Plot ---
p = ggplot(deg_summary, aes(x = DEG_count_signed, y = CellType, fill = Species_Direction)) +
  geom_col(position = position_dodge(width = 0), width = 0.7) +
  facet_wrap(~Tissue, nrow = 1, scales = "free_y") +
  geom_vline(xintercept = 0, color = "black", size = 0.7) +  # vertical zero line
  scale_fill_manual(values = species_fill_map) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    text = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  ) +
  labs(
    x = "Number of DEGs (? Down | Up ?)",
    y = "Cell Type",
    title = "Human-Specific DEGs across Cell Types and Tissues"
  )
# save
ggsave(plot = p, filename = paste0(file_dir, "HS_number_across_CellType_Tissue.pdf"),   width = 12,
  height = 8,
  dpi = 300
)
###################################################
#################
### GO term analysis - 
setwd("/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/DEG_analysis/DESeq2/GO_term/")
library(gprofiler2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

###check the ones that are HS_up, HS_down, CS_up,a nd CS_down
# 1. Create a flat list of comparisons
df_human_up = as.data.frame(df_human_up)
df_human_down = as.data.frame(df_human_down)

# Loop over cell types and tissues
for (Tissue in names(all_species_results)){
  for (CellType in names(all_species_results[[Tissue]])){
    all_genes <- all_genes_tested[[Tissue]][[CellType]]
  	
    HS_up_genes = df_human_up[df_human_up$Tissue == Tissue & df_human_up$CellType == CellType,]$Gene$HC_gene
    HS_down_genes = df_human_down[df_human_down$Tissue == Tissue & df_human_down$CellType == CellType,]$Gene$HC_gene

    # --- Run enrichment if nonempty
    if (length(HS_up_genes) > 0) {
      HS_ego_up <- GOenrich(gns = HS_up_genes, uni = all_genes)
      write.csv(HS_ego_up, paste0("Human_specific_UpGenes_log2FC_0_5_EnrichGO_", Tissue, "_", CellType, ".csv"), row.names = FALSE)
      message("Finished HS up: ", Tissue, "_", CellType)
    } else message("No HS UP genes for ", Tissue, "_", CellType)

    if (length(HS_down_genes) > 0) {
      HS_ego_down <- GOenrich(gns = HS_down_genes, uni = all_genes)
      write.csv(HS_ego_down, paste0("Human_specific_DownGenes_log2FC_0_5_EnrichGO_", Tissue, "_", CellType, ".csv"), row.names = FALSE)
      message("Finished HS down: ", Tissue, "_", CellType)
    } else message("No HS DOWN genes for ", Tissue, "_", CellType)
}
}

# Gather all CSV files
files <- list.files(path = paste0(file_dir, "GO_term/"), pattern =  "EnrichGO_.*\\.csv$", full.names = TRUE)

GO_up = list()
GO_down = list()
GO_up$category = list(
			Human_Putamen_MOL = c("glutamatergic synapse", "regulation of neuron projection", "dendrite development"), 
			Human_Caudate_MOL = c(""),
			Chimp_Putamen_MOL = c("asymmetric synapse", "negative regulation of cell cycle process", "catalytic activity, acting on DNA"),
			Chimp_Caudate_MOL = c("RNA splicing", "postsynaptic specialization", "neuron to neuron synapse"))
GO_down$category = list(
			Human_Putamen_MOL = c("rRNA metabolic process", "helicase activity", "nuclear speck"), 
			Human_Caudate_MOL = c(""),
			Chimp_Putamen_MOL = c("cellular response to epinephrine stimulus", "lipid modification"),
			Chimp_Caudate_MOL = c("regulation of neuron projection development", "histone modification", "methyltransferase complex"))


# Extract condition info from filenames
df_all <- df_all %>%
  mutate(
    source = gsub("EnrichGO_|\\.csv", "", basename(files[as.numeric(source)])),
    Condition = gsub("_Genes_.*", "", source),
    Direction = ifelse(grepl("UpGenes", source), "Up", "Down"),
    Species = ifelse(grepl("Human", source), "Human", "Chimp")
  )

# Choose a significance metric (e.g., -log10(p.adjust))
heat_df <- df_all %>%
  filter(p.adjust < 0.05) %>%
  group_by(Species, Direction, ID, Description) %>%
  summarize(score = -log10(min(p.adjust)), .groups = "drop")

# Spread to wide format
heat_mat <- heat_df %>%
  pivot_wider(
    names_from = c(Species, Direction),
    values_from = score,
    values_fill = 0
  ) %>%
  column_to_rownames("Description")

# Heatmap
pheatmap(
  heat_mat,
  scale = "row",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize = 9,
  fontsize_row = 7,
  main = "GO term enrichment (Human vs. Chimp Up/Down)",
  border_color = NA
)

### dotplot of GO-terms

# Top 5 terms per comparison
top_go <- go_df %>%
  group_by(Tissue, CellType, Comparison) %>%
  slice_min(order_by = p.adjust, n = 5, with_ties = FALSE) %>% 
  ungroup() %>%
  mutate(
    Description_short = str_trunc(Description, width = 50, side = "right"),
    Tissue_CellType = paste(Tissue, CellType, sep = "_")
  )

top_go <- go_df %>%
  mutate(
    Description_short = str_trunc(Description, width = 50, side = "right"),
    Tissue_CellType = paste(Tissue, CellType, sep = "_")
  ) %>%
  group_by(Tissue, CellType, Comparison, Description_short) %>%
  slice_min(order_by = p.adjust, n = 1) %>%  # keep only 1 term per truncated label
  ungroup() %>%
  group_by(Tissue, CellType, Comparison) %>%
  slice_min(order_by = p.adjust, n = 5, with_ties = FALSE) %>%
  ungroup()

# Plot
p <- ggplot(top_go, aes(x = -log10(p.adjust), y = fct_rev(Description_short))) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  scale_color_gradient(low = "skyblue", high = "darkblue") +
  scale_size(range = c(2, 6)) +
  facet_grid(Tissue_CellType ~ Comparison, scales = "free_y", space = "free") +
  labs(
    x = expression(-log[10](FDR)),
    y = NULL,
    size = "Gene Count",
    color = expression(-log[10](FDR)),
    title = "Significant GO Terms per Tissue/CellType/Comparison"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 8, face = "bold"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    panel.spacing = unit(0.3, "lines")
  )

# Save compact figure
ggsave(
  filename = paste0(file_dir, "GO_term/GO_BP_dotplot_compact.pdf"),
  plot = p,
  width = 10,
  height = 16
)


