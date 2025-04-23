# load packages
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
library(DESeq2)
library(edgeR)
library(EnsDb.Hsapiens.v86)
source('~/SCRIPTS/utility_functions.R')
set.seed(1234)

# read data (SPN seurat object)
filtered_spn <- readRDS("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/ANNOTATED/GB_AllTissues_Primates_bat_mouse_ferret_filtered_Integrated_SPN_ncbi_human_pr_coding_orthologs_ANNOTATED.RDS")

ortho_genes <- read.table("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/orthologs_in_7_species_human_pr_codingOnly_2.csv",  header = TRUE)$symbol #15055 pr coding ortho genes
####################################################################
####
## LOAD NEURON DATASETS
####

###extract the count matrices with orthologs for each Species
# human
human <- readRDS("/home2/gkonop/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN/ANNOTATED/human_integrated_cuadate_putamen_ANNOTATED.RDS")
# extract_metadata
human_meta <- human[[]]
# extract count matrix
human_mat <- human@assays$RNA@counts
# subset count matrix with orthologs
human_mat_orth <- human_mat[rownames(human_mat) %in% ortho_genes,]
# create new seurat object
human_orth <- CreateSeuratObject(counts = human_mat_orth, meta.data = human_meta, assay = "RNA")
# Normalize the counts and set the identity in new seurat obj
DefaultAssay(human_orth) <- "RNA"
human_orth <- NormalizeData(human_orth)
Idents(human_orth) <- human$newannot


# chimp
chimp <- readRDS("/home2/gkonop/workdir/01_CHIMP/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN/ANNOTATED/chimp_integrated_cuadate_putamen_ANNOTATED.RDS")
chimp$Species <- "Chimp"
# extract_metadata
chimp_meta <- chimp[[]]

# extract count matrix
chimp_mat <- chimp@assays$RNA@counts

# subset count matrix with orthologs
chimp_mat_orth <- chimp_mat[rownames(chimp_mat) %in% ortho_genes,]

# create new seurat object
chimp_orth <- CreateSeuratObject(counts = chimp_mat_orth, meta.data = chimp_meta, assay = "RNA")

# Normalize the counts and set the identity in new seurat obj
DefaultAssay(chimp_orth) <- "RNA"
chimp_orth <- NormalizeData(chimp_orth)
Idents(chimp_orth) <- chimp_orth$newannot


# macaque
macaque <- readRDS("/home2/gkonop/workdir/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/Caudate_Putamen/ANNOTATED/macaque_integrated_cuadate_putamen_ANNOTATED.RDS")
# extract_metadata
macaque_meta <- macaque[[]]

# extract count matrix
macaque_mat <- macaque@assays$RNA@counts

# subset count matrix with orthologs
macaque_mat_orth <- macaque_mat[rownames(macaque_mat) %in% ortho_genes,]

# create new seurat object
macaque_orth <- CreateSeuratObject(counts = macaque_mat_orth, meta.data = macaque_meta, assay = "RNA")

# Normalize the counts and set the identity in new seurat obj
DefaultAssay(macaque_orth) <- "RNA"
macaque_orth <- NormalizeData(macaque_orth)
Idents(macaque_orth) <- macaque_orth$newannot


# marmoset
marmoset <- readRDS("/home2/gkonop/workdir/01_MARMOSET/MARMOSET_INTEGRATION/Caudate_putamen/ANNOTATED/marmoset_integrated_cuadate_putamen_ANNOTATED.RDS")
# extract_metadata
marmoset_meta <- marmoset[[]]

# extract count matrix
marmoset_mat <- marmoset@assays$RNA@counts

# subset count matrix with orthologs
marmoset_mat_orth <- marmoset_mat[rownames(marmoset_mat) %in% ortho_genes,]

# create new seurat object
marmoset_orth <- CreateSeuratObject(counts = marmoset_mat_orth, meta.data = marmoset_meta, assay = "RNA")

# Normalize the counts and set the identity in new seurat obj
DefaultAssay(marmoset_orth) <- "RNA"
marmoset_orth <- NormalizeData(marmoset_orth)
Idents(marmoset_orth) <- marmoset_orth$newannot

# mouse
mouse_caud_put <- readRDS("/home2/gkonop/workdir/01_MOUSE/03_CLUSTER_ANNOTATE/Mouse_Caudate_Annotated_FINAL.RDS")
mouse_caud_put$Species <- "Mouse"
table(mouse_caud_put$newannot)
mouse_caud_put$Tissue <- "Caudoputamen"
# extract_metadata
mouse_meta <- mouse_caud_put[[]]

# extract count matrix
mouse_mat <- mouse_caud_put@assays$RNA@counts
# convert gene names to upper case
rownames(mouse_mat) <- toupper(rownames(mouse_mat))

# subset count matrix with orthologs
mouse_mat_orth <- mouse_mat[rownames(mouse_mat) %in% ortho_genes,]

# create new seurat object
mouse_orth <- CreateSeuratObject(counts = mouse_mat_orth, meta.data = mouse_meta, assay = "RNA")

# Normalize the counts and set the identity in new seurat obj
DefaultAssay(mouse_orth) <- "RNA"
mouse_orth <- NormalizeData(mouse_orth)
Idents(mouse_orth) <- mouse_orth$newannot


# bat
bat <- readRDS("/home2/gkonop/workdir/01_BAT/03_CLUSTER_ANNOTATE/Caudate_Putamen/ANNOTATED/bat_integrated_cuadate_putamen_ANNOTATED.RDS")
bat_spns <-subset(bat, subset = newannot == "SPN")
# extract_metadata
bat_meta <- bat[[]]

# extract count matrix
bat_mat <- bat@assays$RNA@counts

# subset count matrix with orthologs
bat_mat_orth <- bat_mat[rownames(bat_mat) %in% ortho_genes,]

# create new seurat object
bat_orth <- CreateSeuratObject(counts = bat_mat_orth, meta.data = bat_meta, assay = "RNA")

# Normalize the counts and set the identity in new seurat obj
DefaultAssay(bat_orth) <- "RNA"
bat_orth <- NormalizeData(bat_orth)
Idents(bat_orth) <- bat_orth$newannot

# ferret
ferret <- readRDS("/home2/gkonop/workdir/01_FERRET/03_CLUSTER_ANNOTATE/ANNOTATED/Ferret_Caudate_Krienen_ANNOTATED.RDS")
# create "id" metadata
ferret$id <- ferret$orig.ident
ferret$newannot <- ferret$broad_annot
# extract_metadata
ferret_meta <- ferret[[]]

# extract count matrix
ferret_mat <- ferret@assays$RNA@counts

# subset count matrix with orthologs
ferret_mat_orth <- ferret_mat[rownames(ferret_mat) %in% ortho_genes,]

# create new seurat object
ferret_orth <- CreateSeuratObject(counts = ferret_mat_orth, meta.data = ferret_meta, assay = "RNA")

# Normalize the counts and set the identity in new seurat obj
DefaultAssay(ferret_orth) <- "RNA"
ferret_orth <- NormalizeData(ferret_orth)
Idents(ferret_orth) <- ferret_orth$newannot
####################################################################
### eSPN  markers
## extract eSPN cells for each Species
# change the identity to iSPN,dSPN,eSPN
Idents(filtered_spn) <- filtered_spn$newannot_2

# Initialize an empty list to store the results for each species
results_list <- list()
results_neuron_list <- list()
species_orth_list <-  list()

# List of species names to process
species_list <- c("Human", "Chimp", "Macaque", "Marmoset", "Mouse", "Bat", "Ferret")
# Create a list with species names as keys
species_orth_list <- list(
  "Human" = human_orth,
  "Chimp" = chimp_orth,
  "Macaque" = macaque_orth,
  "Marmoset" = marmoset_orth,
  "Mouse" = mouse_orth,
  "Bat" = bat_orth,
  "Ferret" = ferret_orth
)

# Iterate over each species
for (species in species_list) {
  
  # Step 1: Extract cells for the current species
  species_spns <- subset(filtered_spn, subset = Species == species)
  
  # Step 2: Extract cell barcodes and create a new column for cellbarc_id
  species_spns$cellbarc <- str_split_i(colnames(species_spns@assays$RNA@counts), "-", 1)
  
  # Step 3: Subset to get the eSPN cells
  species_espns <- subset(species_spns, subset = newannot_2 == "eSPN")
  
  # Step 4: Create the cellbarc_id by combining cell barcode and sample information for eSPNs
  species_espns$cellbarc_id <- paste0(str_split_i(colnames(species_espns@assays$RNA@counts), "-", 1), "_", str_split_i(colnames(species_espns@assays$RNA@counts), "-", 2))
  
  # Step 5: Process the original seurat object similarly
  species_orth <- species_orth_list[[species]]
  species_orth$cellbarc_id <- paste0(str_split_i(colnames(species_orth@assays$RNA@counts), "-", 1), "_", str_split_i(colnames(species_orth@assays$RNA@counts), "-", 2))
  
  # Step 6: Assign eSPN labels based on cellbarc_id matches
  species_orth$newannot <- ifelse(names(species_orth$cellbarc_id) %in% names(species_espns$cellbarc_id), 'eSPN', as.character(species_orth$newannot))
  Idents(species_orth) <- species_orth$newannot
  
	  ## Step 7-1: Find markers between eSPN and non-eSPN
	  species_espn_markers <- FindMarkers(species_orth, ident.1 = "eSPN", only.pos = T, logfc.threshold = 0.5)
	  species_espn_markers <- species_espn_markers[species_espn_markers$p_val_adj <= 0.05,]
	  
	  # Add a new column for pct.1 / pct.2
	  species_espn_markers$pct_ratio <- species_espn_markers$pct.1 / species_espn_markers$pct.2
	  
	  # Order by pct_ratio and p_val_adj
	  species_espn_markers_ordered <- species_espn_markers %>% arrange(desc(pct_ratio), p_val_adj)
	  
	  # Add a new column indicating the species
	  species_espn_markers_ordered$Species <- species
	  
	  # Append the current species markers to the results list
	  results_list[[species]] <- species_espn_markers_ordered
	  
	  ## Step 7-2: Find markers between eSPN and other neurons
	  species_espn_neuron_markers <- FindMarkers(species_orth, ident.1 = "eSPN", ident.2 = c("SPN", "Non_SPN"), only.pos = T, logfc.threshold = 0.5)
	  species_espn_neuron_markers <- species_espn_neuron_markers[species_espn_neuron_markers$p_val_adj <= 0.05,]
	  
	  # Add a new column for pct.1 / pct.2
	  species_espn_neuron_markers$pct_ratio <- species_espn_neuron_markers$pct.1 / species_espn_neuron_markers$pct.2
	  
	  # Order by pct_ratio and p_val_adj
	  species_espn_neuron_markers_ordered <- species_espn_neuron_markers %>% arrange(desc(pct_ratio), p_val_adj)
	  
	  # Add a new column indicating the species
	  species_espn_neuron_markers_ordered$Species <- species
	  
	  # Append the current species markers to the results list
	  results_neuron_list[[species]] <- species_espn_neuron_markers_ordered

	# Step 8: Directly update the species_orth object in the species_orth_list
  	species_orth_list[[species]] <- species_orth
}

# Combine all species' results into one dataframe
final_results <- do.call(rbind, results_list)
final_results$gene_name <- str_split_i(rownames(final_results), "\\.", 2)

# save
write.csv(final_results, file = "/home2/gkonop/workdir/03_INTEGRATE_ALL/SPN/ANNOTATED/eSPN_DEGs.csv", quote = F)
saveRDS(final_results, "/home2/gkonop/workdir/03_INTEGRATE_ALL/SPN/ANNOTATED/eSPN_DEGs.RDS")
final_results <- readRDS("/home2/gkonop/workdir/03_INTEGRATE_ALL/SPN/ANNOTATED/eSPN_DEGs.RDS")

# Combine all species' results into one dataframe
final_results_neurons <- do.call(rbind, results_neuron_list)
final_results_neurons$gene_name <- str_split_i(rownames(final_results_neurons), "\\.", 2)

# save
write.csv(final_results_neurons, file = "/home2/gkonop/workdir/03_INTEGRATE_ALL/SPN/ANNOTATED/eSPN_vs_neurons_DEGs.csv", quote = F)
saveRDS(final_results_neurons, "/home2/gkonop/workdir/03_INTEGRATE_ALL/SPN/ANNOTATED/eSPN_vs_neurons_DEGs.RDS")
final_results_neurons <- readRDS("/home2/gkonop/workdir/03_INTEGRATE_ALL/SPN/ANNOTATED/eSPN_vs_neurons_DEGs.RDS")

# save the updated species_ortho
saveRDS(species_orth_list, "/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/ANNOTATED/species_orth_list.RDS")
species_orth_list <- readRDS("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/ANNOTATED/species_orth_list.RDS")
##################################################################################
## find the eSPN DEGs across species
# Load the necessary library
#BiocManager::install("UpSetR")
library(UpSetR)
library(patchwork)

# Extract the gene names for positive markers from each species (inetersection of eSPN vs neuron DEGs and eSPN vs all cells DEGs
# Human
human_genes <- intersect(final_results_neurons[final_results_neurons$Species =="Human",]$gene_name, final_results[final_results$Species =="Human",]$gene_name)
# Chimp
chimp_genes <- intersect(final_results_neurons[final_results_neurons$Species =="Chimp",]$gene_name, final_results[final_results$Species =="Chimp",]$gene_name)
# Macaque
macaque_genes <- intersect(final_results_neurons[final_results_neurons$Species =="Macaque",]$gene_name, final_results[final_results$Species =="Macaque",]$gene_name)
# Marmoset
marmoset_genes <- intersect(final_results_neurons[final_results_neurons$Species == "Marmoset", ]$gene_name, 
                            final_results[final_results$Species == "Marmoset", ]$gene_name)

# Mouse
mouse_genes <- intersect(final_results_neurons[final_results_neurons$Species == "Mouse", ]$gene_name, 
                         final_results[final_results$Species == "Mouse", ]$gene_name)

# Bat
bat_genes <- intersect(final_results_neurons[final_results_neurons$Species == "Bat", ]$gene_name, 
                       final_results[final_results$Species == "Bat", ]$gene_name)

# Ferret
ferret_genes <- intersect(final_results_neurons[final_results_neurons$Species == "Ferret", ]$gene_name, 
                          final_results[final_results$Species == "Ferret", ]$gene_name)

genes_list <- list("Human" =  human_genes,
		   "Chimp" = chimp_genes,
		   "Macaque" = macaque_genes,
		   "Marmoset" = marmoset_genes,
		   "Mouse" = mouse_genes,
		   "Bat" = bat_genes,
		   "Ferret" = ferret_genes)

dotout_combined <- list()
dotout_filtered_list <- list()
dotout_id_list <- list()

# Iterate over each species
for (species in species_list) {
	DefaultAssay(species_orth_list[[species]]) <- "RNA"
	species_orth_list[[species]] <- NormalizeData(species_orth_list[[species]])
	dotout = DotPlot(species_orth_list[[species]], features = genes_list[[species]], group.by = 'newannot')
	dotout_data = dotout$data
	dotout_data$Species = species
        
	# filter the genes with at lest average scaled expression of  1 for eSPNs
	keep_genes <- dotout_data[dotout_data$id == "eSPN" & dotout_data$avg.exp.scaled > 1, ]$features.plot
	keep_genes <- droplevels(keep_genes)
	# split the list accroding to the annotations (id)
	split_dotout_data <- split(dotout_data, dotout_data$id)
	
	# extract the keep_genes from each dataframe
	for (id in dotout_data$id) {
	dotout_filtered <- split_dotout_data[[id]][sort(split_dotout_data[[id]]$features.plot) %in% sort(keep_genes), ]
	dotout_id_list[[id]] <- dotout_filtered
	}

	# convert the list to dataframe
	final_filtered_resuls <- do.call(rbind, dotout_id_list)
	final_filtered_resuls$id <- paste0(final_filtered_resuls$id, "_",species)
        # save
        dotout_combined[[species]] <- final_filtered_resuls
		
}
a <- bind_rows(dotout_combined)

# Extract and clean the 'CellType' based on the pattern in 'id'
a$CellType = gsub('_Human|_Chimp|_Macaque|_Marmoset|_Mouse|_Bat|_Ferret', '', a$id)

# remove COPs
a <-  a[a$CellType != "COP",]
a$CellType = factor(a$CellType, levels = c("eSPN", "SPN", "Non_SPN","MOL", "OPC", "Astrocyte", "Microglia"))
a$Species <- factor(a$Species, levels = c("Human", "Chimp", "Macaque","Marmoset", "Mouse","Bat", "Ferret"))
################################################################
# Create a list to represent the sets
a_list <- split(a, a$Species)

# Find all unique genes across all species
all_genes <- unique(unlist(lapply(a_list, function(species_data) species_data$features.plot)))


# Create a binary presence/absence matrix for each species
gene_presence <- lapply(a_list, function(species_data) {
  # Create a binary vector of presence (1) or absence (0) for each gene across all species
  gene_presence_vector <- as.integer(all_genes %in% species_data$features.plot)
  names(gene_presence_vector) <- all_genes  # Name the vector with all unique genes
  return(gene_presence_vector)
})

# Bind all species' gene presence information into one data frame
upset_data <- as.data.frame(gene_presence)
upset_data[is.na(upset_data)] <- 0  # Convert NA to 0 (absence of gene)
# Convert colnames(upset_data) to a factor with specific levels
colnames(upset_data) <- factor(colnames(upset_data), levels = c('Human', 'Chimp', 'Macaque', 'Marmoset', 'Mouse', 'Ferret', 'Bat'))


#all tissues
pdf("/home2/gkonop/workdir/03_INTEGRATE_ALL/SPN/ANNOTATED/UpsetPlot_Species_ordered_eSPN_pos_markers_log2FC0_5_NCBI_human_pr_coding_orthologs_eSPN_vsAllCells_intersectedwith_eSPN_vsOtherNeurons.pdf", width = 23, height = 15)
upset(upset_data, 
      sets = colnames(upset_data),  # Use the factor levels for sets
      keep.order = TRUE, order.by = "freq",             # Sort by frequency of intersections
      main.bar.color = "deepskyblue2",    # Color for the bars
      sets.bar.color = "orange",     # Color for the sets
      matrix.color = "springgreen3",   # Color for the matrix
      point.size = 3.5, line.size = 1, 
      set_size.show = TRUE, 
      text.scale = c(4, 3, 3, 3, 3, 3))  
dev.off()

# 9 genes are eSPN DEGs within all of the 7 species
# extract and perform gene ontology search
(common_espn_degs <- rownames(upset_data)[apply(upset_data, 1, function(row) all(row == 1))])
# [1] "CASZ1"   "FOXP2"   "TSHZ1"   "RUNX1T1"  "DCLK1" 
##################################################################################
# Combine all data into one dataframe for combined plot
# markers
select_markers = c("CASZ1",   "FOXP2",   "TSHZ1",   "DCLK1", "RUNX1T1")

dotout_combined <- list()

# Iterate over each species
for (species in species_list) {
	DefaultAssay(species_orth_list[[species]]) <- "RNA"
	species_orth_list[[species]] <- NormalizeData(species_orth_list[[species]])
	dotout = DotPlot(species_orth_list[[species]], features = select_markers, group.by = 'newannot')
	dotout_data = dotout$data
	dotout_data$Species = species
        
	# save
        dotout_combined[[species]] <- dotout_data
		
}

# Assuming dotout_combined is a list where each species' data is a data frame
dotout_combined_df <- do.call(rbind, dotout_combined)

# Extract and clean the 'CellType' based on the pattern in 'id'
dotout_combined_df$CellType = gsub('_Human|_Chimp|_Macaque|_Marmoset|_Mouse|_Bat|_Ferret', '', dotout_combined_df$id)
# remove COPs
dotout_combined_df <-  dotout_combined_df[dotout_combined_df$CellType != "COP",]
dotout_combined_df$CellType = factor(dotout_combined_df$CellType, levels = c("eSPN", "SPN", "Non_SPN","MOL", "OPC", "Astrocyte", "Microglia"))
dotout_combined_df$Species <- factor(dotout_combined_df$Species, levels = c("Human", "Chimp", "Macaque","Marmoset", "Mouse","Bat", "Ferret"))


pdf(paste0("eSPN_markers_forAllSpecies_AllTissues.pdf"), width = 30, height = 10)
ggscatter(dotout_combined_df, y = 'features.plot', x = 'CellType', color = 'avg.exp.scaled',
			size = 'pct.exp') +
scale_size(range = c(-0.5,6)) +
xlab('') +
ylab('') +
theme(axis.text.x=element_text(size=20, face = 'bold'),
	axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
theme(text=element_text(size=20, face = 'bold')) +
scale_color_gradient2(low = 'blue', high = 'red', midpoint = 0) +
facet_wrap(~Species, nrow = 1) +
rotate_x_text(45)
dev.off()
