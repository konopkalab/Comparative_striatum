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
library(cowplot)
source('~/SCRIPTS/utility_functions.R')
set.seed(1234)

# set variables
marksToPlot = c('SLC17A7', 'SATB2', 'RBFOX3', 'NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA')

### load the data
human <- readRDS("/home2/gkonop/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN/ANNOTATED/human_integrated_cuadate_putamen_ANNOTATED.RDS")
# featureplot
#pdf("/home2/gkonop/workdir/06_QC_plots/human_FEATUREPLOT_without_COP.pdf", width = 15, height = 35)
#FeaturePlot(human, features = marksToPlot, sort = T, raster = T, pt.size = 2)
#dev.off()

# plot Intronic read ratios across samples
#pdf("/home2/gkonop/workdir/06_QC_plots/human_NuclearFraction_across_samples.pdf")
#ggboxplot(human[[]], x = 'orig.ident', y = 'intronRat') +
#rotate_x_text(90) +
#theme(text=element_text(size=20, face = 'bold'))
#dev.off()


chimp <- readRDS("/home2/gkonop/workdir/01_CHIMP/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN/ANNOTATED/chimp_integrated_cuadate_putamen_ANNOTATED.RDS")
# place COP cells unders OPCs for ease of plotting
a <- gsub("COP", "OPC", chimp$newannot) 
chimp$newannot <- a
# replot the UMAP
#pdf("/home2/gkonop/workdir/01_CHIMP/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN/ANNOTATED/chimp_caudate_putamen_UMAP_without_COPs.pdf")
#DimPlot(chimp, group.by = 'newannot', raster = T, label = T)
#dev.off()

# replot featureplot
#pdf("/home2/gkonop/workdir/06_QC_plots/chimp_FEATUREPLOT_without_COP.pdf", width = 15, height = 35)
#FeaturePlot(chimp, features = marksToPlot, sort = T, raster = T, pt.size = 2)
#dev.off()

# replot nuclear fractions
#pdf("/home2/gkonop/workdir/06_QC_plots/chimp_wihout_COPs_NuclearFraction.pdf")
#ggboxplot(chimp[[]], x = 'newannot', y = 'intronRat') +
#rotate_x_text(90) +
#theme(text=element_text(size=20, face = 'bold'))
#dev.off()


# plot Intronic read ratios across samples
#pdf("/home2/gkonop/workdir/06_QC_plots/chimp_NuclearFraction_across_samples.pdf")
#ggboxplot(chimp[[]], x = 'orig.ident', y = 'intronRat') +
#rotate_x_text(90) +
#theme(text=element_text(size=20, face = 'bold'))
#dev.off()

macaque <- readRDS("/home2/gkonop/workdir/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/Caudate_Putamen/ANNOTATED/macaque_integrated_cuadate_putamen_ANNOTATED.RDS")
# place COP cells unders OPCs for ease of plotting
a <- gsub("COP", "OPC", macaque$newannot) 
macaque$newannot <- a
# replot the UMAP
#pdf("/home2/gkonop/workdir/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/Caudate_Putamen/ANNOTATED/macaque_caudate_putamen_UMAP_without_COPs.pdf")
#DimPlot(macaque, group.by = 'newannot', raster = T, label = T)
#dev.off()
# replot featureplot
#pdf("/home2/gkonop/workdir/06_QC_plots/macaque_FEATUREPLOT_without_COP.pdf", width = 15, height = 35)
#FeaturePlot(macaque, features = marksToPlot, sort = T, raster = T, pt.size = 2)
#dev.off()

# plot Intronic read ratios across samples
pdf("/home2/gkonop/workdir/06_QC_plots/macaque_NuclearFraction_across_samples.pdf")
ggboxplot(macaque[[]], x = 'orig.ident', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# marmoset
marmoset <- readRDS("/home2/gkonop/workdir/01_MARMOSET/MARMOSET_INTEGRATION/Caudate_putamen/ANNOTATED/marmoset_integrated_cuadate_putamen_ANNOTATED.RDS")

# plot Intronic read ratios across samples
pdf("/home2/gkonop/workdir/06_QC_plots/marmoset_NuclearFraction_across_samples.pdf")
ggboxplot(marmoset[[]], x = 'orig.ident', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# mouse
mouse_caud_put <- readRDS("/home2/gkonop/workdir/01_MOUSE/03_CLUSTER_ANNOTATE/Mouse_Caudate_Annotated_FINAL.RDS")
mouse_caud_put$Species <- rep("Mouse", nrow(mouse_caud_put[[]]))
table(mouse_caud_put$newannot)
mouse_caud_put$Tissue <- rep("Caudoputamen", nrow(mouse_caud_put[[]]))
# place COP cells unders OPCs for ease of plotting
a <- gsub("COP", "OPC", mouse_caud_put$newannot) 
mouse_caud_put$newannot <- a
mouse_caud_put$Tissue <- mouse_caud_put$tissue
# replot the UMAP
#pdf("/home2/gkonop/workdir/01_MOUSE/03_CLUSTER_ANNOTATE/ANNOTATED/mouse_caudoputamen_UMAP_without_COPs.pdf")
#DimPlot(mouse_caud_put, group.by = 'newannot', raster = T, label = T)
#dev.off()

# marker genes for mouse
#marksToPlot<- c('Slc17a7', 'Satb2', 'Rbfox3', 'Nrgn', 'Syt1', 'Snap25', 'Gad1', 'Gad2', 'Drd1', 'Drd2', 'Foxp2', 'Tac1', 'Penk', 'Sst', 'Npy', 'Pvalb', 'Mog', 'Plp1', 'Mbp', 'Mobp', 'Slc1a2', 'Slc1a3', 'Aqp4', 'Csf1r', 'Apbb1ip', 'Casz1', 'Ptprz1', 'Pcdh15', 'Sox6', 'Fyn', 'Bcas1', 'Enpp6', 'Gpr17', 'Flt1', 'Dusp1', 'Cobll1', 'Chat', 'Adarb2', 'Pde3a', 'Ncald', 'Tac3', 'Vip', 'Htr1a', 'Th', 'Cck', 'Ttc34', 'Crocc', 'Fhad1', 'Sox4', 'Pdgfrb', 'Pdgfra')
# replot featureplot
#pdf("/home2/gkonop/workdir/06_QC_plots/mouse_FEATUREPLOT_without_COP.pdf", width = 15, height = 35)
#FeaturePlot(mouse_caud_put, features = marksToPlot, sort = T, raster = T, pt.size = 2)
#dev.off()

# bat
bat <- readRDS("/home2/gkonop/workdir/01_BAT/03_CLUSTER_ANNOTATE/Caudate_Putamen/ANNOTATED/bat_integrated_cuadate_putamen_ANNOTATED.RDS")
# place COP cells unders OPCs for ease of plotting
a <- gsub("COP", "OPC", bat$newannot) 
bat$newannot <- a
# replot the UMAP
#pdf("/home2/gkonop/workdir/01_BAT/03_CLUSTER_ANNOTATE/Caudate_Putamen/ANNOTATED/bat_caudate_putamen_UMAP_without_COPs.pdf")
#DimPlot(bat, group.by = 'newannot', raster = T, label = T)
#dev.off()
# set variables
#marksToPlot = c('SLC17A7', 'SATB2', 'RBFOX3', 'NRGN', 'SYT1', 'SNAP25', 'GAD1', 'GAD2', 'DRD1', 'DRD2', 'FOXP2', 'TAC1', 'PENK', 'SST', 'NPY', 'PVALB', 'MOG', 'PLP1', 'MBP', 'MOBP', 'SLC1A2', 'SLC1A3', 'AQP4', 'CSF1R', 'APBB1IP', 'CASZ1', 'PTPRZ1', 'PCDH15', 'SOX6', 'FYN', 'BCAS1', 'ENPP6', 'GPR17', 'FLT1', 'DUSP1', 'COBLL1', 'CHAT', 'ADARB2', 'PDE3A', 'NCALD', 'TAC3', 'VIP', 'HTR1A', 'TH', 'CCK', 'TTC34', 'CROCC', 'FHAD1', 'SOX4', 'PDGFRB', 'PDGFRA')

# replot featureplot
#pdf("/home2/gkonop/workdir/06_QC_plots/bat_FEATUREPLOT_without_COP.pdf", width = 15, height = 35)
#FeaturePlot(bat, features = marksToPlot, sort = T, raster = T, pt.size = 2)
#dev.off()

# plot Intronic read ratios across samples
pdf("/home2/gkonop/workdir/06_QC_plots/bat_NuclearFraction_across_samples.pdf")
ggboxplot(bat[[]], x = 'orig.ident', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()


# ferret
ferret <- readRDS("/home2/gkonop/workdir/01_FERRET/03_CLUSTER_ANNOTATE/ANNOTATED/Ferret_Caudate_Krienen_ANNOTATED.RDS") 
#  featureplot
#pdf("/home2/gkonop/workdir/06_QC_plots/ferret_FEATUREPLOT_without_COP.pdf", width = 15, height = 35)
#FeaturePlot(ferret, features = marksToPlot, sort = T, raster = T, pt.size = 2)
#dev.off()
# plot Intronic read ratios across samples
pdf("/home2/gkonop/workdir/06_QC_plots/ferret_NuclearFraction_across_samples.pdf")
ggboxplot(ferret[[]], x = 'orig.ident', y = 'intronRat') +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

###############################################################
###
## Plot sample-centric plots
###
#### Total UMIs
# extract the metadata
human_meta <- human[[]]
chimp_meta <- chimp[[]]
macaque_meta <- macaque[[]]
marmoset_meta <- marmoset[[]]
bat_meta <- bat[[]]
mouse_meta <- mouse_caud_put[[]]

## merge all the metadata
# Define the columns to extract
columns_to_extract <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "Species", "Tissue", "newannot")

# Extract the specified columns from each metadata dataframe
human_meta_subset <- human_meta[,columns_to_extract]
chimp_meta_subset <- chimp_meta[, columns_to_extract]
macaque_meta_subset <- macaque_meta[, columns_to_extract]
marmoset_meta_subset <- marmoset_meta[, columns_to_extract]
bat_meta_subset <- bat_meta[, columns_to_extract]
mouse_meta_subset <- mouse_meta[, columns_to_extract]

# Combine all subsets into one dataframe
meta_merged <- rbind(human_meta_subset, chimp_meta_subset, macaque_meta_subset, 
                     marmoset_meta_subset, mouse_meta_subset, bat_meta_subset)

# Calculate summary statistics for the bar plots
df_summary <- meta_merged %>%
  group_by(Species, orig.ident, Tissue) %>%
  summarize(
    mean_total_UMI = mean(log10(nCount_RNA), na.rm = TRUE),
    sd_total_UMI = sd(log10(nCount_RNA), na.rm = TRUE),
    n = sum(!is.na(log10(nCount_RNA))),
    sem_total_UMI = sd_total_UMI / sqrt(n),
    .groups = 'drop'
  )
  
# set the order of Species
df_summary$Species <- factor(df_summary$Species, levels = c("Human", "chimp", "Macaque", "Marmoset", "Mouse", "Bat"))

# Create the plot
pdf("/home2/gkonop/workdir/06_QC_plots/Depth_UMI_across_samples_caudate_putamen_all_species_and_tissues.pdf", width= 20, height = 30)
print(
  ggplot(data = df_summary, aes(x = mean_total_UMI, y = interaction(orig.ident, Species, Tissue), fill = factor(Species))) +
  geom_col(position = position_dodge()) +
  geom_errorbar(aes(xmin = mean_total_UMI - sem_total_UMI, xmax = mean_total_UMI + sem_total_UMI),
                position = position_dodge(0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = 'log10(Total UMI)',
    fill = 'Species'
  ) +
  theme(
    legend.position = "right",
    text = element_text(size = 40, face = 'bold'),
    axis.text.x = element_text(size = 30, face = 'bold', hjust = 1, vjust = 0),
    axis.text.y = element_text(size = 30, face = 'bold'),
    axis.title.x = element_text(size = 30, face = 'bold'),
    axis.title.y = element_text(size = 30, face = 'bold')
  )
)
dev.off()
##### Total UMIs of eSPNs ######################################################
## import the SPN seu obj
spn_seu_obj <- readRDS("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/ANNOTATED/GB_SPN_primates_bat_mouse_and_ferret_striosome_matrix_ANNOTATED.RDS")

## extract the SPN labels
# generate unique id for each cell
spn_seu_obj$unique_id <- gsub("-", "_", rownames(spn_seu_obj[[]]))

# extract eSPNs
espns <- subset(spn_seu_obj, newannot_2 == "eSPN")
espn_unique_ids <- espns$unique_id

## merge all the metadata
# Define the columns to extract
columns_to_extract <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "Species", "Tissue", "newannot")

# Extract the specified columns from each metadata dataframe
human_meta_subset <- human_meta[,columns_to_extract]
chimp_meta_subset <- chimp_meta[, columns_to_extract]
macaque_meta_subset <- macaque_meta[, columns_to_extract]
marmoset_meta_subset <- marmoset_meta[, columns_to_extract]
bat_meta_subset <- bat_meta[, columns_to_extract]
mouse_meta_subset <- mouse_meta[, columns_to_extract]

# Combine all subsets into one dataframe
meta_merged <- rbind(human_meta_subset, chimp_meta_subset, macaque_meta_subset, 
                     marmoset_meta_subset, mouse_meta_subset, bat_meta_subset)

# generate unique ids for each cell
meta_merged$unique_id <- gsub("-", "_", rownames(meta_merged))

# put the eSPN label in each species' metadata
a <- ifelse(meta_merged$unique_id %in% espn_unique_ids, "eSPN", meta_merged$newannot)
meta_merged$newannot3 <- a

meta_merged_eSPN <- meta_merged[meta_merged$newannot3 == "eSPN",]

# Calculate summary statistics for the bar plots
df_summary <- meta_merged_eSPN %>%
  group_by(Species, orig.ident) %>%
  summarize(
    mean_total_UMI = mean(log10(nCount_RNA), na.rm = TRUE),
    sd_total_UMI = sd(log10(nCount_RNA), na.rm = TRUE),
    n = sum(!is.na(log10(nCount_RNA))),
    sem_total_UMI = sd_total_UMI / sqrt(n),
    .groups = 'drop'
  )
  
# set the order of Species
df_summary$Species <- factor(df_summary$Species, levels = c("Human", "chimp", "Macaque", "Marmoset", "mouse", "Bat"))
df_summary$orig.ident <- factor(df_summary$orig.ident, levels = unique(df_summary$orig.ident))


# Create the plot
pdf("/home2/gkonop/workdir/06_QC_plots/Depth_UMI_eSPN_across_samples_caudate_putamen_all_species.pdf", width= 20, height = 30)
print(
  ggplot(data = df_summary, aes(x = mean_total_UMI, y = interaction(Species, orig.ident), fill = factor(Species))) +
  geom_col(position = position_dodge()) +
  geom_errorbar(aes(xmin = mean_total_UMI - sem_total_UMI, xmax = mean_total_UMI + sem_total_UMI),
                position = position_dodge(0.9), width = 0.25) +
  theme_minimal() +
  labs(
    x = NULL,
    y = 'log10(Total UMI)',
    fill = 'Species'
  ) +
  theme(
    legend.position = "right",
    text = element_text(size = 40, face = 'bold'),
    axis.text.x = element_text(size = 30, face = 'bold', hjust = 1, vjust = 0),
    axis.text.y = element_text(size = 30, face = 'bold'),
    axis.title.x = element_text(size = 30, face = 'bold'),
    axis.title.y = element_text(size = 30, face = 'bold')
  )
)
dev.off()
######################################################
#### Cell type composition of samples
### Human
stackedbarplot(human[[]], groupx = 'orig.ident', groupfill = 'newannot', "/endosome/work/Neuroinformatics_Core/gkonop/06_QC_plots/Human_integrated_Stacked_AnnotatedSamples.pdf")

stackedbarplot(human[[]], groupx = 'newannot', groupfill = 'orig.ident', "/endosome/work/Neuroinformatics_Core/gkonop/06_QC_plots/Human_integrated_Stacked_AnnotatedSamples2.pdf")

### Chimp
stackedbarplot(chimp[[]], groupx = 'newannot', groupfill = 'orig.ident', "/endosome/work/Neuroinformatics_Core/gkonop/06_QC_plots/Chimp_Stacked_AnnotatedSamples.pdf")

stackedbarplot(chimp[[]], groupx = 'orig.ident', groupfill = 'newannot', "/endosome/work/Neuroinformatics_Core/gkonop/06_QC_plots/Chimp_Stacked_AnnotatedSamples2.pdf")


### Macaque
stackedbarplot(macaque[[]], groupx = 'newannot', groupfill = 'orig.ident', "/endosome/work/Neuroinformatics_Core/gkonop/06_QC_plots/Macaque_Stacked_AnnotatedSamples.pdf")

stackedbarplot(macaque[[]], groupx = 'orig.ident', groupfill = 'newannot', "/endosome/work/Neuroinformatics_Core/gkonop/06_QC_plots/Macaque_Stacked_AnnotatedSamples2.pdf")

### Marmoset
stackedbarplot(marmoset[[]], groupx = 'newannot', groupfill = 'orig.ident', "/endosome/work/Neuroinformatics_Core/gkonop/06_QC_plots/Marmoset_integrated_Stacked_AnnotatedSamples2.pdf")


stackedbarplot(marmoset[[]], groupx = 'orig.ident', groupfill = 'newannot', "/endosome/work/Neuroinformatics_Core/gkonop/06_QC_plots/Marmoset_Stacked_AnnotatedSamples.pdf")

### Mouse
stackedbarplot(mouse_caud_put[[]], groupx = 'newannot', groupfill = 'orig.ident', "/endosome/work/Neuroinformatics_Core/gkonop/06_QC_plots/Mouse_Stacked_AnnotatedSamples.pdf")

stackedbarplot(mouse_caud_put[[]], groupx = 'orig.ident', groupfill = 'newannot', "/endosome/work/Neuroinformatics_Core/gkonop/06_QC_plots/Mouse_Stacked_AnnotatedSamples2.pdf")


### Bat
stackedbarplot(bat[[]], groupx = 'newannot', groupfill = 'orig.ident', "/endosome/work/Neuroinformatics_Core/gkonop/06_QC_plots/Bat_Stacked_AnnotatedSamples.pdf")


stackedbarplot(bat[[]], groupx = 'orig.ident', groupfill = 'newannot', "/endosome/work/Neuroinformatics_Core/gkonop/06_QC_plots/Bat_Stacked_AnnotatedSamples2.pdf")
#################################################################

#### 1. violin plot for each cell type and its best marker(s)
###############################################################
## for loop for each species for specific markers
# Load necessary libraries
library(dplyr)
library(writexl)
library(ggplot2)
library(ggpubr)

# Define a list of Seurat objects, one for each species
species_list <- list(
  human = human,
  chimp = chimp,
  macaque = macaque,
  marmoset = marmoset,
  mouse = mouse_caud_put,
  bat = bat 
)

# Create an empty list to store file paths for each species
file_paths <- list()

# Process each species
for (species_name in names(species_list)) {
  
  # Get the Seurat object for the current species
  seurat_object <- species_list[[species_name]]
  
  # normalize gene expression
  DefaultAssay(seurat_object) <- "RNA"
  seurat_object = NormalizeData(seurat_object)

  # order the cell type names
  # cell annotation order for the plot
  cell_order <- c("OPC","MOL", "Microglia", "Astrocyte", "SPN", "Non_SPN")
  # order the annotations
  seurat_object$newannot <- factor(seurat_object@meta.data$newannot, levels = rev(cell_order))
  # set idents
  Idents(seurat_object) <- seurat_object$newannot

  if(species_name == "mouse") {
  	markers <- c('Pdgfra','Mog','Arhgap15','Apbb1ip','Csf1r','Aqp4','Gja1','Gfap','Ppp1r1b','Tac1', 'Drd1', 'Drd2','Penk','Adora2a','Il1rapl2','Chrm2','Elavl4')
}
  else
  	markers <- c('PDGFRA','MOG','ARHGAP15','APBB1IP','CSF1R','AQP4','GJA1','GFAP','PPP1R1B','TAC1', 'DRD1', 'DRD2','PENK','ADORA2A','IL1RAPL2','CHRM2','ELAVL4')

  # Define file paths for saving plots and summary stats
  plot_dir <- "/home2/gkonop/workdir/06_QC_plots/"
  file_paths[[species_name]] <- list(
    violin_plot_new_markers = paste0(plot_dir, species_name, "_Violin_plot_new_markers.pdf")
  )
    
  # Create and save the total UMI across cell types boxplot
  # Create and save the total UMI across cell types boxplot
  pdf(file_paths[[species_name]]$violin_plot_new_markers)
  # Create the violin plot and explicitly print it
  plot <- VlnPlot(seurat_object, markers, stack = TRUE) +
          theme(legend.position = "none") + 
          ggtitle(paste0(species_name, " Caudate and Putamen")) 
  print(plot)  # Explicitly print the plot
  dev.off()
}
###############################################################
#### 2. heatmap showing each cluster is different in terms of transcriptome 
###############################################################
### Human
# Perform differential expression analysis or extract normalized data
cluster_markers <- FindAllMarkers(human)

# save
saveRDS(cluster_markers, file = "/home2/gkonop/workdir/01_HUMAN/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN/ANNOTATED/human_caudate_putamen_cluster_markers.RDS")

# subset only the positive and min.pct = 0.25, logfc.threshold = 0.25 ones
a <- cluster_markers[(cluster_markers$pct.1 >= 0.25 | cluster_markers$pct.2 >= 0.25 )& cluster_markers$avg_log2FC >= 0.25,]

# get top 10 markers based on log2FC for all the cell types
top_markers <- a %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) # Select top markers per cluster

library(RColorBrewer)
DefaultAssay(human) <- "SCT"
# save
pdf("/home2/gkonop/workdir/06_QC_plots/human_top10_markers_heatmap.pdf")
DoHeatmap(subset(human, downsample = 200),features = top_markers$gene, disp.min = -2, disp.max = 2,group.colors = brewer.pal(10,"Paired"))  +scale_fill_gradientn(colors = brewer.pal(9,"PuOr"))
dev.off()
###############################################################
### Chimp
# Perform differential expression analysis or extract normalized data
chimp_cluster_markers <- FindAllMarkers(chimp,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save
saveRDS(chimp_cluster_markers, file = "/home2/gkonop/workdir/01_CHIMP/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN/ANNOTATED/chimp_caudate_putamen_cluster_markers_filtered.RDS")
chimp_cluster_markers <- readRDS("/home2/gkonop/workdir/01_CHIMP/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN/ANNOTATED/chimp_caudate_putamen_cluster_markers_filtered.RDS")

# get top 10 markers based on log2FC for all the cell types
chimp_top_markers <- chimp_cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) # Select top markers per cluster

library(RColorBrewer)
DefaultAssay(chimp) <- "RNA"
chimp <- ScaleData(chimp)

# save
pdf("/home2/gkonop/workdir/06_QC_plots/chimp_top10_markers_heatmap.pdf")
DoHeatmap(subset(chimp, downsample = 200),features = chimp_top_markers$gene, disp.min = -2, disp.max = 2,group.colors = brewer.pal(10,"Paired"))  +scale_fill_gradientn(colors = brewer.pal(9,"PuOr"))
dev.off()
###############################################################
### Macaque
# Perform differential expression analysis or extract normalized data
macaque_cluster_markers <- FindAllMarkers(macaque,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save
saveRDS(macaque_cluster_markers, file = "/home2/gkonop/workdir/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/macaque_caudate_putamen_cluster_markers_filtered.RDS")
macaque_cluster_markers <- readRDS("/home2/gkonop/workdir/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/macaque_caudate_putamen_cluster_markers_filtered.RDS")

# get top 10 markers based on log2FC for all the cell types
macaque_top_markers <- macaque_cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) # Select top markers per cluster

library(RColorBrewer)
DefaultAssay(macaque) <- "RNA"
macaque <- ScaleData(macaque)

# save
pdf("/home2/gkonop/workdir/06_QC_plots/macaque_top10_markers_heatmap.pdf")
DoHeatmap(subset(macaque, downsample = 200),features = macaque_top_markers$gene, disp.min = -2, disp.max = 2,group.colors = brewer.pal(10,"Paired"))  +scale_fill_gradientn(colors = brewer.pal(9,"PuOr"))
dev.off()
###############################################################
### Marmoset
# Perform differential expression analysis or extract normalized data
marmoset_cluster_markers <- FindAllMarkers(marmoset,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save
saveRDS(marmoset_cluster_markers, file = "/home2/gkonop/workdir/01_MARMOSET/MARMOSET_INTEGRATION/Caudate_putamen/ANNOTATED/marmoset_caudate_putamen_cluster_markers_filtered.RDS")

# get top 10 markers based on log2FC for all the cell types
marmoset_top_markers <- marmoset_cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) # Select top markers per cluster

library(RColorBrewer)
DefaultAssay(marmoset) <- "RNA"
marmoset <- ScaleData(marmoset)

# save
pdf("/home2/gkonop/workdir/06_QC_plots/marmoset_top10_markers_heatmap.pdf")
DoHeatmap(subset(marmoset, downsample = 200),features = marmoset_top_markers$gene, disp.min = -2, disp.max = 2,group.colors = brewer.pal(10,"Paired"))  +scale_fill_gradientn(colors = brewer.pal(9,"PuOr"))
dev.off()
###############################################################
### Mouse
# Perform differential expression analysis or extract normalized data
mouse_cluster_markers <- FindAllMarkers(mouse_caud_put,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save
saveRDS(mouse_caud_put_cluster_markers, file = "/home2/gkonop/workdir/01_MOUSE/03_CLUSTER_ANNOTATE/ANNOTATED/mouse_caudate_putamen_cluster_markers_filtered.RDS")

# get top 10 markers based on log2FC for all the cell types
mouse_top_markers <- mouse_cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) # Select top markers per cluster

library(RColorBrewer)
DefaultAssay(mouse_caud_put) <- "RNA"
mouse_caud_put <- ScaleData(mouse_caud_put)

# save
pdf("/home2/gkonop/workdir/06_QC_plots/mouse_top10_markers_heatmap.pdf")
DoHeatmap(subset(mouse_caud_put, downsample = 200),features = mouse_top_markers$gene, disp.min = -2, disp.max = 2,group.colors = brewer.pal(10,"Paired"))  +scale_fill_gradientn(colors = brewer.pal(9,"PuOr"))
dev.off()
###############################################################
### Bat
# Perform differential expression analysis or extract normalized data
bat_cluster_markers <- FindAllMarkers(bat,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save
saveRDS(bat_cluster_markers, file = "/home2/gkonop/workdir/01_BAT/03_CLUSTER_ANNOTATE/Caudate_Putamen/ANNOTATED/bat_caudate_putamen_cluster_markers_filtered.RDS")

# get top 10 markers based on log2FC for all the cell types
bat_top_markers <- bat_cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) # Select top markers per cluster

library(RColorBrewer)
DefaultAssay(bat) <- "RNA"
bat <- ScaleData(bat)

# save
pdf("/home2/gkonop/workdir/06_QC_plots/bat_top10_markers_heatmap.pdf")
DoHeatmap(subset(bat, downsample = 200),features = bat_top_markers$gene, disp.min = -2, disp.max = 2,group.colors = brewer.pal(10,"Paired"))  +scale_fill_gradientn(colors = brewer.pal(9,"PuOr"))
dev.off()
###############################################################
### Ferret
# Perform differential expression analysis or extract normalized data
# set Idents
Idents(ferret) <- ferret$newannot
ferret_cluster_markers <- FindAllMarkers(ferret,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save
saveRDS(ferret_cluster_markers, file = "/home2/gkonop/workdir/01_FERRET/03_CLUSTER_ANNOTATE/Caudate_Putamen/ANNOTATED/ferret_caudate_cluster_markers_filtered.RDS")

# get top 10 markers based on log2FC for all the cell types
ferret_top_markers <- ferret_cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) # Select top markers per cluster

library(RColorBrewer)
DefaultAssay(ferret) <- "RNA"
ferret <- ScaleData(ferret)

# save
pdf("/home2/gkonop/workdir/06_QC_plots/ferret_top10_markers_heatmap.pdf")
DoHeatmap(subset(ferret, downsample = 200),features = ferret_top_markers$gene, disp.min = -2, disp.max = 2,group.colors = brewer.pal(10,"Paired"))  +scale_fill_gradientn(colors = brewer.pal(9,"PuOr"))
dev.off()

###############################################################
#### 3. Total UMI across samples, tissues and cell types
###############################################################
### Human
# Convert Seurat metadata to a data.frame
human_meta <- human[[]]

# generate a vector with new terms
a <- ifelse(human_meta$newannot %in% c("Astrocyte", "MOL","OPC", "COP", "Microglia"), "Glia", "Neuron")
human$Glia_or_neu <- a


# plot for the total UMI across cell types
pdf("/home2/gkonop/workdir/06_QC_plots/human_caudate_putamen_Depth_UMI_across_tissue_CellType.pdf")
ggboxplot(human_meta, x = 'newannot', y = 'nCount_RNA', ylim = c(0,50000), facet.by = "Tissue") +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# total UMI across samples
pdf("/home2/gkonop/workdir/06_QC_plots/human_caudate_putamen_Depth_UMI_across_samples.pdf")
ggboxplot(human_meta, x = 'orig.ident', y = 'nCount_RNA', ylim = c(0,20000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()
###############################################################
### Chimp
## add a new metadata col
# extract metadata
chimp_meta <- chimp[[]]

# generate a vector with new terms
a <- ifelse(chimp_meta$newannot %in% c("Astrocyte", "MOL","OPC", "COP", "Microglia"), "Glia", "Neuron")
chimp$Glia_or_neu <- a

# Convert Seurat metadata to a data.frame
chimp_meta <- chimp[[]]

# total UMI across tissue
pdf("/home2/gkonop/workdir/06_QC_plots/chimp_caudate_putamen_Depth_UMI_across_tissue.pdf")
ggboxplot(chimp_meta, x = 'Tissue', y = 'nCount_RNA', ylim = c(0,20000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# plot for the total UMI across cell types
pdf("/home2/gkonop/workdir/06_QC_plots/chimp_caudate_putamen_Depth_UMI_across_tissue_CellType.pdf")
ggboxplot(chimp_meta, x = 'newannot', y = 'nCount_RNA', ylim = c(0,50000), facet.by = "Tissue") +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# total UMI across samples
pdf("/home2/gkonop/workdir/06_QC_plots/chimp_caudate_putamen_Depth_UMI_across_samples.pdf")
ggboxplot(chimp_meta, x = 'orig.ident', y = 'nCount_RNA', ylim = c(0,20000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()
###############################################################
### Macaque
## add a new metadata col
# extract metadata
macaque_meta <- macaque[[]]

# generate a vector with new terms
a <- ifelse(macaque_meta$newannot %in% c("Astrocyte", "MOL","OPC", "COP", "Microglia"), "Glia", "Neuron")
macaque$Glia_or_neu <- a

# Convert Seurat metadata to a data.frame
macaque_meta <- macaque[[]]

# total UMI across tissue
pdf("/home2/gkonop/workdir/06_QC_plots/macaque_caudate_putamen_Depth_UMI_across_tissue.pdf")
ggboxplot(macaque_meta, x = 'Tissue', y = 'nCount_RNA', ylim = c(0,20000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# plot for the total UMI across cell types
pdf("/home2/gkonop/workdir/06_QC_plots/macaque_caudate_putamen_Depth_UMI_across_tissue_CellType.pdf")
ggboxplot(macaque_meta, x = 'newannot', y = 'nCount_RNA', ylim = c(0,50000), facet.by = "Tissue") +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# total UMI across samples
pdf("/home2/gkonop/workdir/06_QC_plots/macaque_caudate_putamen_Depth_UMI_across_samples.pdf")
ggboxplot(macaque_meta, x = 'orig.ident', y = 'nCount_RNA', ylim = c(0,20000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

###############################################################
### Marmoset
## add a new metadata col
# extract metadata
marmoset_meta <- marmoset[[]]

# generate a vector with new terms
a <- ifelse(marmoset_meta$newannot %in% c("Astrocyte", "MOL","OPC", "COP", "Microglia"), "Glia", "Neuron")
marmoset$Glia_or_neu <- a

# Convert Seurat metadata to a data.frame
marmoset_meta <- marmoset[[]]

# total UMI across samples
pdf("/home2/gkonop/workdir/06_QC_plots/marmoset_caudate_putamen_Depth_UMI_across_samples.pdf")
ggboxplot(marmoset_meta, x = 'orig.ident', y = 'nCount_RNA', ylim = c(0,20000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()


# total UMI across tissue
pdf("/home2/gkonop/workdir/06_QC_plots/marmoset_caudate_putamen_Depth_UMI_across_tissue.pdf")
ggboxplot(marmoset_meta, x = 'Tissue', y = 'nCount_RNA', ylim = c(0,20000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# plot for the total UMI across cell types
pdf("/home2/gkonop/workdir/06_QC_plots/marmoset_caudate_putamen_Depth_UMI_across_tissue_CellType.pdf")
ggboxplot(marmoset_meta, x = 'newannot', y = 'nCount_RNA', ylim = c(0,50000), facet.by = "Tissue") +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()
###############################################################
### Mouse
## add a new metadata col
# generate a vector with new terms
a <- ifelse(mouse_caud_put$newannot %in% c("Astrocyte", "MOL","OPC", "COP", "Microglia"), "Glia", "Neuron")
mouse_caud_put$Glia_or_neu <- a

# Convert Seurat metadata to a data.frame
mouse_meta <- mouse_caud_put[[]]

# total UMI across samples
pdf("/home2/gkonop/workdir/06_QC_plots/mouse_caudoputamen_Depth_UMI_across_samples.pdf")
ggboxplot(mouse_meta, x = 'orig.ident', y = 'nCount_RNA', ylim = c(0,20000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# total UMI across tissue
pdf("/home2/gkonop/workdir/06_QC_plots/mouse_caudoput_Depth_UMI_across_tissue.pdf")
ggboxplot(mouse_meta, x = 'Tissue', y = 'nCount_RNA', ylim = c(0,20000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# plot for the total UMI across cell types
pdf("/home2/gkonop/workdir/06_QC_plots/mouse_caudoput_Depth_UMI_across_tissue_CellType.pdf")
ggboxplot(mouse_meta, x = 'newannot', y = 'nCount_RNA', ylim = c(0,30000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()
###############################################################
### Bat
## add a new metadata col
# extract metadata
bat_meta <- bat[[]]

# generate a vector with new terms
a <- ifelse(bat_meta$newannot %in% c("Astrocyte", "MOL","OPC", "COP", "Microglia"), "Glia", "Neuron")
bat$Glia_or_neu <- a

# Convert Seurat metadata to a data.frame
bat_meta <- bat[[]]

# total UMI across tissue
pdf("/home2/gkonop/workdir/06_QC_plots/bat_caudate_putamen_Depth_UMI_across_tissue.pdf")
ggboxplot(bat_meta, x = 'Tissue', y = 'nCount_RNA', ylim = c(0,20000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# plot for the total UMI across cell types
pdf("/home2/gkonop/workdir/06_QC_plots/bat_caudate_putamen_Depth_UMI_across_tissue_CellType.pdf")
ggboxplot(bat_meta, x = 'newannot', y = 'nCount_RNA', ylim = c(0,30000), facet.by = "Tissue") +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()
###############################################################
### Ferret
## add a new metadata col
# extract metadata
ferret_meta <- ferret[[]]

# generate a vector with new terms
a <- ifelse(ferret_meta$newannot %in% c("Astrocyte", "MOL","OPC", "COP", "Microglia"), "Glia", "Neuron")
ferret$Glia_or_neu <- a

# Convert Seurat metadata to a data.frame
ferret_meta <- ferret[[]]

# total UMI across tissue
pdf("/home2/gkonop/workdir/06_QC_plots/ferret_caudate_putamen_Depth_UMI_across_tissue.pdf")
ggboxplot(ferret_meta, x = 'Tissue', y = 'nCount_RNA', ylim = c(0,20000)) +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# plot for the total UMI across cell types
pdf("/home2/gkonop/workdir/06_QC_plots/ferret_caudate_putamen_Depth_UMI_across_tissue_CellType.pdf")
ggboxplot(ferret_meta, x = 'newannot', y = 'nCount_RNA', ylim = c(0,50000), facet.by = "Tissue") +
rotate_x_text(90) +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

###############################################################
#### 4. Neuron and  glia plots
# Subset the human_meta data frame to include specific columns
subset_meta <- human_meta %>%
  select(nFeature_RNA, nCount_RNA, orig.ident, newannot, Tissue, Glia_or_neu)

# Group by Tissue, orig.ident, and newannot, then calculate mean and standard error of nFeature_RNA
(summary_stats <- subset_meta %>%
  group_by(Tissue, Glia_or_neu) %>%
  summarise(
    mean_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
    se_nFeature_RNA = sd(nFeature_RNA, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ))

# A tibble: 4 × 4
  Tissue  Glia_or_neu mean_nFeature_RNA se_nFeature_RNA
  <chr>   <chr>                   <dbl>           <dbl>
1 Caudate Glia                    1760.            3.98
2 Caudate Neuron                  5064.           19.2 
3 Putamen Glia                    1402.            3.58
4 Putamen Neuron                  3817.           13.2 


# A tibble: 12 × 4
   Tissue  newannot  mean_nFeature_RNA se_nFeature_RNA
   <chr>   <chr>                 <dbl>           <dbl>
 1 Caudate Astrocyte             2408.           19.2 
 2 Caudate MOL                   1691.            3.94
 3 Caudate Microglia             1308.           10.1 
 4 Caudate Non_SPN               4758.           51.3 
 5 Caudate OPC                   2229.           17.4 
 6 Caudate SPN                   5126.           20.6 
 7 Putamen Astrocyte             2147.           17.1 
 8 Putamen MOL                   1286.            3.44
 9 Putamen Microglia             1102.            9.84
10 Putamen Non_SPN               5002.           43.3 
11 Putamen OPC                   2092.           16.4 
12 Putamen SPN                   3618.           13.2 

# repeat for total UMI
(summary_stats <- subset_meta %>%
  group_by(Tissue, Glia_or_neu) %>%
  summarise(
    mean_nFeature_RNA = mean(nCount_RNA, na.rm = TRUE),
    se_nFeature_RNA = sd(nCount_RNA, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ))
# A tibble: 4 × 4
  Tissue  Glia_or_neu mean_nFeature_RNA se_nFeature_RNA
  <chr>   <chr>                   <dbl>           <dbl>
1 Caudate Glia                    4009.            16.3
2 Caudate Neuron                 23758.           172. 
3 Putamen Glia                    2959.            11.9
4 Putamen Neuron                 14714.            91.7

# A tibble: 12 × 4
  Tissue  newannot  mean_nFeature_RNA se_nFeature_RNA
   <chr>   <chr>                 <dbl>           <dbl>
 1 Caudate Astrocyte             6301.            84.1
 2 Caudate MOL                   3777.            16.5
 3 Caudate Microglia             2486.            31.0
 4 Caudate Non_SPN              20035.           455. 
 5 Caudate OPC                   5441.            73.0
 6 Caudate SPN                  24523.           185. 
 7 Putamen Astrocyte             5144.            63.3
 8 Putamen MOL                   2628.            11.1
 9 Putamen Microglia             1919.            25.9
10 Putamen Non_SPN              24185.           394. 
11 Putamen OPC                   5074.            64.3
12 Putamen SPN                  13123.            79.4

#############################
# Load necessary libraries
library(dplyr)
library(writexl)
library(ggplot2)
library(ggpubr)

# Define a list of Seurat objects, one for each species
species_list <- list(
  human = human,
  chimp = chimp,
  macaque = macaque
)

# Create an empty list to store summary statistics for all species
all_summary_stats_genes <- list()
all_summary_stats_UMI <- list()

# Create an empty list to store file paths for each species
file_paths <- list()

# Process each species
for (species_name in names(species_list)) {
  
  # Get the Seurat object for the current species
  seurat_object <- species_list[[species_name]]
  
  # Generate a vector with new terms
  a <- ifelse(seurat_object@meta.data$newannot %in% c("Astrocyte", "MOL", "OPC", "COP", "Microglia"), "Glia", "Neuron")
  seurat_object$Glia_or_neu <- a
  
  # Convert Seurat metadata to a data.frame
  meta_data <- seurat_object[[]]
  
  # Define file paths for saving plots and summary stats
  plot_dir <- "/home2/gkonop/workdir/06_QC_plots/"
  file_paths[[species_name]] <- list(
    UMI_across_tissue = paste0(plot_dir, species_name, "_Depth_UMI_across_tissue.pdf"),
    UMI_across_celltype = paste0(plot_dir, species_name, "_Depth_UMI_across_tissue_CellType.pdf"),
    summary_stats = paste0(plot_dir, species_name, "_summary_stats.xlsx")
  )
  
  # Create and save the total UMI across tissue boxplot
  pdf(file_paths[[species_name]]$UMI_across_tissue)
  ggboxplot(meta_data, x = 'Tissue', y = 'nCount_RNA', ylim = c(0, 20000)) +
    rotate_x_text(90) +
    theme(text = element_text(size = 20, face = 'bold'))
  dev.off()
  
  # Create and save the total UMI across cell types boxplot
  pdf(file_paths[[species_name]]$UMI_across_celltype)
  ggboxplot(meta_data, x = 'newannot', y = 'nCount_RNA', ylim = c(0, 50000), facet.by = "Tissue") +
    rotate_x_text(90) +
    theme(text = element_text(size = 20, face = 'bold'))
  dev.off()
  
  # Subset the meta_data to include specific columns
  subset_meta <- meta_data %>%
    select(nFeature_RNA, nCount_RNA, orig.ident, newannot, Tissue, Glia_or_neu)
  
  # Calculate summary statistics for total number of genes (nFeature_RNA)
  summary_stats_genes <- subset_meta %>%
    group_by(Tissue, Glia_or_neu) %>%
    summarise(
      mean_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
      se_nFeature_RNA = sd(nFeature_RNA, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
  
  # Calculate summary statistics for total number of UMIs (nCount_RNA)
  summary_stats_UMI <- subset_meta %>%
    group_by(Tissue, Glia_or_neu) %>%
    summarise(
      mean_nCount_RNA = mean(nCount_RNA, na.rm = TRUE),
      se_nCount_RNA = sd(nCount_RNA, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
  
  # Store summary statistics in the list
  all_summary_stats_genes[[species_name]] <- summary_stats_genes
  all_summary_stats_UMI[[species_name]] <- summary_stats_UMI
}

# Save all summary statistics to separate Excel files for each species
for (species_name in names(all_summary_stats_genes)) {
  write_xlsx(
    list(
      "Total_Genes" = all_summary_stats_genes[[species_name]],
      "Total_UMI" = all_summary_stats_UMI[[species_name]]
    ),
    path = file_paths[[species_name]]$summary_stats
  )
}


#### 5. % of each cell type within the species per Tissue
#########################################################################
### human
# Group by species and cell type and calculate percentages
cell_type_percentage <- human@meta.data %>%
  group_by(Tissue, newannot) %>%
  tally() %>%
  group_by(Tissue) %>%
  mutate(percentage = n / sum(n) * 100)


pdf("/home2/gkonop/workdir/06_QC_plots/human_cell_type_percentage.pdf")
# Create the stacked bar plot
ggplot(cell_type_percentage, aes(x = Tissue, y = percentage, fill = newannot)) +
  geom_bar(stat = "identity", color = "black", position = "stack") +
  theme_minimal() +
  labs(title = "Cell Type Percentages in Caudate and Putamen",
       x = "Tissue",
       y = "Percentage",
       fill = "Cell Type")  
dev.off()
#########################################################################
### chimp
# Group by species and cell type and calculate percentages
cell_type_percentage <- chimp@meta.data %>%
  group_by(Tissue, newannot) %>%
  tally() %>%
  group_by(Tissue) %>%
  mutate(percentage = n / sum(n) * 100)

# Save the plot as a PDF
pdf("/home2/gkonop/workdir/06_QC_plots/chimp_cell_type_percentage.pdf")
# Create the stacked bar plot
ggplot(cell_type_percentage, aes(x = Tissue, y = percentage, fill = newannot)) +
  geom_bar(stat = "identity", color = "black", position = "stack") +
  theme_minimal() +
  labs(title = "Cell Type Percentages in Caudate and Putamen",
       x = "Tissue",
       y = "Percentage",
       fill = "Cell Type")
dev.off()
#########################################################################
### macaque
# Group by species and cell type and calculate percentages
cell_type_percentage <- macaque@meta.data %>%
  group_by(Tissue, newannot) %>%
  tally() %>%
  group_by(Tissue) %>%
  mutate(percentage = n / sum(n) * 100)

# Save the plot as a PDF
pdf("/home2/gkonop/workdir/06_QC_plots/macaque_cell_type_percentage.pdf")
# Create the stacked bar plot
ggplot(cell_type_percentage, aes(x = Tissue, y = percentage, fill = newannot)) +
  geom_bar(stat = "identity", color = "black", position = "stack") +
  theme_minimal() +
  labs(title = "Cell Type Percentages in Caudate and Putamen",
       x = "Tissue",
       y = "Percentage",
       fill = "Cell Type")
dev.off()
#########################################################################
### marmoset
# Group by species and cell type and calculate percentages
cell_type_percentage <- marmoset@meta.data %>%
  group_by(Tissue, newannot) %>%
  tally() %>%
  group_by(Tissue) %>%
  mutate(percentage = n / sum(n) * 100)

# Save the plot as a PDF
pdf("/home2/gkonop/workdir/06_QC_plots/marmoset_cell_type_percentage.pdf")
# Create the stacked bar plot
ggplot(cell_type_percentage, aes(x = Tissue, y = percentage, fill = newannot)) +
  geom_bar(stat = "identity", color = "black", position = "stack") +
  theme_minimal() +
  labs(title = "Cell Type Percentages in Caudate and Putamen",
       x = "Tissue",
       y = "Percentage",
       fill = "Cell Type")
dev.off()
#########################################################################
### mouse
# Group by species and cell type and calculate percentages
cell_type_percentage <- mouse_caud_put@meta.data %>%
  group_by(Tissue, newannot) %>%
  tally() %>%
  group_by(Tissue) %>%
  mutate(percentage = n / sum(n) * 100)

# Save the plot as a PDF
pdf("/home2/gkonop/workdir/06_QC_plots/mouse_cell_type_percentage.pdf")
# Create the stacked bar plot
ggplot(cell_type_percentage, aes(x = Tissue, y = percentage, fill = newannot)) +
  geom_bar(stat = "identity", color = "black", position = "stack") +
  theme_minimal() +
  labs(title = "Cell Type Percentages in Caudoputamen",
       x = "Tissue",
       y = "Percentage",
       fill = "Cell Type")
dev.off()
#################################################################################
### bat
# Group by species and cell type and calculate percentages
cell_type_percentage <- bat@meta.data %>%
  group_by(Tissue, newannot) %>%
  tally() %>%
  group_by(Tissue) %>%
  mutate(percentage = n / sum(n) * 100)

# Save the plot as a PDF
pdf("/home2/gkonop/workdir/06_QC_plots/bat_cell_type_percentage.pdf")
# Create the stacked bar plot
ggplot(cell_type_percentage, aes(x = Tissue, y = percentage, fill = newannot)) +
  geom_bar(stat = "identity", color = "black", position = "stack") +
  theme_minimal() +
  labs(title = "Cell Type Percentages in Caudate and Putamen",
       x = "Tissue",
       y = "Percentage",
       fill = "Cell Type")
dev.off()
