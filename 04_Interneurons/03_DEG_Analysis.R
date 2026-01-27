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
table(bat_filt$newannot3)

         Astrocyte           CCK_VIP-           CCK_VIP+               CHAT 
             13610                302                 95                537 
               COP         FOXP2_EYA2    FOXP2_TSHZ2_BAT           LMO3_BAT 
               269                306                229                832 
         Microglia                MOL                OPC              PDGFD 
              3295              27952               4627                168 
PDGFD_PTHLH_PVALB- PDGFD_PTHLH_PVALB+                SPN            SST_NPY 
               293                526              63799               1025 
              TAC3                 TH 
               355                122 
###########################################################################
#### find markers of LMO3_Bat within bat putamen compared to other interneurons
# subset Putamen
sub <- subset(bat_filt, Tissue == "Putamen")
# subset the interneurons
sub = subset(sub, subset = newannot3 %in% c("MOL", "OPC", "COP", "Astrocyte", "Microglia", "SPN"), invert = T)

# set identity
Idents(sub) <- sub$newannot3

## filter the genes
# Step 1: Retrieve the count data from the Seurat object
counts_data <- GetAssayData(sub, assay = "RNA", slot = "counts")

# Step 2: Define the smallest group size for filtering
smallestGroupSize <- 4

# Step 3: Apply the filtering logic: Keep genes with at least 5 counts in at least 'smallestGroupSize' samples
keep <- rowSums(counts_data >= 5) >= smallestGroupSize
sum(keep) #[1] 2370

# Step 4: Subset the Seurat object to keep only the filtered genes
sub_filtered <- subset(sub, features = rownames(counts_data)[keep])

# Now you can perform FindMarkers between LMO3_BAT and the other interneurons
# Define your list of interneuron types
interneuron_types <- c("CCK_VIPneg", "CCK_VIPpos", "CHAT", "FOXP2_EYA2", "FOXP2_TSHZ2_BAT", 
                       "SST_NPY", "TAC3", "TH", "PDGFD", "PDGFD_PTHLH_PVALBneg", "PDGFD_PTHLH_PVALBpos")

# Initialize an empty list to store the results
all_markers_filtered_lmo3 <- list()

Idents(sub_filtered) <- sub_filtered$newannot3

# Step 5: Perform FindMarkers between LMO3_BAT and each interneuron type
for (interneuron in interneuron_types) {
  markers <- FindMarkers(sub_filtered, ident.1 = "LMO3_BAT", ident.2 = interneuron, logfc.threshold = 0)
  
  # Add pct_ratio column
  markers$pct_ratio <- markers$pct.1 / markers$pct.2
  
  # Store the results in the list
  all_markers_filtered_lmo3[[interneuron]] <- markers
}

# Save the filtered results
saveRDS(all_markers_filtered_lmo3, "batPutamen_LMO3vsInterneuronFindMarkersDEGs_Filtered.RDS")

# save in an excel (each list element in a spreadsheet)
# Install and load the openxlsx package if not already installed
#if (!require("openxlsx")) install.packages("openxlsx", dependencies = TRUE)
library(openxlsx)

# Define the file path where you want to save the Excel file
file_path <- "batPutamen_LMO3vsInterneurons_FindMarkersDEGs_Filtered.xlsx"

# Create a new workbook
wb <- createWorkbook()

# Loop through each element in combined_markers_filtered and add it to a separate sheet
for (interneuron in names(all_markers_filtered_lmo3)) {
  # Get the markers for the current interneuron
  markers <- all_markers_filtered_lmo3[[interneuron]]
    
  # Add the data frame as a new sheet in the workbook
  addWorksheet(wb, sheetName = interneuron)
  writeData(wb, sheet = interneuron, markers, rowNames = TRUE)
}

# Save the workbook to the file
saveWorkbook(wb, file_path, overwrite = TRUE)

# Message to confirm the file has been saved
cat("Excel file saved at:", file_path, "\n")
####################################################
sub_filtered$newannot4 <- ifelse(sub_filtered$newannot3 == "LMO3_BAT", "LMO3_BAT", "Other")
Idents(sub_filtered) <- sub_filtered$newannot4
table(sub_filtered$newannot4)
#LMO3_BAT    Other 
#     832     1784 

# Perform FindMarkers for LMO3 vs all to generate p vals in a single plo
a <- FindMarkers(sub_filtered, ident.1 = "LMO3_BAT", ident.2 = "Other", logfc.threshold = 0)

# Add pct_ratio column
a$pct_ratio <- a$pct.1 / a$pct.2

# Add a column for DEG status based on log2FC and p-value adjusted
a$diffexpressed <- ifelse(a$avg_log2FC > 0.6 & a$p_val_adj < 0.05, "UP",
                                   ifelse(a$avg_log2FC < -0.6 & a$p_val_adj < 0.05, "DOWN", "NO"))
  
# Label the significantly differentially expressed genes
labeled <- a[abs(a$avg_log2FC) > 0.7 & a$p_val_adj < 0.05, ]

# Ensure "LMO3" is labeled by adding it to a delabel column or directly referring to it
a$delabel <- ifelse(rownames(a) == "LMO3", "LMO3", rownames(a))
 
# save
write.csv(a, "batPutamen_LMO3vsOtherInterneurons_FindMarkersDEGs_Filtered.csv", quote = F, row.names = T) 


# Create the volcano plot
plot <- ggplot(
  a,
  aes(
    x = avg_log2FC,
    y = -log10(p_val_adj),
    col = diffexpressed,
    label = delabel
  )
) +
  theme_bw() +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "grey", linetype = "longdash") +
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "longdash") +
  scale_colour_manual(name = "DEG", values = c("blue", "black", "darkgreen")) +

  geom_text_repel(
    data = labeled,
    aes(label = rownames(labeled)),
    nudge_y = 0.2,
    box.padding = 0.3,
    point.padding = 0.3
  ) +
  geom_text_repel(
    data = subset(a, delabel == "LMO3"),
    aes(label = delabel),
    nudge_y = 0.2,
    box.padding = 0.3,
    point.padding = 0.3
  ) +

  scale_x_continuous(
    limits = c(-5.5, 5.5),
    breaks = c(-5, -2.5, 0, 2.5, 5)
  ) +
  scale_y_continuous(
    limits = c(0, 250),
    breaks = seq(0, 250, by = 50)
  ) +

  ggtitle("Volcano Plot: LMO3 vs All Other Interneurons")

# Save the plot for each comparison to a PDF file
pdf(paste0("VolcanoPlot_LMO3_vs_AllOtherInterneurons_withinBat_FindMarkers_DEGs.pdf"))
print(plot)
dev.off()
####################################################
#### find markers of FOXP2_TSHZ2_BAT within bat putamen compared to other interneurons
# Now you can perform FindMarkers between FOXP2_TSHZ2_BAT and the other interneurons
# Define your list of interneuron types
interneuron_types <- c("CCK_VIPneg", "CCK_VIPpos", "CHAT", "FOXP2_EYA2", "LMO3_BAT", 
                       "SST_NPY", "TAC3", "TH", "PDGFD", "PDGFD_PTHLH_PVALBneg", "PDGFD_PTHLH_PVALBpos")

# Initialize an empty list to store the results
all_markers_filtered_foxp2_tshz2 <- list()

Idents(sub_filtered) <- sub_filtered$newannot3

# Step 5: Perform FindMarkers between FOXP2_TSHZ2_BAT and each interneuron type
for (interneuron in interneuron_types) {
  markers <- FindMarkers(sub_filtered, ident.1 = "FOXP2_TSHZ2_BAT", ident.2 = interneuron, logfc.threshold = 0)
  
  # Add pct_ratio column
  markers$pct_ratio <- markers$pct.1 / markers$pct.2
  
  # Store the results in the list
  all_markers_filtered_foxp2_tshz2[[interneuron]] <- markers
}

# Save the filtered results
saveRDS(all_markers_filtered_foxp2_tshz2, "batPutamen_FOXP2_TSHZ2vsInterneuronFindMarkersDEGs_Filtered.RDS")

# save in an excel (each list element in a spreadsheet)
# Install and load the openxlsx package if not already installed
#if (!require("openxlsx")) install.packages("openxlsx", dependencies = TRUE)
library(openxlsx)

# Define the file path where you want to save the Excel file
file_path <- "batPutamen_FOXP2_TSHZ2vsInterneurons_FindMarkersDEGs_Filtered.xlsx"

# Create a new workbook
wb <- createWorkbook()

# Loop through each element in combined_markers_filtered and add it to a separate sheet
for (interneuron in names(all_markers_filtered_foxp2_tshz2)) {
  # Get the markers for the current interneuron
  markers <- all_markers_filtered_foxp2_tshz2[[interneuron]]
    
  # Add the data frame as a new sheet in the workbook
  addWorksheet(wb, sheetName = interneuron)
  writeData(wb, sheet = interneuron, markers, rowNames = TRUE)
}

# Save the workbook to the file
saveWorkbook(wb, file_path, overwrite = TRUE)

# Message to confirm the file has been saved
cat("Excel file saved at:", file_path, "\n")
####################################################
sub_filtered$newannot4 <- ifelse(sub_filtered$newannot3 == "FOXP2_TSHZ2_BAT", "FOXP2_TSHZ2_BAT", "Other")
Idents(sub_filtered) <- sub_filtered$newannot4
table(sub_filtered$newannot4)
#FOXP2_TSHZ2_BAT    Other 
#     227     2389 

# Perform FindMarkers for FOXP2_TSHZ2 vs all to generate p vals in a single plot
a <- FindMarkers(sub_filtered, ident.1 = "FOXP2_TSHZ2_BAT", ident.2 = "Other", logfc.threshold = 0)

# Add pct_ratio column
a$pct_ratio <- a$pct.1 / a$pct.2

# Add a column for DEG status based on log2FC and p-value adjusted
a$diffexpressed <- ifelse(a$avg_log2FC > 0.6 & a$p_val_adj < 0.05, "UP",
                                   ifelse(a$avg_log2FC < -0.6 & a$p_val_adj < 0.05, "DOWN", "NO"))
  
# Label the significantly differentially expressed genes
labeled <- a[abs(a$avg_log2FC) > 0.7 & a$p_val_adj < 0.05, ]


# Ensure "FOXP2_TSHZ2" is labeled by adding it to a delabel column or directly referring to it
a$delabel <- ifelse(rownames(a) %in% c("FOXP2", "TSHZ2"), rownames(a), NA

a_foxp2_tshz2 <- a
# save
write.csv(a_foxp2_tshz2, "batPutamen_FOXP2_TSHZ2vsOtherInterneurons_FindMarkersDEGs_Filtered.csv", quote = F, row.names = T) 


# Create the volcano plot
plot <- ggplot(
  a_foxp2_tshz2,
  aes(
    x = avg_log2FC,
    y = -log10(p_val_adj),
    col = diffexpressed,
    label = delabel
  )
) +
  theme_bw() +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "grey", linetype = "longdash") +
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "longdash") +
  scale_colour_manual(name = "DEG", values = c("blue", "black", "darkgreen")) +

geom_text_repel(data = labeled, aes(label = rownames(labeled)), nudge_y = 0.2, box.padding = 0.3, point.padding = 0.3) +
  # Label significant genes
  geom_text_repel(data = subset(a_foxp2_tshz2, delabel == "TSHZ2"), aes(label = delabel), nudge_y = 0.2, box.padding = 0.3, point.padding = 0.3) +  # Ensure TSHZ2 is labeled

  scale_x_continuous(
    limits = c(-5.5, 5.5),
    breaks = c(-5, -2.5, 0, 2.5, 5)
  ) +
  scale_y_continuous(
    limits = c(0, 250),
    breaks = seq(0, 250, by = 50)
  ) +

  ggtitle("Volcano Plot: FOXP2_TSHZ2 vs All Other Interneurons")


# Save the plot for each comparison to a PDF file
pdf(paste0("VolcanoPlot_FOXP2_TSHZ2_vs_AllOtherInterneurons_withinBat_FindMarkers_DEGs.pdf"))
print(plot)
dev.off()
###########################################################################
### GO term analysis - FOXP2_TSHZ2_BAT
library(gprofiler2)
library(clusterProfiler)
#install.packages("DT")
library(DT)

#### gprofiler2
# get upreg genes
foxp2_tshz2_upreg_genes <- a_foxp2_tshz2[a_foxp2_tshz2$diffexpressed == "UP",]

# Perform GO enrichment analysis using gprofiler2
foxp2_tshz2_gprofiler2_results <- gost(query = rownames(foxp2_tshz2_upreg_genes), custom_bg = rownames(counts_data), organism = "hsapiens", user_threshold = 0.05, sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG"))

# View the results
summary(foxp2_tshz2_gprofiler2_results)

# To view the GO term enrichment details
head(foxp2_tshz2_gprofiler2_results$result)

###@@@ 
# Convert list-type columns to strings
foxp2_tshz2_gprofiler2_results$result$genes <- sapply(foxp2_tshz2_gprofiler2_results$result$genes, 
                                                      function(x) paste(x, collapse = ", "))

# Save the results to a CSV file
write.csv(foxp2_tshz2_gprofiler2_results$result, 
          file = "foxp2_tshz2_gprofiler2_results.csv", 
          row.names = FALSE)



### EnrichGO
# get the upregulated genes
foxp2_tshz2_upreg_genes <- rownames(a_foxp2_tshz2[a_foxp2_tshz2$diffexpressed == "UP",])

# perform GO-term enrichment
go1_foxp2_tshz2 <- enrichGO(foxp2_tshz2_upreg_genes, OrgDb = "org.Hs.eg.db", keyType= 'SYMBOL', ont = "BP", maxGSSize = 5000, universe = rownames(counts_data))

#write
go1_foxp2_tshz2_summary <- data.frame(go1_foxp2_tshz2) %>%  write.table(file="EnrichGO_FOXP2_TSHZ2_withinBat_AcrossInterneuron_DEG_up.csv", sep=",", quote=F, row.names = F)

# get the downregulated genes
foxp2_tshz2_downreg_genes <- rownames(a_foxp2_tshz2[a_foxp2_tshz2$diffexpressed == "DOWN",])

# perform GO-term enrichment
go2_foxp2_tshz2 <- enrichGO(foxp2_tshz2_downreg_genes, OrgDb = "org.Hs.eg.db", keyType= 'SYMBOL', ont = "BP", maxGSSize = 5000, universe = rownames(counts_data))

#write
go2_foxp2_tshz2_summary <- data.frame(go2_foxp2_tshz2) %>%  write.table(file="EnrichGO_FOXP2_TSHZ2_withinBat_AcrossInterneuron_DEG_down.csv", sep=",", quote=F, row.names = F)

### GO term analysis - LMO3_BAT
library(gprofiler2)
#install.packages("DT")
library(DT)

# get the upregulated genes
lmo3_upreg_genes <- rownames(a[a$diffexpressed == "UP",])

# perform GO-term enrichment
go1 <- enrichGO(lmo3_upreg_genes, OrgDb = "org.Hs.eg.db", keyType= 'SYMBOL', ont = "BP", maxGSSize = 5000)

#write
go1_summary <- data.frame(go1) %>%  write.table(file="LMO3_withinBat_AcrossInterneuron_DEG_up.csv", sep=",", quote=F)

# get the downregulated genes
lmo3_downreg_genes <- rownames(bat_lmo3_vs_interneu_markers[bat_lmo3_vs_interneu_markers$diffexpressed == "DOWN",])

# perform GO-term enrichment
go2 <- enrichGO(lmo3_downreg_genes, OrgDb = "org.Hs.eg.db", keyType= 'SYMBOL', ont = "BP", maxGSSize = 5000)

#write
go2_summary <- data.frame(go2) %>%  write.table(file="LMO3_withinBat_AcrossInterneuron_DEG_down.csv", sep=",", quote=F)
#############################################################################
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
 [1] ggrepel_0.9.6               glmGamPoi_1.10.2           
 [3] edgeR_3.40.2                limma_3.54.2               
 [5] ashr_2.2-63                 apeglm_1.20.0              
 [7] SingleCellExperiment_1.20.1 DESeq2_1.38.3              
 [9] SummarizedExperiment_1.28.0 Biobase_2.58.0             
[11] MatrixGenerics_1.10.0       matrixStats_0.63.0         
[13] GenomicRanges_1.50.0        GenomeInfoDb_1.34.8        
[15] IRanges_2.32.0              S4Vectors_0.36.2           
[17] BiocGenerics_0.44.0         GeneOverlap_1.34.0         
[19] harmony_1.0.1               Rcpp_1.0.10                
[21] data.table_1.14.8           rio_1.0.1                  
[23] reshape2_1.4.4              ggpubr_0.6.0               
[25] lubridate_1.9.2             forcats_1.0.0              
[27] stringr_1.5.0               dplyr_1.1.4                
[29] purrr_1.0.1                 readr_2.1.4                
[31] tidyr_1.3.0                 tibble_3.2.0               
[33] tidyverse_2.0.0             plyr_1.8.8                 
[35] ggplot2_3.4.4               Matrix_1.6-1               
[37] SeuratObject_4.1.4          Seurat_4.4.0               
[39] openxlsx_4.2.7.1           

loaded via a namespace (and not attached):
  [1] utf8_1.2.4             spatstat.explore_3.0-6 reticulate_1.26       
  [4] R.utils_2.12.2         tidyselect_1.2.0       RSQLite_2.3.1         
  [7] AnnotationDbi_1.60.2   htmlwidgets_1.6.2      grid_4.2.3            
 [10] BiocParallel_1.32.5    Rtsne_0.16             munsell_0.5.0         
 [13] codetools_0.2-19       ica_1.0-3              future_1.33.1         
 [16] miniUI_0.1.1.1         withr_2.5.2            spatstat.random_3.1-3 
 [19] colorspace_2.1-0       progressr_0.13.0       ROCR_1.0-11           
 [22] ggsignif_0.6.4         tensor_1.5             listenv_0.9.0         
 [25] labeling_0.4.3         bbmle_1.0.25.1         GenomeInfoDbData_1.2.9
 [28] mixsqp_0.3-54          polyclip_1.10-4        farver_2.1.1          
 [31] bit64_4.0.5            coda_0.19-4            parallelly_1.35.0     
 [34] vctrs_0.6.5            generics_0.1.3         timechange_0.2.0      
 [37] R6_2.5.1               invgamma_1.1           locfit_1.5-9.8        
 [40] bitops_1.0-7           spatstat.utils_3.0-1   cachem_1.0.7          
 [43] DelayedArray_0.24.0    promises_1.2.0.1       scales_1.2.1          
 [46] gtable_0.3.3           globals_0.16.2         goftest_1.2-3         
 [49] rlang_1.1.0            splines_4.2.3          rstatix_0.7.2         
 [52] lazyeval_0.2.2         spatstat.geom_3.0-6    broom_1.0.3           
 [55] abind_1.4-5            backports_1.4.1        httpuv_1.6.9          
 [58] tools_4.2.3            tcltk_4.2.3            ellipsis_0.3.2        
 [61] gplots_3.1.3           RColorBrewer_1.1-3     ggridges_0.5.4        
 [64] zlibbioc_1.44.0        RCurl_1.98-1.12        deldir_1.0-6          
 [67] pbapply_1.7-0          cowplot_1.1.1          zoo_1.8-11            
 [70] cluster_2.1.4          magrittr_2.0.3         scattermore_1.2       
 [73] lmtest_0.9-40          RANN_2.6.1             truncnorm_1.0-9       
 [76] SQUAREM_2021.1         mvtnorm_1.2-3          fitdistrplus_1.1-8    
 [79] hms_1.1.2              patchwork_1.1.2        mime_0.12             
 [82] xtable_1.8-4           XML_3.99-0.14          emdbook_1.3.13        
 [85] gridExtra_2.3          bdsmatrix_1.3-7        compiler_4.2.3        
 [88] KernSmooth_2.23-20     crayon_1.5.2           R.oo_1.25.0           
 [91] htmltools_0.5.4        later_1.3.0            tzdb_0.3.0            
 [94] geneplotter_1.76.0     DBI_1.1.3              MASS_7.3-58.3         
 [97] car_3.1-2              cli_3.6.2              R.methodsS3_1.8.2     
[100] parallel_4.2.3         igraph_1.5.1           pkgconfig_2.0.3       
[103] numDeriv_2016.8-1.1    sp_1.6-0               plotly_4.10.1         
[106] spatstat.sparse_3.0-0  annotate_1.76.0        XVector_0.38.0        
[109] digest_0.6.31          sctransform_0.4.0      RcppAnnoy_0.0.20      
[112] spatstat.data_3.0-0    Biostrings_2.66.0      leiden_0.4.3          
[115] uwot_0.1.14            shiny_1.7.4            gtools_3.9.4          
[118] lifecycle_1.0.3        nlme_3.1-162           jsonlite_1.8.4        
[121] carData_3.0-5          viridisLite_0.4.1      fansi_1.0.6           
[124] pillar_1.9.0           lattice_0.21-8         KEGGREST_1.38.0       
[127] fastmap_1.1.1          httr_1.4.5             survival_3.5-3        
[130] glue_1.6.2             zip_2.2.2              png_0.1-8             
[133] bit_4.0.5              stringi_1.7.12         blob_1.2.3            
[136] caTools_1.18.2         memoise_2.0.1          irlba_2.3.5.1         
[139] future.apply_1.10.0  
