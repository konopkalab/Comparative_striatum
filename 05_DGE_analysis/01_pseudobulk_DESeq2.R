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
library(WGCNA)
library(flashClust)

# path to save files
outdir = "/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/DEG_analysis/DESeq2/"

# load data
file_outdir = "/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/seu_objs/"
merged_pb_meta = readRDS(file = paste0(file_outdir,"primates_merged_pb_meta_PostgliaCleanup.RDS"))
merged_sub_seu_pb_count = readRDS(file = paste0(file_outdir,"primates_merged_sub_seu_pb_count_PostgliaCleanup.RDS"))
merged_meta_df = readRDS(file = paste0(file_outdir,"primates_merged_meta_df_PostgliaCleanup.RDS"))
  
keep_genes_list = list()
### run Variance analysis and DESeq2 for each Tissue-cellType
for (t in names(merged_pb_meta)){
    for (ctype in names(merged_pb_meta[[t]])){
      ## filter the genes which have at least 3 copies in >= 75% of the samples
      # count mat
      mat = merged_sub_seu_pb_count[[t]][[ctype]]
      # Extract species and tissue from columns
      parts <- strsplit(colnames(mat), "_")
      species <- sapply(parts, function(x) tail(x, 1))
      tissue  <- sapply(parts, function(x) x[length(x)-1])

      group_ids <- paste(species, tissue, sep = "_")

      # Split sample names by species–tissue
      groups <- split(colnames(mat), group_ids)

      # For each group, compute whether genes have >=3 counts in >=75% of samples
      keep_genes_list <- lapply(groups, function(samps) {
          submat <- mat[, samps, drop = FALSE]

          gene_ok <- submat >= 3
          prop_ok <- rowMeans(gene_ok)

          names(prop_ok)[prop_ok >= 0.75]
      })
      keep_genes <- Reduce(intersect, keep_genes_list)

       # normalize the counts
      norm_aggr_counts = t(edgeR::cpm(merged_sub_seu_pb_count[[t]][[ctype]][keep_genes,], log = TRUE))
      
      counts_for_deseq2 = merged_sub_seu_pb_count[[t]][[ctype]][keep_genes,]
      # update the rownames
      rownames(merged_pb_meta[[t]][[ctype]]) = merged_pb_meta[[t]][[ctype]]$Sample
      ####### look how var.s are explained:
      form = ~ Sex + log10_sum_ncount_RNA + Humanized_age + Species

      varPart <- fitExtractVarPartModel(t(norm_aggr_counts), form, merged_pb_meta[[t]][[ctype]])

      # sort variables (i.e. columns) by median fraction of variance explained
      vp <- sortCols(varPart)

      # Bar plot of variance fractions for the first 10 genes
      plotPercentBars(vp[1:10, ])
      # violin plot of contribution of each variable to total variance
      pdf(paste0(outdir, t, "_", ctype, "_vars_explained.pdf"))
      print(plotVarPart(vp))
      dev.off()

      # Create DESeq2 object        
      dds <- DESeqDataSetFromMatrix(counts_for_deseq2, 
		              colData = merged_pb_meta[[t]][[ctype]], 
		              design = ~ Sex + log10_sum_ncount_RNA + Humanized_age + Species)

      ##### QC
      # Transform counts for data visualization
      rld <- rlog(dds, blind=TRUE)
      #### if you'd like to check the other PCs
      rld_mat <- assay(rld)
      pca <- prcomp(t(rld_mat))

       # Create data frame with metadata and PC3 and PC4 values for input to ggplot
       #df <- cbind(merged_pb_meta[[t]][[ctype]], pca$x)
       #ggplot(df) + geom_point(aes(x=PC3, y=PC4, color = Humanized_age))
	# Plot PCA
	p <- DESeq2::plotPCA(rld, ntop = 500, intgroup = "Species")

	# save to file
	ggsave(paste0("/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/DEG_analysis/DESeq2/", t, "_", ctype, "_Species_DESeq2_PCA_Humanized_age_asCoVar.png"), plot = p, width = 6, height = 5, dpi = 300)
	# Save the plot with sample names
	pcaplot = p + geom_text_repel(aes(label = rownames(colData(rld))), size = 3)
	ggsave(paste0("/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/DEG_analysis/DESeq2/", t, "_", ctype, "_Species_DESeq2_PCA_withSampleNames_Humanized_age_asCoVar.png"), plot = pcaplot, width = 20, height = 8, dpi = 300)
	##### Run DESeq2 differential expression analysis
	dds <- DESeq(dds)

	# Plot dispersion estimates
	plotDispEsts(dds)

	# Check the coefficients for the comparison
	resultsNames(dds)
	# [1] "Intercept"                 "Sex_Male_vs_Female"       
#[3] "log10_sum_ncount_RNA"      "Humanized_age"            
#[5] "Species_Chimp_vs_Human"    "Species_Macaque_vs_Human" 
#[7] "Species_Marmoset_vs_Human"

	# drop the species that did not survive the cutoff
	merged_pb_meta[[t]][[ctype]]$Species = factor(merged_pb_meta[[t]][[ctype]]$Species)
	merged_pb_meta[[t]][[ctype]]$Species = droplevels(merged_pb_meta[[t]][[ctype]]$Species)
	# Get all species levels
	species_levels <- levels(merged_pb_meta[[t]][[ctype]]$Species)

	# Generate all unique species pairs
	pairs <- combn(species_levels, 2, simplify = FALSE)
	
	# Directory to save results
	output_dir <- "/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/DEG_analysis/DESeq2/log10_sum_ncount_RNA_Humanized_age_covar_"

# Loop over all pairs
	for (pair in pairs) {	
	  sp1 <- pair[1]
	  sp2 <- pair[2]

	  message("Running comparison: ", sp1, " vs ", sp2)
	  
	  # Get DESeq2 results for this pair (Sex is still covariate)
	  res <- results(dds, contrast = c("Species", sp1, sp2), alpha = 0.05)


	  # Shrink log2 fold changes with apeglm
	  # interpret results as: log2FoldChange = log2(expression_in_sp1 / expression_in_sp2) 
	  res <- lfcShrink(dds, contrast = c("Species", sp1, sp2), res = res, type = "ashr")
	  
	  # Remove rows with NA p-values
	  res <- res[!is.na(res$padj), ]
	  
	  # Convert to tibble
	  res_tbl <- res %>%
	    data.frame() %>%
	    rownames_to_column(var = "gene") %>%
	    as_tibble() %>%
	    arrange(padj)
	  
	  # Add DEG status
	  res_tbl <- res_tbl %>%
	    mutate(diffexpressed = case_when(
	      log2FoldChange > 0.6 & padj < 0.05  ~ "UP",
	      log2FoldChange < -0.6 & padj < 0.05 ~ "DOWN",
	      TRUE ~ "NO"
	    ))
	  
	  # Significant genes for volcano labels
	  labeled <- filter(res_tbl, abs(log2FoldChange) > 0.7 & padj < 0.05)
	  
	  # Save CSV
	  write.csv(res_tbl,
		    file = paste0(output_dir, t, "_", ctype, "_", sp1, "_vs_", sp2, "_all_genes.csv"),
		    quote = FALSE, row.names = FALSE)
	  
	  # Volcano plot
	  plot <- ggplot(res_tbl, aes(x = log2FoldChange, y = -log10(padj), 
			              col = diffexpressed, label = gene)) +
	    theme_bw() +
	    geom_point(alpha = 0.5) +
	    geom_vline(xintercept = c(-0.6, 0.6), col = "grey", linetype = "longdash") +
	    geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "longdash") +
	    scale_colour_manual(name = "DEG", values = c("blue", "black", "darkgreen")) +
	    geom_text_repel(data = labeled, aes(label = gene), nudge_y = 0.2,
			    box.padding = 0.3, point.padding = 0.3) +
	    ggtitle(paste("Volcano Plot:", t, "_", ctype, sp1, "vs", sp2))
	  
	  # Save PDF
	  pdf(paste0(output_dir, "Volcano_", t, "_",ctype, "_", sp1, "_vs_", sp2, ".pdf"))
	  print(plot)
	  dev.off()
	}
}
}

############################################################
############################################################
sessionInfo()
R version 4.2.3 (2023-03-15)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux 9.4 (Plow)

Matrix products: default
BLAS/LAPACK: /project/backup_modules/apps/rstudio-desktop/2022.12.0/lib/libopenblasp-r0.3.21.so

locale:
 [1] LC_CTYPE=C.UTF-8           LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] ComplexHeatmap_2.14.0       ComplexUpset_1.3.6         
 [3] DT_0.31                     clusterProfiler_4.6.2      
 [5] gprofiler2_0.2.3            UpSetR_1.4.0               
 [7] flashClust_1.01-2           WGCNA_1.73                 
 [9] fastcluster_1.3.0           dynamicTreeCut_1.63-1      
[11] curl_5.2.2                  Matrix.utils_0.9.8         
[13] variancePartition_1.28.9    BiocParallel_1.32.5        
[15] ashr_2.2-63                 edgeR_3.40.2               
[17] limma_3.54.0                ggrepel_0.9.3              
[19] DESeq2_1.38.3               harmony_1.1.0              
[21] Rcpp_1.0.10                 data.table_1.14.8          
[23] rio_1.0.1                   reshape2_1.4.4             
[25] ggpubr_0.6.0                lubridate_1.9.3            
[27] forcats_1.0.0               stringr_1.5.0              
[29] purrr_1.0.1                 readr_2.1.4                
[31] tidyr_1.3.0                 tibble_3.2.1               
[33] tidyverse_2.0.0             plyr_1.8.9                 
[35] ggplot2_3.4.4               Matrix_1.6-5               
[37] dplyr_1.1.4                 DropletQC_0.0.0.9000       
[39] DropletUtils_1.18.1         SingleCellExperiment_1.20.0
[41] SummarizedExperiment_1.28.0 Biobase_2.58.0             
[43] GenomicRanges_1.50.0        GenomeInfoDb_1.34.9        
[45] IRanges_2.32.0              S4Vectors_0.36.0           
[47] BiocGenerics_0.44.0         MatrixGenerics_1.10.0      
[49] matrixStats_0.63.0          rhdf5_2.42.1               
[51] BPCells_0.3.1               Seurat_5.3.0               
[53] SeuratObject_5.1.0          sp_1.6-0                   
[55] patchwork_1.3.0.9000       

loaded via a namespace (and not attached):
  [1] Hmisc_5.1-1               ica_1.0-3                
  [3] foreach_1.5.2             lmtest_0.9-40            
  [5] crayon_1.5.2              rbibutils_2.2.16         
  [7] MASS_7.3-58.3             rhdf5filters_1.10.1      
  [9] nlme_3.1-162              backports_1.4.1          
 [11] impute_1.72.3             GOSemSim_2.24.0          
 [13] rlang_1.1.4               XVector_0.38.0           
 [15] HDO.db_0.99.1             ROCR_1.0-11              
 [17] irlba_2.3.5.1             nloptr_2.0.3             
 [19] rjson_0.2.21              bit64_4.0.5              
 [21] glue_1.6.2                mixsqp_0.3-54            
 [23] sctransform_0.4.1         pbkrtest_0.5.2           
 [25] parallel_4.2.3            spatstat.sparse_3.0-0    
 [27] AnnotationDbi_1.60.2      dotCall64_1.1-0          
 [29] DOSE_3.24.2               spatstat.geom_3.0-6      
 [31] tidyselect_1.2.0          fitdistrplus_1.1-8       
 [33] XML_3.99-0.14             zoo_1.8-11               
 [35] xtable_1.8-4              RcppHNSW_0.4.1           
 [37] magrittr_2.0.3            evaluate_0.20            
 [39] Rdpack_2.6                scuttle_1.8.0            
 [41] cli_3.6.2                 zlibbioc_1.44.0          
 [43] rstudioapi_0.14           miniUI_0.1.1.1           
 [45] rpart_4.1.19              fastmatch_1.1-4          
 [47] aod_1.3.3                 EnvStats_2.8.1           
 [49] treeio_1.22.0             fastDummies_1.7.3        
 [51] shiny_1.7.4               xfun_0.47                
 [53] clue_0.3-64               gson_0.1.0               
 [55] cluster_2.1.4             caTools_1.18.2           
 [57] tidygraph_1.2.3           KEGGREST_1.38.0          
 [59] ape_5.7-1                 listenv_0.9.0            
 [61] Biostrings_2.66.0         png_0.1-8                
 [63] future_1.33.1             withr_2.5.2              
 [65] bitops_1.0-7              ggforce_0.4.1            
 [67] dqrng_0.3.0               pillar_1.9.0             
 [69] GlobalOptions_0.1.2       gplots_3.1.3             
 [71] cachem_1.0.7              GetoptLong_1.0.5         
 [73] DelayedMatrixStats_1.20.0 vctrs_0.6.5              
 [75] ellipsis_0.3.2            generics_0.1.3           
 [77] tools_4.2.3               foreign_0.8-85           
 [79] remaCor_0.0.16            munsell_0.5.0            
 [81] tweenr_2.0.2              fgsea_1.24.0             
 [83] DelayedArray_0.24.0       fastmap_1.1.1            
 [85] compiler_4.2.3            abind_1.4-5              
 [87] httpuv_1.6.9              plotly_4.10.1            
 [89] GenomeInfoDbData_1.2.9    gridExtra_2.3            
 [91] lattice_0.21-8            deldir_1.0-6             
 [93] utf8_1.2.4                later_1.3.0              
 [95] jsonlite_1.8.4            scales_1.3.0             
 [97] tidytree_0.4.5            pbapply_1.7-0            
 [99] carData_3.0-5             sparseMatrixStats_1.10.0 
[101] lazyeval_0.2.2            promises_1.2.0.1         
[103] car_3.1-2                 doParallel_1.0.17        
[105] R.utils_2.12.2            goftest_1.2-3            
[107] spatstat.utils_3.1-4      reticulate_1.42.0        
[109] checkmate_2.3.0           rmarkdown_2.28           
[111] cowplot_1.1.1             Rtsne_0.16               
[113] downloader_0.4.1          uwot_0.1.14              
[115] igraph_1.4.2              HDF5Array_1.26.0         
[117] survival_3.5-3            SQUAREM_2021.1           
[119] htmltools_0.5.8.1         memoise_2.0.1            
[121] locfit_1.5-9.8            graphlayouts_1.0.0       
[123] viridisLite_0.4.1         digest_0.6.31            
[125] RhpcBLASctl_0.23-42       mime_0.12                
[127] spam_2.10-0               RSQLite_2.3.1            
[129] yulab.utils_0.0.8         future.apply_1.10.0      
[131] blob_1.2.3                R.oo_1.25.0              
[133] preprocessCore_1.60.2     splines_4.2.3            
[135] Formula_1.2-5             Rhdf5lib_1.20.0          
[137] RCurl_1.98-1.12           broom_1.0.3              
[139] hms_1.1.2                 colorspace_2.1-0         
[141] base64enc_0.1-3           shape_1.4.6              
[143] aplot_0.2.0               nnet_7.3-18              
[145] mclust_6.0.0              RANN_2.6.1               
[147] circlize_0.4.15           mvtnorm_1.2-3            
[149] enrichplot_1.18.4         fansi_1.0.6              
[151] tzdb_0.3.0                truncnorm_1.0-9          
[153] parallelly_1.35.0         R6_2.5.1                 
[155] ggridges_0.5.4            lifecycle_1.0.3          
[157] ggsignif_0.6.4            minqa_1.2.6              
[159] qvalue_2.30.0             RcppAnnoy_0.0.20         
[161] RColorBrewer_1.1-3        iterators_1.0.14         
[163] spatstat.explore_3.0-6    htmlwidgets_1.6.2        
[165] beachmat_2.14.0           polyclip_1.10-4          
[167] shadowtext_0.1.2          gridGraphics_0.5-1       
[169] timechange_0.2.0          globals_0.16.2           
[171] htmlTable_2.4.2           spatstat.random_3.1-3    
[173] progressr_0.13.0          codetools_0.2-19         
[175] invgamma_1.1              GO.db_3.16.0             
[177] gtools_3.9.4              prettyunits_1.1.1        
[179] RSpectra_0.16-1           R.methodsS3_1.8.2        
[181] gtable_0.3.3              DBI_1.1.3                
[183] ggfun_0.1.2               tensor_1.5               
[185] httr_1.4.5                KernSmooth_2.23-20       
[187] stringi_1.7.12            progress_1.2.2           
[189] farver_2.1.1              annotate_1.76.0          
[191] viridis_0.6.2             ggtree_3.6.2             
[193] boot_1.3-28.1             grr_0.9.5                
[195] lme4_1.1-35.1             geneplotter_1.76.0       
[197] ggplotify_0.1.2           scattermore_1.2          
[199] bit_4.0.5                 scatterpie_0.2.1         
[201] spatstat.data_3.0-0       ggraph_2.1.0             
[203] pkgconfig_2.0.3           rstatix_0.7.2            
[205] knitr_1.48  
