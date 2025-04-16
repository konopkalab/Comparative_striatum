# load the modules
source ~/load_modules.sh
module load hdf5/gcc/1.12.0

# change the dir
cd /home2/gkonop/project/02_MATRIX_FOR_CELLBENDER/KRIENEN_FERRET

# start R session 
R
##########################################################
### Remove feature-cell barcode count matrix from
# count matrix (cellranger count output)

# packages
library(rhdf5)
library(Seurat)
library(Matrix)


# load data
samps = c("SRR11921037", "SRR11921038")
folddir = c('/home2/gkonop/project/00_BAM_DOWNLOADED/KRIENEN_2020/FERRET')
folder_name = c('repaired_krienen_ferret_rxn1','repaired_krienen_ferret_rxn2')
outdir = '/home2/gkonop/project/02_MATRIX_FOR_CELLBENDER/KRIENEN_FERRET'


####
## Loop through the samples
####

for(i in 1:length(samps)){
	# Get the Gene expression matrix
	# read the matrix (Gene Expression)
	h <- Read10X_h5(paste0(folddir, '/', samps[i],'/','cellranger_count/', folder_name[i], '/', 'outs/raw_feature_bc_matrix.h5'))
	 # Select top 80,000 cells by expression
  # Get the top 80,000 cells which express the most amount of genes
  a <- sort(colSums(h), decreasing = TRUE)[1:80000]
  
  # Extract the count matrix for the selected cells
  b <- h[, names(a)]
  
  # Extract barcodes corresponding to the selected cells
  barcodes <- names(a)
  # Extract genes and barcodes
  genes <- rownames(b)
  
  # create the output dir
  dir.create(paste0(outdir, '/', samps[i]))
  dir.create(paste0(outdir, '/', samps[i], '/snRNA_count_data_top80000'))
 
  # Write genes and barcodes to files
  write.table(genes, file = paste0(outdir, '/', samps[i], '/snRNA_count_data_top80000/genes.tsv'), sep = '\t', col.names = FALSE, row.names = genes, quote = FALSE)
  write.table(barcodes, file = paste0(outdir, '/', samps[i], '/snRNA_count_data_top80000/barcodes.tsv'), sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)

  # Convert matrix to Matrix Market format and write to file
  b_sparse <- as(b, "sparseMatrix")
  writeMM(b_sparse, file = paste0(outdir, '/', samps[i], '/snRNA_count_data_top80000/matrix.mtx'))
}

###### Packages loaded
# sessionInfo()
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
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] Matrix_1.6-1       SeuratObject_4.1.4 Seurat_4.4.0       rhdf5_2.42.1      

loaded via a namespace (and not attached):
  [1] nlme_3.1-162           matrixStats_0.63.0     spatstat.sparse_3.0-0 
  [4] bit64_4.0.5            RcppAnnoy_0.0.20       RColorBrewer_1.1-3    
  [7] httr_1.4.5             sctransform_0.4.0      tools_4.2.3           
 [10] utf8_1.2.4             R6_2.5.1               irlba_2.3.5.1         
 [13] KernSmooth_2.23-20     uwot_0.1.14            lazyeval_0.2.2        
 [16] colorspace_2.1-0       rhdf5filters_1.10.1    sp_1.6-0              
 [19] tidyselect_1.2.0       gridExtra_2.3          bit_4.0.5             
 [22] compiler_4.2.3         progressr_0.13.0       cli_3.6.2             
 [25] hdf5r_1.3.8            spatstat.explore_3.0-6 plotly_4.10.1         
 [28] scales_1.2.1           lmtest_0.9-40          spatstat.data_3.0-0   
 [31] ggridges_0.5.4         pbapply_1.7-0          goftest_1.2-3         
 [34] stringr_1.5.0          digest_0.6.31          spatstat.utils_3.0-1  
 [37] pkgconfig_2.0.3        htmltools_0.5.4        parallelly_1.35.0     
 [40] fastmap_1.1.1          htmlwidgets_1.6.2      rlang_1.1.0           
 [43] shiny_1.7.4            generics_0.1.3         zoo_1.8-11            
 [46] jsonlite_1.8.4         spatstat.random_3.1-3  ica_1.0-3             
 [49] dplyr_1.1.4            magrittr_2.0.3         patchwork_1.1.2       
 [52] Rcpp_1.0.10            munsell_0.5.0          Rhdf5lib_1.20.0       
 [55] fansi_1.0.6            abind_1.4-5            reticulate_1.26       
 [58] lifecycle_1.0.3        stringi_1.7.12         MASS_7.3-58.3         
 [61] Rtsne_0.16             plyr_1.8.8             grid_4.2.3            
 [64] parallel_4.2.3         listenv_0.9.0          promises_1.2.0.1      
 [67] ggrepel_0.9.3          deldir_1.0-6           miniUI_0.1.1.1        
 [70] lattice_0.21-8         cowplot_1.1.1          splines_4.2.3         
 [73] tensor_1.5             pillar_1.9.0           igraph_1.5.1          
 [76] spatstat.geom_3.0-6    future.apply_1.10.0    reshape2_1.4.4        
 [79] codetools_0.2-19       leiden_0.4.3           glue_1.6.2            
 [82] BiocManager_1.30.20    data.table_1.14.8      png_0.1-8             
 [85] vctrs_0.6.5            httpuv_1.6.9           polyclip_1.10-4       
 [88] gtable_0.3.3           RANN_2.6.1             purrr_1.0.1           
 [91] tidyr_1.3.0            scattermore_1.2        future_1.33.1         
 [94] ggplot2_3.4.4          mime_0.12              xtable_1.8-4          
 [97] later_1.3.0            survival_3.5-3         viridisLite_0.4.1     
[100] tibble_3.2.0           cluster_2.1.4          globals_0.16.2        
[103] fitdistrplus_1.1-8     ellipsis_0.3.2         ROCR_1.0-11   
