# load packages
library(ggplot2)
library(Seurat)
library(dplyr)
library(Matrix)
library(plyr)
library(tidyr)
library(reshape2)
library(data.table)
library(ggrepel)
library(ggpubr)
library(WGCNA)
library(DESeq2)
library(fastDummies)
set.seed(1234)
source("/bmapfs/archive/konopkalab1/SCRIPTS_pr/SCRIPTS/utility_functions.R")
library(curl)
#conda install bioconda::r-wgcna
#BiocManager::install('WGCNA')
#install.packages("fastDummies")
library(scater)
options(warn = 0)

# path to save files
outdir = "/bmapfs/archive/konopkalab1/Gozde_workdir/projects/comparative_striatum/revision/UTSW_Bat_neurons/"

## load the data
str_seu = readRDS("/bmapfs/archive/konopkalab1/Gozde_workdir/projects/comparative_striatum/seu_objs/merged_clean_AllCellTypes_AllTissues_seu.RDS")

# subset bat interneurons
bat_interneu = subset(str_seu, subset = Species == "Bat" & CellType == "Non_SPN")

# bat putamen subset
bat_putamen_interneu = subset(bat_interneu, subset = Tissue == "Putamen")

# merge some of the inteneuron subtypes: PDGFD interneurons --> PDGFD_merged; TAC3 & TH --> TAC3_TH, CCK_VIP- & CCK_VIP+ --> CCK
bat_putamen_interneu$wgcna_annot = bat_putamen_interneu$detailed_annot
bat_putamen_interneu$wgcna_annot = gsub(pattern = "^(PDGFD|PDGFD_PTHLH_PVALB-|PDGFD_PTHLH_PVALB\\+)$",   replacement = "PDGFD_merged", x = bat_putamen_interneu$wgcna_annot) 
bat_putamen_interneu$wgcna_annot = gsub(pattern = "^(TAC3|TH)$",  replacement = "TAC3_TH", x = bat_putamen_interneu$wgcna_annot) 
bat_putamen_interneu$wgcna_annot = gsub(pattern = "^(CCK_VIP-|CCK_VIP\\+)$",   replacement = "CCK", x = bat_putamen_interneu$wgcna_annot)

#check 
table(bat_putamen_interneu$wgcna_annot)
## generate all the matrices with Sample_CellType instead of Sample bcs we're interested in comparing all these celltypes with each other (and not Species within each Sample)
## generate pseudobulk count mat by filtering genes: at least 1 copy in at least 3 Sample_`wgcna_annot` for each seu object (i.e. bat putamen interneuron seu)
# generate pseudobulk count matrix at the Sample_`wgcna_annot` granularity
bat_putamen_interneu_pb_count <- pseudobulk_sample_cellType_v5(
seurObj = bat_putamen_interneu,
sample_col = "wgcna_annot",
n_copy = 1,
sample_number = 3)

#save
saveRDS(bat_putamen_interneu_pb_count, file = "/bmapfs/archive/konopkalab1/Gozde_workdir/projects/comparative_striatum/seu_objs/BMAP_pseudobulk_count_matrices_Bat_Putamen_Samples_interneuron_subTypes.RDS")
#####################################################
# normalize the counts and transpose the mat
dge <- DGEList(counts = bat_putamen_interneu_pb_count)
cpm_norm <- edgeR::cpm(dge, log = FALSE)   # library-size normalized counts (not log yet)
# apply variance stabilizing transformation
norm_aggr_counts <- t(log2(cpm_norm + 1))

# find the mean UMI counts for each sample in the species-wgcna_annot pair
bat_putamen_interneu$Sample_wgcna_annot = paste0(bat_putamen_interneu$Sample, "_", bat_putamen_interneu$wgcna_annot)

meta_Sample_wgcna_annot <- bat_putamen_interneu@meta.data %>%
dplyr::group_by(Sample_wgcna_annot) %>%
dplyr::summarise(
log10_sum_ncount_RNA = log10(sum(nCount_RNA, na.rm = TRUE)),
n = dplyr::n(),
.groups = "drop")
meta_Sample_wgcna_annot = as.data.frame(meta_Sample_wgcna_annot)

# set rownames
rownames(meta_Sample_wgcna_annot) = meta_Sample_wgcna_annot$Sample_wgcna_annot

# check if the rownames in metadata matches the rownames of the aggregated count matrix for each celltype
(match_ok <- all(rownames(norm_aggr_counts) == rownames(meta_Sample_wgcna_annot)))
####################################################
# in aggr_counts data.frame, the columns should represent genes and rows represent samples
gsg <-goodSamplesGenes(norm_aggr_counts)

summary(gsg)
gsg$allOK

# Genes that failed the check
badGenes <- colnames(norm_aggr_counts)[!gsg$goodGenes]
length(badGenes)
head(badGenes)

if (!gsg$allOK)
{
if (sum(!gsg$goodGenes)>0) 
printFlush(paste("Removing genes:", paste(names(norm_aggr_counts)[!gsg$goodGenes], collapse = ", "))); #Identifies and prints outlier genes
if (sum(!gsg$goodSamples)>0)
printFlush(paste("Removing samples:", paste(rownames(norm_aggr_counts)[!gsg$goodSamples], collapse = ", "))); #Identifies and prints oulier samples
norm_aggr_counts <- norm_aggr_counts[gsg$goodSamples == TRUE, gsg$goodGenes == TRUE] # Removes the offending genes and samples from the data
}

sampleTree <- hclust(dist(norm_aggr_counts), method = "average") #Clustering samples based on distance 

#Plotting the cluster dendrogram
pdf(paste0(outdir, "BMAP_BatPutamen_interneurons_sampleclustering_aggregated_interneu.pdf"), width = 30, height = 10)  

# Set graphical parameters
par(cex = 0.6)
par(mar = c(0,4,2,0))

# Plot dendrogram
plot(sampleTree,
     main = "Sample clustering to detect outliers",
     sub = "",
     xlab = "",
     cex.lab = 2,
     cex.axis = 2,
     cex.main = 2,
     cex = 2)

dev.off()

################ MODULE construction ######################
################  signed network  #########################
## Power analysis
datExpr = norm_aggr_counts
powers = c(seq(2,30,2))

sft=pickSoftThreshold(datExpr,powerVector=powers,verbose = 5, blockSize= 14000, networkType = "signed",RsquaredCut = 0.85) 

# plot the thresholds
pdf(paste0(outdir,"BMAP_SoftThresholdingPower_signed_batPutamen_interneurons.pdf"))
par(mfrow = c(1,2), mar=c(5.1,5.1,4.1,2.1));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red");
abline(h=0.5,col="red"); abline(h=0.8,col="blue");abline(h=0.9,col="black")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

(PWR=sft$powerEstimate) #14

# calculate module eigengenes block-wise from all genes
net = blockwiseModules(datExpr,
corType="bicor",
maxBlockSize = 15000,
networkType="signed",
minCoreKME = 0.4, 
minKMEtoStay = 0.5,
power=PWR, 
checkMissingData = TRUE,
minModuleSize=25,
nThreads=15,
TOMType = "signed",
TOMDenom = "mean",
deepSplit=0.5,
verbose=1,
mergeCutHeight=0.10,
reassignThreshold = 1e-10,
numericLabels=TRUE)

moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
MEsAutomatic=net$MEs
unique(moduleColorsAutomatic)
table(moduleColorsAutomatic)
write.table(moduleColorsAutomatic, paste0(outdir, "BMAP_batPutamen_interneurons_WGCNA_module_colors.txt",sep="\t",quote=F))
save(net,file=paste0(outdir,"BMAP_batPutamen_interneurons_WGCNA_net.RData"))

# Calculate module membership for each gene
KMEs<-signedKME(datExpr, net$MEs,corFnc = "bicor") # bicor: Calculate biweight midcorrelation efficiently for matrices.
kme=data.frame(colnames(datExpr), moduleColorsAutomatic, KMEs)
colnames(kme)[1]="Symbol"
rownames(kme)=NULL
write.table(kme,paste0(outdir,"BMAP_Gene-module_membership_batPutamen_interneurons.txt"),sep="\t",quote=F)
save(kme,file=paste0(outdir,"BMAP_batPutamen_interneurons_WGCNA_kme.RData"))
load(paste0(outdir,"BMAP_batPutamen_interneurons_WGCNA_kme.RData")) # loaded as kme

# Dat Trait Heatmap
nGenes = ncol(datExpr) 
nSamples = nrow(datExpr)

### find the mean UMI counts for each sample in the species-wgcna_annot pair
bat_putamen_interneu$Sample_wgcna_annot = paste0(bat_putamen_interneu$Sample, "_", bat_putamen_interneu$wgcna_annot)

# calc the log10 of mean UMI counts for each sample
meta_Sample_wgcna_annot <- bat_putamen_interneu@meta.data %>%
dplyr::group_by(Sample_wgcna_annot) %>%
dplyr::summarise(
log10_sum_ncount_RNA = log10(sum(nCount_RNA, na.rm = TRUE)),
n = dplyr::n(),
.groups = "drop")
meta_Sample_wgcna_annot = as.data.frame(meta_Sample_wgcna_annot)

# set rownames
rownames(meta_Sample_wgcna_annot) = meta_Sample_wgcna_annot$Sample_wgcna_annot

# set datTraits dataframe 
Samples <- rownames(norm_aggr_counts)
traitRows <- match(Samples, meta_Sample_wgcna_annot$Sample_wgcna_annot)
datTraits <- meta_Sample_wgcna_annot[traitRows,]
rownames(datTraits) <- meta_Sample_wgcna_annot[traitRows, 1]
datTraits$wgcna_annot = str_split_i(datTraits$Sample_wgcna_annot, "_Putamen_", 2)
datTraits$Sample = str_split_i(datTraits$Sample_wgcna_annot, "_Putamen_", 1)
# remove unnecessary cols
datTraits$Sample_wgcna_annot = NULL
datTraits_numeric <- fastDummies::dummy_cols(datTraits, remove_selected_columns = TRUE)
rownames(datTraits_numeric) = rownames(datTraits)

MEs0 = moduleEigengenes(datExpr,moduleLabelsAutomatic)$eigengenes
MEsTemp = MEs0
MEsTemp$MEgrey=NULL
# use determen correlation (spearman) between data traits and modules 
modTraitCor= cor(MEsTemp, datTraits_numeric,method="spearman")
write.table(modTraitCor,paste0(outdir,"BMAP_modTraitCor_batPutamen_interneurons.txt"),sep="\t",quote=F)
modTraitP = corPvalueStudent(modTraitCor, nSamples)
write.table(modTraitP,paste0(outdir,"BMAP_modTraitP_batPutamen_interneurons.txt"),sep="\t",quote=F)

# find the markers of all interneurons
Idents(bat_putamen_interneu) = bat_putamen_interneu$wgcna_annot
markers_bat_putamen_interneu = FindAllMarkers(bat_putamen_interneu, logfc.threshold = 0.5)

signif_pos_markers_bat_putamen_interneu = markers_bat_putamen_interneu[markers_bat_putamen_interneu$p_val_adj < 0.05& markers_bat_putamen_interneu$avg_log2FC > 0,]

table(signif_pos_markers_bat_putamen_interneu$cluster)

# generate a list of these markers
markersL = split(signif_pos_markers_bat_putamen_interneu, signif_pos_markers_bat_putamen_interneu$cluster) %>% lapply(., function(x){x$gene})
# generate a list of genes belonging to each module
kmeL = split(kme, kme$moduleColorsAutomatic) %>% lapply(., function(x){x$Symbol})

### run Fisher's exact test between the marker of each celltype and eigengenes of modules
geneOvEnr(gnL1 = markersL, gnL2 = kmeL, bcg = nrow(kme), cutoff = 0.01, plot = T, hg = 18, wd = 14, fn = paste0(outdir, "BMAP_Batputamen_interneurons_marker_module_gene_enrichment")) # plot can be set to false (F) to return a data frame instead of the plot.

# extract all genes
all_genes = kme$Symbol #[1] 11513 genes

# get the appropriate module genes
blue_genes = kmeL[["blue"]]

lmo3_genes_df = signif_pos_markers_bat_putamen_interneu[signif_pos_markers_bat_putamen_interneu$cluster == "LMO3_BAT",]
lmo3_genes_df_ordered =lmo3_genes_df[order(lmo3_genes_df$p_val_adj),] 

# check
pdf(paste0(outdir, file = "BMAP_DotPlot_WGCNA_BlueModulegenes.pdf"), width = 35, height =10)
DotPlot(bat_putamen_interneu, features = blue_genes, group.by = "wgcna_annot") + xlab('') +
ylab('') +
theme(axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()

# get the appropriate module genes
cyan_genes =kmeL[["cyan"]]

tshz2_genes_df = signif_pos_markers_bat_putamen_interneu[signif_pos_markers_bat_putamen_interneu$cluster == "FOXP2_TSHZ2_BAT",]
tshz2_genes_df_ordered =tshz2_genes_df[order(tshz2_genes_df$p_val_adj),] 

# check
pdf(paste0(outdir, file = "BMAP_DotPlot_WGCNA_CyanModulegenes.pdf"), width = 25, height =10)
DotPlot(bat_putamen_interneu, features = cyan_genes, group.by = "wgcna_annot") + xlab('') +
ylab('') +
theme(axis.text.y=element_text(size=20, face = 'bold'),
	legend.text = element_text(size=20, face = 'bold'),
	legend.title = element_text(size=20, face = 'bold')) +
rotate_x_text(45)
dev.off()
#############################################
### GO-term
lmo3_up = GOenrich(gns = blue_genes, uni = all_genes)

write.csv(lmo3_up, paste0(outdir,"BMAP_LMO3_UpGenes_WGCNA_within_Bat_Putamen_compared_to_other_interneurons_EnrichGO.csv"), row.names = FALSE)

### dotplot of GO-terms
# Top terms per comparison
top_go <- lmo3_up %>%
  group_by(GO) %>%
  mutate(minuslog10_fdr = -log10(p.adjust)) %>% # Calculate -log10(p.adjust) for each group
  slice_max(order_by = minuslog10_fdr, n = 20, with_ties = FALSE) %>% 
  ungroup() %>%
  mutate(
    Description_short = str_trunc(Description, width = 50, side = "right"),
    Description_short = fct_reorder(Description_short, minuslog10_fdr, .desc = TRUE)  # Reorder based on -log10(FDR)
  )


# Plot
p <- ggplot(top_go, aes(x = minuslog10_fdr, y = fct_rev(Description_short))) +
  geom_point(aes(size = Count, color = minuslog10_fdr)) +
  scale_color_gradient(low = "skyblue", high = "darkblue") +
  scale_size(range = c(2, 6)) +
  labs(
    x = expression(-log[10](FDR)),
    y = NULL,
    size = "Gene Count",
    color = expression(-log[10](FDR)),
    title = "Significant GO Terms LMO3_BAT interneuron"
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
  filename = paste0(outdir, "BMAP_GO_dotplot_LMO3_WGCNA_upGenes_BatPutamen_Interneurons.pdf"),
  plot = p,
  width = 7,
  height = 5
)

#### perform GO-term on green module
lmo3_up = GOenrich(gns = green_genes, uni = all_genes)

write.csv(lmo3_up, paste0(outdir,"GreenModule_LMO3_UpGenes_WGCNA_within_Bat_Putamen_compared_to_other_interneurons_EnrichGO.csv"), row.names = FALSE)

### dotplot of GO-terms
# Top terms per comparison
top_go <- lmo3_up %>%
  group_by(GO) %>%
  mutate(minuslog10_fdr = -log10(p.adjust)) %>% # Calculate -log10(p.adjust) for each group
  slice_max(order_by = minuslog10_fdr, n = 20, with_ties = FALSE) %>% 
  ungroup() %>%
  mutate(
    Description_short = str_trunc(Description, width = 50, side = "right"),
    Description_short = fct_reorder(Description_short, minuslog10_fdr, .desc = TRUE)  # Reorder based on -log10(FDR)
  )


# Plot
p <- ggplot(top_go, aes(x = minuslog10_fdr, y = fct_rev(Description_short))) +
  geom_point(aes(size = Count, color = minuslog10_fdr)) +
  scale_color_gradient(low = "skyblue", high = "darkblue") +
  scale_size(range = c(2, 6)) +
  labs(
    x = expression(-log[10](FDR)),
    y = NULL,
    size = "Gene Count",
    color = expression(-log[10](FDR)),
    title = "GreenModule Significant GO Terms LMO3_BAT interneuron"
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
  filename = paste0(outdir, "GreenModule_GO_dotplot_LMO3_WGCNA_upGenes_BatPutamen_Interneurons.pdf"),
  plot = p,
  width = 7,
  height = 5
)

#####
## GO-term for FOXP2_TSHZ2 interneurons
tshz2_up = GOenrich(gns = cyan_genes, uni = all_genes)

write.csv(tshz2_up, paste0(outdir,"FOXP2_TSHZ2_UpGenes_WGCNA_within_Bat_Putamen_compared_to_other_interneurons_EnrichGO.csv"), row.names = FALSE, quote = T)

### dotplot of GO-terms
# Top terms per comparison
top_go <- tshz2_up[tshz2_up$GO == "CC",] %>%
  group_by(GO) %>%
  mutate(minuslog10_fdr = -log10(p.adjust)) %>% # Calculate -log10(p.adjust) for each group
  slice_max(order_by = minuslog10_fdr, n = 20, with_ties = FALSE) %>% 
  ungroup() %>%
  mutate(
    Description_short = str_trunc(Description, width = 50, side = "right"),
    Description_short = fct_reorder(Description_short, minuslog10_fdr, .desc = TRUE)  # Reorder based on -log10(FDR)
  )


# Plot
p <- ggplot(top_go, aes(x = minuslog10_fdr, y = fct_rev(Description_short))) +
  geom_point(aes(size = Count, color = minuslog10_fdr)) +
  scale_color_gradient(low = "skyblue", high = "darkblue") +
  scale_size(range = c(2, 6)) +
  labs(
    x = expression(-log[10](FDR)),
    y = NULL,
    size = "Gene Count",
    color = expression(-log[10](FDR)),
    title = "Significant GO Terms FOXP2_TSHZ2 interneuron"
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
  filename = paste0(outdir, "GO_dotplot_FOXP2_TSHZ2_WGCNA_upGenes_BatPutamen_Interneurons.pdf"),
  plot = p,
  width = 7,
  height = 5
)

#################################################

############################################################
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
  
