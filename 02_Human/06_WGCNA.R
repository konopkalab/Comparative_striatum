# load packages
library(patchwork)
library(Seurat)
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
library(variancePartition)
library(fastDummies)
#BiocManager::install("scuttle")
library(scuttle)
library(clusterProfiler)
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
set.seed(1234)
source("/bmapfs/archive/konopkalab1/SCRIPTS_pr/SCRIPTS/utility_functions.R")
library(curl)
#conda install bioconda::r-wgcna
#BiocManager::install('WGCNA')
library(WGCNA)
library(flashClust)
options(warn = 0)

# path to save files
outdir = "/bmapfs/archive/konopkalab1/Gozde_workdir/projects/comparative_striatum/revision/UTSW_human_wgcna/"

## load the data
str_seu = readRDS("/bmapfs/archive/konopkalab1/Gozde_workdir/projects/comparative_striatum/seu_objs/merged_clean_AllCellTypes_AllTissues_seu.RDS")
# str_seu = readRDS("/u/project/gkonopka/gozdeb/gozde_biohpc/workdir_pr/s422071/projects/comparative_striatum/seu_objs/merged_clean_AllCellTypes_AllTissues_seu.RDS")


# subset human interneurons
human_str = subset(str_seu, subset = Species == "Human" & newannot_2 != "Non_SPN")
# put COP together with OPC
human_str$detailed_annot = gsub("COP", "OPC", human_str$detailed_annot)

human_str_pb_count <- pseudobulk_sample_cellType_v5(
seurObj = human_str,
sample_col = "detailed_annot",
n_copy = 1,
sample_number = 3) #13674 genes

#####################################################
# normalize the counts and transpose the mat
dge <- DGEList(counts = human_str_pb_count)
cpm_norm <- edgeR::cpm(dge, log = FALSE)   # library-size normalized counts (not log yet)
# apply variance stabilizing transformation
norm_aggr_counts <- t(log2(cpm_norm + 1))

# find the mean UMI counts for each sample in the species-detailed_annot pair
human_str$Sample_detailed_annot = paste0(human_str$Sample, "_", human_str$detailed_annot)

meta_Sample_detailed_annot <- human_str@meta.data %>%
dplyr::group_by(Sample_detailed_annot) %>%
dplyr::summarise(
log10_sum_ncount_RNA = log10(sum(nCount_RNA, na.rm = TRUE)),
n = dplyr::n(),
.groups = "drop")
meta_Sample_detailed_annot = as.data.frame(meta_Sample_detailed_annot)

# set rownames
rownames(meta_Sample_detailed_annot) = meta_Sample_detailed_annot$Sample_detailed_annot

# check if the rownames in metadata matches the rownames of the aggregated count matrix for each celltype
(match_ok <- all(rownames(norm_aggr_counts) == rownames(meta_Sample_detailed_annot)))


# If FALSE, correct the metadata order
#meta_Sample_detailed_annot = meta_Sample_detailed_annot[match(rownames(norm_aggr_counts), rownames(meta_Sample_detailed_annot)),]
# check if the rownames in metadata matches the colnames of the aggregated count matrix for each celltype

# check again
#(match_ok <- all(rownames(norm_aggr_counts) ==rownames(meta_Sample_detailed_annot)))
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
pdf(paste0(outdir, "BMAP_human_AllCellTypes_sampleclustering_aggregated.pdf"), width = 30, height = 10)  

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
pdf(paste0(outdir,"SoftThresholdingPower_signed_human_AllCellTypes.pdf"))
par(mfrow = c(1,2), mar=c(5.1,5.1,4.1,2.1));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red");
abline(h=0.5,col="red"); abline(h=0.8,col="blue");abline(h=0.9,col="black")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

(PWR=sft$powerEstimate) #16

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
deepSplit=0.1,
verbose=1,
mergeCutHeight=0.10,
reassignThreshold = 1e-10,
numericLabels=TRUE)


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
MEsAutomatic=net$MEs
unique(moduleColorsAutomatic)
table(moduleColorsAutomatic)
write.table(moduleColorsAutomatic, paste0(outdir, "human_AllCellTypes_WGCNA_module_colors.txt",sep="\t",quote=F))
save(net,file=paste0(outdir,"human_AllCellTypes_WGCNA_net.RData"))

# Calculate module membership for each gene
KMEs<-signedKME(datExpr, net$MEs,corFnc = "bicor") # bicor: Calculate biweight midcorrelation efficiently for matrices.
kme=data.frame(colnames(datExpr), moduleColorsAutomatic, KMEs)
colnames(kme)[1]="Symbol"
rownames(kme)=NULL
write.table(kme,"Gene-module_membership_human_AllCellTypes.txt",sep="\t",quote=F)
save(kme,file=paste0(outdir,"human_AllCellTypes_WGCNA_kme.RData"))
kme = get(load(paste0(outdir,"human_AllCellTypes_WGCNA_kme.RData")))

# generate a list of genes belonging to each module
kmeL = split(kme, kme$moduleColorsAutomatic) %>% lapply(., function(x){x$Symbol})


# Dat Trait Heatmap
nGenes = ncol(datExpr) 
nSamples = nrow(datExpr)

### find the mean UMI counts for each sample in the species-detailed_annot pair
human_str$Sample_detailed_annot = paste0(human_str$Sample, "_", human_str$detailed_annot)

# calc the log10 of mean UMI counts for each sample
meta_Sample_detailed_annot <- human_str@meta.data %>%
dplyr::group_by(Sample_detailed_annot) %>%
dplyr::summarise(
log10_sum_ncount_RNA = log10(sum(nCount_RNA, na.rm = TRUE)),
n = dplyr::n(),
.groups = "drop")
meta_Sample_detailed_annot = as.data.frame(meta_Sample_detailed_annot)

# set rownames
rownames(meta_Sample_detailed_annot) = meta_Sample_detailed_annot$Sample_detailed_annot

# set datTraits dataframe 
Samples <- rownames(norm_aggr_counts)
traitRows <- match(Samples, meta_Sample_detailed_annot$Sample_detailed_annot)
datTraits <- meta_Sample_detailed_annot[traitRows,]
rownames(datTraits) <- meta_Sample_detailed_annot[traitRows, 1]
datTraits$detailed_annot <- stringr::str_extract(
  datTraits$Sample_detailed_annot,
  "(?<=_(Putamen|Caudate)_).*"
)
datTraits$Sample <- stringr::str_extract(
  datTraits$Sample_detailed_annot,
  "^.*(?=_(Putamen|Caudate)_)"
)
# remove unnecessary cols
datTraits$Sample_detailed_annot = NULL
datTraits_numeric <- fastDummies::dummy_cols(datTraits, remove_selected_columns = TRUE)
rownames(datTraits_numeric) = rownames(datTraits)

MEs0 = moduleEigengenes(datExpr,moduleLabelsAutomatic)$eigengenes
MEsTemp = MEs0
MEsTemp$MEgrey=NULL
# use determen correlation (spearman) between data traits and modules 
modTraitCor= cor(MEsTemp, datTraits_numeric,method="spearman")
write.table(modTraitCor,paste0(outdir, "modTraitCor_human_AllCellTypes.txt"),sep="\t",quote=F)
modTraitP = corPvalueStudent(modTraitCor, nSamples)
write.table(modTraitP,paste0(outdir,"modTraitP_human_AllCellTypes.txt"),sep="\t",quote=F)



# extract all genes
all_genes = kme$Symbol #[1] 13674 genes
##### pick celltype-specific modules: pink = SPN, cyan = MOL, magenta = Ast, midnightblue = OPC, yellow = microglia

#### run Fisher's exact test between the HS DEGs of each celltype and eigengenes of modules
# read HS DEGs up
hs_degs_up_df = read.csv("/bmapfs/archive/konopkalab1/Gozde_workdir/projects/comparative_striatum/revision/DEG_analysis/DESeq2/Human_specific_Upgenes_logFC_0_5.csv")
hs_degs_down_df = read.csv("/bmapfs/archive/konopkalab1/Gozde_workdir/projects/comparative_striatum/revision/DEG_analysis/DESeq2/Human_specific_Downgenes_logFC_0_5.csv")

hs_degs_up_list = split(hs_degs_up_df, hs_degs_up_df$CellType)
hs_degs_up <- lapply(hs_degs_up_list, function(df) df$`Gene.HC_gene`)

hs_degs_down_list = split(hs_degs_down_df, hs_degs_down_df$CellType)
hs_degs_down <- lapply(hs_degs_down_list, function(df) df$`Gene.HC_gene`)


geneOvEnr(gnL1 = hs_degs_up, gnL2 = kmeL, bcg = nrow(kme), plot = T, hg = 18, wd = 14, cutoff = 0.01, fn = paste0(outdir, "HS_DEG_UpGenes_module_gene_enrichment")) # plot can be set to false (F) to return a data frame instead of the plot.

geneOvEnr(gnL1 = hs_degs_down, gnL2 = kmeL, bcg = nrow(kme), plot = T, hg = 18, wd = 14, cutoff = 0.01, fn = paste0(outdir, "HS_DEG_DownGenes_module_gene_enrichment")) # plot can be set to false (F) to return a data frame instead of the plot.
# for upregulated HS genes light cyan and blue modules are dSPN and iSPN specific. red module in addition to dSPN and iSPN, is also specific to the eSPN-HS upgenes.
###########################
### GO-term
for (module in names(kmeL)) {
	go_res = GOenrich(gns = kmeL[[module]], uni = all_genes)

	write.csv(go_res, paste0(outdir, module, "_genes_WGCNA_human_AllTissues_AllCellTypes_EnrichGO.csv"), row.names = FALSE)

	
	# if there are no results, skip to next module
  	if (is.null(go_res)) {
   	 next
  	}
	### dotplot of GO-terms
	# Top terms per comparison
	top_go <- go_res %>%
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
	    title = "Significant GO Terms"
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
	  filename = paste0(outdir, "GO_dotplot_", module, "_WGCNA_genes_human_AllTissues_AllCellTypes.pdf"),
	  plot = p,
	  width = 7,
	  height = 5
	)
}

###### violin plots of important genes
## these genes are found in red and blue modules (SPN-specific ones) and based on the log2FC difference between human vs chimp DESEq2 DEG analysis
# set dir
setwd("/bmapfs/archive/konopkalab1/Gozde_workdir/projects/comparative_striatum/revision/UTSW_human_wgcna/")

# subset to only SPNs
SPN_seu = subset(str_seu, subset = newannot_2 %in% c("dSPN", "iSPN", "eSPN") & Species %in% c("Human", "chimp", "Macaque", "Marmoset"))

DefaultAssay(SPN_seu) = "RNA"
SPN_seu <- NormalizeData(SPN_seu)

#### get the genes on interest
## blue module genes
blue_module_df = read.csv("/bmapfs/archive/konopkalab1/Gozde_workdir/projects/comparative_striatum/revision/UTSW_human_wgcna/blue_genes_WGCNA_human_AllTissues_AllCellTypes_EnrichGO.csv")

# get the signif ones and sort (lowest p.adjust to highest)
blue_module_df = blue_module_df[blue_module_df$p.adjust < 0.01,]
#blue_module_df_sorted = arrange(blue_module_df, p.adjust)

# get top BP
blue_top_BP_cat_df =  blue_module_df[blue_module_df$Description %in% c("regulation of ion transmembrane transport", "regulation of membrane potential", "synaptic vesicle exocytosis", "neurotransmitter secretion", "synapse organization"),]
 
# Combine gene IDs into a single vector, separated by "/"
blue_top_BP_cat_genes <- paste(blue_top_BP_cat_df$geneID, collapse = "/")
# get the genes in a vector

# Split the string into individual gene IDs
blue_top_BP_cat_genes <- strsplit(blue_top_BP_cat_genes, "/")[[1]]

## red module genes
red_module_df = read.csv("/bmapfs/archive/konopkalab1/Gozde_workdir/projects/comparative_striatum/revision/UTSW_human_wgcna/red_genes_WGCNA_human_AllTissues_AllCellTypes_EnrichGO.csv")
# Get the significant ones and sort (lowest p.adjust to highest)
red_module_df = red_module_df[red_module_df$p.adjust < 0.01,]

# Filter for GO category "BP" and take the top
red_top_BP_cat_df = red_module_df[red_module_df$Description %in% c("modulation of chemical synaptic transmission", "regulation of synaptic plasticity", "signal release from synapse", "synapse assembly", "glutamate receptor signaling pathway"),]


### plot only the top categories
# compute global x-axis limits
x_limits <- range(
  -log10(blue_top_BP_cat_df$p.adjust),
  -log10(red_top_BP_cat_df$p.adjust),
  na.rm = TRUE
)

# add padding
x_limits <- c(0, ceiling(x_limits[2]))

# Top terms per comparison
blue_top_go <- blue_top_BP_cat_df %>%
  mutate(minuslog10_fdr = -log10(p.adjust)) %>% # Calculate -log10(p.adjust) for each group
  ungroup() %>%
  mutate(
    Description_short = str_trunc(Description, width = 50, side = "right"),
    Description_short = fct_reorder(Description_short, minuslog10_fdr, .desc = TRUE)  # Reorder based on -log10(FDR)
  )
     
# Plot
p <- ggplot(blue_top_go, aes(x = minuslog10_fdr, y = fct_rev(Description_short))) +
	  geom_point(aes(size = Count, color = minuslog10_fdr)) +
	  scale_color_gradient(low = "skyblue", high = "darkblue") +
	  scale_size(range = c(2, 6)) +
    coord_cartesian(xlim = x_limits) + 
	  labs(
	    x = expression(-log[10](FDR)),
	    y = NULL,
	    size = "Gene Count",
	    color = expression(-log[10](FDR)),
	    title = "Significant GO Terms"
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
  filename = paste0(outdir, "Top_BP_GO_dotplot_bluemodule_WGCNA_genes_human_AllTissues_AllCellTypes.pdf"),
  plot = p,
  width = 7,
  height = 5
)

### red module plot
# Top terms per comparison
red_top_go <- red_top_BP_cat_df %>%
  mutate(minuslog10_fdr = -log10(p.adjust)) %>% # Calculate -log10(p.adjust) for each group
  ungroup() %>%
  mutate(
    Description_short = str_trunc(Description, width = 50, side = "right"),
    Description_short = fct_reorder(Description_short, minuslog10_fdr, .desc = TRUE)  # Reorder based on -log10(FDR)
  )
    
# Plot
p <- ggplot(red_top_go, aes(x = minuslog10_fdr, y = fct_rev(Description_short))) +
	  geom_point(aes(size = Count, color = minuslog10_fdr)) +
	  scale_color_gradient(low = "skyblue", high = "darkblue") +
	  scale_size(range = c(2, 6)) +
    coord_cartesian(xlim = x_limits) +
	  labs(
	    x = expression(-log[10](FDR)),
	    y = NULL,
	    size = "Gene Count",
	    color = expression(-log[10](FDR)),
	    title = "Significant GO Terms"
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
  filename = paste0(outdir, "Top_BP_GO_dotplot_redmodule_WGCNA_genes_human_AllTissues_AllCellTypes.pdf"),
  plot = p,
  width = 7,
  height = 5
)


# Combine gene IDs into a single vector, separated by "/"
red_top_BP_cat_genes <- paste(red_top_BP_cat_df$geneID, collapse = "/")

# Get the genes in a vector
# Split the string into individual gene IDs
red_top_BP_cat_genes <- strsplit(red_top_BP_cat_genes, "/")[[1]]

HS_genes = read.csv("/bmapfs/archive/konopkalab1/Gozde_workdir/projects/comparative_striatum/revision/DEG_analysis/DESeq2/Human_specific_Upgenes_logFC_0_5.csv")
SPN_HS_genes = HS_genes[HS_genes$CellType %in% c("dSPN", "iSPN", "eSPN"),]
SPN_HS_genes = SPN_HS_genes[order(SPN_HS_genes$Gene.HC_log2FoldChange), ]$Gene.HC_gene

(genes_of_interest = union(intersect(unique(SPN_HS_genes), unique(blue_top_BP_cat_genes)), intersect(unique(SPN_HS_genes), unique(red_top_BP_cat_genes)))) # "PKD1L1" "CALB1"  "TPBG"   "KCNIP2" "HRH3"   "UNC80"  "GRIN1"  "SNAP91" "SYN2"   "SNAP25"

plot_gene_violin <- function(seu, gene, group.by = "Species") {
  
  # fetch expression + metadata ONLY for this gene
  df <- FetchData(
    object = seu,
    vars = c(gene, group.by)
  )
  
  colnames(df) <- c("expression", "group")
  
  pdf(paste0(gene, "_violinPlot_blue_red_module_primates_SPNs.pdf"))
  p <- ggplot(df, aes(x = group, y = expression)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.3) +
    theme_classic() +
    labs(
      title = gene,
      x = group.by,
      y = "Expression"
    )
    print(p)
    dev.off()
}

for (g in genes_of_interest) {
  print(plot_gene_violin(SPN_seu, g, group.by = "Species"))
}

SPN_seu$Species = factor(SPN_seu$Species, levels = c("Human", "chimp", "Macaque", "Marmoset"))
### try seurat vlnplot with conda activate striatum_R423_seurat_safe
p <- VlnPlot(
  SPN_seu,
  features = genes_of_interest,
  group.by = "Species",
  pt.size = 0
)

# standardize the y axis scale
p <- p & scale_y_continuous(limits = c(0, 5), breaks = 0:5) 

pdf("BP_VlnPlot_blue_red_module_primates_SPNs.pdf")
p
dev.off()

pdf("BP_VlnPlot_blue_red_module_primates_SPNs_separated.pdf")
VlnPlot(SPN_seu, features = genes_of_interest, group.by = "newannot_2",
  split.by = "Species", pt.size = 0)
dev.off()

pdf("BP_VlnPlot_blue_red_module_primates_SPNs_single_plot.pdf")
VlnPlot(
  SPN_seu,
  features = genes_of_interest,
  split.by = "Species",
  pt.size = 0,
  stack = TRUE,
  flip = TRUE
)
dev.off()
