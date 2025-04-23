### convert Seurat object to python AnnData object
# Method 1
# Tips: 
#1. set default assay to RNA before covert to h5ad.
#2. if raw read count need to be imported to anndata, you should only contain counts slot in your seurat object before convertion (https://zqfang.github.io/2020-04-28-seurat2scanpy/)
library(Seurat)
library(SeuratDisk)

allseur_integrated <- readRDS("/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/ANNOTATED/GB_AllTissues_Primates_bat_mouse_ferret_Integrated_SPN_ncbi_human_pr_coding_orthologs_res1_5_res_0_5_ANNOTATED.RDS")

# step 1: Slim down a Seurat object. So you get raw counts, lognorm counts

seu = DietSeurat(
  allseur_integrated,
  counts = TRUE, # so, raw counts save to adata.raw.X 
  data = TRUE, # so, log1p counts save to adata.X
  scale.data = FALSE, # set to false, or else will save to adata.X
  features = rownames(allseur_integrated), # export all genes, not just top highly variable genes
  assays = "RNA",
  dimreducs = c("pca","umap"),
  graphs = c("RNA_nn", "RNA_snn"), # to RNA_nn -> distances, RNA_snn -> connectivities
  misc = TRUE
)

# step 2: factor to character, or else your factor will be number in adata 
i <- sapply(seu@meta.data, is.factor)
seu@meta.data[i] <- lapply(seu@meta.data[i], as.character)

# step 3: convert 
SaveH5Seurat(seu, filename = "/endosome/work/Neuroinformatics_Core/gkonop/03_INTEGRATE_ALL/SPN/ANNOTATED/GB_SPN_AllSpecies_human_pr_coding_orth_integrated_annotated.h5seurat", overwrite = TRUE)
Convert("/endosome/work/Neuroinformatics_Core/gkonop/03_INTEGRATE_ALL/SPN/ANNOTATED/GB_SPN_AllSpecies_human_pr_coding_orth_integrated_annotated.h5seurat", "/endosome/work/Neuroinformatics_Core/gkonop/03_INTEGRATE_ALL/SPN/ANNOTATED/GB_SPN_AllSpecies_human_pr_coding_orth_integrated_annotated.h5ad", assay="RNA", overwrite = TRUE)
################################################################
