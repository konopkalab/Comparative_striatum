
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
source('~/SCRIPTS/utility_functions.R')
set.seed(1234)

####
## PREPARE DATASETS
####

# Load broadly annotated objects
human_caud_put <- readRDS('/endosome/work/Neuroinformatics_Core/gkonop/01_HUMAN/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN/ANNOTATED/human_integrated_cuadate_putamen_ANNOTATED.RDS')

chimp_caud_put <- readRDS('/endosome/work/Neuroinformatics_Core/gkonop/01_CHIMP/03_CLUSTER_ANNOTATE/CAUDATE_PUTAMEN/ANNOTATED/chimp_integrated_cuadate_putamen_ANNOTATED.RDS')

macaque_caud_put <- readRDS("/endosome/work/Neuroinformatics_Core/gkonop/01_MACAQUE/INTEGRATED/03_CLUSTER_ANNOTATE/Caudate_Putamen/ANNOTATED/macaque_integrated_cuadate_putamen_ANNOTATED.RDS")

marmoset_caud_put <- readRDS("/endosome/work/Neuroinformatics_Core/gkonop/01_MARMOSET/MARMOSET_INTEGRATION/Caudate_putamen/ANNOTATED/marmoset_integrated_cuadate_putamen_ANNOTATED.RDS")

mouse = read_rds('~/workdir/01_MOUSE/03_CLUSTER_ANNOTATE/Mouse_Caudate_Annotated_FINAL.RDS')

bat_caud_put <- readRDS('/endosome/work/Neuroinformatics_Core/gkonop/01_BAT/03_CLUSTER_ANNOTATE/Caudate_Putamen/ANNOTATED/bat_integrated_cuadate_putamen_ANNOTATED.RDS')

# check
dim(human_caud_put)
table(human_caud_put$Tissue)

dim(chimp_caud_put)
table(chimp_caud_put$Tissue)

dim(macaque_caud_put)
table(macaque_caud_put$Tissue)

dim(marmoset_caud_put)
table(marmoset_caud_put$Tissue)

dim(mouse)

dim(bat_caud_put)
table(bat_caud_put$Tissue)

# Reshape and combine metadata
mouse$Tissue = 'CaudoPutamen'
mouse$Species = 'Mouse'


human_meta = human_caud_put@meta.data[, c('orig.ident', 'newannot', 'Tissue')]
human_meta$Species = 'Human'

chimp_meta = chimp_caud_put@meta.data[, c('orig.ident', 'newannot', 'Tissue')]
chimp_meta$Species = 'Chimp'

macaque_meta = macaque_caud_put@meta.data[, c('orig.ident', 'newannot', 'Tissue')]
macaque_meta$Species = 'Macaque'

marmoset_meta = marmoset_caud_put@meta.data[, c('orig.ident', 'newannot', 'Tissue')]
marmoset_meta$Species = 'Marmoset'

bat_meta = bat_caud_put@meta.data[, c('orig.ident', 'newannot', 'Tissue')]
bat_meta$Species = 'Bat'

# Combine putamen metadata
all_meta = do.call(rbind, list(human_meta[, c('orig.ident', 'newannot', 'Tissue', 'Species')],
				mouse@meta.data[, c('orig.ident', 'newannot', 'Tissue', 'Species')],
				chimp_meta, macaque_meta, marmoset_meta, bat_meta))

write_rds(all_meta, '~/workdir/02_PROPORTIONAL_COMPARISONS/GB_Combined_metadata.RDS')


####
## GLIA / NEURON
####

# Assign glia or neuron
all_meta$GorN = ifelse(all_meta$newannot %in% c('SPN', 'Non_SPN'), 'Neuron', 'Glia')


# Calculate glia / neuron ratio for tissues
df = all_meta %>% group_by(orig.ident,Tissue,Species) %>% group_keys %>% as.data.frame
df$totsize = all_meta %>% group_by(orig.ident,Tissue,Species) %>% group_size
df$gliasize = all_meta %>% filter(GorN == 'Glia') %>% group_by(orig.ident,Tissue,Species) %>% group_size
df$neuronsize = all_meta %>% filter(GorN == 'Neuron') %>% group_by(orig.ident,Tissue,Species) %>% group_size
df$grn = df$gliasize / df$neuronsize
df$Species = factor(df$Species, levels = c('Human', 'Chimp', 'Macaque', 'Marmoset', 'Bat', 'Mouse'))

df_no_mouse = df[df$Species != 'Mouse',]

# Compare and plot the glia / neuron ratio between caudate and putament per species
comps = list(c('Caudate', 'Putamen'))

pdf("GB_Glia_to_Neuron_Ratio_Caudate_Putamen.pdf", width = 10, height = 6)
ggboxplot(df_no_mouse, x = 'Tissue', y = 'grn', color = 'Tissue', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('Glia / Neuron Ratio') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
facet_wrap(~Species, nrow = 1) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
dev.off()

### correct for multi-test
compare_means(neuronRat ~ Tissue, data = df_no_mouse, group.by = "Species", method = "t.test", p.adjust.method = "BH")

# Compare and plot the glia / neuron ratio across species
comps = list(c('Marmoset', 'Macaque'), c('Marmoset', 'Human'), c('Macaque', 'Human'), c('Chimp', 'Human'), c('Human', 'Bat'))

df_no_mouse_caud = df_no_mouse[df_no_mouse$Tissue == 'Caudate',]

pdf("GB_Glia_to_Neuron_Ratio_Caudate_Across_Species.pdf", width = 6, height = 6)
ggboxplot(df_no_mouse_caud, x = 'Species', y = 'grn', color = 'Species', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('Glia / Neuron Ratio in Caudate') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
dev.off()

compare_means(neuronRat ~ Species, data = df_no_mouse_caud,  method = "t.test", p.adjust.method = "BH")

df_no_mouse_put = df_no_mouse[df_no_mouse$Tissue == 'Putamen',]

pdf("GB_Glia_to_Neuron_Ratio_Putamen_Across_Species.pdf", width = 6, height = 6)
ggboxplot(df_no_mouse_put, x = 'Species', y = 'grn', color = 'Species', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('Glia / Neuron Ratio in Putamen') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
dev.off()

compare_means(neuronRat ~ Species, data = df_no_mouse_put,  method = "t.test", p.adjust.method = "BH")

# Comparison with mouse (caudate+putamen)
comps = list(c('Marmoset', 'Macaque'), c('Marmoset', 'Human'), c('Macaque', 'Human'), c('Chimp', 'Human'), c('Human', 'Bat'), c('Human', 'Mouse'), c('Bat', 'Mouse'))

pdf("GB_Glia_to_Neuron_Ratio_DorsalStriatum_Across_Species.pdf", width = 6, height = 6)
ggboxplot(df, x = 'Species', y = 'grn', color = 'Species', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('Glia / Neuron Ratio') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
dev.off()

compare_means(neuronRat ~ Species, data = df,  method = "t.test", p.adjust.method = "BH")
####
## OPC / NEURON
####

# Assign glia or neuron
all_meta$OorN = ifelse(all_meta$newannot %in% c('SPN', 'Non_SPN'), 'Neuron', ifelse(all_meta$newannot %in% c('OPC', 'COP'), 'OPC', 'Glia'))

# Calculate the ratio for tissues
df = all_meta %>% group_by(orig.ident,Tissue,Species) %>% group_keys %>% as.data.frame
df$totsize = all_meta %>% group_by(orig.ident,Tissue,Species) %>% group_size
df$OPCsize = all_meta %>% filter(OorN == 'OPC') %>% group_by(orig.ident,Tissue,Species) %>% group_size
df$neuronsize = all_meta %>% filter(OorN == 'Neuron') %>% group_by(orig.ident,Tissue,Species) %>% group_size
df$OPC_to_N = df$OPCsize / df$neuronsize
df$Species = factor(df$Species, levels = c('Human', 'Chimp', 'Macaque', 'Marmoset', 'Bat', 'Mouse'))

df_no_mouse = df[df$Species != 'Mouse',]

# Compare and plot the glia / neuron ratio between caudate and putament per species
comps = list(c('Caudate', 'Putamen'))

pdf("GB_OPC_to_Neuron_Ratio_Caudate_Putamen.pdf", width = 10, height = 6)
ggboxplot(df_no_mouse, x = 'Tissue', y = 'OPC_to_N', color = 'Tissue', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('OPC / Neuron Ratio') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
facet_wrap(~Species, nrow = 1) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
dev.off()

# multi-test correction
compare_means(N_to_OPC ~ Tissue, group.by = "Species", data = df,  method = "t.test", p.adjust.method = "BH")

# Compare and plot the ratio across species
comps = list(c('Marmoset', 'Macaque'), c('Marmoset', 'Human'), c('Macaque', 'Human'), c('Chimp', 'Human'), c('Human', 'Bat'))


df_caud = df[df$Tissue != 'Putamen',]

pdf(paste0(outdir,"Neuron_to_OPC_Ratio_Caudate_Across_Species.pdf"), width = 6, height = 6)
ggboxplot(df_caud, x = 'Species', y = 'N_to_OPC', color = 'Species', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('Neuron / OPC Ratio in Caudate') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
dev.off()

# multi-test correction
compare_means(N_to_OPC ~ Species, data = df_caud,  method = "t.test", p.adjust.method = "BH")

####
## MOL / NEURON
####

# Assign glia or neuron
all_meta$MorN = ifelse(all_meta$newannot %in% c('SPN', 'Non_SPN'), 'Neuron', ifelse(all_meta$newannot %in% c('MOL'), 'MOL', 'Glia'))


# Calculate the ratio for tissues
df = all_meta %>% group_by(orig.ident,Tissue,Species) %>% group_keys %>% as.data.frame
df$totsize = all_meta %>% group_by(orig.ident,Tissue,Species) %>% group_size
df$MOLsize = all_meta %>% filter(MorN == 'MOL') %>% group_by(orig.ident,Tissue,Species) %>% group_size
df$neuronsize = all_meta %>% filter(MorN == 'Neuron') %>% group_by(orig.ident,Tissue,Species) %>% group_size
df$MOL_to_N = df$MOLsize / df$neuronsize
df$Species = factor(df$Species, levels = c('Human', 'Chimp', 'Macaque', 'Marmoset', 'Bat', 'Mouse'))

df_no_mouse = df[df$Species != 'Mouse',]

# Compare and plot the glia / neuron ratio between caudate and putament per species
comps = list(c('Caudate', 'Putamen'))

pdf("GB_MOL_to_Neuron_Ratio_Caudate_Putamen.pdf", width = 10, height = 6)
ggboxplot(df_no_mouse, x = 'Tissue', y = 'MOL_to_N', color = 'Tissue', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('MOL / Neuron Ratio') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
facet_wrap(~Species, nrow = 1) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) 
dev.off()

# multi-test correction
compare_means(N_to_MOL ~ Tissue, group.by = "Species", data = df_no_mouse,  method = "t.test", p.adjust.method = "BH")

# Compare and plot the ratio across species
comps = list(c('Marmoset', 'Macaque'), c('Marmoset', 'Human'), c('Macaque', 'Human'), c('Chimp', 'Human'), c('Human', 'Bat'))


df_caud = df[df$Tissue != 'Putamen',]

pdf(paste0(outdir,"Neuron_to_MOL_Ratio_Caudate_Across_Species.pdf"), width = 6, height = 6)
ggboxplot(df_caud, x = 'Species', y = 'N_to_MOL', color = 'Species', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('Neuron / MOL Ratio in Caudate') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
dev.off()

# multi-test correction
compare_means(N_to_MOL ~ Species, data = df_caud,  method = "t.test", p.adjust.method = "BH")

####
## Astrocyte / NEURON
####

# Assign glia or neuron
all_meta$AorN = ifelse(all_meta$newannot %in% c('SPN', 'Non_SPN'), 'Neuron', ifelse(all_meta$newannot == 'Astrocyte', 'Astrocyte', 'Glia'))

# Calculate the ratio for tissues
df = all_meta %>% group_by(orig.ident,Tissue,Species) %>% group_keys %>% as.data.frame
df$totsize = all_meta %>% group_by(orig.ident,Tissue,Species) %>% group_size
df$Astsize = all_meta %>% filter(AorN == 'Astrocyte') %>% group_by(orig.ident,Tissue,Species) %>% group_size
df$neuronsize = all_meta %>% filter(AorN == 'Neuron') %>% group_by(orig.ident,Tissue,Species) %>% group_size
df$Ast_to_N = df$Astsize / df$neuronsize
df$Species = factor(df$Species, levels = c('Human', 'Chimp', 'Macaque', 'Marmoset', 'Bat', 'Mouse'))

df_no_mouse = df[df$Species != 'Mouse',]

# Compare and plot the glia / neuron ratio between caudate and putament per species
comps = list(c('Caudate', 'Putamen'))

pdf("GB_Ast_to_Neuron_Ratio_Caudate_Putamen.pdf", width = 10, height = 6)
ggboxplot(df_no_mouse, x = 'Tissue', y = 'Ast_to_N', color = 'Tissue', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('Ast / Neuron Ratio') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
facet_wrap(~Species, nrow = 1) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
dev.off()

# multi-test correction
compare_means(N_to_Ast ~ Tissue, group.by = "Species", data = df_no_mouse,  method = "t.test", p.adjust.method = "BH")

# Compare and plot the ratio across species
comps = list(c('Marmoset', 'Macaque'), c('Marmoset', 'Human'), c('Macaque', 'Human'), c('Chimp', 'Human'), c('Human', 'Bat'))


df_caud = df[df$Tissue != 'Putamen',]

pdf(paste0(outdir,"Neuron_to_Ast_Ratio_Caudate_Across_Species.pdf"), width = 6, height = 6)
ggboxplot(df_caud, x = 'Species', y = 'N_to_Ast', color = 'Species', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('Neuron / Ast Ratio in Caudate') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
dev.off()

# multi-test correction
compare_means(N_to_Ast ~ Species, data = df_caud,  method = "t.test", p.adjust.method = "BH")

####
## Microglia / NEURON
####

# Assign glia or neuron
all_meta$MiorN = ifelse(all_meta$newannot %in% c('SPN', 'Non_SPN'), 'Neuron', ifelse(all_meta$newannot == "Microglia", 'Microglia', 'Glia'))


# Calculate the ratio for tissues
df = all_meta %>% group_by(orig.ident,Tissue,Species) %>% group_keys %>% as.data.frame
df$totsize = all_meta %>% group_by(orig.ident,Tissue,Species) %>% group_size
df$Migsize = all_meta %>% filter(MiorN == 'Microglia') %>% group_by(orig.ident,Tissue,Species) %>% group_size
df$neuronsize = all_meta %>% filter(MiorN == 'Neuron') %>% group_by(orig.ident,Tissue,Species) %>% group_size
df$Mig_to_N = df$Migsize / df$neuronsize
df$Species = factor(df$Species, levels = c('Human', 'Chimp', 'Macaque', 'Marmoset', 'Bat', 'Mouse'))

df_no_mouse = df[df$Species != 'Mouse',]

# Compare and plot the glia / neuron ratio between caudate and putament per species
comps = list(c('Caudate', 'Putamen'))

pdf("GB_Microglia_to_Neuron_Ratio_Caudate_Putamen.pdf", width = 10, height = 6)
ggboxplot(df_no_mouse, x = 'Tissue', y = 'Mig_to_N', color = 'Tissue', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('Microglia / Neuron Ratio') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
facet_wrap(~Species, nrow = 1) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
dev.off()

# multi-test correction
compare_means(N_to_Mig ~ Tissue, group.by = "Species", data = df_no_mouse,  method = "t.test", p.adjust.method = "BH")

df_caud = df[df$Tissue != 'Putamen',]

pdf(paste0(outdir,"Neuron_to_Ast_Ratio_Caudate_Across_Species.pdf"), width = 6, height = 6)
ggboxplot(df_caud, x = 'Species', y = 'N_to_Mig', color = 'Species', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('Neuron to Microglia Ratio in Caudate') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
dev.off()

# multi-test correction
compare_means(N_to_Mig ~ Species, data = df_caud,  method = "t.test", p.adjust.method = "BH")

### convert seurat obj to AnnData for scCODA
library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)
library(stringr)

# path to save files
outdir = "/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/revision/prop_analysis/"

sub_obj_new_list = readRDS("/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/seu_objs/sub_obj_new_list_with_allsamples.RDS")

for (sp in names(sub_obj_new_list)) {
  sub_obj_new_list[[sp]][["integrated"]] <- NULL
  sub_obj_new_list[[sp]][["SCT"]] <- NULL
}

# merge seu obj
all_species_seu = Reduce(merge, sub_obj_new_list)

all_species_seu$Species = gsub("chimp","Chimp", all_species_seu$Species)
all_species_seu$broadannot = gsub("COP", "OPC", all_species_seu$newannot)

# Primate vs mouse, bat, and ferret
meta = all_species_seu@meta.data
orig_meta = meta

meta$Species_Broad = 'Primate'
meta[meta$Species == 'Bat', 'Species_Broad'] = 'Bat'
meta[meta$Species == 'Ferret', 'Species_Broad'] = 'Ferret'
meta[meta$Species == 'Mouse', 'Species_Broad'] = 'Mouse'
meta$id = paste0(meta$orig.ident, "_", meta$Tissue)

## adjust the metadata
# read species traits taken from https://genomics.senescence.info/species/entry.php?species=Mus_musculus, https://genomics.senescence.info/species/entry.php?species=Mustela_nigripes, https://genomics.senescence.info/species/entry.php?species=Phyllostomus_hastatus, https://animaldiversity.org/accounts/Phyllostomus_hastatus/, https://genomics.senescence.info/species/entry.php?species=Phyllostomus_discolor
species_traits = read.csv(paste0(outdir, "primate_lifeTraits.csv"))
species_traits$Chimp = species_traits$chimp 

sample_metadata = read.table(paste0(outdir, "comparative_striatum_sample_metadata.csv"), sep = "\t", header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
sample_metadata$id = as.character(paste0(sample_metadata$`Barcode/ID`, "_", sample_metadata$`Brain Area`))

meta$id = as.character(meta$id)
meta$cellbarc = str_split_i(rownames(meta), "-", 1)
meta$cellID = paste0(meta$cellbarc, "_", meta$id)

# add the ages in years to the single cell-level metadata
meta_df <- merge(meta, sample_metadata, by = "id", all.x = TRUE)
meta_df$Species <- meta_df$Species.x
meta_df$Age <- meta_df$Age.y 
meta_df$Sex <- meta_df$Sex.y

# check 
unique(meta_df[is.na(meta_df$Sex),"id"])
unique(meta_df[is.na(meta_df$Age),"id"])
unique(meta_df[is.na(meta_df$Species),"id"])

meta_df <- meta_df[, !grepl("\\.x$|\\.y$", names(meta_df))]
meta_df$Brain.Area = NULL

# Add orig.ident only if rownames are not already "cellbarc-orig.ident"
rownames(meta_df) <- ifelse(
  grepl("[-_]", meta_df$cellbarc),  # check for presence of "-" or "_"
  meta_df$cellbarc,                  # if present, keep as is
  paste0(meta_df$cellbarc, "-", meta_df$orig.ident)  # else, append orig.ident
)

# restore original cell order (this prevents scrambling)
meta_df <- meta_df[rownames(meta), ]

# Check if same set but different order
setequal(rownames(meta_df), rownames(meta))


humanize_all_ages <- function(meta_df, species_traits) {
  # Initialize vector for humanized ages (same length and order)
  humanized_ages <- numeric(nrow(meta_df))
  
  # Get unique species in the metadata
  species_list <- unique(meta_df$Species)
  
  for (species in species_list) {
    # Convert ages to numerical
    species_ages <- meta_df$Age_in_years[meta_df$Species == species]
    species_ages_num <- as.numeric(species_ages)
    if (species %in% c("Human")) {
      # For humans, humanized age = original age
      humanized_ages[meta_df$Species == species] <- species_ages_num
    } else {
      # Check if species exists in lifetraits
      if (species %in% colnames(species_traits)) {
        # Build linear model: species ~ Human (lifetraits)
        model <- lm(species_traits[[species]] ~ species_traits$Human)
        
        # Calculate humanized ages using the model coefficients
        humanized <- (species_ages_num - model$coefficients[1]) / model$coefficients[2]
        
        # Assign back to correct rows in humanized_ages vector
        humanized_ages[meta_df$Species == species] <- humanized
      } else {
        warning(paste("Species", species, "not found in lifetraits. Assigning NA"))
        humanized_ages[meta_df$Species == species] <- NA
      }
    }
  }
  
  return(humanized_ages)
}


meta_df$Humanized_age = humanize_all_ages(meta_df, species_traits)

# update the metadata of the seu obj
# Join by cell ID (ensure cell names match exactly)
all_species_seu@meta.data = meta_df

# save seu obj
saveRDS(all_species_seu, file = "/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/seu_objs/all_species_seu.RDS")
all_species_seu = readRDS("/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/seu_objs/all_species_seu.RDS")

## convert seu obj to .h5ad
# 1. extract raw counts matrix
counts_mat <- all_species_seu[["RNA"]]$counts

# 2. metadata (convert factors → chars)
meta_df <- data.frame(lapply(all_species_seu@meta.data, as.character),
                      stringsAsFactors = FALSE)
# check
nrow(meta_df) == ncol(counts_mat)

# 3. create SCE object
sce <- SingleCellExperiment(
assays = list(counts = counts_mat)
)
colData(sce) <- S4Vectors::DataFrame(meta_df)

# 4. write h5ad directly
out_path <- paste0(
"/project/Neuroinformatics_Core/Konopka_lab/s422071/workdir_pr/s422071/projects/comparative_striatum/seu_objs/str_all_species_withallSamples_annotated_seu_obj.h5ad")

zellkonverter::writeH5AD(sce, file = out_path, X_name = "counts")
