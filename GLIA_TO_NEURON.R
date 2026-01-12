
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
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
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
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
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
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
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
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
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
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
dev.off()

# multi-test correction
compare_means(N_to_OPC ~ Tissue, group.by = "Species", data = df,  method = "t.test", p.adjust.method = "BH")

# Compare and plot the ratio across species
comps = list(c('Marmoset', 'Macaque'), c('Marmoset', 'Human'), c('Macaque', 'Human'), c('Chimp', 'Human'), c('Human', 'Bat'))

df_no_mouse_caud = df_no_mouse[df_no_mouse$Tissue == 'Caudate',]

pdf("GB_OPC_to_Neuron_Ratio_Caudate_Across_Species.pdf", width = 6, height = 6)
ggboxplot(df_no_mouse_caud, x = 'Species', y = 'OPC_to_N', color = 'Species', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('OPC / Neuron Ratio in Caudate') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
dev.off()


df_no_mouse_put = df_no_mouse[df_no_mouse$Tissue == 'Putamen',]

pdf("GB_OPC_to_Neuron_Ratio_Putamen_Across_Species.pdf", width = 6, height = 6)
ggboxplot(df_no_mouse_put, x = 'Species', y = 'OPC_to_N', color = 'Species', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('OPC / Neuron Ratio in Putamen') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
dev.off()

# Comparison with mouse (caudate+putamen)
comps = list(c('Marmoset', 'Macaque'), c('Marmoset', 'Human'), c('Macaque', 'Human'), c('Chimp', 'Human'), c('Human', 'Bat'), c('Human', 'Mouse'), c('Bat', 'Mouse'))

pdf("GB_OPC_to_Neuron_Ratio_DorsalStriatum_Across_Species.pdf", width = 6, height = 6)
ggboxplot(df, x = 'Species', y = 'OPC_to_N', color = 'Species', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('OPC / Neuron Ratio') +
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
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
dev.off()

# Compare and plot the ratio across species
comps = list(c('Marmoset', 'Macaque'), c('Marmoset', 'Human'), c('Macaque', 'Human'), c('Chimp', 'Human'), c('Human', 'Bat'))

df_no_mouse_caud = df_no_mouse[df_no_mouse$Tissue == 'Caudate',]

pdf("GB_MOL_to_Neuron_Ratio_Caudate_Across_Species.pdf", width = 6, height = 6)
ggboxplot(df_no_mouse_caud, x = 'Species', y = 'MOL_to_N', color = 'Species', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('MOL / Neuron Ratio in Caudate') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
dev.off()

df_no_mouse_put = df_no_mouse[df_no_mouse$Tissue == 'Putamen',]

pdf("GB_MOL_to_Neuron_Ratio_Putamen_Across_Species.pdf", width = 6, height = 6)
ggboxplot(df_no_mouse_put, x = 'Species', y = 'MOL_to_N', color = 'Species', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('MOL / Neuron Ratio in Putamen') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
dev.off()

# Comparison with mouse (caudate+putamen)
comps = list(c('Marmoset', 'Macaque'), c('Marmoset', 'Human'), c('Macaque', 'Human'), c('Chimp', 'Human'), c('Human', 'Bat'), c('Human', 'Mouse'), c('Bat', 'Mouse'))

pdf("GB_MOL_to_Neuron_Ratio_DorsalStriatum_Across_Species.pdf", width = 6, height = 6)
ggboxplot(df, x = 'Species', y = 'MOL_to_N', color = 'Species', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('MOL / Neuron Ratio') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
dev.off()

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

# Compare and plot the ratio across species
comps = list(c('Marmoset', 'Macaque'), c('Marmoset', 'Human'), c('Macaque', 'Human'), c('Chimp', 'Human'), c('Human', 'Bat'))

df_no_mouse_caud = df_no_mouse[df_no_mouse$Tissue == 'Caudate',]

pdf("GB_Ast_to_Neuron_Ratio_Caudate_Across_Species.pdf", width = 6, height = 6)
ggboxplot(df_no_mouse_caud, x = 'Species', y = 'Ast_to_N', color = 'Species', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('Ast / Neuron Ratio in Caudate') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
dev.off()

df_no_mouse_put = df_no_mouse[df_no_mouse$Tissue == 'Putamen',]

pdf("GB_Ast_to_Neuron_Ratio_Putamen_Across_Species.pdf", width = 6, height = 6)
ggboxplot(df_no_mouse_put, x = 'Species', y = 'Ast_to_N', color = 'Species', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('Ast / Neuron Ratio in Putamen') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
dev.off()

# Comparison with mouse (caudate+putamen)
comps = list(c('Marmoset', 'Macaque'), c('Marmoset', 'Human'), c('Macaque', 'Human'), c('Chimp', 'Human'), c('Human', 'Bat'), c('Human', 'Mouse'), c('Bat', 'Mouse'))

pdf("GB_Ast_to_Neuron_Ratio_DorsalStriatum_Across_Species.pdf", width = 6, height = 6)
ggboxplot(df, x = 'Species', y = 'Ast_to_N', color = 'Species', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('Ast / Neuron Ratio') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
dev.off()

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

# Compare and plot the ratio across species
comps = list(c('Marmoset', 'Macaque'), c('Marmoset', 'Human'), c('Macaque', 'Human'), c('Chimp', 'Human'), c('Human', 'Bat'))

df_no_mouse_caud = df_no_mouse[df_no_mouse$Tissue == 'Caudate',]

pdf("GB_Microglia_to_Neuron_Ratio_Caudate_Across_Species.pdf", width = 6, height = 6)
ggboxplot(df_no_mouse_caud, x = 'Species', y = 'Mig_to_N', color = 'Species', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('Microglia / Neuron Ratio in Caudate') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
dev.off()

df_no_mouse_put = df_no_mouse[df_no_mouse$Tissue == 'Putamen',]

pdf("GB_Microglia_to_Neuron_Ratio_Putamen_Across_Species.pdf", width = 6, height = 6)
ggboxplot(df_no_mouse_put, x = 'Species', y = 'Mig_to_N', color = 'Species', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('Microglia / Neuron Ratio in Putamen') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
dev.off()

# Comparison with mouse (caudate+putamen)
comps = list(c('Marmoset', 'Macaque'), c('Marmoset', 'Human'), c('Macaque', 'Human'), c('Chimp', 'Human'), c('Human', 'Bat'), c('Human', 'Mouse'), c('Bat', 'Mouse'))

pdf("GB_Microglia_to_Neuron_Ratio_DorsalStriatum_Across_Species.pdf", width = 6, height = 6)
ggboxplot(df, x = 'Species', y = 'Mig_to_N', color = 'Species', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('Microglia / Neuron Ratio') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
dev.off()


####
## NON_SPN / SPN
####

# Assign glia or neuron
neuron_meta = all_meta[all_meta$newannot %in% c('SPN', 'Non_SPN'),]

# Calculate glia / neuron ratio for tissues
gr1 = 'Non_SPN'
gr2 = 'SPN'
df2 = neuron_meta %>% group_by(orig.ident,Tissue,Species) %>% group_keys %>% as.data.frame
df2$totsize = neuron_meta %>% group_by(orig.ident,Tissue,Species) %>% group_size
df2$group_1 = neuron_meta %>% filter(newannot == gr1) %>% group_by(orig.ident,Tissue,Species) %>% group_size
df2$group_2 = neuron_meta %>% filter(newannot == gr2) %>% group_by(orig.ident,Tissue,Species) %>% group_size
df2$ratio = df2$group_1 / df2$group_2
df2$Species = factor(df$Species, levels = c('Human', 'Chimp', 'Macaque', 'Marmoset', 'Bat', 'Mouse'))

df2_no_mouse = df2[df2$Species != 'Mouse',]

# Compare and plot the ratio between caudate and putament per species
comps = list(c('Caudate', 'Putamen'))

pdf("GB_Non_SPN_To_SPN_Ratio_Caudate_Putamen.pdf", width = 10, height = 6)
ggboxplot(df2_no_mouse, x = 'Tissue', y = 'ratio', color = 'Tissue', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('Non_SPN / SPN Ratio') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
facet_wrap(~Species, nrow = 1) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
dev.off()

# Compare and plot the ratio across species
comps = list(c('Marmoset', 'Macaque'), c('Marmoset', 'Human'), c('Macaque', 'Human'), c('Chimp', 'Human'), c('Human', 'Bat'))

df2_no_mouse_caud = df2_no_mouse[df2_no_mouse$Tissue == 'Caudate',]

pdf("GB_Non_SPN_To_SPN_Ratio_Caudate_Across_Species.pdf", width = 7, height = 7)
ggboxplot(df2_no_mouse_caud, x = 'Species', y = 'ratio', color = 'Species', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('Non_SPN / SPN Ratio in Caudate') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
stat_compare_means(comparisons = comps, method = 't.test')
dev.off()

df2_no_mouse_put = df2_no_mouse[df2_no_mouse$Tissue == 'Putamen',]

pdf("GB_Non_SPN_To_SPN_Ratio_Putamen_Across_Species.pdf", width = 6, height = 7)
ggboxplot(df2_no_mouse_put, x = 'Species', y = 'ratio', color = 'Species', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('Non_SPN / SPN Ratio in Putamen') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
dev.off()

# Comparison with mouse (caudate+putamen)
comps = list(c('Marmoset', 'Macaque'), c('Marmoset', 'Human'), c('Macaque', 'Human'), c('Chimp', 'Human'), c('Human', 'Bat'), c('Human', 'Mouse'), c('Bat', 'Mouse'))

pdf("GB_Non_SPN_To_SPN_Ratio_DorsalStriatum_Across_Species.pdf", width = 6, height = 6)
ggboxplot(df2, x = 'Species', y = 'ratio', color = 'Species', add = 'dotplot', size = 1) +
NoLegend() +
xlab('') +
ylab('Non_SPN / SPN Ratio') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
stat_compare_means(comparisons = comps, method = 't.test')
dev.off()
