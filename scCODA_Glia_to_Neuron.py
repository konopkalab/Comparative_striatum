#!/usr/bin/env python3

# start a python session
#python
# Imports
import warnings
import scanpy as sc
import rpy2
import os
import sys
import io
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from rpy2 import robjects
from rpy2.robjects import pandas2ri
# Activate automatic conversion between R and pandas DataFrames
#pandas2ri.activate()
import rpy2.robjects as robjects
# scCODA Imports
import pandas as pd
import anndata as ad
import sccoda
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz
from sccoda.util import comp_ana as mod

import sccoda.datasets as scd

warnings.filterwarnings("ignore")
import sys
#sys.path.append("/project/Neuroinformatics_Core/Konopka_lab/s422071/SCRIPTS_pr/SCRIPTS/comparative_striatum_revision")

sys.path.append("/bmapfs/archive/konopkalab1/SCRIPTS_pr/SCRIPTS/comparative_striatum_revision")

from python_functions import aggregate_cellType_pairs_in_STR
from python_functions import process_species_comparison
from python_functions import process_tissue_comparison

# load the data
data_scanpy = sc.read_h5ad("/bmapfs/archive/konopkalab1/Gozde_workdir/projects/comparative_striatum/seu_objs/str_all_species_withallSamples_annotated_seu_obj.h5ad")
####################
### scCODA - extract the count matrix which scCODA requires
######################################
# Swap COP and OPC labels
tmp = data_scanpy.obs['newannot'].copy()

data_scanpy.obs['newannot'] = (
    tmp.replace({
        'COP': 'OPC',
    })
)
#comps = 'Tissue'
#n = "Neuron"
#d = "Glia"
 
# define the comparisons
comparisons = ('Species', 'Tissue')
# define which ratio (cell type) to calculate
numerators = ('Neuron',)
denominators = ('Glia', 'MOL', 'OPC', 'Astrocyte','Microglia')


### Compare across each 'comparison' type (Species, Tissue)
for comp in comparisons:
    for d in denominators:
        for n in numerators:
            # For each comparison type, append the tuple (comparison, numerator, denominator)
            
            aggregated_anndata = aggregate_cellType_pairs_in_STR(n, d, data_scanpy)

            #if comp == "Species":
                # All tissues
                process_species_comparison(aggregated_anndata, n, d)

                # Caudate + Caudoputamen
                process_species_comparison(aggregated_anndata, n, d, tissues=['Caudate', 'Caudoputamen'])

                # Putamen + Caudoputamen
                process_species_comparison(aggregated_anndata, n, d, tissues=['Putamen', 'Caudoputamen'])

                
            if comp == "Tissue":
                process_tissue_comparison(aggregated_anndata, n, d)
                    
