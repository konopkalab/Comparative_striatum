#!/usr/bin/env python3
# conda activate sccoda_env
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
import numpy as np
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
from python_functions import run_interneuron_subtype_sccoda

# load the data
Non_SPN_scanpy  = sc.read_h5ad("/bmapfs/archive/konopkalab1/Gozde_workdir/projects/comparative_striatum/Gena_workdir_08252025/03_INTEGRATE_ALL/Non_SPN/ANNOTATION/GB_Non_SPN_AllSpecies_AllTissues_human_pr_coding_orthologs_ANNOTATED.h5ad")

import os

os.chdir("/ifshome/gbuyukka/konopkalab1/Gozde_workdir/projects/comparative_striatum/revision/prop_analysis/Fig_5A_scCODA/")
####################
### scCODA - extract the count matrix which scCODA requires
######################################
# n = 'TAC3'
interneuron_types = (
    Non_SPN_scanpy.obs["newannot"]
    .dropna()
    .unique()
    .tolist()
)

for n in interneuron_types:
    print(f"Running scCODA for interneuron subtype: {n}")

    # Caudate + Caudoputamen
    run_interneuron_subtype_sccoda(
        Non_SPN_scanpy,
        target_subtype=n,
        tissues=["Caudate", "Caudoputamen"],
        out_prefix=f"scCODA_{n}_vs_OtherInterneurons_Primate_vs_Non_Primate_Caud_Caudoput"
    )

    # Putamen + Caudoputamen
    run_interneuron_subtype_sccoda(
        Non_SPN_scanpy,
        target_subtype=n,
        tissues=["Putamen", "Caudoputamen"],
        out_prefix=f"scCODA_{n}_vs_OtherInterneurons_Primate_vs_Non_Primate_Put_Caudoput"
    )
