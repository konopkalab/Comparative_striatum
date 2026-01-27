## python_functions.py

import matplotlib.pyplot as plt
import warnings
import sccoda
import scanpy as sc
import rpy2
import os
import sys
import io
import matplotlib
from rpy2 import robjects
from rpy2.robjects import pandas2ri
# Activate automatic conversion between R and pandas DataFrames
import rpy2.robjects as robjects
# scCODA Imports
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
import numpy as np
import warnings

from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz
from sccoda.util import comp_ana as mod
import sccoda.datasets as scd

warnings.filterwarnings("ignore")
import sys

# ------------------------------------------------------------------
# Function 1: Aggregate counts based on input cellType pairs
# ------------------------------------------------------------------
def aggregate_cellType_pairs_in_STR(cellType1, cellType2, anndata_obj):
    anndata_obj = anndata_obj.copy()
    metadata_df = anndata_obj.obs[['id', 'newannot']].copy()
    metadata_df['aggregated_cellType'] = None

    # Condition 1: Neuron + Glia
    if ('Neuron' in [cellType1, cellType2]) and ('Glia' in [cellType1, cellType2]):
        metadata_df.loc[metadata_df['newannot'].isin(['SPN', 'Non_SPN']), 'aggregated_cellType'] = 'Neuron'
        glia_types = ['Astrocyte', 'MOL', 'Microglia', 'OPC', 'COP']
        metadata_df.loc[metadata_df['newannot'].isin(glia_types), 'aggregated_cellType'] = 'Glia'

    # Condition 2: Neuron + other
    elif 'Neuron' in [cellType1, cellType2]:
        metadata_df.loc[metadata_df['newannot'].isin(['SPN', 'Non_SPN']), 'aggregated_cellType'] = 'Neuron'
        other_cell_type = cellType1 if cellType2 == 'Neuron' else cellType2
        metadata_df.loc[metadata_df['newannot'] == other_cell_type, 'aggregated_cellType'] = other_cell_type

    # Condition 3: any other pair
    else:
        metadata_df['aggregated_cellType'] = metadata_df['newannot']

    anndata_obj.obs['aggregated_cellType'] = metadata_df['aggregated_cellType']
    return anndata_obj


# ------------------------------------------------------------------
# Function 2: Process species comparison
# ------------------------------------------------------------------
"""
    Run scCODA species comparison for numerator 'n' vs denominator 'd'.
    
    Parameters
    ----------
    aggregated_anndata : AnnData
        Annotated single-cell object with 'aggregated_cellType'.
    n : str
        Numerator cell type (e.g., 'Neuron').
    d : str
        Denominator cell type (reference for scCODA, e.g., 'Glia').
    tissues : list of str, optional
        If provided, only include these tissues. Default is None (all tissues).
 """
# ------------------------------------------------------------------
def process_species_comparison(aggregated_anndata, n, d, tissues=None):
    # 1. Subset metadata
    metadata_df = aggregated_anndata.obs[['id', 'Species', 'aggregated_cellType', 'Tissue']].copy()
    if tissues is not None:
        metadata_df = metadata_df[metadata_df['Tissue'].isin(tissues)].copy()

    tissue_str = "_".join(tissues) if tissues else "AllTissues"

    # 2. Aggregate counts
    count_df = metadata_df.groupby(['id', 'Species', 'aggregated_cellType']).size().reset_index(name='count')
    cell_counts = count_df.pivot_table(
        index=['id', 'Species'],
        columns='aggregated_cellType',
        values='count',
        fill_value=0
    ).reset_index()

    sample_species_lookup = metadata_df.set_index('id')['Species'].to_dict()
    cell_counts_filtered = cell_counts[
        cell_counts.apply(lambda row: row['Species'] == sample_species_lookup.get(row['id'], None), axis=1)
    ]

    counts_only = cell_counts_filtered[['id', 'Species', n, d]]
    sccoda_data = dat.from_pandas(counts_only, covariate_columns=['id', 'Species'])

    # 3. Order species and concatenate
    species_order_all = ["Human", "Chimp", "Macaque", "Marmoset", "Mouse", "Bat", "Ferret"]
    data_sorted_all_list = [sccoda_data[sccoda_data.obs['Species'] == s] for s in species_order_all if s in sccoda_data.obs['Species'].values]

    if not data_sorted_all_list:
        print("No matching species found in data!")
        return

    data_sorted_all = data_sorted_all_list[0]
    for adata in data_sorted_all_list[1:]:
        data_sorted_all = data_sorted_all.concatenate(adata)

    # 4. Visualization
    viz.stacked_barplot(data_sorted_all, feature_name='id', figsize=(25, 20))
    plt.xticks(fontsize=20, rotation=45, ha='right')
    plt.yticks(fontsize=20)
    plt.xlabel('Sample ID', fontsize=20)
    plt.ylabel('Frequency', fontsize=20)
    plt.title(f'Species_{d}_{n}_Stacked_Barplot for {tissue_str} (Ordered by Species)', fontsize=24)
    plt.subplots_adjust(bottom=0.2)
    plt.savefig(f"{n}_{d}_stacked_barplot_{tissue_str}_species_ordered.svg", format='svg')
    #plt.show()

    viz.rel_abundance_dispersion_plot(data=sccoda_data)

    # 5. Run scCODA model
    model_allTissues = mod.CompositionalAnalysis(
        data_sorted_all,
        formula=f"C(Species, Treatment('Human'))",
        reference_cell_type=d
    )
    allTissues_results = model_allTissues.sample_hmc()
    allTissues_results.set_fdr(0.1)

    # 6. Save results
    human_allTissues_raw_results = allTissues_results.credible_effects(data_sorted_all)
    allTissues_df = human_allTissues_raw_results.reset_index()
    allTissues_df.columns = ['Covariate', 'Cell Type', 'Final Parameter']
    allTissues_df.to_csv(f'scCODA_Species_{d}_{n}_{tissue_str}_HumanvsRest_results_fdr0_1.csv', index=False)

    # Save summary text
    old_stdout = sys.stdout
    sys.stdout = io.StringIO()
    allTissues_results.summary_extended()
    summary_text = sys.stdout.getvalue()
    sys.stdout = old_stdout

    with open(f"scCODA_Species_{n}_vs_{d}_{tissue_str}_HumanvsRest_summary_extended.txt", "w") as file:
        file.write(summary_text)

    print(f"Finished Species comparison for numerator={n}, denominator={d}, tissues={tissue_str}")


# ------------------------------------------------------------------
# Function 3: Process tissue comparison
# ------------------------------------------------------------------
def process_tissue_comparison(aggregated_anndata, n, d):

    # ------------------------------------------------------------------
    # 1. Extract metadata
    # ------------------------------------------------------------------
    metadata_df = aggregated_anndata.obs[
        ["id", "Species", "Tissue", "aggregated_cellType"]
    ].copy()

    species_list = ["Human", "Chimp", "Macaque", "Marmoset", "Bat"]

    # ------------------------------------------------------------------
    # 2. Loop over species (SEPARATE MODELS)
    # ------------------------------------------------------------------
    for species in species_list:

        species_df = metadata_df[metadata_df["Species"] == species].copy()

        if species_df.empty:
            print(f"Skipping {species}: no samples")
            continue

        # ------------------------------------------------------------------
        # 3. Count cells per (sample, tissue, celltype)
        # ------------------------------------------------------------------
        count_df = (
            species_df
            .groupby(["id", "Tissue", "aggregated_cellType"])
            .size()
            .reset_index(name="count")
        )

        # ------------------------------------------------------------------
        # 4. Pivot one row per sample
        # ------------------------------------------------------------------
        cell_counts = (
            count_df
            .pivot_table(
                index=["id", "Tissue"],
                columns="aggregated_cellType",
                values="count",
                fill_value=0
            )
            .reset_index()
        )

        # ------------------------------------------------------------------
        # 5. Keep ONLY the two cell types of interest
        # ------------------------------------------------------------------
        for ct in [n, d]:
            if ct not in cell_counts.columns:
                cell_counts[ct] = 0

        cell_counts = cell_counts[["id", "Tissue", n, d]]

        # Drop samples with zero total counts
        cell_counts = cell_counts[(cell_counts[n] + cell_counts[d]) > 0]

        if cell_counts.empty:
            print(f"Skipping {species}: no valid samples after filtering")
            continue

        # ------------------------------------------------------------------
        # 6. Build scCODA object
        # ------------------------------------------------------------------
        sccoda_data = dat.from_pandas(
            cell_counts,
            covariate_columns=["id", "Tissue"]
        )

        # ------------------------------------------------------------------
        # 7. Visualization
        # ------------------------------------------------------------------
        viz.stacked_barplot(
            sccoda_data,
            feature_name="id",
            figsize=(14, 6)
        )
        plt.title(f"{species}: {n} vs {d} (Caudate vs Putamen)")
        plt.tight_layout()
        plt.savefig(f"{species}_{n}_vs_{d}_stacked_barplot.svg")
        plt.close()

        # ------------------------------------------------------------------
        # 8. scCODA model
        # ------------------------------------------------------------------
        model = mod.CompositionalAnalysis(
            sccoda_data,
            formula="C(Tissue, Treatment('Putamen'))",
            reference_cell_type=d
        )

        results = model.sample_hmc()
        results.set_fdr(0.1)

        # ------------------------------------------------------------------
        # 9. Save results (credible effects + extended summary)
        # ------------------------------------------------------------------
        # Extract credible effects (full table)
        raw_results = results.credible_effects(sccoda_data)
        results_df = raw_results.reset_index()  # keep all columns
        results_df.columns = ['Covariate', 'Cell Type', 'Final Parameter']
        
        # Save CSV with all columns
        results_df.to_csv(f'scCODA_{species}_{n}_vs_{d}_stacked_barplot_results_fdr0_1.csv', index=False)

        # Save extended summary text
        old_stdout = sys.stdout
        sys.stdout = io.StringIO()
        with pd.option_context('display.max_columns', None, 'display.width', 1000):
          results.summary_extended()
        summary_text = sys.stdout.getvalue()
        sys.stdout = old_stdout

        with open(f'scCODA_{species}_{n}_vs_{d}_stacked_barplot_summary_extended.txt', "w") as file:
            file.write(summary_text)
      

        print(f"Finished Species comparison for numerator={n}, denominator={d}, species={species}")
        
# Function 5: 
def run_interneuron_subtype_sccoda(
    adata,
    target_subtype,
    tissues,
    out_prefix
):
    """
    Compare one interneuron subtype vs all other interneurons
    between Primates and Non-primates within selected tissues
    """
    PRIMATES = {"Human", "Chimp", "Macaque", "Marmoset"}

    # ----------------------------
    # Subset metadata
    # ----------------------------
    meta = adata.obs[
        ["id", "Species", "Tissue", "newannot"]
    ].copy()

    meta = meta[meta["Tissue"].isin(tissues)]

    # ----------------------------
    # Define group (Primate vs Non)
    # ----------------------------
    meta["Group"] = meta["Species"].apply(
        lambda x: "Primate" if x in PRIMATES else "Non_primate"
    )

    # ----------------------------
    # Collapse cell types
    # ----------------------------
    meta["CellType"] = meta["newannot"].apply(
        lambda x: target_subtype if x == target_subtype else "Other_Interneurons"
    )

    # ----------------------------
    # Aggregate counts per sample
    # ----------------------------
    count_df = (
        meta
        .groupby(["id", "Group", "CellType"])
        .size()
        .reset_index(name="count")
    )

    cell_counts = (
        count_df
        .pivot_table(
            index=["id", "Group"],
            columns="CellType",
            values="count",
            fill_value=0
        )
        .reset_index()
    )

    # drop samples with zero total
    cell_counts = cell_counts[
        (cell_counts[target_subtype] +
         cell_counts["Other_Interneurons"]) > 0
    ]

    # ----------------------------
    # Create scCODA object
    # ----------------------------
    sccoda_data = dat.from_pandas(
        cell_counts,
        covariate_columns=["id", "Group"]
    )

    # ----------------------------
    # Visualization
    # ----------------------------
    viz.stacked_barplot(
        sccoda_data,
        feature_name="id",
        figsize=(18, 10)
    )
    plt.savefig(f"{out_prefix}_stacked_barplot.pdf")
    plt.close()

    # ----------------------------
    # Model
    # ----------------------------
    model = mod.CompositionalAnalysis(
        sccoda_data,
        formula="C(Group, Treatment('Primate'))",
        reference_cell_type="Other_Interneurons"
    )

    results = model.sample_hmc()
    results.set_fdr(0.1)

    # ----------------------------
    # Save full summary
    # ----------------------------
     # Extract credible effects (full table)
    raw_results = results.credible_effects(sccoda_data)
    results_df = raw_results.reset_index()  # keep all columns
    results_df.columns = ['Covariate', 'Cell Type', 'Final Parameter']
    
    # Save CSV with all columns
    results_df.to_csv(f'out_prefix_stacked_barplot_results_fdr0_1.csv', index=False)
            
    # Save extended summary text
    old_stdout = sys.stdout
    sys.stdout = io.StringIO()
    with pd.option_context('display.max_columns', None, 'display.width', 1000):
      results.summary_extended()
      summary_text = sys.stdout.getvalue()
      sys.stdout = old_stdout

    with open(f'{out_prefix}_summary_extended.txt', "w") as file:
           file.write(summary_text)

 
