# qcfilter.py
# -------------------------------

import scanpy as sc
import numpy as np
import pandas as pd

# -------------------------------

def qc_filter(adatas, args):
    """
    Performs iterative quality control (QC) filtering on single-cell AnnData objects.

    The filtering loop ensures convergence between gene- and cell-level filters, addressing
    potential instability from sequential thresholding. Cells are filtered by gene count, UMI count,
    and mitochondrial gene content. Genes are filtered based on the number of expressing cells.

    Args:
        adatas (dict): Dictionary of {sample_id: AnnData} objects representing individual samples.
        args (Namespace): Parsed command-line arguments containing QC threshold parameters.

    Returns:
        qc_adatas (dict): Dictionary of filtered AnnData objects post-QC.
        qc_summary_df (pd.DataFrame): Summary DataFrame with the number of cells and genes retained.
    """

    # Extract QC thresholds from args
    min_cells = args.min_cell
    min_genes_per_cell = args.min_genes_per_cell
    max_genes_per_cell = np.inf if args.max_genes_per_cell is None else args.max_genes_per_cell
    min_umi_per_cell = args.min_umi_per_cell
    max_umi_per_cell = np.inf if args.max_umi_per_cell is None else args.max_umi_per_cell
    max_mt_percent = args.max_mt_percent
    species = args.species

    qc_adatas = {}
    qc_summary = []

    # Work with a copied list of adatas to avoid side effects
    adatas_list = [(sample_id, adata.copy()) for sample_id, adata in adatas.items()]

    for i in range(len(adatas_list)):  # Access using index
        sample_id, adata = adatas_list[i]  # Unpack values
        sample_name = adata.obs['sample'].unique()[0]
        print(f"*** 🔄 Processing QC for: {sample_id}")

        # Pre-QC stats
        initial_cells = (adata.X.sum(axis=1) > 0).sum()  # Cells with non-zero UMI counts
        initial_genes = (adata.X.sum(axis=0) > 0).sum()  # Genes with non-zero UMI counts
        print(f"*** 📊 Initial: {initial_cells:,} cells, {initial_genes:,} genes with non-zero UMI counts")

        prev_shape = None  # Store previous shape to detect convergence
        iteration = 0
        initial_cell_count = adata.n_obs  # Total cells before QC

        while True:
            # -------------------------------
            # Apply gene filtering
            sc.pp.filter_genes(adata, min_cells=min_cells)

            # -------------------------------
            # Apply cell filtering based on gene count
            sc.pp.filter_cells(adata, min_genes=min_genes_per_cell)
            if max_genes_per_cell < np.inf:
                sc.pp.filter_cells(adata, max_genes=max_genes_per_cell)

            # -------------------------------
            # Apply cell filtering based on UMI count (total_counts)
            sc.pp.filter_cells(adata, min_counts=min_umi_per_cell)
            if max_umi_per_cell < np.inf:
                sc.pp.filter_cells(adata, max_counts=max_umi_per_cell)

            # -------------------------------
            # Calculate mitochondrial % per cell
            if species == "hs":
                adata.var["mito"] = adata.var_names.str.startswith("MT-")
            elif species == "mm":
                adata.var["mito"] = adata.var_names.str.startswith("mt-")

            sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], percent_top=None, log1p=False, inplace=True)

            # -------------------------------
             # Filter cells with high mitochondrial content
            adata = adata[adata.obs["pct_counts_mito"] <= max_mt_percent, :].copy()   # copy() ensures adata is a full copy before modifying it           

            # -------------------------------
            # Convergence check (no change in shape)
            if prev_shape == adata.shape:
                print(f"*** ✅ Convergence reached after {iteration} iterations.")
                break

            prev_shape = adata.shape  # Update previous shape for next iteration
            iteration += 1
            
            # Get the number of non-zero UMI genes and cells **before filtering**
            post_filter_cells = (adata.X.sum(axis=1) > 0).sum()
            post_filter_genes = (adata.X.sum(axis=0) > 0).sum()
            print(f"*** 🔄 QC iteration {iteration}: {post_filter_cells:,} cells, {post_filter_genes:,} genes with non-zero UMI counts")

        # -------------------------------
        # Report removed cells
        removed_cells = initial_cell_count - adata.n_obs 
        print(f"*** 📊 Removed {removed_cells:,} cells in sample '{sample_id}'")

        # -------------------------------
        # Final counts after QC
        final_cells = (adata.X.sum(axis=1) > 0).sum()  # Cells with non-zero UMI counts
        final_genes = (adata.X.sum(axis=0) > 0).sum()  # Genes with non-zero UMI counts
        print(f"*** 📊 Final: {final_cells:,} cells, {final_genes:,} genes with non-zero UMI counts")

        # -------------------------------
        # Retain only non-empty AnnData objects
        if adata.shape[0] > 0 and adata.shape[1] > 0:
            qc_adatas[sample_id] = adata

            qc_summary.append({
                "sample": sample_name,
                "sample_id": sample_id,
                "type": "post-qc",
                "genes": (adata.X.sum(axis=0) > 0).sum(), # non-zero UMI Number of genes
                "cells": (adata.X.sum(axis=1) > 0).sum()  # non-zero UMI Number of cells
            })
        else:
            print(f"*** ⚠️ WARNING: {sample_name} resulted in an empty dataset post-QC. Skipping.")

    # -------------------------------
    # Compile QC summary
    qc_summary_df = pd.DataFrame(qc_summary)

    print("*** 📊 Summary Table: Number of Genes & Cells per Sample.")
    print(qc_summary_df.to_string(index=False))  # Print table without row index
    print(f"*** ✅ QC filtering completed. {len(qc_adatas)} samples retained.")

    return qc_adatas, qc_summary_df