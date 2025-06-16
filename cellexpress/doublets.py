# doublets.py
# -------------------------------

import scrublet as scr
import scanpy as sc
import pandas as pd

# -------------------------------

def doublets_id(adatas, args):
    """
    Identifies and optionally removes doublets using the specified method (currently supports Scrublet).

    Doublets are detected on a per-sample basis. Doublet scores and predictions are stored in 
    `adata.obs`. Cells with doublet scores above a specified threshold are removed. 
    A QC summary is printed and returned along with Scrublet scores.

    Args:
        adatas (dict): Dictionary of {sample_id: AnnData} objects (one per sample).
        args (Namespace): Parsed command-line arguments, must include:
            - doublet_method (str): Currently only 'scrublet' is supported.
            - scrublet_cutoff (float): Score threshold above which doublets are removed.

    Returns:
        qc_adatas (dict): Dictionary of filtered AnnData objects post-doublet removal.
        qc_summary_df (pd.DataFrame): QC summary table showing genes/cells per sample.
        scrublet_scores (dict): Per-sample Scrublet score distributions for observed and simulated cells.
    """

    qc_summary = []         # List to collect QC stats
    qc_adatas = {}          # Output dictionary of filtered AnnData
    scrublet_scores = {}    # Store Scrublet scores per sample

    # Convert adatas to list for incremental access
    adatas_list = list(adatas.items()) # For index-based access

    for i in range(len(adatas_list)): 
        sample_id, adata = adatas_list[i]
        sample_name = adata.obs['sample'].unique()[0]
        print(f"*** 🔄 Running {args.doublet_method} doublet detection for: {sample_id}")

        # --------- Scrublet Method ---------
        if args.doublet_method == "scrublet":
            # Initialize Scrublet
            scrub = scr.Scrublet(adata.X, random_state=123)

            # Run doublet detection
            doublet_scores, predicted_doublets = scrub.scrub_doublets()

            # Store Scrublet outputs in AnnData
            adata.obs['doublet_scores_obs'] = doublet_scores
            adata.obs['predicted_doublet'] = predicted_doublets
            # print(adata.obs.columns)  

            # Store Scrublet scores externally (for export or plotting)
            scrublet_scores[sample_id] = {
                "doublet_scores_obs": doublet_scores.tolist(),
                "doublet_scores_sim": scrub.doublet_scores_sim_.tolist()
            }

            # Filter cells with doublet scores above the cutoff
            initial_cell_count = adata.n_obs  # Total cells before QC
            # Filter doublets based on user-defined cutoff
            adata = adata[adata.obs['doublet_scores_obs'] <= args.scrublet_cutoff].copy()
            removed_cells = initial_cell_count - adata.n_obs 
            print(f"*** 📊 Removed {removed_cells:,} cells with doublet scores above {args.scrublet_cutoff} in sample '{sample_name}'")

            # Store filtered AnnData object
            qc_adatas[sample_id] = adata  

            # Report final stats
            final_cells = (adata.X.sum(axis=1) > 0).sum()  # Cells with non-zero UMI counts
            final_genes = (adata.X.sum(axis=0) > 0).sum()  # Genes with non-zero UMI counts
            print(f"*** 📊 Final: {final_cells:,} cells, {final_genes:,} genes with non-zero UMI counts")

            # Save summary for this sample
            qc_summary.append({
                "sample": sample_name,
                "sample_id": sample_id,
                "type": "post-qc",
                "genes": (adata.X.sum(axis=0) > 0).sum(),  # non-zero UMI Number of genes
                "cells": (adata.X.sum(axis=1) > 0).sum()   # non-zero UMI Number of cells
            })

        else:
            raise ValueError(f"🚨 Unsupported doublet identification method: {args.doublet_method}")

    # Convert QC stats to DataFrame
    qc_summary_df = pd.DataFrame(qc_summary)
    
    print("*** 📊 Summary Table: Number of Genes & Cells per Sample.")
    print(qc_summary_df.to_string(index=False))  # Print table without row index
    print(f"*** ✅ Doublet identification & filtering completed using {args.doublet_method}. {len(qc_adatas)} samples retained.")

    return qc_adatas, qc_summary_df, scrublet_scores