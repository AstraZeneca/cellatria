# readin.py
# -------------------------------

import os
import sys
import argparse
import pandas as pd
import scanpy as sc
from helper import check_gene_names_format, fix_gene_names, try_load_sample_from_path, fix_file_format

# -------------------------------

def read_in(args):
    """
    Loads single-cell gene expression data for all samples described in a metadata file.

    This function supports loading 10x Genomics formatted data:
        - HDF5 format (`filtered_feature_bc_matrix.h5`)
        - Matrix Market format (`matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`)
        - Preprocessed `.h5ad` files (only one per sample directory)

    Args:
        args (Namespace): Parsed command-line arguments. Key attributes include:
            - args.input (str): Path to directory containing sample subfolders.
            - args.metadata (str): Metadata type (e.g., "sample_based").
            - args.fix_gene_names (str): Column name in `.var` to use for renaming gene symbols.
            - args.species (str): Species ID (e.g., "hs" or "mm").

    Returns:
        dict: Dictionary of AnnData objects indexed by `sample_id`.
        pd.DataFrame: Full metadata table.
        pd.DataFrame: Summary table with number of genes and cells per sample.
    """

    # -------------------------------
    # Load and standardize metadata
    metadata_df = pd.read_csv(os.path.join(args.input, "metadata.csv"))
    # Normalize column names: lowercase + replace spaces with underscores
    metadata_df.columns = metadata_df.columns.str.lower().str.replace(r"\s+", "_", regex=True)
    # Add sample_id column directly to metadata_df
    metadata_df["sample_id"] = [f"S{i+1}" for i in range(len(metadata_df))]

    # -------------------------------
    # Prepare containers for output
    adatas = {}
    summary_data = []

    # -------------------------------
    # Iterate through each sample in metadata
    for i in range(len(metadata_df)):
        sample_name = metadata_df.loc[i, "sample"]
        sample_id = metadata_df.loc[i, "sample_id"]  # Pull sample_id from updated metadata_df
        sample_path = os.path.join(args.input, sample_name)  # Adjusted to use args.input

        if not os.path.exists(sample_path):
            raise FileNotFoundError(f"🚨 Sample directory not found: {sample_path}")

        # this function makes sure data is compatible
        fix_log = fix_file_format(sample_path)
        print(fix_log)

        # 1️⃣ Try flat structure first
        adata = try_load_sample_from_path(sample_path, sample_id)

        # 2️⃣ Try nested directory if flat failed
        if adata is None:
            subdirs = [os.path.join(sample_path, d) for d in os.listdir(sample_path)
                    if os.path.isdir(os.path.join(sample_path, d))]
            if len(subdirs) == 1:
                adata = try_load_sample_from_path(subdirs[0], sample_id)

        # 3️⃣ Fail if nothing worked
        if adata is None:
            raise FileNotFoundError(f"🚨 No valid files found for sample {sample_name} in {sample_path}.")

        # -------------------------------
        # Report basic stats
        cells = (adata.X.sum(axis=1) > 0).sum()  # Cells with non-zero UMI counts
        genes = (adata.X.sum(axis=0) > 0).sum()  # Genes with non-zero UMI counts
        print(f"*** 📊 Loaded {sample_id} with {cells:,} cells, {genes:,} genes with non-zero UMI counts.")

        # -------------------------------
        # Assign sample metadata to adata.obs
        adata.obs['sample_id'] = sample_id  # This is critical for merging later
        for col in metadata_df.columns:
            adata.obs[col] = metadata_df.loc[i, col]
        
        # Convert low-cardinality columns to categorical
        categorical_cols = metadata_df.columns[(metadata_df.nunique() < 30)].tolist()
        for col in categorical_cols:
            adata.obs[col] = adata.obs[col].astype("category")

        # -------------------------------
        # Optional gene name fixing
        if args.fix_gene_names != "no":
            adata = fix_gene_names(adata, column_name=args.fix_gene_names)
        else: 
            check_gene_names_format(adata)

        # -------------------------------
        # Annotate mitochondrial genes and calculate QC metrics
        if args.species == "hs":
            adata.var["mito"] = adata.var_names.str.startswith("MT-")
        elif args.species == "mm":
            adata.var["mito"] = adata.var_names.str.startswith("mt-")
        # Calculate QC metrics, including pct_counts_mito
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], percent_top=None, log1p=False, inplace=True)
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], percent_top=None, log1p=False, inplace=True)
        # -------------------------------
        # Store loaded object
        adatas[sample_id] = adata
        
        # Add sample-level summary
        summary_data.append({
            "sample": sample_name, 
            "sample_id": sample_id,
            "type": "pre-qc",  # Indicate that this summary is before QC filtering
            "genes": (adata.X.sum(axis=0) > 0).sum(),  # non-zero UMI Number of genes
            "cells": (adata.X.sum(axis=1) > 0).sum()  # non-zero UMI Number of cells
            })

    # -------------------------------
    # Summarize and report
    summary_df = pd.DataFrame(summary_data)
    print("*** 📊 Summary Table: Number of Genes & Cells per Sample.")
    print(summary_df.to_string(index=False))  # Print table without row index
    print(f"*** ✅ Successfully loaded {len(adatas)} samples.")

    return adatas, metadata_df, summary_df
