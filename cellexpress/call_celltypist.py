# annot_celltypist.py
# -------------------------------

import os
import sys
import pandas as pd
import scanpy as sc
import celltypist
from celltypist import models
from helper import ontology_map, summary_by_abundance

# -------------------------------

def run_celltypist(raw_counts, adata, args):
    """
    Perform CellTypist-based cell type annotation.

    This function:
      - Loads a CellTypist model from the specified path.
      - Normalizes the raw count data (log1p-normalized to 10,000 counts per cell).
      - Annotates cells using the CellTypist model.
      - Stores cell state predictions and maps them to higher-level cell types using an ontology.
      - Adds results to `adata.obs["cellstate_celltypist"]` and `adata.obs["celltype_celltypist"]`.

    Args:
        raw_counts (AnnData): AnnData object containing raw UMI counts.
        adata (AnnData): AnnData object to store celltype predictions (shared obs).
        args (Namespace): Command-line arguments with:
            - cty_model_path (str): Path to the model directory.
            - cty_model_name (str): Name of the CellTypist model (without `.pkl` extension).

    Returns:
        AnnData: The `adata` object with CellTypist annotations added.
    """

    # -------------------------------
    # Construct full path to model file
    model_file = os.path.join(args.cty_model_path, args.cty_model_name + ".pkl")
    
    # Load the CellTypist model
    model = models.Model.load(model = model_file)
    
    print(f"*** 🔍 Using CellTypist model: '{args.cty_model_name}'")
    print(f"*** 📂 Model path: {model_file}")
    print(f"*** 📂 Model details:\n{model}")

    # -------------------------------
    # Create new AnnData object with raw counts, same obs and var (metadata)    
    adata_tmp = raw_counts.copy()

    # Ensure CellTypist expectation to get log1p normalized expression to 10000 counts per cell;
    # Normalize data (log1p-normalized to 10,000 UMI per cell)
    print("*** 🔄 Normalizing data (log1p norm to 1e4 per cell)...")
    # Normalize each cell to the same total UMI count (default 10,000 UMIs per cell)
    # This is to address CellTypist warning: 
    # 👀 Invalid expression matrix in `.X`, expect log1p normalized expression to 10000 counts per cell; will use `.raw.X` instead
    sc.pp.normalize_total(adata_tmp, target_sum=1e4)
    # Apply natural log (log1p) transformation to stabilize variance
    sc.pp.log1p(adata_tmp)

    # -------------------------------
    # Run CellTypist annotation
    predictions = celltypist.annotate(adata_tmp, model=model_file)

    # -------------------------------
    # Store cell state predictions
    adata.obs["cellstate_celltypist"] = predictions.predicted_labels

    # -------------------------------
    # Map predicted states to higher-level cell types
    unmapped_states = set(adata.obs["cellstate_celltypist"].unique()) - set(ontology_map.keys())
    # Print a warning if there are unmapped cell states
    if unmapped_states:
        print(f"*** ⚠️  WARNING: The following cell states do not exist in the ontology map and will remain unchanged: {', '.join(unmapped_states)}")

    adata.obs["celltype_celltypist"] = adata.obs["cellstate_celltypist"].map(ontology_map).fillna(adata.obs["cellstate_celltypist"])    

    # -------------------------------
    # Convert predictions to categorical for efficiency
    adata.obs["cellstate_celltypist"] = adata.obs["cellstate_celltypist"].astype("category")
    adata.obs["celltype_celltypist"] = adata.obs["celltype_celltypist"].astype("category")

    # -------------------------------    
    # Summarize annotation results
    summary_by_abundance(adata, "celltype_celltypist")
    summary_by_abundance(adata, "cellstate_celltypist")

    # -------------------------------
    print("*** ✅ Celltypist annotation completed successfully.")
    return adata