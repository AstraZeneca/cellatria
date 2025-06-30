# annot_scimilarity.py
# -------------------------------

import os
import sys
import pandas as pd
import scanpy as sc
from scimilarity.utils import lognorm_counts, align_dataset
from scimilarity import CellAnnotation
from helper import ontology_map, summary_by_abundance

# -------------------------------

def run_scimilarity(raw_counts, adata, args):
    """
    Perform cell type annotation using SCimilarity.

    This function:
      - Initializes the SCimilarity model using the provided model path.
      - Aligns the raw count data to the model's gene order.
      - Applies SCimilarity-specific normalization (log-normalized counts).
      - Computes SCimilarity embeddings.
      - Predicts cell states using k-nearest neighbors.
      - Maps predicted cell states to higher-level cell types using a predefined ontology.
      - Stores results in `adata.obs["cellstate_scimilarity"]` and `adata.obs["celltype_scimilarity"]`.

    Args:
        raw_counts (AnnData): Unmodified, raw-count AnnData object.
        adata (AnnData): Processed AnnData object used for downstream integration.
        args (Namespace): Parsed command-line arguments, including `sci_model_path`.

    Returns:
        AnnData: Updated `adata` object with SCimilarity-based annotations.
    """

    # -------------------------------
    print("*** 🔄 Initializing SCimilarity for cell annotation...")

    # -------------------------------
    # Load the SCimilarity model
    ca = CellAnnotation(model_path=args.sci_model_path)
    overlap = len(set(adata.var_names) & set(ca.gene_order))
    print(f"*** 🔔 Gene overlap with SCimilarity reference: {overlap}")

    # -------------------------------
    # Assign raw counts to a new 'counts' layer for SCimilarity compatibility
    adata_tmp = raw_counts.copy()
    adata_tmp.layers["counts"] = adata_tmp.X.copy()

    # -------------------------------
    # Align the gene order to the SCimilarity model
    print("*** 🔄 Aligning dataset to SCimilarity gene order...")
    adata_tmp = align_dataset(adata_tmp, ca.gene_order)

    # -------------------------------
    # Normalize counts using SCimilarity's method
    print("*** 🔄 Applying SCimilarity normalization (lognorm tp10k)...")
    adata_tmp = lognorm_counts(adata_tmp)

    # -------------------------------
    # Compute SCimilarity embeddings
    print("*** 🔄 Computing SCimilarity embeddings...")
    adata_tmp.obsm["X_scimilarity"] = ca.get_embeddings(adata_tmp.X)

    # -------------------------------
    # Perform cell annotation (unconstrained)
    print("*** 🔄 Running unconstrained cell type annotation...")
    predictions, nn_idxs, nn_dists, nn_stats = ca.get_predictions_knn(adata_tmp.obsm["X_scimilarity"])

    # -------------------------------
    # Store predictions in adata.obs (ensure correct index alignment)
    print("*** 🔄 Storing predictions in adata.obs['predictions_unconstrained']...")
    adata.obs["cellstate_scimilarity"] = pd.Series(predictions.values, index=adata.obs.index)

    # -------------------------------
    # Map to ontology-based cell type using a dictionary
    unmapped_states = set(adata.obs["cellstate_scimilarity"].unique()) - set(ontology_map.keys())
    if unmapped_states:
        print(f"*** ⚠️  WARNING: The following cell states do not exist in the ontology map and will remain unchanged: {', '.join(unmapped_states)}")

    adata.obs["celltype_scimilarity"] = adata.obs["cellstate_scimilarity"].map(ontology_map).fillna(adata.obs["cellstate_scimilarity"])

    # -------------------------------
    # Cast annotations to categorical for efficiency and clarity
    adata.obs["cellstate_scimilarity"] = adata.obs["cellstate_scimilarity"].astype("category")
    adata.obs["celltype_scimilarity"] = adata.obs["celltype_scimilarity"].astype("category")

    # -------------------------------
    # Display abundance summaries
    summary_by_abundance(adata, "celltype_scimilarity")
    summary_by_abundance(adata, "cellstate_scimilarity")

    # -------------------------------
    print("*** ✅ SCimilarity annotation completed successfully.")
    return adata
