# dea_analysis.py
# -------------------------------

import warnings
import scanpy as sc
import pandas as pd
import pandas as pd
import numpy as np

# -------------------------------

def compute_degs(adata, groupby, dea_method="wilcoxon", 
                 pval_threshold=0.05, logfc_threshold=0.25, pts_threshold=0.1,
                 pts=True, n_genes=100):
    """
    Computes differentially expressed genes (DEGs) across groups in a single-cell dataset.

    This function filters and returns significant DEGs for each group using the specified method
    (e.g., Wilcoxon), with thresholds for p-value, log fold-change, and percentage of cells expressing the gene.

    Args:
        adata (AnnData): Processed AnnData object with normalized/log-transformed data.
        groupby (str): Column name in `adata.obs` used to define groups for comparison.
        dea_method (str): Statistical test to use ("wilcoxon", "t-test", etc.). Default is "wilcoxon".
        pval_threshold (float): Adjusted p-value threshold for significance. Default is 0.05.
        logfc_threshold (float): Minimum log fold-change for DEGs. Default is 0.25.
        pts_threshold (float): Minimum percentage of expressing cells. Default is 0.1.
        pts (bool): Whether to compute percentage of cells expressing each gene. Default is True.
        n_genes (int): Maximum number of DEGs to keep per group. Default is 100.

    Returns:
        AnnData: The input `adata` object with filtered DEGs stored in `adata.uns['degs_filtered_<groupby>']`.
    """
    # -------------------------------

    if groupby not in adata.obs.columns:
        raise ValueError(f"*** 🚨 Grouping column '{groupby}' not found in adata.obs.")

    # -------------------------------
    # Count cells per group and exclude groups with <2 cells
    group_counts = adata.obs[groupby].value_counts()
    small_groups = group_counts[group_counts < 2].index.tolist()

    # Determine if need to subset
    if small_groups:
        print(f"*** ⚠️  Skipping DE analysis for groups with < 2 cells: {', '.join(small_groups)}")
        adata_used = adata[adata.obs[groupby].isin(group_counts[group_counts >= 2].index), :].copy()
    else:
        # No filtering needed, use original adata
        adata_used = adata

    # -------------------------------

    print(f"*** 🔄 Running DEA using {dea_method} test for {groupby} groups...")
    key_name = f"degs_filtered_{groupby}"

    # -------------------------------

    # Suppress pandas PerformanceWarnings globally: "PerformanceWarning: DataFrame is highly fragmented."
    warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

    # Perform DE analysis with top n_genes
    sc.tl.rank_genes_groups(
        adata_used,
        groupby=groupby,
        method=dea_method,
        pts=pts,
        key_added=key_name  # Store results under unique key
    )

    # Re-enable warnings after DE computation
    warnings.simplefilter(action="default", category=pd.errors.PerformanceWarning)

    # -------------------------------
    # Extract DE results from AnnData
    degs = []
    gene_names = adata_used.uns[key_name]["names"]
    logFC_vals = adata_used.uns[key_name]["logfoldchanges"]
    pvals = adata_used.uns[key_name]["pvals"]
    pvals_adj = adata_used.uns[key_name]["pvals_adj"]
    pts_vals = adata_used.uns[key_name]["pts"] if pts else None

    for group in gene_names.dtype.names:
        num_genes = len(gene_names[group])  # Ensure uniform length
        group_degs = pd.DataFrame({
            "gene": list(gene_names[group]),
            "logFC": list(logFC_vals[group])[:num_genes],  
            "pval": list(pvals[group])[:num_genes],  
            "pval_adj": list(pvals_adj[group])[:num_genes],  
            "pts": list(pts_vals[group])[:num_genes] if pts_vals is not None else None,  
            "group": group
        })
        
        degs.append(group_degs)

    # Concatenate all results **before filtering**
    degs = pd.concat(degs, ignore_index=True)

    # -------------------------------
    # Apply filtering thresholds and retain top N genes per group
    degs_filtered = (
        degs[
            (degs["pval_adj"] < pval_threshold) &
            (degs["logFC"] > logfc_threshold) &
            (degs["pts"] > pts_threshold)
        ]
        .sort_values(by=["group", "logFC"], ascending=[True, False])  # Sort within each group
        .groupby("group")  # Apply per-group filtering
        .head(n_genes)  # Keep only top N DEGs per group
    )

    # Warn if some groups have no significant DEGs
    all_groups = set(degs["group"].unique())
    filtered_groups = set(degs_filtered["group"].unique())
    missing_groups = all_groups - filtered_groups
    if missing_groups:
        print(f"*** ⚠️  Warning: No significant DEGs found for {len(missing_groups)} group(s): {', '.join(missing_groups)}")

    # -------------------------------
    # Store filtered DEGs in AnnData for report export
    adata.uns[key_name] = degs_filtered

    # -------------------------------
    print(f"*** 📁 Filtered DEGs stored in `adata.uns['{key_name}']`.")
    print(f"*** ✅ DE analysis complete. Identified {len(degs_filtered):,} significant DEGs.")    
    return adata
