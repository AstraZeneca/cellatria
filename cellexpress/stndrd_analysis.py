# stndrd_analysis.py
# -------------------------------

import os
import sys
sys.path.append("/mnt/work/projects/cellatria")

# -------------------------------

import scanpy as sc
import harmonypy as hm
from cellexpress.helper import parse_vars, graph_pipeline

# -------------------------------

def run_analysis(adata, args):
    """
    Runs the standard Scanpy single-cell analysis pipeline, including normalization,
    variable gene selection, scaling, PCA, optional regression and batch correction,
    clustering, and dimensionality reduction (UMAP, optionally TSNE).

    The function also preserves:
      - raw counts (for annotation),
      - normalized counts (for visualization & DEA),
      - scaled data (for PCA/clustering),
      - harmonized and non-harmonized versions (if batch correction is enabled).

    Args:
        adata (AnnData): Combined AnnData object after sample merging.
        args (Namespace): Parsed command-line arguments including normalization, 
                          batch correction, PCA, clustering, and regression settings.

    Returns:
        tuple:
            - raw_counts (AnnData): Copy of unmodified input counts (for cell annotation).
            - adata (AnnData): Processed and optionally batch-corrected dataset.
            - adata_nohm (AnnData or None): Non-harmonized version (if batch correction applied).
    """

    # -------------------------------
    # 0. Store Raw Data (preprocessing snapshot)
    # This preserves the raw expression for celltype annotation later
    raw_counts = adata.copy()  

    # -------------------------------
    # 1. Normalization (log1p)
    print("*** 🔄 Normalizing data (log1p)...")
    sc.pp.normalize_total(adata, target_sum=args.norm_target_sum) # (default 10,000 UMIs per cell)
    sc.pp.log1p(adata)

    # -------------------------------
    # 2. HVG Selection
    print("*** 🔄 Selecting highly variable genes...")
    sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=args.n_top_genes)

    # -------------------------------
    # 3. Preserve Log-Normalized Data
    adata.raw = adata.copy() # Log1p-normalized counts saved for visualization/DEA

    # -------------------------------
    # 4. Filter dataset to only keep highly variable genes (can be skipped if desired)
    adata = adata[:, adata.var['highly_variable']]

    # -------------------------------
    # 5. Regress Out Technical Effects (optional)
    # Regress out total counts and mitochondrial percentage to remove technical bias.
    # It requires dense data, so it can be expensive on large datasets.
    if args.regress_out.lower() == "yes":
        print("*** 🔄 Regressing out total counts and pct_counts_mito...")
        sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mito'])
    else:
        print("*** 🚫 Skipping regression step (per user input).")

    # -------------------------------
    # 6. Scaling
    print(f"*** 🔄 Scaling data with max_value={args.scale_max_value}...")    
    sc.pp.scale(adata, max_value=args.scale_max_value)

    # -------------------------------
    # 7. Principal Component Analysis (PCA)
    print(f"*** 🔄 Running PCA with {args.n_pcs} components...")
    sc.tl.pca(adata, svd_solver='arpack', n_comps=args.n_pcs)

    # -------------------------------

    # 8. Batch Correction
    # Initialize adata_nohm to None — this will be returned if no batch correction happens
    adata_nohm = None
    if args.batch_correction and args.batch_correction.lower() == "harmony":        
        btchvrs = parse_vars(args.batch_vars) # Parse and validate vars_batch (comma-separated columns).        
        missing_vars = [var for var in btchvrs if var not in adata.obs.columns] # Ensure all specified columns exist in adata.obs.
        if missing_vars:
            raise ValueError(f"*** 🚨 The following columns specified in --batch_vars are missing from adata.obs: {', '.join(missing_vars)}")        

        adata_nohm = adata.copy() # Store a copy of adata before batch correction

        # Perform Harmony batch correction
        print(f"*** 🔄 Performing batch correction using Harmony on: {' and '.join(btchvrs)} ...")
        harmony_out = hm.run_harmony(adata.obsm['X_pca'], adata.obs, btchvrs)        
        adata.obsm['X_pca'] = harmony_out.Z_corr.T # Replace PCA embeddings in adata with Harmony-corrected version
        print("*** ✅ Harmony correction complete.")

    else:
        print("*** 🚫 No batch correction applied (per user input).")

    # -------------------------------
    # 9. Neighborhood graph, UMAP, TSNE, and clustering
    graph_pipeline(adata, args)

    # Run for non-harmonized data if available
    if adata_nohm is not None:
        print("*** 🔄 Processing non-harmonized (original) data...")
        graph_pipeline(adata_nohm, args, label = "no batch corrected")

    # -------------------------------
    print("*** ✅ Standard analysis pipeline completed.")
    return raw_counts, adata, adata_nohm