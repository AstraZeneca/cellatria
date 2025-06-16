# merge_samples.py
# -------------------------------

import scanpy as sc
import anndata as ad

# -------------------------------

def merge_samples(adatas, args):
    """
    Merges a dictionary of AnnData objects into a single unified AnnData object.
    
    This function ensures all cell barcodes are globally unique by appending `sample_id` as a suffix.
    It then concatenates the individual AnnData objects using an outer join on gene features,
    preserving all genes across samples and aligning missing values with zeros. After merging,
    it recomputes QC metrics including mitochondrial content.

    Args:
        adatas (dict): Dictionary of {sample_id: AnnData} objects, post-QC.
        args (Namespace): Command-line arguments, must include `args.species` to infer mitochondrial genes.

    Returns:
        AnnData: A merged AnnData object with updated QC metrics and unified sample identifiers.
    """

    # -------------------------------   
    # Step 1: Make cell barcodes globally unique using sample_id suffix
    adata_list = []
    sample_ids = []

    for sample_name, adata in adatas.items():
        adata = adata.copy() # decouple from the original, so qc_adatas barcodes remains original and not suffix-ed
        sample_id = adata.obs['sample_id'].unique()[0]
        adata.obs_names = [f"{barcode}_{sample_id}" for barcode in adata.obs_names]
        adata_list.append(adata)
        sample_ids.append(sample_id)

    # ------------------------------- 
    # Step 2: Merge all samples (skip if only one)
    if len(adata_list) == 1:
        adata = adata_list[0]
        print("*** ✅ Only one sample found — skipping merge.")
    else:
        # Merge using anndata.concat
        adata = ad.concat(
            adata_list,
            join="outer",  # Ensures all genes across datasets are retained (union)
            label="sample_id",  # This creates obs column "sample_id"
            keys=sample_ids,  # Use sample_id as the key
            fill_value=0  # Fill missing genes in some samples with zeros
        )

    # -------------------------------      
    # Step 3: Recompute basic QC metrics including mitochondrial content
    print("*** 🔄 recompute basic qc metrics...")    
    if args.species == "hs":
        adata.var["mito"] = adata.var_names.str.startswith("MT-")
    elif args.species == "mm":
        adata.var["mito"] = adata.var_names.str.startswith("mt-")

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], percent_top=None, log1p=False, inplace=True)    

    # -------------------------------   
    # Step 4: Print QC summary by sample
    print(f"*** 📊 Combined dataset contains {adata.n_obs:,} cells and {adata.n_vars:,} genes across {len(adatas)} samples.")

    # Generate and print sanity check table: number of cells per sample_id
    count_table = (
        adata.obs['sample_id']
        .value_counts()
        .reset_index()
        .sort_values(by='sample_id')
    )
    print(count_table.to_string(index=False))
    print(f"*** ✅ Concatenation completed.")

    return adata
