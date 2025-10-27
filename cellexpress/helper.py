# helper.py
# -------------------------------

import os
import sys
import importlib
import random
import string
import platform
import pkg_resources
import re
import requests
import urllib3
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import json
import gzip
import shutil
import glob

# -------------------------------

def fix_duplicate_gene_names(adata):
    """
    Resolve duplicate gene names in an AnnData object.

    Checks whether `adata.var_names` contains duplicate gene symbols and 
    ensures uniqueness by appending suffixes (e.g., '-1', '-2') using 
    `AnnData.var_names_make_unique()`.

    Args:
        adata (AnnData): The annotated data matrix containing gene expression profiles.

    Returns:
        AnnData: A modified object where `var_names` are guaranteed to be unique.

    Notes:
        - Only the `var_names` index is changed; underlying data and metadata remain intact.
        - A summary of detected duplicates (up to 10 examples) will be printed if any are found.
    """

    # Identify duplicate gene names
    duplicate_genes = adata.var_names[adata.var_names.duplicated()]

    if len(duplicate_genes) > 0:
        print(f"*** ⚠️  Warning: Found {len(duplicate_genes)} duplicate gene names!")
        print("*** ⚠️  Duplicate Genes:", ", ".join(duplicate_genes[:10]))  # Print first 10 duplicates
        if len(duplicate_genes) > 10:
            print(f"... and {len(duplicate_genes) - 10} more.")

        # Make gene names unique
        adata.var_names_make_unique()
        print("*** ✅ Gene names have been made unique.")

    return adata  # Return fixed AnnData object

# -------------------------------

def create_unique_ids(n, char_len=5, seed_no=None):
    """
    Generate a list of unique alphanumeric IDs.

    Each ID consists of randomly selected characters from [A-Za-z0-9]. Uniqueness is guaranteed
    by storing the results in a set until the desired number is reached.

    Args:
        n (int): Number of unique IDs to generate.
        char_len (int, optional): Length of each ID. Default is 5.
        seed_no (int, optional): Random seed to ensure reproducibility. Default is None.

    Returns:
        list of str: List containing `n` unique alphanumeric ID strings.

    Raises:
        ValueError: If `n` is negative or too large relative to `char_len`.
    """

    if n <= 0:
        raise ValueError("*** 🚨 Number of IDs (`n`) must be greater than zero.")
    if char_len <= 0:
        raise ValueError("*** 🚨 `char_len` must be greater than zero.")

    if seed_no is not None:
        random.seed(seed_no)  # Set seed for reproducibility

    pool = string.ascii_letters + string.digits  # A-Z, a-z, 0-9
    unique_ids = set()  # Use a set to avoid duplicates

    # In theory, uniqueness can't be guaranteed beyond len(pool)**char_len
    max_possible = len(pool) ** char_len
    if n > max_possible:
        raise ValueError(f"Cannot generate {n} unique IDs of length {char_len} — maximum possible is {max_possible}.")

    while len(unique_ids) < n:
        new_id = "".join(random.choices(pool, k=char_len))
        unique_ids.add(new_id)  # Ensures uniqueness automatically

    return list(unique_ids)  # Convert set to list before returning

# -------------------------------

def convert_none(value):
    """
    Convert a string 'None' (case-insensitive) to Python None.

    This utility is useful when parsing user-provided input or configuration
    values where the string "None" should be interpreted as a Python None object.

    Args:
        value (Any): The input value to check.

    Returns:
        Any: Returns None if the input is a case-insensitive string "None"; 
             otherwise, returns the original input value unchanged.
    """
    return None if isinstance(value, str) and value.strip().lower() == "none" else value

# -------------------------------

def none_or_int(value):
    value = convert_none(value)
    return None if value is None else int(value)

def none_or_float(value):
    value = convert_none(value)
    return None if value is None else float(value)

def none_or_str(value):
    value = convert_none(value)
    return None if value is None else str(value)

# -------------------------------

def total_unique_genes(adatas):
    """
    Calculate the total number of unique genes across multiple AnnData objects.

    Args:
        adatas (dict): Dictionary mapping sample names to AnnData objects.

    Returns:
        int: Total count of unique gene names across all samples.
    """
    unique_genes = set()
    for adata in adatas.values():
        unique_genes.update(adata.var_names)
    return len(unique_genes)

# -------------------------------

def parse_vars(vars_str):
    """
    Parse a comma-separated string into a list of clean column names.

    Args:
        vars_str (str): A comma-separated string of column names.

    Returns:
        list: A list of column names, stripped of leading/trailing whitespace.
              Returns an empty list if the input is None or empty.
    """
    if vars_str is None or vars_str.strip() == "":
        return []
    return [col.strip() for col in vars_str.split(",") if col.strip()]

# -------------------------------

def graph_pipeline(adata, args, label=""):
    """
    Construct neighborhood graph, run clustering, and compute embeddings (UMAP/t-SNE) on an AnnData object.

    Args:
        adata (AnnData): The AnnData object to process.
        args (Namespace): Parsed command-line arguments containing parameters like neighbors, resolution, etc.
        label (str, optional): Label to distinguish outputs (e.g., "nohm" for non-harmonized). Default is "".
    """
    suffix = f" ({label})" if label else ""

    # Step 1: Build k-NN graph
    print(f"*** 🔄 Computing neighborhood graph{suffix}...")
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs)

    # Step 2: Perform Leiden clustering
    print(f"*** 🔄 Running Leiden clustering with resolution {args.resolution}{suffix} ...")
    sc.tl.leiden(adata, resolution=args.resolution, flavor="igraph", n_iterations=2, directed=False)

    # Step 3: Reorder clusters by size (optional utility function)
    print(f"*** 🔄 reordering clusters by size...")
    reorder_clusters(adata)

    # rename the "leiden" column in adata.obs to "leiden_cluster"
    adata.obs.rename(columns={"leiden": "leiden_cluster"}, inplace=True)

    # Step 5: Report clustering summary
    # Get cluster sizes
    cluster_counts = adata.obs["leiden_cluster"].value_counts().sort_values(ascending=False)
    print(f"*** 📊 Identified {len(cluster_counts)} clusters.")
    print(f"*** 📊 Largest cluster: {cluster_counts.idxmax()} ({cluster_counts.max():,} cells)")
    print(f"*** 📊 Smallest cluster: {cluster_counts.idxmin()} ({cluster_counts.min():,} cells)")

    # Step 6: Compute UMAP embedding
    print(f"*** 🔄 Computing UMAP{suffix}...")
    sc.tl.umap(adata)   

    # Step 7: Optionally compute t-SNE
    if args.compute_tsne.lower() == "yes":
        print(f"*** 🔄 Computing TSNE{suffix}...")
        sc.tl.tsne(adata, n_pcs=args.n_pcs)
    else:
        print(f"*** 🚫 Skipping TSNE embedding (per user input){suffix}.")

    # Final Step: Report pipeline completion
    print(f"*** ✅ Graph construction, clustering, and embeddings complete{suffix}.")

# -------------------------------

def summary_by_abundance(adata, column):
    """ 
    Summarizes the abundance of annotations in a specified column of adata.obs.

    Args:
        adata (AnnData): AnnData object containing single-cell expression data.
        column (str): Column in adata.obs to summarize (e.g., 'cell_type', 'leiden_cluster').

    Returns:
        pd.DataFrame: DataFrame with counts of each category, sorted in descending order.
    """

    if column not in adata.obs.columns:
        raise ValueError(f"🚨 Column '{column}' not found in adata.obs!")

    # Step 1: Count occurrences of each category in the specified column
    category_counts = adata.obs[column].value_counts()

    # Step 2: Convert to DataFrame for better visualization
    summary_df = category_counts.reset_index()
    summary_df.columns = ["Cell Type", "Count"]

    # Step 3: Display the table nicely
    print(f"*** 📊 Annotated Cells ({column})")
    print(summary_df.to_string(index=False))  # Display without index

# -------------------------------

def compute_qc_stats_objs(metadata_df, adatas):
    """
    Computes per-sample quality control (QC) summaries including:
    - UMI counts per cell
    - Number of genes per cell
    - Mitochondrial percentage per cell

    Args:
        metadata_df (pd.DataFrame): DataFrame containing metadata for each sample, 
            including at least 'sample' and 'sample_id' columns.
        adatas (dict): Dictionary mapping sample IDs to AnnData objects (after QC).

    Returns:
        tuple:
            - List of dictionaries for UMI per cell summaries (count_by_cell_ls)
            - List of dictionaries for gene per cell summaries (gene_by_cell_ls)
            - List of dictionaries for mitochondrial percentage summaries (pct_mito_ls)

            Each list is structured for JSON serialization and includes per-sample stats
            like min, max, quantiles, etc.
    """

    count_by_cell = {}
    gene_by_cell = {}
    pct_mito = {}

    for i in range(len(metadata_df)):
        sample_name = metadata_df.loc[i, "sample"]
        sample_id = metadata_df.loc[i, "sample_id"]
        if sample_id not in adatas:
            print(f"*** ⚠️ Sample '{sample_id}' was removed during QC and will be skipped.")
            continue
        adata = adatas[sample_id]

        # Extract QC summary statistics
        count_by_cell[sample_id], gene_by_cell[sample_id], pct_mito[sample_id] = qc_quantiles(adata, sample_name = sample_name, sample_id = sample_id)

    # Convert to list of records for JSON compatibility
    count_by_cell_ls = [{"sample": k, **v} for k, v in count_by_cell.items()]
    gene_by_cell_ls = [{"sample": k, **v} for k, v in gene_by_cell.items()]
    pct_mito_ls = [{"sample": k, **v} for k, v in pct_mito.items()]

    return count_by_cell_ls, gene_by_cell_ls, pct_mito_ls

# -------------------------------

def compute_qc_stats_obj(adata):
    """
    Computes QC summary statistics for a single AnnData object, including:
    - Total UMI counts per cell
    - Number of genes detected per cell
    - Mitochondrial gene expression percentage per cell

    Args:
        adata (AnnData): Merged and QC-filtered AnnData object.

    Returns:
        tuple:
            - List of dictionaries summarizing UMI counts per cell.
            - List of dictionaries summarizing gene counts per cell.
            - List of dictionaries summarizing mitochondrial percentage per cell.
            
            Each list contains records in the form: {"metric": ..., "value": ...},
            making them ready for JSON serialization.
    """

    count_by_cell, gene_by_cell, pct_mito = qc_quantiles(adata)

    # Format summaries as list-of-dictionaries for JSON compatibility
    def to_list(summary_dict):
        return [{"metric": k, "value": v} for k, v in summary_dict.items()]

    return to_list(count_by_cell), to_list(gene_by_cell), to_list(pct_mito)

# -------------------------------

def qc_quantiles(adata, sample_name = None, sample_id = None):
    """
    Computes quantile statistics for QC metrics of an AnnData object:
    - Total UMI counts per cell
    - Number of genes per cell
    - Mitochondrial gene percentage per cell

    Args:
        adata (AnnData): The AnnData object to compute statistics on.
        sample_name (str, optional): Descriptive name of the sample (for tracking/reporting).
        sample_id (str, optional): Unique sample identifier (for tracking/reporting).

    Returns:
        tuple:
            - count_by_cell (dict): Summary statistics of total UMI counts per cell.
            - gene_by_cell (dict): Summary statistics of genes detected per cell.
            - pct_mito (dict): Summary statistics of mitochondrial gene percentage per cell.

        Each dictionary contains:
            - sample and sample_id (if provided)
            - min, q0, q25, q50, q75, q100, max values
    """

    # UMI per cell
    count_by_cell = np.array(adata.obs["total_counts"]) # The total number of raw counts (UMIs or reads)
    count_by_cell = {
        "sample": sample_name,
        "sample_id": sample_id,
        "min": float(np.min(count_by_cell)),
        "q0": float(np.quantile(count_by_cell, 0.0)),
        "q25": float(np.quantile(count_by_cell, 0.25)),
        "q50": float(np.quantile(count_by_cell, 0.5)),
        "q75": float(np.quantile(count_by_cell, 0.75)),
        "q100": float(np.quantile(count_by_cell, 1.0)),
        "max": float(np.max(count_by_cell))
    }

    # Gene per cell
    gene_by_cell = np.array(adata.obs["n_genes_by_counts"]) # The number of genes with non-zero counts per cell
    gene_by_cell = {
        "sample": sample_name,
        "sample_id": sample_id,
        "min": float(np.min(gene_by_cell)),
        "q0": float(np.quantile(gene_by_cell, 0.0)),
        "q25": float(np.quantile(gene_by_cell, 0.25)),
        "q50": float(np.quantile(gene_by_cell, 0.5)),
        "q75": float(np.quantile(gene_by_cell, 0.75)),
        "q100": float(np.quantile(gene_by_cell, 1.0)),
        "max": float(np.max(gene_by_cell))
    }

    # Mitochondrial (%) per cell
    pct_mito = np.array(adata.obs["pct_counts_mito"])
    pct_mito = {
        "sample": sample_name,
        "sample_id": sample_id,
        "min": float(np.min(pct_mito)),
        "q0": float(np.quantile(pct_mito, 0.0)),
        "q25": float(np.quantile(pct_mito, 0.25)),
        "q50": float(np.quantile(pct_mito, 0.5)),
        "q75": float(np.quantile(pct_mito, 0.75)),
        "q100": float(np.quantile(pct_mito, 1.0)),
        "max": float(np.max(pct_mito))
    }

    return count_by_cell, gene_by_cell, pct_mito

# -------------------------------

def generate_qc_plots_and_filters(adatas, metadata_df, args):
    """
    Generates long-format QC data for plotting and defines threshold-based filters.

    Args:
        adatas (dict): Dictionary of {sample_id: AnnData} containing pre-QC data.
        metadata_df (pd.DataFrame): Metadata DataFrame with at least 'sample' and 'sample_id' columns.
        args (Namespace): Argument object containing QC parameters such as:
            - min_umi_per_cell, max_umi_per_cell
            - min_genes_per_cell, max_genes_per_cell
            - max_mt_percent

    Returns:
        tuple: A tuple containing:
            - db_plots (list of dict): Long-format data for QC plotting (JSON-serializable).
            - filters (list of dict): QC threshold filters for plotting overlays (JSON-serializable).
    """

    # Step 1: Create long-format QC data for all samples
    db_plots = []

    # Iterate through all original (pre-QC) samples
    for i in range(len(metadata_df)):
        sample_name = metadata_df.loc[i, "sample"]
        sample_id = metadata_df.loc[i, "sample_id"]
        adata = adatas[sample_id]  # Use pre-QC adata

        # Extract per-cell metrics into a wide-format DataFrame
        sample_data = pd.DataFrame({
            "barcode": adata.obs.index,
            "nCount_RNA": adata.obs["total_counts"],
            "nFeature_RNA": adata.obs["n_genes_by_counts"],
            "percent_mt": adata.obs["pct_counts_mito"],
            "sample_id": sample_id
            })

        # Melt to long-format (type = metric, feature = value)
        sample_data = sample_data.melt(
            id_vars=["barcode", "sample_id"],
            value_vars=["nCount_RNA", "nFeature_RNA", "percent_mt"],
            var_name="type",
            value_name="feature"
            )

        # Append to combined db_plot
        db_plots.append(sample_data)

    # Combine all samples into a single DataFrame
    db_plots = pd.concat(db_plots, ignore_index=True)

    # Sanitize NaNs BEFORE converting to records (JSON cannot handle NaN)
    db_plots = db_plots.replace({np.nan: None})

    # Convert to JSON-serializable list
    db_plots = db_plots.to_dict(orient="records")

    # Step 2: Create QC filter thresholds
    filters = []
    for sample_id in metadata_df["sample_id"]:
        filters.extend([
            {"sample_id": sample_id, "type": "nCount_RNA", "threshold": float(args.min_umi_per_cell[sample_id])},
            {"sample_id": sample_id, "type": "nCount_RNA", "threshold": float(args.max_umi_per_cell[sample_id])},
            {"sample_id": sample_id, "type": "nFeature_RNA", "threshold": float(args.min_genes_per_cell[sample_id])},
            {"sample_id": sample_id, "type": "nFeature_RNA", "threshold": float(args.max_genes_per_cell[sample_id])},
            {"sample_id": sample_id, "type": "percent_mt", "threshold": float(args.max_mt_percent[sample_id])}
        ])

    return db_plots, filters

# -------------------------------

def compute_barcode_overlap_matrices(adatas):
    """
    Computes barcode overlap metrics between all pairs of samples.

    Args:
        adatas (dict): Dictionary of {sample_name: AnnData} objects.

    Returns:
        - bc_qc_raw (list of dict): Raw barcode overlap counts (JSON-serializable).
        - bc_qc_ji (list of dict): Jaccard index (% overlap) between samples (JSON-serializable).
    """

    samples = list(adatas.keys())
    n_sample = len(samples)

    # Initialize empty matrices
    bc_qc_raw = pd.DataFrame(np.zeros((n_sample, n_sample), dtype=int), 
                             index=samples, columns=samples)

    bc_qc_ji = pd.DataFrame(np.zeros((n_sample, n_sample), dtype=float), 
                            index=samples, columns=samples)

    # Compute pairwise barcode overlaps and Jaccard indices
    for i, sample_i in enumerate(samples):
        for j, sample_j in enumerate(samples):
            if i <= j:  # Only need upper triangle (saves some time)
                barcodes_i = set(adatas[sample_i].obs_names)
                barcodes_j = set(adatas[sample_j].obs_names)

                # Intersection size
                intersect_count = len(barcodes_i.intersection(barcodes_j))

                # Fill symmetric matrices
                bc_qc_raw.loc[sample_i, sample_j] = intersect_count
                bc_qc_raw.loc[sample_j, sample_i] = intersect_count

                # Jaccard index relative to average sample size
                mean_cells = np.mean([adatas[sample_i].n_obs, adatas[sample_j].n_obs])
                jaccard_index = (intersect_count / mean_cells) * 100  # as percentage

                bc_qc_ji.loc[sample_i, sample_j] = jaccard_index
                bc_qc_ji.loc[sample_j, sample_i] = jaccard_index

    # Optional filtering: remove samples with zero overlap
    if n_sample > 1:
        non_zero_rows = (bc_qc_raw.sum(axis=1) > 0).values
        non_zero_cols = (bc_qc_raw.sum(axis=0) > 0).values

        bc_qc_raw = bc_qc_raw.loc[non_zero_rows, non_zero_cols]
        bc_qc_ji = bc_qc_ji.loc[non_zero_rows, non_zero_cols]
        
    # Round Jaccard index values for clean display
    bc_qc_ji = bc_qc_ji.round(2)

    # Convert to JSON-serializable list of records
    bc_qc_raw = bc_qc_raw.reset_index().to_dict(orient="records")
    bc_qc_ji = bc_qc_ji.reset_index().to_dict(orient="records")

    return bc_qc_raw, bc_qc_ji

# -------------------------------

def scrublet_data(scrublet_scores):
    """
    Prepares long-format scrublet doublet scores for both observed and simulated distributions.

    Args:
        scrublet_scores (dict): Dictionary where keys are sample IDs and values are dictionaries with:
            - "doublet_scores_obs": List of observed doublet scores.
            - "doublet_scores_sim": List of simulated doublet scores.

    Returns:
        Two lists of dictionaries:
            - scrublet_data_obs: Observed scores with keys {"sample_id", "score"}.
            - scrublet_data_sim: Simulated scores with keys {"sample_id", "score"}.
    """

    # Initialize containers for long-format score data
    scrublet_data_obs = []
    scrublet_data_sim = []

    # Iterate over all sample entries in the input dictionary
    for sample_id, scores in scrublet_scores.items():
        obs_scores = scores["doublet_scores_obs"] # Observed scores for the sample
        sim_scores = scores["doublet_scores_sim"] # Simulated scores for the sample

        # Format each observed score as a dict and append to the list
        scrublet_data_obs.extend([
            {"sample_id": sample_id, "score": float(score)} for score in obs_scores
        ])
        # Format each simulated score as a dict and append to the list
        scrublet_data_sim.extend([
            {"sample_id": sample_id, "score": float(score)} for score in sim_scores
        ])

    # Return both observed and simulated scores in long-format
    return scrublet_data_obs, scrublet_data_sim

# -------------------------------

def qc_impact_data(adatas):
    """
    Computes summary statistics on UMI and gene detection counts per cell for each sample.

    Args:
        adatas (dict): Dictionary where keys are sample IDs and values are AnnData objects.

    Returns:
        list: List of dictionaries (one per sample), each containing:
            - sample_id (str): Sample identifier.
            - sample_name (str): Sample name from adata.obs['sample'].
            - avg_count_g (float): Average total UMI counts per cell.
            - sd_count_g (float): Standard deviation of UMI counts per cell.
            - avg_count_c (float): Average number of genes detected per cell.
            - sd_count_c (float): Standard deviation of gene detection counts.
    """

    # Convert adatas to list for incremental access
    adatas_list = list(adatas.items())
    db_plot = []

    # Iterate over each sample and compute metrics
    for i in range(len(adatas_list)):  
        # Extract sample name from metadata
        sample_id, adata = adatas_list[i]
        sample_name = adata.obs['sample'].unique()[0]

        # UMI counts per cell (sum of expression values per row)
        gene_counts = np.array(adata.X.sum(axis=1)).flatten()

        # Number of genes detected per cell (non-zero gene counts)
        detected_genes_counts = np.array((adata.X > 0).sum(axis=1)).flatten()

        # Compute statistics
        avg_count_g = float(np.mean(gene_counts))  
        sd_count_g = float(np.std(gene_counts))   
        avg_count_c = float(np.mean(detected_genes_counts))
        sd_count_c = float(np.std(detected_genes_counts))  

        # Append structured result
        db_plot.append({
            "sample_id": sample_id,
            "sample_name": sample_name,
            "avg_count_g": avg_count_g,
            "sd_count_g": sd_count_g,
            "avg_count_c": avg_count_c,
            "sd_count_c": sd_count_c
        })

    return db_plot  

# -------------------------------   

def prepare_visNetwork(adata, ontology, from_cell, to_cell):
    """
    Prepares network data for visualizing cell state → cell type mappings.

    This function constructs a visNetwork-compatible structure showing connections
    between cell states (fine-grained annotations) and cell types (coarse annotations)
    based on an ontology.

    Args:
        adata (AnnData): Annotated data matrix with `obs` containing cell labels.
        ontology (dict): Dictionary mapping each cell state to its parent cell type.
        from_cell (str): Column in `adata.obs` containing cell states (source).
        to_cell (str): Column in `adata.obs` containing cell types (target).

    Returns:
        dict: A dictionary with two keys:
            - "nodes": list of node dictionaries (id, label, group, value, title)
            - "edges": list of edge dictionaries (from, to)
    """

    # Subset and drop missing entries
    df = adata.obs[[to_cell, from_cell]].dropna()

    # Total number of cells (used to calculate proportions)
    total_cells = adata.shape[0]

    # Combine both columns to count node frequencies (both cell types and states)
    node_counts = (
        pd.concat([
            df[to_cell],
            df[from_cell]
        ])
        .value_counts()
        .reset_index()
        .rename(columns={"index": "label", 0: "count"})
    )

    # Ensure columns are strings to avoid category mismatch
    df[to_cell] = df[to_cell].astype(str)
    df[from_cell] = df[from_cell].astype(str)

    # Identify self-loops (same label appears as both state and type)
    self_occurrences = df[df[to_cell] == df[from_cell]][to_cell].value_counts()

    # Adjust node counts by removing self-loop duplications
    node_counts["count"] -= node_counts["label"].map(self_occurrences).fillna(0).astype(int)

    # Construct visNetwork nodes
    nodes = []
    for _, row in node_counts.iterrows():
        label = row["label"]
        count = int(row["count"])
        proportion = (count / total_cells) * 100
        group = "celltype" if label in df[to_cell].values else "cellstate"

        nodes.append({
            "id": label,
            "label": label,
            # "title": f"{label}<br>Count: {count}<br>Proportion: {proportion:.2f}%",
            "title": f"{label}<br>Count: {count:,}<br>Proportion: {proportion:.2f}%",
            "value": count,
            "group": group
        })

    # Construct visNetwork edges using ontology (cellstate → celltype)
    edges = []
    for _, row in df.drop_duplicates().iterrows():
        cellstate = row[from_cell]
        celltype = ontology.get(cellstate)
        if celltype and celltype != cellstate:  # avoid self-loop edges
            edges.append({"from": cellstate, "to": celltype})

    return {"nodes": nodes, "edges": edges}

# -------------------------------   

def get_python_environment():
    """
    Collects information about the current Python runtime environment.

    Returns:
        dict: A dictionary containing:
            - 'python_version': Current Python version as a string (e.g., "3.9.16")
            - 'packages': List of installed packages with their names and versions,
                          each represented as a dictionary with 'Package' and 'Version' keys.
    """

    packages = [
        {"Package": d.project_name, "Version": d.version}
        for d in pkg_resources.working_set
    ]

    return {
        "python_version": platform.python_version(),
        "packages": packages
    }
    
# -------------------------------

def computed_metadata_description(adata):
    """
    Generate structured descriptions for metadata fields present in `adata.obs`.

    This function cross-references known computed fields against `adata.obs` and returns
    a list of field-description mappings. Unknown fields are labeled with a default note
    and reported via console warning.

    Args:
        adata (AnnData): Annotated AnnData object containing metadata in `adata.obs`.

    Returns:
        list: A list of dictionaries, each with:
              - "field": the name of the metadata column.
              - "description": human-readable explanation of the field.
    """

    # Predefined dictionary of known computed metadata and their descriptions
    known_descriptions = {        
        "sample": "Sample name provided by the user, matching the input file or experimental label.",
        "sample_id": "Unique identifier for each sample, used to track each sample.",
        "n_genes": "Number of detected genes per cell",
        "n_counts": "Total number of UMIs per cell",
        "n_genes_by_counts": "Number of detected genes per cell",
        "total_counts": "Total number of UMIs per cell",
        "total_counts_mito": "Total UMIs from mitochondrial genes per cell",
        "pct_counts_mito": "Percentage of mitochondrial UMIs per cell",
        "doublet_scores_obs": "Scrublet doublet score for each observed cell",
        "predicted_doublet": "Doublet prediction label (True/False)",
        "leiden_cluster": "Leiden cluster assignment",
        "cellstate_scimilarity": "Predicted cell state from scimilarity",
        "celltype_scimilarity": "Curated cell type using pipeline ontology",
        "cellstate_celltypist": "Predicted cell state from celltypist",
        "celltype_celltypist": "Curated cell type using pipeline ontology",
        "scevan_class": "Indicates the tumor classification assigned by SCEVAN",
        "min_cell": "Minimum number of cells in which a gene must be detected to be retained.",
        "min_umi_per_cell": "Minimum UMI counts required per cell.",
        "max_umi_per_cell": "Maximum UMI counts allowed per cell.",
        "min_genes_per_cell": "Minimum number of detected genes per cell.",
        "max_genes_per_cell": "Maximum number of detected genes per cell.",
        "max_mt_percent": "Maximum allowed percentage of mitochondrial content per cell."
    }

    # Iterate over metadata columns and annotate with known descriptions
    mtdta = []
    for column in adata.obs.columns:
        description = known_descriptions.get(column, "No description available.")
        if column not in known_descriptions:
            print(f"*** ⚠️  Skipping unknown computed metadata field: '{column}' (no impact on the analysis).")
        mtdta.append({"field": column, "description": description})

    return mtdta

# -------------------------------

# Granular -> Higher-level cell type mapping
ontology_map = {
    # B cells
    "B cell": "B CELL",
    "B cells": "B CELL",
    "b_naive": "B CELL",
    "bmem_switched": "B CELL",
    "bmem_unswitched": "B CELL",
    "Age-associated B cells": "B CELL",
    "class switched memory B cell": "B CELL",
    "follicular B cell": "B CELL",
    "Follicular B cells": "B CELL",
    "germinal center B cell": "B CELL",
    "Germinal center B cells": "B CELL",
    "immature B cell": "B CELL",
    "memory B cell": "B CELL",
    "Memory B cells": "B CELL",
    "naive B cell": "B CELL",
    "Naive B cells": "B CELL",
    "Cycling B cells": "B CELL",
    "plasmablast": "B CELL",
    "Plasmablasts": "B CELL",
    "plasma cell": "B CELL",
    "Plasma cells": "B CELL",
    "IgA plasma cell": "B CELL",
    "IgG plasma cell": "B CELL",
    "IgM plasma cell": "B CELL",
    "plasma_IgG": "B CELL",
    "plasma_IgA": "B CELL",
    "precursor B cell": "B CELL",
    "pro-B cell": "B CELL",
    "Pro-B cells": "B CELL",
    "Pre-pro-B cells": "B CELL",
    "Large pre-B cells": "B CELL",
    "Small pre-B cells": "B CELL",
    "Proliferative germinal center B cells": "B CELL",
    "Transitional B cells": "B CELL",

    # T cells
    "T cell": "T CELL",
    "T cells": "T CELL",
    "CD4-Treg": "T CELL",
    "CD4-Tem": "T CELL",
    "CD4-positive helper T cell": "T CELL",
    "CD4-positive, alpha-beta T cell": "T CELL",
    "CD4-activated": "T CELL",
    "CD8-positive, alpha-beta T cell": "T CELL",
    "CD8a/a": "T CELL",
    "CD8a/b(entry)": "T CELL",
    "CD8-activated": "T CELL",
    "CD4-positive, alpha-beta cytotoxic T cell": "T CELL",
    "CD8-positive, alpha-beta cytotoxic T cell": "T CELL",
    "CD4-positive, alpha-beta memory T cell": "T CELL",
    "CD8-positive, alpha-beta memory T cell": "T CELL",
    "Memory CD4+ cytotoxic T cells": "T CELL",
    "T-helper 1 cell": "T CELL",
    "Type 1 helper T cells": "T CELL",
    "T-helper 17 cell": "T CELL",
    "Type 17 helper T cells": "T CELL",
    "T follicular helper cell": "T CELL",
    "Follicular helper T cells": "T CELL",
    "gamma-delta T cell": "T CELL",
    "CRTAM+ gamma-delta T cells": "T CELL",
    "gamma-delta T cells": "T CELL",
    "mucosal invariant T cell": "T CELL",
    "MAIT cells": "T CELL",
    "naive T cell": "T CELL",
    "memory T cell": "T CELL",
    "naive thymus-derived CD4-positive, alpha-beta T cell": "T CELL",
    "naive thymus-derived CD8-positive, alpha-beta T cell": "T CELL",
    "regulatory T cell": "T CELL",
    "Regulatory T cells": "T CELL",
    "Treg(diff)": "T CELL",
    "CD4-positive, CD25-positive, alpha-beta regulatory T cell": "T CELL",
    "naive regulatory T cell": "T CELL",
    "mature NK T cell": "T CELL",
    "NKT cells": "T CELL",
    "alpha-beta T cell": "T CELL",
    "activated CD4-positive, alpha-beta T cell": "T CELL",
    "activated CD8-positive, alpha-beta T cell": "T CELL",
    "central memory CD4-positive, alpha-beta T cell": "T CELL",
    "central memory CD8-positive, alpha-beta T cell": "T CELL",
    "effector memory CD4-positive, alpha-beta T cell": "T CELL",
    "effector memory CD8-positive, alpha-beta T cell": "T CELL",
    "effector memory CD8-positive, alpha-beta T cell, terminally differentiated": "T CELL",
    "Tcm/Naive cytotoxic T cells": "T CELL",
    "Tcm/Naive helper T cells": "T CELL",
    "Tem/Effector helper T cells": "T CELL",
    "Tem/Effector helper T cells PD1+": "T CELL",
    "Tem/Temra cytotoxic T cells": "T CELL",
    "Tem/Trm cytotoxic T cells": "T CELL",
    "Trm cytotoxic T cells": "T CELL",
    "T(agonist)": "T CELL",
    "Early lymphoid/T lymphoid": "T CELL",
    "Cycling T cells": "T CELL",
    "Double-positive thymocytes": "T CELL",
    "double-positive, alpha-beta thymocyte": "T CELL",
    "double negative thymocyte": "T CELL",
    "Double-negative thymocytes": "T CELL",
    "thymocyte": "T CELL",
    "CD4 T cells": "T CELL",
    "effector CD8-positive, alpha-beta T cell": "T CELL",
    "effector CD4-positive, alpha-beta T cell": "T CELL",
    "CD8 T cells": "T CELL",
    "T cells proliferating": "T CELL",
    "T_prol": "T CELL", # Proliferating T cell
    "GD": "T CELL", # Gamma delta T cell

    # Monocytes & Macrophages
    "monocyte": "MYELOIDS",
    "Monocytes": "MYELOIDS",
    "classical monocyte": "MYELOIDS",
    "Classical monocytes": "MYELOIDS",
    "non-classical monocyte": "MYELOIDS",
    "Non-classical monocytes": "MYELOIDS",
    "intermediate monocyte": "MYELOIDS",
    "Intermediate macrophages": "MYELOIDS",
    "inflammatory macrophage": "MYELOIDS",
    "macrophage": "MYELOIDS",
    "Macrophages": "MYELOIDS",
    "alveolar macrophage": "MYELOIDS",
    "Alveolar macrophages": "MYELOIDS",
    "Kupffer cell": "MYELOIDS",
    "Kupffer cells": "MYELOIDS",
    "Kidney-resident macrophages": "MYELOIDS",
    "Erythrophagocytic macrophages": "MYELOIDS",
    "Intestinal macrophages": "MYELOIDS",   
    "CD14-low, CD16-positive monocyte": "MYELOIDS",
    "CD14-positive monocyte": "MYELOIDS",
    "CD14-positive, CD16-positive monocyte": "MYELOIDS",
    "myeloid cell": "MYELOIDS",
    "Interstitial Mph perivascular": "MYELOIDS",
    "Alveolar Mph CCL3+": "MYELOIDS",
    "Monocyte-derived Mph": "MYELOIDS",
    "Alveolar Mph proliferating": "MYELOIDS",
    "lung macrophage": "MYELOIDS",
    "Mono-mac": "MYELOIDS",
    "Monocyte precursor": "MYELOIDS",
    "Langerhans cells": "MYELOIDS",
    "Langerhans cell": "MYELOIDS",
    "Alveolar Mph MT-positive": "MYELOIDS",
    "Macro-lipo": "MYELOIDS",
    "Macro-m2-CXCL": "MYELOIDS",
    "Mono-non-classical": "MYELOIDS",
    "mye-prol": "MYELOIDS", # Proliferating myeloid progenitor
    "Macro-m1-CCL": "MYELOIDS",
    "Macro-m2": "MYELOIDS",
    "Macro-m1": "MYELOIDS",

    # Dendritic Cells
    "dendritic cell": "DENDRITIC CELL",
    "DC": "DENDRITIC CELL",
    "DC1": "DENDRITIC CELL",
    "DC2": "DENDRITIC CELL",
    "DC3": "DENDRITIC CELL",
    "cDC1": "DENDRITIC CELL",
    "cDC2": "DENDRITIC CELL",
    "Migratory DCs": "DENDRITIC CELL",
    "DC precursor": "DENDRITIC CELL",
    "Transitional DC": "DENDRITIC CELL",
    "pDC": "DENDRITIC CELL",
    "pDC precursor": "DENDRITIC CELL",
    "conventional dendritic cell": "DENDRITIC CELL",
    "plasmacytoid dendritic cell": "DENDRITIC CELL",
    "CD1c-positive myeloid dendritic cell": "DENDRITIC CELL",
    "CD141-positive myeloid dendritic cell": "DENDRITIC CELL",
    "plasmacytoid dendritic cell, human": "DENDRITIC CELL",
    "myeloid dendritic cell": "DENDRITIC CELL",
    "myeloid dendritic cell, human": "DENDRITIC CELL",
    "dendritic cell, human": "DENDRITIC CELL",
    "follicular dendritic cell": "DENDRITIC CELL",
    "Plasmacytoid DCs": "DENDRITIC CELL",

    # NK cells
    "natural killer cell": "NK CELL",
    "NK cells": "NK CELL",
    "NK": "NK CELL",
    "CD16+ NK cells": "NK CELL",
    "CD16- NK cells": "NK CELL",
    "Transitional NK": "NK CELL",
    "Cycling NK cells": "NK CELL",
    "CD16-negative, CD56-bright natural killer cell, human": "NK CELL",
    "CD16-positive, CD56-dim natural killer cell, human": "NK CELL",

    # Granulocyte Cells
    "neutrophil": "GRANULOCYTE",
    "Neutrophils": "GRANULOCYTE",
    "Neutrophil": "GRANULOCYTE",
    "Neutrophil-myeloid progenitor": "GRANULOCYTE",
    "Granulocytes": "GRANULOCYTE",
    "mast cell": "GRANULOCYTE",
    "Mast": "GRANULOCYTE",
    "Mast cells": "GRANULOCYTE",
    "Myelocytes": "GRANULOCYTE",
    "granulocyte": "GRANULOCYTE",

    # Erythroid Lineage
    "erythrocyte": "ERYTHROID",
    "ERYTHROID": "ERYTHROID",
    "Erythrocytes": "ERYTHROID",
    "Early erythroid": "ERYTHROID",
    "Late erythroid": "ERYTHROID",
    "Mid erythroid": "ERYTHROID",
    "erythroid progenitor cell": "ERYTHROID",
    "erythroid lineage cell": "ERYTHROID",
    "Erythroid": "ERYTHROID",

    # Platelets
    "platelet": "PLATELET",
    "Platelet": "PLATELET",
    "megakaryocyte": "PLATELET",
    "Megakaryocytes/platelets": "PLATELET",
    "Megakaryocyte precursor": "PLATELET",
    "Megakaryocyte-erythroid-mast cell progenitor": "PLATELET",
    "Early MK": "PLATELET",
    "MEMP": "PLATELET", # Megakaryocyte-Erythroid-Mast cell
    
    # Epithelial Cells
    "epithelial cell": "EPITHELIAL",
    "ciliated cell": "EPITHELIAL",
    "goblet cell": "EPITHELIAL",
    "club cell": "EPITHELIAL",
    "Club (nasal)": "EPITHELIAL",
    "cholangiocyte": "EPITHELIAL",
    "kidney loop of Henle thin descending limb epithelial cell": "EPITHELIAL",
    "luminal epithelial cell of mammary gland": "EPITHELIAL",
    "type I pneumocyte": "EPITHELIAL",
    "type II pneumocyte": "EPITHELIAL",
    "respiratory epithelial cell": "EPITHELIAL",
    "basal cell": "EPITHELIAL",
    "basal": "EPITHELIAL",
    "respiratory basal cell": "EPITHELIAL",
    "lung ciliated cell": "EPITHELIAL",
    "mucus secreting cell": "EPITHELIAL",
    "secretory cell": "EPITHELIAL",
    "pulmonary ionocyte": "EPITHELIAL",
    "tracheal goblet cell": "EPITHELIAL",
    "basal cell of prostate epithelium": "EPITHELIAL",
    "lung secretory cell": "EPITHELIAL",
    "enterocyte": "EPITHELIAL",
    "intestinal tuft cell": "EPITHELIAL",
    "kidney collecting duct principal cell": "EPITHELIAL",
    "SMG mucous": "EPITHELIAL", # Submucosal Gland (SMG) Mucous Cells
    "SMG duct": "EPITHELIAL", # Submucosal Gland Duct
    "SMG serous (nasal)": "EPITHELIAL",
    "AT0": "EPITHELIAL", # (Alveolar Type 0 Epithelial Cell)
    "AT1": "EPITHELIAL", # (Alveolar Type 1 Epithelial Cell)
    "AT2": "EPITHELIAL", # (Alveolar Type 2 Epithelial Cell)
    "AT2 proliferating": "EPITHELIAL",
    "AT3": "EPITHELIAL", # (Alveolar Type 3 Epithelial Cell)
    "Multiciliated (nasal)": "EPITHELIAL",
    "Neuroendocrine": "EPITHELIAL",
    "neuroendocrine cell": "EPITHELIAL",
    "Ionocyte": "EPITHELIAL",
    "ionocyte": "EPITHELIAL",
    "Tuft": "EPITHELIAL",
    "Goblet (subsegmental)": "EPITHELIAL",
    "Goblet (bronchial)": "EPITHELIAL",
    "Goblet (nasal)": "EPITHELIAL",
    "Hillock-like": "EPITHELIAL",
    "Deuterosomal": "EPITHELIAL",
    "Club (non-nasal)": "EPITHELIAL",
    "Basal resting": "EPITHELIAL",
    "Suprabasal": "EPITHELIAL",
    "pre-TB secretory": "EPITHELIAL", # terminal bronchiolar (TB)
    "Multiciliated (non-nasal)": "EPITHELIAL",
    "squamous epithelial cell": "EPITHELIAL",
    "Epithelial cells": "EPITHELIAL",
    "epithelial cell of proximal tubule": "EPITHELIAL",
    "kidney loop of Henle thick ascending limb epithelial cell": "EPITHELIAL",
    "melanocyte": "EPITHELIAL",
    "lung neuroendocrine cell": "EPITHELIAL",
    "parietal epithelial cell": "EPITHELIAL",
    "kidney connecting tubule epithelial cell": "EPITHELIAL",
    "mesothelial cell": "EPITHELIAL",
    "duct epithelial cell": "EPITHELIAL",
    "pancreatic D cell": "EPITHELIAL",
    "kidney epithelial cell": "EPITHELIAL",
    "enteroendocrine cell": "EPITHELIAL",
    "luminal cell of prostate epithelium": "EPITHELIAL",
    "kidney proximal convoluted tubule epithelial cell": "EPITHELIAL",
    "type B pancreatic cell": "EPITHELIAL",
    "kidney collecting duct intercalated cell": "EPITHELIAL",
    "medullary thymic epithelial cell": "EPITHELIAL",
    "kidney distal convoluted tubule epithelial cell": "EPITHELIAL",
    "keratinocyte": "EPITHELIAL",
    "glandular epithelial cell": "EPITHELIAL",
    "Mesothelium": "EPITHELIAL",
    "SMG serous (bronchial)": "EPITHELIAL",
    "Lumsec-HLA": "EPITHELIAL",
    "Lumsec-basal": "EPITHELIAL",
    "Lumsec-prol": "EPITHELIAL", # Proliferating luminal secretory
    "LummHR-major": "EPITHELIAL", # Major luminal hormone-responsive cell
    "LummHR-SCGB": "EPITHELIAL", # Luminal secretory epithelial
    "LummHR-active": "EPITHELIAL",
    "Lumsec-KIT": "EPITHELIAL", # KIT+ luminal secretory epithelial cell
    "colon epithelial cell": "EPITHELIAL",

    # Fibroblasts
    "Fibroblasts": "FIBROBLAST",
    "fibroblast": "FIBROBLAST",
    "Fibro-matrix": "FIBROBLAST",
    "myofibroblast cell": "FIBROBLAST",
    "fibroblast of cardiac tissue": "FIBROBLAST",
    "fibroblast of lung": "FIBROBLAST",
    "kidney interstitial fibroblast": "FIBROBLAST",
    "Myofibroblasts": "FIBROBLAST",
    "Alveolar fibroblasts": "FIBROBLAST",
    "Subpleural fibroblasts": "FIBROBLAST",
    "Adventitial fibroblasts": "FIBROBLAST",
    "Peribronchial fibroblasts": "FIBROBLAST",
    "Fibro-SFRP4": "FIBROBLAST",

    # Endothelial Cells
    "Endothelial cells": "ENDOTHELIAL",
    "endothelial cell": "ENDOTHELIAL",
    "blood vessel endothelial cell": "ENDOTHELIAL",
    "capillary endothelial cell": "ENDOTHELIAL",
    "lung endothelial cell": "ENDOTHELIAL",
    "endothelial cell of artery": "ENDOTHELIAL",
    "endothelial cell of vascular tree": "ENDOTHELIAL",
    "EC venous pulmonary": "ENDOTHELIAL",
    "vein endothelial cell": "ENDOTHELIAL",
    "endothelial cell of lymphatic vessel": "ENDOTHELIAL",
    "endothelial cell of hepatic sinusoid": "ENDOTHELIAL",
    "cardiac endothelial cell": "ENDOTHELIAL",
    "Lymphatic EC proliferating": "ENDOTHELIAL",
    "EC aerocyte capillary": "ENDOTHELIAL",
    "EC venous systemic": "ENDOTHELIAL",
    "Lymphatic EC differentiating": "ENDOTHELIAL",
    "EC general capillary": "ENDOTHELIAL",
    "Lymphatic EC mature": "ENDOTHELIAL",
    "EC arterial": "ENDOTHELIAL",
    "Lymph-valve1": "ENDOTHELIAL",
    "Lymph-valve2": "ENDOTHELIAL",
    "Vas-arterial": "ENDOTHELIAL",
    "Vas-capillary": "ENDOTHELIAL",
    "Vas-venous": "ENDOTHELIAL",

    # Glial Cells
    "astrocyte": "GLIAL",
    "microglial cell": "GLIAL",
    "oligodendrocyte": "GLIAL",
    "radial glial cell": "GLIAL",
    "glial cell": "GLIAL",

    # Neurons
    "neuron": "NEURON",
    "cardiac neuron": "NEURON",
    "retinal cone cell": "NEURON",
    "retinal ganglion cell": "NEURON",
    "glutamatergic neuron": "NEURON",
    "neural cell": "NEURON",
    "ON-bipolar cell": "NEURON",
    "amacrine cell": "NEURON",
    "retinal rod cell": "NEURON",
    "retinal bipolar neuron": "NEURON",

    # Stem cell
    "hematopoietic stem cell": "STEM CELL",
    "mesenchymal stem cell": "STEM CELL",
    "progenitor cell": "STEM CELL",
    "common lymphoid progenitor": "STEM CELL",
    "hematopoietic precursor cell": "STEM CELL",
    "stem cell": "STEM CELL",
    "Hematopoietic stem cells": "STEM CELL",
    "HSC/MPP": "STEM CELL",

    # Stromal
    "stromal cell": "STROMAL",
    "mesenchymal cell": "STROMAL",
    "Pericytes": "STROMAL",
    "pericytes": "STROMAL",
    "Smooth muscle": "STROMAL",
    "enteric smooth muscle cell": "STROMAL",
    "pericyte": "STROMAL",
    "fat cell": "STROMAL",
    "stromal cell of ovary": "STROMAL",
    "epicardial adipocyte": "STROMAL",

    # Hepatocyte
    "hepatocyte": "HEPATOCYTE",
    "periportal region hepatocyte": "HEPATOCYTE",

    # Pancreatic
    "pancreatic acinar cell": "PANCREATIC CELL",
    "pancreatic ductal cell": "PANCREATIC CELL",
    "pancreatic A cell": "PANCREATIC CELL",
    "pancreatic stellate cell": "PANCREATIC CELL",

    # Miscellaneous
    "lymphocyte": "LYMPHOCYTE",
    "Lymph-major": "LYMPHOCYTE",
    "leukocyte": "LEUKOCYTE",
    "native cell": "NATIVE CELL",
    "Smooth Muscle FAM83D+": "MESENCHYMAL",
    "Smooth muscle FAM83D+": "MESENCHYMAL",
    "animal cell": "Other",
    "Cycling cells": "Other",

    # Muscle cell
    "cardiac muscle cell": "MUSCLE CELL",
    "fast muscle cell": "MUSCLE CELL",
    "smooth muscle cell of prostate": "MUSCLE CELL",
    "bronchial smooth muscle cell": "MUSCLE CELL",
    "vascular associated smooth muscle cell": "MUSCLE CELL",
    "smooth muscle cell": "MUSCLE CELL",
    "SM activated stress response": "MUSCLE CELL",
    "vsmc": "MUSCLE CELL", # Vascular smooth muscle cell

    # ILCs
    "ILC": "INNATE LYMPHOID CELL",
    "ILC1": "INNATE LYMPHOID CELL",
    "ILC2": "INNATE LYMPHOID CELL",
    "ILC3": "INNATE LYMPHOID CELL",
    "ILC precursor": "INNATE LYMPHOID CELL",
    "innate lymphoid cell": "INNATE LYMPHOID CELL", 
    "group 3 innate lymphoid cell": "INNATE LYMPHOID CELL"
}

# Any unmapped cell types will retain their original labels.

# -------------------------------

pipeline_arguments = {
    "species": "Species type used in the analysis.",
    "tissue": "Tissue type for the dataset.",
    "disease": "Disease name",
    "only_qc": "If enabled, the pipeline only runs quality control steps without further analysis.",
    "input": "Path to the input directory containing raw data for processing.",
    "project": "Project name identifier for the analysis run.",
    "pipe_version": "Version number of the pipeline used for execution.",
    "package_path": "Directory path where the pipeline package and dependencies are stored.",
    "date": "Date when the pipeline was executed.",
    "outputs_path": "Directory path where all pipeline-generated outputs will be stored.",
    "ui": "Unique identifier for the pipeline execution instance.",
    
    # QC parameters
    "min_umi_per_cell": "Minimum UMI counts required per cell.",
    "max_umi_per_cell": "Maximum UMI counts allowed per cell.",
    "min_genes_per_cell": "Minimum number of detected genes per cell.",
    "max_genes_per_cell": "Maximum number of detected genes per cell.",
    "max_mt_percent": "Maximum allowed percentage of mitochondrial content per cell.",
    "doublet_method": "Method used for doublet identification.",
    "scrublet_cutoff": "Threshold for filtering doublets based on Scrublet score.",
    "min_cell": "Minimum number of cells in which a gene must be detected to be retained.",

    # Processing settings
    "norm_target_sum": "Target sum for normalization, defining the total counts per cell.",
    "n_top_genes": "Number of highly variable genes to retain for downstream analysis.",
    "regress_out": "Whether total counts and mitochondrial percentages are regressed out.",
    "scale_max_value": "Maximum value for data scaling.",
    "n_pcs": "Number of principal components used in PCA.",
    "batch_correction": "Batch correction method applied.",
    "batch_vars": "Metadata column(s) used for batch correction.",
    "n_neighbors": "Number of neighbors used for kNN graph construction.",
    "resolution": "Resolution parameter for Leiden clustering.",
    "compute_tsne": "Indicates whether t-SNE embedding is computed.",
    
    # Annotation settings
    "annotation_method": "Method used for cell annotation.",
    "sci_model_path": "Path to the SCimilarity model for annotation.",
    "cty_model_path": "Path to the Celltypist model for annotation.",
    "cty_model_name": "Celltypist model used.",

    # Differentially expressed analysis arguments
    "pval_threshold": "P-value threshold for filtering differentially expressed genes (DEGs). ",
    "logfc_threshold": "Log fold-change (logFC) threshold for filtering DEGs.",
    "dea_method": "Statistical test used for differential expression analysis.",
    "top_n_deg_leidn": "Number of top differentially expressed genes (DEGs) to return per 'leiden' clusters.",
    "top_n_deg_scim": "Number of top differentially expressed genes (DEGs) to return per 'scimilarity' annotated cells.",
    "top_n_deg_cltpst": "Number of top differentially expressed genes (DEGs) to return per 'celltypist' annotated cells.",
    "pts_threshold": "Minimum fraction of cells expressing a gene for it to be considered a DEG.",

    # tumor id
    "tumor_id": "Tumor cell identification method.",

    # extra
    "runtime_minute": "pipeline run time in minute.",
    "doc_url": "URL to the published documentation associated with the analysis.",
    "data_url": "URL to the publicly available dataset or repository.",
    "limit_threads": "Apply thread limits to avoid memory crashes",
    "fix_gene_names": "Fix gene names if Ensembl IDs are detected.",
}

# -------------------------------

def extract_pipeline_arguments(args, argument_dict):
    """
    Extract relevant pipeline arguments and their metadata for reporting purposes.

    This function filters the command-line arguments used in the pipeline,
    matches them with human-readable descriptions, and structures them into
    a list suitable for logging, reporting, or summary tables.

    Args:
        args (Namespace): Parsed arguments from argparse.
        argument_dict (dict): Dictionary mapping argument names to their descriptions.

    Returns:
        list: List of dictionaries, each containing:
              - "Argument": argument name.
              - "Description": brief explanation of the argument.
              - "Value": stringified user-specified value.
    """

    # Arguments to be excluded from the report (e.g., internal references, non-user facing)
    SKIP_ARGS = {
        "project", "metadata", "genesets", "docker", "input", "help", "ui", "pipe_version", "date", "outputs_path",  
        "only_qc", "paper_url", "data_url", "package_path", "runtime_minute", "pre_qc_cells", "post_qc_cells", 
        "num_samples", "database_path"
    }

    used_args = []
    missing_descriptions = []

    for arg_name in vars(args):
        if arg_name in SKIP_ARGS:
            continue  # Skip specified arguments

        # Add argument with its description and actual value
        if arg_name in argument_dict:            
            used_args.append({
                "Argument": arg_name,
                "Description": argument_dict[arg_name],
                "Value": str(getattr(args, arg_name))
            })
        else:
            # Flag argument as undocumented for developer follow-up
            missing_descriptions.append(arg_name)

    # Issue a warning if any argument is missing from pipeline_arguments
    if missing_descriptions:
        print(f"*** ⚠️  The following arguments have no description: {', '.join(missing_descriptions)}, (no impact on the analysis).")

    return used_args

# -------------------------------

def prepare_qc_density_data(adata, include_tsne):
    """
    Prepare metadata and embedding data for QC-based density plots.

    This function extracts cell-level metadata and dimensionality reduction embeddings
    (UMAP and optionally t-SNE) to create a JSON-serializable format suitable for
    downstream visualization (e.g., density plots).

    Args:
        adata (AnnData): Annotated single-cell object containing metadata and embeddings.
        include_tsne (bool): If True, includes t-SNE coordinates if available in adata.obsm.

    Returns:
        list: List of dictionaries (one per cell) containing coordinates and metadata.
    """

    # Extract primary features and UMAP embedding
    df = pd.DataFrame({
        "sample_id": adata.obs["sample_id"].tolist(),
        "nCount_RNA": adata.obs["total_counts"].tolist(),
        "nFeature_RNA": adata.obs["n_genes_by_counts"].tolist(),
        "percent_mt": adata.obs["pct_counts_mito"].tolist(),
        "UMAP1": adata.obsm["X_umap"][:, 0].tolist(),
        "UMAP2": adata.obsm["X_umap"][:, 1].tolist()
    })

    # Optionally include t-SNE coordinates
    if include_tsne:
        df["TSNE1"] = adata.obsm["X_tsne"][:, 0].tolist()
        df["TSNE2"] = adata.obsm["X_tsne"][:, 1].tolist()

    # Return as JSON-compatible list of records (rows as dictionaries)
    return df.to_dict(orient="records")

# -------------------------------

def reorder_clusters(adata):
    """
    Reorders Leiden clusters in descending order of size.

    This function reassigns Leiden cluster labels such that:
    - Cluster 0 corresponds to the largest cluster,
    - Cluster 1 to the second largest,
    - and so on.

    The updated labels are written back to `adata.obs['leiden']`.

    Args:
        adata (AnnData): Annotated data matrix with a 'leiden' column in `adata.obs`.

    Returns:
        None: The function modifies `adata.obs['leiden']` in place.
    """

    # Count the number of cells in each Leiden cluster and sort in descending order
    cluster_counts = adata.obs["leiden"].value_counts().sort_values(ascending=False)

    # Generate a mapping from the original cluster labels to new labels (ordered by size)
    cluster_mapping = {old: str(new) for new, old in enumerate(cluster_counts.index)}

    # Apply the new labels to reorder clusters
    adata.obs["leiden"] = adata.obs["leiden"].map(cluster_mapping)

# -------------------------------

def degs_json(adata, groupby):
    """
    Extracts differentially expressed genes (DEGs) from `adata.uns` and formats them for JSON output.

    This function restructures the DEG table such that each group (cluster, condition, etc.)
    becomes a key and the associated DEGs are stored as a list of gene names.

    Args:
        adata (AnnData): Annotated data object containing DEG results in `adata.uns`.
        groupby (str): The column used for grouping cells during differential expression analysis.

    Returns:
        dict: Dictionary where each key is a group name and the value is a list of DEGs (gene names).
              This structure is JSON-serializable.
    """

    # Construct the expected key in adata.uns
    key_name = f"degs_filtered_{groupby}"

    # Retrieve the filtered DEG DataFrame
    degs_filtered = adata.uns[key_name]

    # Select relevant columns
    degs_table = degs_filtered[["gene", "group"]]

    # Restructure into a dictionary: group → [list of DEGs]
    degs_json = degs_table.groupby("group")["gene"].apply(list).to_dict()

    return degs_json

# -------------------------------

def check_gene_names_format(adata, preview_n=5):
    """
    Check if the `adata.var_names` primarily contain Ensembl gene identifiers and suggest potential replacements.

    This function inspects the gene names and identifies whether a majority follow the Ensembl ID format.
    If so, it searches for alternative columns in `adata.var` that might contain human-readable gene symbols.

    Args:
        adata (AnnData): Annotated single-cell object with gene metadata in `var`.
        threshold (float): Proportion of Ensembl-like names required to trigger a warning (default is 0.5).
        preview_n (int): Number of var_names to preview (default: 10).

    Returns:
        None. Prints diagnostic information and guidance to the user.
    """

    var_names = adata.var_names.astype(str)
 
    # Define Ensembl gene ID pattern (e.g., ENSG00000123456)
    ensembl_pattern = re.compile(r"^ENSG\d{11}$")
    ensembl_count = sum(bool(ensembl_pattern.match(name)) for name in var_names)
 
    total_genes = len(var_names)
    ensembl_ratio = ensembl_count / total_genes
 
    print(f"*** 🧬 Detected {ensembl_count:,} Ensembl-style gene names out of {total_genes:,} total genes "
          f"({ensembl_ratio:.2%}).")

    print("*** ⚠️ Preview of `adata.var_names`:")
    for name in var_names[:preview_n]:
        print(f"    - {name}")
 
    # Search for alternative gene symbol columns
    candidate_cols = [col for col in adata.var.columns if "symbol" in col.lower() or "gene" in col.lower()]
 
    if candidate_cols:
        print("*** 🔍 Candidate columns in `adata.var` that may contain gene symbols:")            
        for col in candidate_cols:
            print(f"***  > '{col}'")            
        print(f"{list(adata.var.columns)}")
    else:
        print("*** ⚠️  No obvious gene symbol columns found.")
        print("*** 📋 Available columns in `adata.var`:")
        print(f"   {list(adata.var.columns)}")
 
    print("*** 💡 To change this automatically, re-run and activate: --fix_gene_names following the candidate gene name column.")
 
# -------------------------------
 
def fix_gene_names(adata, column_name):
    """
    Replaces `adata.var_names` with values from a specified column in `adata.var` and ensures uniqueness.

    This function is used when gene symbols (e.g., HGNC names) are stored in a column other than `adata.var_names`,
    often after loading a dataset with Ensembl IDs or unnamed features.

    Args:
        adata (AnnData): Annotated single-cell expression object.
        column_name (str): Name of the column in `adata.var` to use as the new gene names.

    Returns:
        AnnData: The same object with updated and uniquely renamed `var_names`.
    """

    # Check if the specified column exists
    if column_name not in adata.var.columns:
        raise ValueError(
            f"*** 🚨 Column '{column_name}' not found in `adata.var.columns`.\n"
            f"*** 🚨 Available columns: {list(adata.var.columns)}"
        )
 
    # Replace gene names with the specified column
    adata.var_names = adata.var[column_name].astype(str) # Convert to string Index (avoid pandas Categorical issues)

    # Make sure names are unique (required for AnnData)
    adata.var_names_make_unique()

    print(f"*** ✅ Gene names replaced using `adata.var['{column_name}']`.")
 
    return adata
 
# -------------------------------

def qc_verbose_ls(metadata_df, adatas):
    """
    Prints a terminal-friendly summary of key QC metrics (UMIs, genes, mitochondrial %) for each sample.

    Args:
        metadata_df (pd.DataFrame): Metadata containing at least 'sample' and 'sample_id' columns.
        adatas (dict): Dictionary of AnnData objects keyed by sample_id.

    Returns:
        None: Prints formatted QC summaries to the terminal.
    """

    # Compute QC statistics using helper function
    umi_ls, gene_ls, mt_ls = compute_qc_stats_objs(metadata_df, adatas)

    # Helper to convert a list of dicts into a DataFrame and prefix metric columns
    def _to_df(lst, prefix):
        return pd.DataFrame(lst).rename(columns=lambda c: f"{prefix}_{c}" if c not in ["sample", "sample_id"] else c)

    # Convert each metric set to a DataFrame with prefixed column names
    df_umi = _to_df(umi_ls, "umi")
    df_gene = _to_df(gene_ls, "gene")
    df_mito = _to_df(mt_ls, "mito")

    # Merge all metrics on 'sample' and 'sample_id'
    df_qc = df_umi.merge(df_gene, on=["sample", "sample_id"]).merge(df_mito, on=["sample", "sample_id"])

    # Round numeric QC metrics for clean output
    metric_cols = [col for col in df_qc.columns if any(col.startswith(p) for p in ["umi_", "gene_", "mito_"])]
    df_qc[metric_cols] = df_qc[metric_cols].round(2)

    # Print formatted QC metrics per sample
    print("*** 🧬 QC Summary:")
    for _, row in df_qc.iterrows():
        print(f"*** 🧬 {row['sample']} -> {row['sample_id']}")
        print(f"*** ** UMI counts: Min: {row['umi_min']}, Q25: {row['umi_q25']}, Q50: {row['umi_q50']}, Q75: {row['umi_q75']}, Q100: {row['umi_q100']}")
        print(f"*** ** Gene counts: Min: {row['gene_min']}, Q25: {row['gene_q25']}, Q50: {row['gene_q50']}, Q75: {row['gene_q75']}, Q100: {row['gene_q100']}")
        print(f"*** ** Mitochondrial: Min: {row['mito_min']}%, Q25: {row['mito_q25']}%, Q50: {row['mito_q50']}%, Q75: {row['mito_q75']}%, Q100: {row['mito_q100']}%")

# -------------------------------

def qc_verbose(adata):
    """
    Prints a terminal-friendly QC summary for a single AnnData object,
    including statistics on UMIs per cell, genes per cell, and mitochondrial content.

    Args:
        adata (AnnData): Annotated single-cell object after QC.

    Returns:
        None: Prints QC summary to terminal.
    """

    # Compute QC summaries from helper
    umi_ls, gene_ls, mt_ls = compute_qc_stats_obj(adata)

    # Helper to convert list of {metric, value} to a flat dictionary with prefixed keys
    def to_dict(lst, prefix):
        result = {}
        for d in lst:
            key = f"{prefix}_{d['metric']}"
            val = d["value"]
            result[key] = round(val, 2) if isinstance(val, (int, float)) else val
        return result

    # Format all metrics into prefixed flat dictionaries
    umi_stats = to_dict(umi_ls, "umi")
    gene_stats = to_dict(gene_ls, "gene")
    mito_stats = to_dict(mt_ls, "mito")

    # Merge all stats into one dict for printing
    stats = {**umi_stats, **gene_stats, **mito_stats}

    # Print formatted QC metrics
    print("*** 🧬 QC Summary:")
    print(f"*** ** UMI counts: Min: {stats['umi_min']}, Q25: {stats['umi_q25']}, Q50: {stats['umi_q50']}, Q75: {stats['umi_q75']}, Q100: {stats['umi_q100']}")
    print(f"*** ** Gene counts: Min: {stats['gene_min']}, Q25: {stats['gene_q25']}, Q50: {stats['gene_q50']}, Q75: {stats['gene_q75']}, Q100: {stats['gene_q100']}")
    print(f"*** ** Mitochondrial: Min: {stats['mito_min']}%, Q25: {stats['mito_q25']}%, Q50: {stats['mito_q50']}%, Q75: {stats['mito_q75']}%, Q100: {stats['mito_q100']}%")

# -------------------------------

def summarize_adata_structure(adata):
    print("*** 🔍 AnnData summary:")
    print(f"• Shape (cells x genes): {adata.shape}")
    
    # .obs
    print("*** 🔍 adata.obs (cell metadata):")
    print(f"• Columns: {list(adata.obs.columns)}")
    
    # .var
    print("*** 🔍 adata.var (gene metadata):")
    print(f"• Columns: {list(adata.var.columns)}")

    # .uns
    print("*** 🔍 adata.uns (unstructured info):")
    for k in adata.uns.keys():
        v = adata.uns[k]
        desc = type(v).__name__
        if isinstance(v, dict):
            desc += f" (keys: {list(v.keys())[:5]})"
        print(f"  - {k}: {desc}")

    # .obsm
    print("*** 🔍 adata.obsm (multi-dimensional obs annotations):")
    for k in adata.obsm.keys():
        print(f"  - {k}: shape {adata.obsm[k].shape}")

    # .varm
    print("*** 🔍 adata.varm (multi-dimensional var annotations):")
    for k in adata.varm.keys():
        print(f"  - {k}: shape {adata.varm[k].shape}")

    # .layers
    print("*** 🔍 adata.layers (alternative expression matrices):")
    for k in adata.layers.keys():
        print(f"  - {k}: shape {adata.layers[k].shape}")

    # .raw
    if adata.raw is not None:
        print("*** 🔍 adata.raw available:")
        print(f"  - raw.shape: {adata.raw.shape}")
        print(f"  - raw.var_names: {list(adata.raw.var_names[:5])} ...")
    else:
        print("*** 🔍 adata.raw: None")

    return {
        "n_obs": adata.n_obs,
        "n_vars": adata.n_vars,
        "layers": list(adata.layers.keys()) if adata.layers else [],
        "obs_keys": list(adata.obs_keys()),
        "var_keys": list(adata.var_keys()),
        "uns_keys": list(adata.uns_keys()),
        "obsm_keys": list(adata.obsm_keys()),
        "varm_keys": list(adata.varm_keys()),
        "X_type": str(type(adata.X)),
        "X_shape": adata.X.shape,
    }

# -------------------------------

def read_study_metadata(path_to_json):
    """
    Reads required metadata fields from a CellExpress config JSON file.

    Args:
        path_to_json (str): Path to the JSON configuration file.

    Returns:
        dict: Dictionary containing 'outputs_path', 'ui', 'date', 'project', and 'species'.
    """
    if not os.path.exists(path_to_json):
        raise FileNotFoundError(f"*** 🚨 JSON file not found: {path_to_json}")

    with open(path_to_json, 'r') as f:
        config = json.load(f)

    try:
        return {
            "project": config["settings"]["project"],
            "outputs_path": config["settings"]["outputs_path"],
            "ui": config["execution_summary"]["ui"],
            "date": config["execution_summary"]["date_iso"],
            "species": config["execution_summary"]["species"]
        }
    except KeyError as e:
        raise ValueError(f"*** 🚨 Missing key in JSON: {e}")

# -------------------------------

def clean_value(val, key):
    # Treat all standard invalid entries + recognized nulls
    invalid_values = {"", " ", "NA", "na", "NaN", "nan", "NAN", "N/A", "n/a"}
    
    if pd.isna(val) or str(val).strip().lower() in {"none"}:
        return np.inf if key.startswith("max_") else None

    if isinstance(val, str):
        val = val.strip()
    
    val_str = str(val).strip().lower()

    # Reject ambiguous placeholders
    if val_str in {v.lower() for v in invalid_values}:
        raise ValueError(f"*** 🚨 Invalid QC value for '{key}': '{val}' is not allowed. Provide a positive number or 'None' for upper limits.")

    # Handle explicit Inf
    if val_str in {"inf", "infinity"}:
        return np.inf

    # Try parsing and enforce > 0
    try:
        numeric_val = float(val) if "percent" in key else int(val)
        if numeric_val <= 0:
            raise ValueError(f"*** 🚨 Invalid QC value for '{key}': must be > 0. Got: {numeric_val}")
        return numeric_val
    except Exception:
        raise ValueError(f"*** 🚨 Failed to parse QC value for '{key}': '{val}' is not a valid number.")


# -------------------------------
# Uniform per-sample QC dictionary builder
def build_qc_dicts(args, metadata_df, qc_keys):
    sample_ids = metadata_df["sample_id"].tolist()

    for key in qc_keys:
        if key in metadata_df.columns:
            qc_dict = {
                sample_id: clean_value(val, key)
                for sample_id, val in zip(metadata_df["sample_id"], metadata_df[key])
            }
        else:
            raw_default = getattr(args, key, None)
            if isinstance(raw_default, dict):
                qc_dict = raw_default  # already a valid dict
            else:
                cleaned_default = clean_value(raw_default, key)
                qc_dict = {sample_id: cleaned_default for sample_id in sample_ids}

        setattr(args, key, qc_dict)

# -------------------------------

def r_sanitize(value, key=None):
    if isinstance(value, dict):
        if key == "opt":
            # opt will be handled separately
            return value
        items = []
        for k, v in value.items():
            v_str = (
                "'Inf'" if isinstance(v, float) and np.isinf(v)
                else "NULL" if v is None
                else f"'{v}'" if isinstance(v, str)
                else str(v)
            )
            items.append(f"{k}={v_str}")
        return f"list({', '.join(items)})"
    elif isinstance(value, str):
        return f"'{value}'"
    elif value is None:
        return "NULL"
    elif isinstance(value, float) and np.isinf(value):
        return "'Inf'"
    else:
        return str(value)

# -------------------------------

def sanitize_inf(obj):
    if isinstance(obj, dict):
        return {k: sanitize_inf(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [sanitize_inf(x) for x in obj]
    elif isinstance(obj, float) and np.isinf(obj):
        return "Inf"  # JSON-safe string, will be handled in R
    else:
        return obj

# -------------------------------

def upgrade_features_to_v3(features_path: str) -> bool:

    # Read first few lines
    with gzip.open(features_path, "rt") as f:
        lines = [line.strip().split("\t") for line in f if line.strip()]
        col_count = len(lines[0]) if lines else 0

    if col_count == 3:
        return False  # Already v3-style
    elif col_count == 2:
        print(f"🧬 Detected 2-column v2-style features.tsv.gz — upgrading to v3-compatible format.")
        df = pd.read_csv(features_path, sep="\t", header=None)
        df["feature_type"] = "Gene Expression"
        df.to_csv(features_path, sep="\t", header=False, index=False, compression="gzip")
        return True
    elif col_count == 1:
        print(f"🧬 Detected 1-column legacy file — upgrading to full 3-column features.tsv.gz.")
        gene_ids = [row[0] for row in lines]
        df = pd.DataFrame({
            0: gene_ids,
            1: gene_ids,  # duplicate as gene_symbol
            2: "Gene Expression"
        })
        df.to_csv(features_path, sep="\t", header=False, index=False, compression="gzip")
        return True
    else:
        raise ValueError(f"❌ Unexpected column count in {features_path}: found {col_count} columns")

# -------------------------------

def try_load_sample_from_path(path: str, sample_id: str):
    """Attempt to load sample from path via supported formats."""
    # 10x single h5
    h5_file = glob.glob(os.path.join(path, "*.h5"))
    # processed h5ad
    h5ad_files = glob.glob(os.path.join(path, "*.h5ad"))
    # 10x triple
    matrix_file = os.path.join(path, "matrix.mtx.gz")
    features_file = os.path.join(path, "features.tsv.gz")
    barcodes_file = os.path.join(path, "barcodes.tsv.gz")
    # ParseBio triple 
    parsebio_genes = os.path.join(path, "all_genes.csv")
    parsebio_cells = os.path.join(path, "cell_metadata.csv")
    parsebio_mtx   = os.path.join(path, "count_matrix.mtx")

    # (1) HDF5
    if len(h5_file) == 1 and os.path.exists(h5_file[0]):
        print(f"*** 🔄 Loading sample from H5 file: {h5_file[0]} -> {sample_id}.")
        adata = sc.read_10x_h5(h5_file[0])
        return fix_duplicate_gene_names(adata)

    # (2) 10X-style Matrix
    if all(os.path.exists(f) for f in [matrix_file, features_file, barcodes_file]):
        print(f"*** 🔄 Loading sample from CSV files: {path} -> {sample_id}.")
        upgraded = upgrade_features_to_v3(features_file)
        try:
            # Try sparse first
            adata = sc.read_10x_mtx(path, var_names="gene_symbols", cache=False)
            return fix_duplicate_gene_names(adata)
        except Exception as e:
            print(f"⚠️ `matrix.mtx.gz` is not a Matrix Market file — falling back to dense TSV parsing.")
            # dense matrix load
            df = pd.read_csv(matrix_file, sep="\t", compression="gzip", index_col=0)
            adata = ad.AnnData(X=df.T)
            adata.var_names = df.index.astype(str)
            adata.obs_names = df.columns.astype(str)
            return fix_duplicate_gene_names(adata)

    # (3) ParseBio
    if all(os.path.exists(f) for f in [parsebio_genes, parsebio_cells, parsebio_mtx]):
            print(f"*** 🔄 Loading sample from ParseBio triple: {path} -> {sample_id}.")            
            # Read matrix, genes, and barcodes
            adata = sc.read_mtx(parsebio_mtx)
            gene_data = pd.read_csv(parsebio_genes)
            cell_meta = pd.read_csv(parsebio_cells)
            # Filter genes with NaN names and subset matrix columns
            if "gene_name" not in gene_data.columns:
                raise KeyError("ParseBio genes file must contain a 'gene_name' column.")
            gene_data = gene_data[gene_data.gene_name.notnull()].reset_index(drop=True)
            valid_gene_indices = gene_data.index.to_list()
            adata = adata[:, valid_gene_indices]
            # Assign var (genes)
            adata.var = gene_data
            adata.var.set_index('gene_name', inplace=True)
            adata.var.index.name = None
            adata.var_names_make_unique()
            # Assign obs (cells)
            if "bc_wells" not in cell_meta.columns:
                raise KeyError("ParseBio cell metadata must contain a 'bc_wells' column.")
            adata.obs = cell_meta
            adata.obs.set_index('bc_wells', inplace=True)
            adata.obs.index.name = None
            adata.obs_names_make_unique()
            return adata

    # (4) AnnData (.h5ad)
    if len(h5ad_files) == 1:
        print(f"*** 🔄 Loading sample from h5ad file: {h5ad_files[0]} -> {sample_id}.")
        adata = sc.read_h5ad(h5ad_files[0])
        metadata_cols = list(adata.obs.columns)
        if metadata_cols:
            print(f"*** 📋 Metadata fields available in `.obs`: {', '.join(metadata_cols)}")
        else:
            print("*** ⚠️  No metadata fields found in `.obs`.")
        return adata

    if len(h5ad_files) > 1:
        raise FileExistsError(f"🚨 Multiple `.h5ad` files found in {path}. Cannot determine which one to use.")

    # (5) Plain .txt.gz matrix 
    txtgz_files = glob.glob(os.path.join(path, "*.txt.gz"))
    if len(txtgz_files) == 1:
        print(f"*** 🔄 Loading matrix from plain txt.gz file: {txtgz_files[0]} -> {sample_id}.")
        df = pd.read_csv(txtgz_files[0], sep="\t", index_col=0)
        df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
        adata = ad.AnnData(X=df.T)
        adata.var_names = df.index.astype(str)
        adata.obs_names = df.columns.astype(str)
        return fix_duplicate_gene_names(adata)

    # (6) CSV-style matrix (.csv.gz)
    csvgz_files = glob.glob(os.path.join(path, "*.csv.gz"))
    if len(csvgz_files) == 1:
        print(f"*** 🔄 Loading matrix from csv.gz file: {csvgz_files[0]} -> {sample_id}.")
        df = pd.read_csv(csvgz_files[0], index_col=0)
        adata = ad.AnnData(X=df.T)
        adata.var_names = df.index.astype(str)
        adata.obs_names = df.columns.astype(str)
        return fix_duplicate_gene_names(adata)

    return None

# -------------------------------

import os, re, glob, tarfile, gzip, shutil

# Canonical file expectations
expected_patterns = {
    "matrix.mtx.gz": re.compile(r".*matrix.*\.mtx(\.gz)?$|.*matrix\.tsv$"),
    "features.tsv.gz": re.compile(r".*(features|genes).*\.tsv(\.gz)?$|.*genes\.tsv$"),
    "barcodes.tsv.gz": re.compile(r".*barcodes.*\.tsv(\.gz)?$|.*barcode\.tsv$"),
}

alt_to_canonical = {
    # 🔹 Legacy TSV file patterns
    "count_matrix_sparse.mtx": "matrix.mtx.gz",
    "count_matrix_genes.tsv": "features.tsv.gz",
    "count_matrix_barcodes.tsv": "barcodes.tsv.gz",

    # 🔹 GEO/GSM-style CSV variants
    "genes.csv.gz": "features.tsv.gz",
    "barcode.csv.gz": "barcodes.tsv.gz",
    "counts.mtx.gz": "matrix.mtx.gz",

    # 🔹 Uncompressed 10X-style .tsv variants
    "genes.tsv": "features.tsv.gz",
    "barcode.tsv": "barcodes.tsv.gz",
    "matrix.tsv": "matrix.mtx.gz",

    # 🔹 Common spelling alternates
    "barcodes.tsv": "barcodes.tsv.gz",
    "features.tsv": "features.tsv.gz",
    "matrix.mtx": "matrix.mtx.gz",
}

def extract_tar_gz(file_path: str, extract_dir: str) -> None:
    """Extracts a .tar.gz file to the given directory."""
    with tarfile.open(file_path, "r:gz") as tar:
        tar.extractall(path=extract_dir)

# -------------------------------

def fix_file_format(sample_path: str) -> str:
    """
    Ensures 10X-style format in sample_dir:
    - Extracts .tar.gz if present
    - Converts legacy names to canonical
    - Renames pattern-matched files to canonical names
    - Removes original mismatched files
    """

    report = []

    # 1️⃣ Extract TAR if present
    tar_files = glob.glob(os.path.join(sample_path, "*.tar.gz"))
    if len(tar_files) == 1:
        report.append(f"*** 📦 Extracting TAR file: {os.path.basename(tar_files[0])}")
        extract_tar_gz(tar_files[0], extract_dir=sample_path)
    elif len(tar_files) > 1:
        raise ValueError(f"🚨 Multiple .tar.gz files found in `{sample_path}`")

    # 2️⃣ Check for early exit: if any non-10X format exists, skip
    if glob.glob(os.path.join(sample_path, "*.h5")) \
    or glob.glob(os.path.join(sample_path, "*.h5ad")) \
    or glob.glob(os.path.join(sample_path, "*.txt.gz")) \
    or len(glob.glob(os.path.join(sample_path, "*.csv.gz"))) == 1:
        report.append(f"✅ Non-10X format detected — skipping format fix for `{sample_path}`")
        return "\n".join(report)

    # 3️⃣ Scan all subdirectories
    found_valid_dir = False
    for subdir, _, files in os.walk(sample_path):
        canonical_found = {k: None for k in expected_patterns}
        renamed = []

        # 2a — Legacy count_matrix_* renaming + gzipping
        for alt_name, canonical_name in alt_to_canonical.items():
            alt_path = os.path.join(subdir, alt_name)
            if os.path.exists(alt_path):
                gz_path = os.path.join(subdir, canonical_name)
                with open(alt_path, "rb") as f_in, gzip.open(gz_path, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
                os.remove(alt_path)
                renamed.append(f"*** - Gzipped & renamed `{alt_name}` → `{canonical_name}`")
                canonical_found[canonical_name] = canonical_name

        # 2b — Regex-match nonstandard canonical files and rename
        for f in files:
            full_path = os.path.join(subdir, f)
            for canonical, pattern in expected_patterns.items():
                if canonical_found[canonical] is None and pattern.match(f):
                    new_path = os.path.join(subdir, canonical)
                    if f != canonical:
                        shutil.copy2(full_path, new_path)
                        os.remove(full_path)
                        renamed.append(f"*** - Renamed & removed `{f}` → `{canonical}`")
                    canonical_found[canonical] = canonical

        # 2c — Final validation
        if all(canonical_found.values()):
            if renamed:
                report.append(f"📁 `{subdir}`")
                report.extend(renamed)
            found_valid_dir = True
            break

    # 3️⃣ Fail gracefully if no valid dir
    if not found_valid_dir:
        raise FileNotFoundError(
            f"🚨 Could not locate a directory under `{sample_path}` containing all required 10X files: "
            "`matrix.mtx.gz`, `features.tsv.gz` (or `genes.tsv.gz`), `barcodes.tsv.gz`."
        )

    report.insert(0, f"✅ File format validated and fixed under `{sample_path}`")
    return "\n".join(report)

# -------------------------------