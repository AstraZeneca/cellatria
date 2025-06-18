# report_wrapper.py
# -------------------------------

import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# -------------------------------

import json
import subprocess
import tempfile
import numpy as np
import pandas as pd
from helper import (compute_qc_stats_objs, generate_qc_plots_and_filters, compute_barcode_overlap_matrices, 
                    scrublet_data, computed_metadata_description, qc_impact_data, compute_qc_stats_obj,
                    prepare_visNetwork, ontology_map, get_python_environment, pipeline_arguments, 
                    extract_pipeline_arguments, prepare_qc_density_data, parse_vars, degs_json,
                    r_sanitize, sanitize_inf)

# -------------------------------

def generate_report(adatas, raw_counts, qc_adatas, adata, adata_nohm, metadata_df, 
                    scrublet_scores, summary_df, qc_summary_df, rmd_file, output_file, 
                    disease_label, tissue_label,
                    date, runtime_minute, ui, args):
    """
    Generate a full RMarkdown report by summarizing analysis results from the pipeline.

    This function:
    - Gathers all necessary summary statistics, QC results, annotations, and markers.
    - Builds a structured report object and passes it to RMarkdown via a temporary file.
    - Calls `rmarkdown::render()` to generate an HTML report from the RMarkdown template.

    Parameters:
        adatas, raw_counts, qc_adatas: Input and filtered AnnData objects.
        adata, adata_nohm: Merged and optionally batch-corrected AnnData objects.
        metadata_df, summary_df, qc_summary_df: Sample metadata and summary tables.
        scrublet_scores: Doublet scores if applicable.
        rmd_file: Path to the RMarkdown template.
        output_file: Desired HTML output path.
        disease_label, tissue_label: Ontology metadata.
        date: Date of execution.
        runtime_minute: Total runtime of the pipeline in minutes.
        ui: Unique pipeline ID.
        args: Full command-line argument namespace.
    """
    
    # -------------------------------
    # Compute average UMI and gene counts per cell (before QC)    
    mtx = raw_counts.X # Use raw counts stored 
    # Total UMI per cell 
    count_by_cell = round(np.array(mtx.sum(axis=1)).flatten().mean())
    # Number of genes detected per cell
    gene_by_cell = round(np.array((mtx > 0).sum(axis=1)).flatten().mean())

    # -------------------------------
    # Compute per-sample QC stats
    umi_per_cell_list, gene_per_cell_list, mt_per_cell_list = compute_qc_stats_objs(metadata_df, qc_adatas)
    count_by_cell_f, gene_by_cell_f, pct_mito_f = compute_qc_stats_obj(adata)

    # -------------------------------
    # Prepare plots and threshold filters
    db_plots, filters = generate_qc_plots_and_filters(adatas, metadata_df, args)

    # -------------------------------
    # Compute sample barcode overlap matrices
    bc_raw, bc_ji = compute_barcode_overlap_matrices(adatas)
    bc_raw_qc, bc_ji_qc = compute_barcode_overlap_matrices(qc_adatas)

    # -------------------------------
    # Initialize report object
    metadata_df = metadata_df.fillna("") # in case these is NaN from QuartzBio
    rprt = {
        "fobj_data": {
            "umap": adata.obsm["X_umap"].tolist(),
            "cluster": adata.obs["leiden_cluster"].tolist(),
            "sample": adata.obs["sample"].tolist(),
            "sample_id": adata.obs["sample_id"].tolist(),
            "total_counts": adata.obs["total_counts"].tolist(),
            "n_genes_by_counts": adata.obs["n_genes_by_counts"].tolist(),
            "celltype_scimilarity": adata.obs["celltype_scimilarity"].tolist() if "celltype_scimilarity" in adata.obs.columns else None,
            "cellstate_scimilarity": adata.obs["cellstate_scimilarity"].tolist() if "cellstate_scimilarity" in adata.obs.columns else None,
            "celltype_celltypist": adata.obs["celltype_celltypist"].tolist() if "celltype_celltypist" in adata.obs.columns else None,
            "cellstate_celltypist": adata.obs["cellstate_celltypist"].tolist() if "cellstate_celltypist" in adata.obs.columns else None,
            "scevan_class": adata.obs["scevan_class"].tolist() if "scevan_class" in adata.obs.columns else None
        },
        "metadata_df": metadata_df.to_dict(orient="records"),
        "summary_df": summary_df.to_dict(orient="records"),
        "qc_summary_df": qc_summary_df.to_dict(orient="records"),
        "count_by_cell": count_by_cell,
        "gene_by_cell": gene_by_cell,
        "umi_per_cell": umi_per_cell_list,
        "gene_per_cell": gene_per_cell_list,
        "mt_per_cell": mt_per_cell_list,
        "count_by_cell_f": count_by_cell_f, 
        "gene_by_cell_f": gene_by_cell_f, 
        "pct_mito_f": pct_mito_f,
        "db_plots": db_plots, 
        "filters": filters,
        "bc_raw": bc_raw,
        "bc_ji": bc_ji,
        "bc_raw_qc": bc_raw_qc, 
        "bc_ji_qc": bc_ji_qc
    }

    # -------------------------------
    # Conditionally add tSNE if requested
    if args.compute_tsne == "yes":
        rprt["fobj_data"]["tsne"] = adata.obsm["X_tsne"].tolist()

    # -------------------------------
    # If non-harmonized data is available, include it
    if adata_nohm is not None:
        rprt["non_hm_obj"] = {
            "umap": adata_nohm.obsm["X_umap"].tolist(),
            "cluster": adata_nohm.obs["leiden_cluster"].tolist(),            
            "sample": adata_nohm.obs["sample"].tolist(),
            "sample_id": adata.obs["sample_id"].tolist()
        }
        btchvrs = parse_vars(args.batch_vars)
        # Add batch variable values for each cell
        for var in btchvrs:
            rprt["fobj_data"][var] = adata.obs[var].tolist()
            rprt["non_hm_obj"][var] = adata_nohm.obs[var].tolist()

    # -------------------------------
    # Conditionally add tSNE if requested
    if args.compute_tsne == "yes" and args.batch_correction is not None:
        rprt["non_hm_obj"]["tsne"] = adata_nohm.obsm["X_tsne"].tolist()

    # -------------------------------
    # extract scrublet data
    if args.doublet_method is not None:
        scrublet_data_obs, scrublet_data_sim = scrublet_data(scrublet_scores)

        rprt["scrublet_data_obs"] = scrublet_data_obs
        rprt["scrublet_data_sim"] = scrublet_data_sim

    # -------------------------------
    # extract metadata description
    rprt["metadata_description"] = computed_metadata_description(adata)

    # -------------------------------
    # extract the impact of QC
    rprt["qc_impact"] = qc_impact_data(qc_adatas)

    # -------------------------------
    # extract clsuers markers
    if args.top_n_deg_leidn != 0: 
        rprt["clusters_markers"] = degs_json(adata, groupby="leiden_cluster") 

    # -------------------------------
    # extract celltype annotation outcomes
    if args.annotation_method:
        methods = parse_vars(args.annotation_method) # Allow multiple methods        
        
        if "scimilarity" in methods:
            rprt["scimilarity_network"] = prepare_visNetwork(adata, ontology_map, from_cell="cellstate_scimilarity", to_cell="celltype_scimilarity")
            if args.top_n_deg_scim != 0:
                rprt["scimilarity_markers"] = degs_json(adata, groupby="cellstate_scimilarity")  # extract _scimilarity markers

        if "celltypist" in methods:
            rprt["celltypist_network"] = prepare_visNetwork(adata, ontology_map, from_cell="cellstate_celltypist", to_cell="celltype_celltypist")
            if args.top_n_deg_cltpst != 0:
                rprt["celltypist_markers"] = degs_json(adata, groupby="cellstate_celltypist")  # extract celltypist markers

    # -------------------------------
    # Environment and config metadata
    rprt["python_environment"] = get_python_environment()
    rprt["pipeline_arguments"] = extract_pipeline_arguments(args, pipeline_arguments)

    # -------------------------------
    # Prepare density plot data
    rprt["density_plot_data"] = prepare_qc_density_data(adata, include_tsne="X_tsne" in adata.obsm)
    rprt = sanitize_inf(rprt)

    # -------------------------------
    # Write the JSON to a temp file (safer than passing huge JSON string directly via CLI)
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".json") as tmpfile:
        json.dump(rprt, tmpfile)
        tmpfile_path = tmpfile.name

     # Convert args to dictionary from argparse.Namespace
    args_dict = vars(args)
    opt_str = ", ".join(f"{k}={r_sanitize(v)}" for k, v in sorted(args_dict.items()))

    params_str = (
        f"snapshot_file='{tmpfile_path}', "
        f"opt=list({opt_str}), "
        f"date={r_sanitize(date)}, "
        f"runtime_minute={r_sanitize(runtime_minute)}, "
        f"ui={r_sanitize(ui)}, "
        f"disease_label={r_sanitize(disease_label)}, "
        f"tissue_label={r_sanitize(tissue_label)}"
    )

    # -------------------------------
    # Call RMarkdown to render report
    cmd = [
        "Rscript", "-e",
        (
            f"options(pandoc.stack.size = '4096m'); "  # Set stack size to 4GB
            f"rmarkdown::render("
            f"'{rmd_file}', "
            f"output_file='{output_file}', "
            f"params=list({params_str})"
            f")"
        )
    ]

    try:
        subprocess.run(cmd, check=True)
        print(f"*** ✅ Summary report generated: {output_file}")
        os.remove(tmpfile_path)
    except subprocess.CalledProcessError as e:
        print(f"*** 🚨 Error running report: {e}")
        os.remove(tmpfile_path)
        sys.exit(1)  # Exit with a non-zero status to indicate failure

    # -------------------------------