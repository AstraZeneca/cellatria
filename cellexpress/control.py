# control.py
# -------------------------------

import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# -------------------------------

import time
import json
import argparse
import shutil
from datetime import datetime
from checks import checks_args
from readin import read_in
from qcfilter import qc_filter
from doublets import doublets_id
from merge_samples import merge_samples
from stndrd_analysis import run_analysis
from report_wrapper import generate_report
from dea_analysis import compute_degs
from call_celltypist import run_celltypist
from call_scimilarity import run_scimilarity
from helper import (create_unique_ids, total_unique_genes, parse_vars, qc_verbose, 
                    qc_verbose_ls, summarize_adata_structure)

# -------------------------------
def control_pipe(args):
    """
    Main orchestrator for the CellExpress single-cell RNA-seq analysis pipeline.

    This function coordinates all stages of the pipeline including:
    - Argument validation
    - Data loading and QC filtering
    - Doublet detection
    - Data integration and dimensionality reduction
    - Cell type annotation
    - Tumor identification
    - DEA computation
    - Summary report generation and configuration saving

    Args:
        args (Namespace): Parsed command-line arguments from argparse.
    """

    # -------------------------------
    print("*** 🚀 Starting pipeline...")
    start_time = time.time()  # Record start time

    # -------------------------------
    # Set pipeline constants
    setattr(args, "pipe_version", "1-0-0") # add pipeline version
    setattr(args, "package_path", "/mnt/work/projects/cellatria")  # add package path
    setattr(args, "cellexpress_path", os.path.join(args.package_path, "cellexpress"))  # add package path
    setattr(args, "resources_path", os.path.join(args.package_path, "resources"))

    # -------------------------------
    # Argument validation
    args = checks_args(args)

    # -------------------------------
    # Generate unique run ID
    ui = f"cellexpress_v{args.pipe_version}_{create_unique_ids(1, char_len=7)[0]}" 
    print(f"*** ✅ Created unique ID '{ui}'")

    # -------------------------------
    # Create output directory
    outputs_path = os.path.join(f"{args.input}", f"outputs_{ui}")
    os.makedirs(outputs_path, exist_ok=True)  
    setattr(args, "outputs_path", outputs_path)
    print(f"*** ✅ Created (or found existing) outputs directory: {outputs_path}")

    # -------------------------------
    # Load pre-QC data
    print("*** 🔄 loading in the single-cell data...")
    adatas, metadata_df, summary_df = read_in(args)
    qc_verbose_ls(metadata_df, adatas)

    # -------------------------------
    # Summarize total number of retained cells across all samples
    pre_qc_cells = sum(adata.n_obs for adata in adatas.values())
    genes = total_unique_genes(adatas)
    print(f"*** 🔔 total {pre_qc_cells:,} cells and {genes:,} genes retained across {len(adatas):,} samples.")

    # -------------------------------
     # Apply QC filtering
    print("*** 🔄 Removing low quality cells and genes...")
    qc_adatas, qc_summary_df = qc_filter(adatas, metadata_df, args)
    qc_verbose_ls(metadata_df, qc_adatas)

    # -------------------------------
    # Doublet detection
    scrublet_scores = None
    if args.doublet_method is None:
        print("*** 🚫 Doublet identification skipped (no method specified).")
    else:
        print("*** 🔄 Running doublet identification...")
        qc_adatas, qc_summary_df, scrublet_scores = doublets_id(qc_adatas, args)
        qc_verbose_ls(metadata_df, qc_adatas)

    # -------------------------------
    # Summarize total number of retained cells across all samples
    post_qc_cells = sum(adata.n_obs for adata in qc_adatas.values())
    genes = total_unique_genes(qc_adatas)
    print(f"*** 🔔 total {post_qc_cells:,} cells and {genes:,} genes retained across {len(qc_adatas):,} samples  after QC.")

    # -------------------------------
    ## Merge samples
    print("*** 🔄 Merging all samples into a combined dataset...")
    adata = merge_samples(qc_adatas, args)
    qc_verbose(adata)
    
    # -------------------------------
    # Run standard analysis
    print("*** 🔄 Running standard analysis pipeline...")
    raw_counts, adata, adata_nohm = run_analysis(adata, args)

    # Run DEA for cluster
    if args.top_n_deg_leidn != 0:
        adata = compute_degs(adata, groupby="leiden_cluster", pts=True, dea_method=args.dea_method, n_genes=args.top_n_deg_leidn,
                            pval_threshold=args.pval_threshold, logfc_threshold=args.logfc_threshold, pts_threshold=args.pts_threshold)

    # -------------------------------
    # Cell Annotation (SCimilarity and/or CellTypist)
    if args.annotation_method:
        methods = parse_vars(args.annotation_method) # Allow multiple methods

        if "scimilarity" in methods:
            print("*** 🔄 Running SCimilarity annotation...")
            adata = run_scimilarity(raw_counts, adata, args)  # Run SCimilarity
            # Run DEA for celltype_scimilarity
            if args.top_n_deg_scim != 0:
                adata = compute_degs(adata, groupby="cellstate_scimilarity", pts=True, dea_method=args.dea_method, n_genes=args.top_n_deg_scim,
                                    pval_threshold=args.pval_threshold, logfc_threshold=args.logfc_threshold, pts_threshold=args.pts_threshold)

        if "celltypist" in methods:
            print("*** 🔄 Running CellTypist annotation...")
            adata = run_celltypist(raw_counts, adata, args)  # Run CellTypist
            # Run DEA for celltype_celltypist
            if args.top_n_deg_cltpst != 0:
                adata = compute_degs(adata, groupby="cellstate_celltypist", pts=True,  dea_method=args.dea_method, n_genes=args.top_n_deg_cltpst,
                                    pval_threshold=args.pval_threshold, logfc_threshold=args.logfc_threshold, pts_threshold=args.pts_threshold)

    else:
        print("*** 🚫 No cell annotation performed (annotation_method='none').")

    # -------------------------------
    end_time = time.time()
    elapsed_time = end_time - start_time
    runtime_minute = round(elapsed_time / 60, 2)
    print(f"*** ⏱️  Pipeline completed in {int(elapsed_time // 60)} min {int(elapsed_time % 60)} sec.")

    # -------------------------------
    # Add date
    date = datetime.now().strftime("%d %B, %Y")
    # number of samples processed
    num_samples = adata.obs["sample_id"].nunique()

    # -------------------------------
    # Generate HTML report
    print("*** 🔄 Creating report summary...")
    generate_report(adatas = adatas,
                    raw_counts = raw_counts,
                    qc_adatas = qc_adatas,
                    adata = adata, 
                    adata_nohm = adata_nohm,
                    metadata_df = metadata_df,
                    scrublet_scores = scrublet_scores,
                    summary_df = summary_df,
                    qc_summary_df = qc_summary_df,
                    rmd_file = os.path.join(args.cellexpress_path, "report.Rmd"), 
                    output_file = os.path.join(args.outputs_path, f"report_{ui}_{datetime.today().date().isoformat()}.html"), 
                    disease_label = args.disease,
                    tissue_label = args.tissue,
                    date = date,
                    runtime_minute = runtime_minute,
                    ui = ui,
                    args = args)

    # -------------------------------
    # Store processed data
    print("*** 🔄 Storing data...")

    fileneme = os.path.join(args.outputs_path, f"counts-qced_{ui}_{datetime.today().date().isoformat()}.h5ad")
    raw_counts.write(fileneme)
    print(f"*** ✅ Data stord to: {fileneme}")

    fileneme = os.path.join(args.outputs_path, f"adata_{ui}_{datetime.today().date().isoformat()}.h5ad")
    json_adata = summarize_adata_structure(adata.raw.to_adata())
    adata.raw.to_adata().write(fileneme)    
    print(f"*** ✅ Data stord to: {fileneme}")

    # -------------------------------
    # Save configuration JSON
    json_dict = {
        "execution_summary": {
            "ui": ui,
            "date_iso": datetime.today().date().isoformat(),
            "species": args.species,
            "disease_label": args.disease,
            "tissue_label": args.tissue,
            "num_samples": num_samples,
            "pre_qc_cells": pre_qc_cells,
            "post_qc_cells": post_qc_cells,
            "runtime_minute": runtime_minute
            },
        "sample_metadata": metadata_df.to_dict(orient="records"),  # List of dictionaries, one per sample
        "settings": vars(args),  # all argparse parameters
        "adata_structure": json_adata
        }

    # Construct filename
    json_filename = os.path.join(args.outputs_path, f"config_{ui}_{datetime.today().date().isoformat()}.json")
    
    # Save to JSON
    with open(json_filename, "w") as f:
        json.dump(json_dict, f, indent=4)
    print(f"*** ✅ Configs saved to: {json_filename}")

    # -------------------------------
    print("*** 🎯 Pipeline execution completed successfully.")    
