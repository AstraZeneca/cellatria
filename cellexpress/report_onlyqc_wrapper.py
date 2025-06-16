# report_onlyqc_wrapper.py
# -------------------------------

import os
import sys
sys.path.append("/mnt/work/projects/cellatria")

# -------------------------------

import json
import subprocess
import tempfile
import numpy as np
import pandas as pd
from cellexpress.helper import (compute_qc_stats_objs, generate_qc_plots_and_filters, compute_barcode_overlap_matrices)

# -------------------------------

def generate_onlyqc_report(adatas, metadata_df, summary_df, rmd_file, output_file, 
                            disease_id, disease_label, tissue_id, tissue_label, args):
    """
    Wrapper to generate an RMarkdown report, passing data from two AnnData objects
    (harmonized and non-harmonized) directly into the RMarkdown report via parameters.
    
    Inputs:
        adata - harmonized AnnData object (final object)
        adata_nohm - non-harmonized AnnData object
        rmd_file - path to the RMarkdown file
        output_file - desired output report file (HTML)
    """
    
    # -------------------------------
    # Compute per-sample QC stats
    umi_per_cell_list, gene_per_cell_list, mt_per_cell_list = compute_qc_stats_objs(metadata_df, adatas)

    # -------------------------------
    # Prepare plots and threshold filters
    db_plots, filters = generate_qc_plots_and_filters(adatas, metadata_df, args)

    # -------------------------------
    # Compute sample barcode overlap matrices
    bc_qc_raw, bc_qc_ji = compute_barcode_overlap_matrices(adatas)

    # -------------------------------
    # Initialize report object
    metadata_df = metadata_df.fillna("") # in case these is NaN from QuartzBio
    rprt = {
        "metadata_df": metadata_df.to_dict(orient="records"),
        "summary_df": summary_df.to_dict(orient="records"),
        "umi_per_cell": umi_per_cell_list,
        "gene_per_cell": gene_per_cell_list,
        "mt_per_cell": mt_per_cell_list,
        "db_plots": db_plots, 
        "filters": filters,
        "bc_qc_raw": bc_qc_raw, 
        "bc_qc_ji": bc_qc_ji 
    }

    # -------------------------------
    # Write the JSON to a temp file (safer than passing huge JSON string directly via CLI)
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".json") as tmpfile:
        json.dump(rprt, tmpfile)
        tmpfile_path = tmpfile.name

     # Convert args to dictionary from argparse.Namespace
    args_dict = vars(args)

    # Build the params object for RMarkdown — include `snapshot_file` and `opt`
    params = {
        "snapshot_file": tmpfile_path,
        "opt": args_dict,  # All pipeline args go into opt as a nested structure
        "disease_id": disease_id, 
        "disease_label": disease_label,
        "tissue_id": tissue_id, 
        "tissue_label": tissue_label
    }

    # Function to format params into R-compatible strings for rmarkdown::render
    def format_param(key, value):
        if isinstance(value, dict):
            # Use single quotes around keys/values for R compatibility, but avoid backslashes
            items = [f"{k}='{v}'" for k, v in value.items()]   # No backslashes needed
            return f"{key}=list({', '.join(items)})"
        elif isinstance(value, str):
            return f"{key}='{value}'"   # Use single quotes here too
        else:
            return f"{key}={value}"    # Numbers and booleans directly

    params_str = ", ".join(format_param(k, v) for k, v in params.items())

    # -------------------------------
    # Call RMarkdown with all params
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