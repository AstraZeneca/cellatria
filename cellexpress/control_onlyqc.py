# control.py
# -------------------------------

import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# -------------------------------

import contextlib
import io
import shutil
from checks import checks_args
from readin import read_in
from qcfilter import qc_filter
from report_onlyqc_wrapper import generate_onlyqc_report
from helper import total_unique_genes, qc_verbose_ls
from datetime import datetime

# -------------------------------

def control_onlyqc_pipe(args):
    """
    Executes the QC-only workflow of the CellExpress pipeline.

    This mode performs an initial assessment of the raw single-cell dataset
    without applying any filtering criteria. It computes basic quality control
    (QC) metrics such as UMI counts, gene counts, and mitochondrial content per cell,
    and generates a summary report for inspection.

    Specifically, the function:
        - Loads the metadata and per-sample expression data
        - Computes the total number of cells and genes across samples
        - Forces all QC thresholds to neutral values (i.e., no filtering)
        - Computes QC metrics for each sample, including mitochondrial content
        - Generates an RMarkdown-based report summarizing per-sample metrics

    Args:
        args (Namespace): Parsed command-line arguments.
    """

    # -------------------------------
    # Set pipeline constants
    setattr(args, "pipe_version", "1-0-0") # add pipeline version
    setattr(args, "package_path", "/mnt/work/projects/cellatria")  # add package path
    setattr(args, "cellexpress_path", os.path.join(args.package_path, "cellexpress"))  # add package path

    # -------------------------------
    # Argument validation
    args = checks_args(args)

    # -------------------------------
    # Load pre-QC data
    if hasattr(args, "vault_id"):
        print("*** 🔄 Accessing QuartzBio for downloading data...")
        data_path = run_quartzbio(args)
        setattr(args, "vault_path", args.input)  # retrive new data path
        setattr(args, "input", data_path)  # retrive new data path
        delattr(args, "qb_token")
    
    # -------------------------------
    # Create outputs directory for saving results
    outputs_path = os.path.join(args.input, "outputs")
    os.makedirs(outputs_path, exist_ok=True)  
    setattr(args, "outputs_path", outputs_path)
    print(f"*** ✅ Created (or found existing) outputs directory: {outputs_path}")

    # -------------------------------
    # Load expression matrices and sample metadata
    print("*** 🔄 loading in the single-cell data...")
    adatas, metadata_df, summary_df = read_in(args)
    qc_verbose_ls(metadata_df, adatas)

    # -------------------------------
    # Summarize cell and gene statistics across all samples
    total_cells = sum(adata.n_obs for adata in adatas.values())
    total_genes = total_unique_genes(adatas)
    print(f"*** 🔔 total {total_cells:,} cells and {total_genes:,} genes retained across {len(adatas):,} samples.")
    setattr(args, "pre_qc_cells", total_cells)

    # -------------------------------
    # Override QC args to be very permissive (liberal)
    args.min_cell = 1                   # Keep genes expressed in at least 1 cell
    args.min_genes_per_cell = 1         # Keep cells with ≥1 gene
    args.max_genes_per_cell = None      # No upper limit
    args.min_umi_per_cell = 1           # Keep cells with ≥1 count
    args.max_umi_per_cell = None        # No upper limit
    args.max_mt_percent = 100.0         # Allow up to 100% mitochondrial content

    print("*** 🔄 Removing empty cells and genes...")
    qc_adatas, qc_summary_df = qc_filter(adatas, metadata_df, args)
    qc_verbose_ls(metadata_df, qc_adatas)

    # -------------------------------
    # Generate RMarkdown-based summary report
    print("*** 🔄 Creating report summary...")
    generate_onlyqc_report(adatas = qc_adatas,
                           metadata_df = metadata_df,
                           summary_df = qc_summary_df,
                           rmd_file = os.path.join(args.cellexpress_path, "report_onlyqc.Rmd"), 
                           output_file = os.path.join(args.outputs_path, "report_onlyqc_cellexpress.html"), 
                           disease_label = args.disease,
                           tissue_label = args.tissue,
                           args = args)

    # -------------------------------
    if hasattr(args, "vault_id"):
        print(f"*** 🧹 Removing {len(metadata_df)} samples downloaded from QuartzBio project: {args.vault_project}")
        for sample in metadata_df["sample"]:
            shutil.rmtree(os.path.join(args.input, sample))
            print(f"*** ✅ Removed: {sample}")
    
    # -------------------------------
    print("*** 🎯 QC pipeline execution completed successfully.")  