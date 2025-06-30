#!/usr/bin/env python3
# -------------------------------

import os
import sys
# sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# -------------------------------

# Apply thread limits BEFORE any numerical imports. Improves stability in large scale single cell data
# To prevent memory exhaustion, thread over-subscription, and thread collisions across libraries
if "--limit_threads" in sys.argv:
    idx = sys.argv.index("--limit_threads")
    try:
        n_threads = int(sys.argv[idx + 1])
    except (IndexError, ValueError):
        print("*** 🚨 Error: --limit_threads must be followed by an integer (e.g. --limit_threads 1)")
        sys.exit(1)
else:
    n_threads = 1  # Default if not specified

os.environ["OPENBLAS_NUM_THREADS"] = str(n_threads) # Controls the number of threads used by OpenBLAS, which backs many NumPy and SciPy operations (Harmony related).
os.environ["OMP_NUM_THREADS"] = str(n_threads) # Controls threads used by OpenMP, which is leveraged by many C/C++/Fortran numerical libraries (Scikit-learn related).
os.environ["MKL_NUM_THREADS"] = str(n_threads) # Controls threads for Intel MKL, the math kernel library used by NumPy when compiled with Intel compilers.
os.environ["NUMEXPR_NUM_THREADS"] = str(n_threads) # Controls threading in numexpr, a library that speeds up array expressions

print(f"*** 🔒 Thread limits applied: {n_threads} thread(s) per library.")

# -------------------------------

import argparse
from helper import none_or_int, none_or_float, none_or_str

# -------------------------------

# Argument parser setup
parser = argparse.ArgumentParser(description="Pipeline for single-cell data processing")

# General arguments
parser.add_argument("--input",    
                    type=none_or_str, 
                    required=True, 
                    default=None, 
                    help=("Path to input data directory. Must meet the following criteria:\n"
                        "• Absolute path (e.g., /domino/datasets/local/CTPM-CellExpress/projects/my_project)\n"
                        "• No trailing slash (/)\n"
                        "• Directory must exist on the file system"))
parser.add_argument("--project",  
                    type=none_or_str, 
                    required=True, 
                    default=None, 
                    help="Project name")
parser.add_argument("--species",  
                    type=none_or_str, 
                    required=True, 
                    default=None, 
                    help=("Options: 'hs' or 'mm'. Use 'hs' for human and 'mm' for mouse." 
                        "Important to claculate the mitochondrial gene percentage."))
parser.add_argument("--tissue",   
                    type=none_or_str, 
                    required=True, 
                    default=None, 
                    help=("Tissue name."))
parser.add_argument("--disease",   
                    type=none_or_str, 
                    required=True, 
                    default=None, 
                    help=("Disease name."))

# Quality control arguments
parser.add_argument("--min_umi_per_cell", 
                    type=none_or_int, 
                    default=750, 
                    help="Minimum UMI counts per cell")
parser.add_argument("--max_umi_per_cell", 
                    type=none_or_int, 
                    default=None, 
                    help="Maximum UMI counts per cell")
parser.add_argument("--min_genes_per_cell", 
                    type=none_or_int, 
                    default=250, 
                    help="Minimum genes per cell")
parser.add_argument("--max_genes_per_cell", 
                    type=none_or_int, 
                    default=None, 
                    help="Maximum genes per cell")
parser.add_argument("--min_cell", 
                    type=none_or_int, 
                    default=3, 
                    help="Genes detected in at least this many cells")
parser.add_argument("--max_mt_percent", 
                    type=none_or_float, 
                    default=15, 
                    help="Max mitochondrial gene percentage")
parser.add_argument("--doublet_method", 
                    type=none_or_str, 
                    default=None, 
                    help="Method for doublet identification. Options: 'scrublet'. Default is None")
parser.add_argument("--scrublet_cutoff", 
                    type=none_or_float, 
                    default=0.25, 
                    help="Maximum allowed doublet score (applies only if doublet_method='scrublet')")

# Analysis arguments
parser.add_argument("--norm_target_sum", 
                    type=none_or_float, 
                    default=1e4,
                    help="Target sum for total counts per cell during normalization (default: 1e4)")
parser.add_argument("--n_top_genes", 
                    type=none_or_int, 
                    default=2000,
                    help="Number of top highly variable genes to retain (default: 2000)")
parser.add_argument("--regress_out",
                    type=none_or_str,
                    default="no",
                    help=("Options: 'yes' or 'no'.Regress out total counts and mitochondrial percentage before scaling. "
                        "This is recommended for smaller datasets (<20k cells) where technical noise could distort biological signals. "
                        "Consider enabling this if your data comes from different platforms, batches, or protocols with variable sequencing depth. "
                        "For larger datasets, skipping this step significantly improves speed and reduces memory usage. (default: no)"))
parser.add_argument("--scale_max_value", 
                    type=none_or_float, 
                    default=10,
                    help="Clip (truncate) to this value after scaling. (default: 10)")
parser.add_argument("--n_pcs",
                    type=none_or_int, 
                    default=30,
                    help="Number of principal components to compute (default: 30)")
parser.add_argument("--batch_correction",
                    type=none_or_str,
                    default=None,
                    help="Choose batch correction method: 'harmony'. Default is None.")
parser.add_argument("--batch_vars", 
                    type=none_or_str, 
                    default=None,
                    help="Column(s) in adata.obs to use for Harmony batch correction "
                         "(comma-separated if multiple, e.g., 'sample_id,time_point').")
parser.add_argument("--n_neighbors", 
                    type=none_or_int, 
                    default=15,
                    help="Number of neighbors for k-NN graph (default: 15)")
parser.add_argument("--resolution", 
                    type=none_or_float, 
                    default=0.6,
                    help="Resolution parameter for Leiden clustering (default: 0.6)")
parser.add_argument("--compute_tsne",
                    type=none_or_str,
                    default="no",
                    help="Options: 'yes' or 'no'. Option to compute t-SNE embedding in addition to UMAP. Time consuming in large datasets. (default: no)")

# Annotation method argument
parser.add_argument("--annotation_method",
                    type=none_or_str,
                    default=None,
                    help=("Annotation method(s) to use: 'scimilarity', 'celltypist', or both. "
                        "For multiple methods, separate with commas (e.g., 'scimilarity,celltypist'). "
                        "Default: None.")
    )
# SCimilarity model directory argument
parser.add_argument("--sci_model_path",
                    type=str,
                    default=None,
                    help=("Path to the SCimilarity model directory. "
                        "Required if --annotation_method='scimilarity'. "
                        "Ensure this path contains the trained SCimilarity model files."))
# CellTypist model directory argument
parser.add_argument("--cty_model_path",
                    type=str,
                    default=None,
                    help=("Path to the directory containing the CellTypist model file. "
                          "Required if using CellTypist for annotation."))
# CellTypist model name argument
parser.add_argument("--cty_model_name",
                    type=none_or_str,
                    default=None,
                    help=("Name of the CellTypist model to use for annotation (without '.pkl' extension). "
                          "Users can browse available models at: https://www.celltypist.org/models. "
                          "Example: 'Immune_All_High'.")) 

# Differentially expressed analysis arguments
parser.add_argument("--pval_threshold",
                    type=none_or_float,
                    default=0.05,
                    help="P-value threshold for filtering marker genes (DEGs). "
                         "Genes with adjusted p-values below this threshold are considered significant.")
parser.add_argument("--logfc_threshold",
                    type=none_or_float,
                    default=0.25,
                    help="Log fold-change (logFC) threshold for filtering marker genes. "
                         "Genes with absolute logFC above this value are retained.")
parser.add_argument("--dea_method",
                    type=none_or_str,
                    default="wilcoxon",
                    help="Statistical test used for marker genes identification. "
                         "Options: 'wilcoxon' (default), 't-test', 't-test_overestim_var', or 'logreg'.")
parser.add_argument("--top_n_deg_leidn",
                    type=none_or_int,
                    default=100,  # Default to 100 genes per group
                    help="Number of top marker genes to return per 'leiden' groups. "
                        "If set to 0, skip.")
parser.add_argument("--top_n_deg_scim",
                    type=none_or_int,
                    default=100,  # Default to 100 genes per group
                    help="Number of top marker genes to return per 'scimilarity' annotated cells. "
                        "If set to 0, skip.")
parser.add_argument("--top_n_deg_cltpst",
                    type=none_or_int,
                    default=100,  # Default to 100 genes per group
                    help="Number of top marker genes to return per 'celltypist' annotated cells. "
                        "If set to 0, skip.")
parser.add_argument("--pts_threshold",
                    type=none_or_float,
                    default=0.1,  # Example default, adjust as needed
                    help="Minimum fraction of cells expressing a gene for it to be considered a marker. (default: 0.1)")

# references
parser.add_argument("--doc_url",
                    type=none_or_str,
                    default=None,
                    help="URL to the documentation associated with the data (optional).")
parser.add_argument("--data_url",
                    type=none_or_str,
                    default=None,
                    help="URL to the publicly available dataset or repository (optional).")

# only QC
parser.add_argument("--only_qc",
                    type=none_or_str,
                    default="no",
                    help=("Options: 'yes' or 'no'. Set to 'yes' to generate an overview of quality control (QC) metrics"
                           "prior to running the main pipeline (default: 'no')."))

parser.add_argument("--fix_gene_names",
                    type=none_or_str,
                    default="no",
                    help=("Fix gene names in `adata.var_names` if Ensembl IDs are detected."
                         "Specify column name in adata.var that may contain gene symbols. (default: 'no')."))
parser.add_argument("--limit_threads",
                    type=none_or_int,
                    default=1,
                    help="Apply thread limits to avoid OpenBLAS/OMP crashes (recommended for large datasets (>500K)).")

# -------------------------------

# Parse arguments
args = parser.parse_args()

# -------------------------------

from control import control_pipe
from control_onlyqc import control_onlyqc_pipe

# -------------------------------

# Call control pipeline
if args.only_qc.lower() == "no":
    control_pipe(args)
else:
    control_onlyqc_pipe(args)
 
# -------------------------------