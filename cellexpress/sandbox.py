import os
import sys

os.environ["OPENBLAS_NUM_THREADS"] = "1" # Controls the number of threads used by OpenBLAS, which backs many NumPy and SciPy operations (Harmony related).
os.environ["OMP_NUM_THREADS"] = "1" # Controls threads used by OpenMP, which is leveraged by many C/C++/Fortran numerical libraries (Scikit-learn related).
os.environ["MKL_NUM_THREADS"] = "1" # Controls threads for Intel MKL, the math kernel library used by NumPy when compiled with Intel compilers.
os.environ["NUMEXPR_NUM_THREADS"] = "1" # Controls threading in numexpr, a library that speeds up array expressions
print("*** 🔒 Thread limits applied to avoid memory issues.")

import argparse

# Manually create an argparse Namespace object
args = argparse.Namespace(
    input="/mnt/work/projects/cellatria/data/my_projects/infants_vaccine",
    project="Newborn_PBMC_Immunization",
    species="hs",
    tissue="blood",
    disease="influenza",
    metadata="sample_based",
    min_umi_per_cell=750,
    max_umi_per_cell=None,
    max_mt_percent=10.0,
    min_genes_per_cell=250,
    max_genes_per_cell=None,
    min_cell=3,
    pipe_version="1.0.0",
    doublet_method=None,
    scrublet_cutoff=0.25,
    norm_target_sum=1e4,
    n_top_genes=2000,
    regress_out="no",
    scale_max_value=10,
    n_pcs=30,
    batch_correction=None, #"harmony",
    batch_vars="sample",
    n_neighbors=15,
    resolution=0.6,
    compute_tsne="no",
    annotation_method = "scimilarity,celltypist",
    sci_model_path = "/mnt/work/projects/cellatria/data/scimilarity/model_v1.1",
    cty_model_path = "/mnt/work/projects/cellatria/data/celltypist/model_v1.6.3",
    cty_model_name= "Immune_All_Low",
    pval_threshold=0.05,
    logfc_threshold=0.25,
    dea_method="wilcoxon",
    top_n_deg_leidn=0,
    top_n_deg_scim=0,
    top_n_deg_cltpst=0,
    pts_threshold=0.1, 
    doc_url="https://www.nature.com/articles/s41467-023-43758-2",
    data_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204716",
    # tumor_id= "scevan",
    only_qc="yes",
    qb_token=None,
    bt_token=None,
    fix_gene_names="no", # "original_gene_symbols"
    env="ws"
)