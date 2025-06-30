# checks.py
# -------------------------------

import os
import sys
import re
import pandas as pd
from helper import (convert_none, parse_vars)

# -------------------------------

def checks_args(args):
    """
    Validates and normalizes user-provided arguments to ensure pipeline consistency.

    Performs the following:
    - Converts "None" strings to actual `None` objects.
    - Verifies existence of the input directory and metadata file.
    - Checks required metadata fields (e.g., 'sample').
    - Checks forbidden metadata fields (e.g., 'sample_id').
    - Validates species, doublet detection method, batch correction settings.
    - Ensures SCimilarity and CellTypist models are configured if selected.
    - Confirms valid TSNE setting.

    Args:
        args (Namespace): Parsed command-line arguments.

    Returns:
        tuple: (args)

    Raises:
        FileNotFoundError: If input or metadata file is missing.
        ValueError: For invalid species, missing metadata columns, or improper settings.
    """
    
    # -------------------------------
    # Convert all 'None' string values to real Python None
    for arg_name, arg_value in vars(args).items():
        setattr(args, arg_name, convert_none(arg_value))  # Apply conversion

    # -------------------------------

    # Ensure input path is absolute path
    if not os.path.isabs(args.input):
        raise ValueError("*** 🚨 Please provide an absolute path for --input (e.g., /path-to/project), not a relative one like ./project")

    # Ensure no trailing slash
    if args.input.endswith("/") or args.input.endswith("\\"):
        raise ValueError("*** 🚨 The input path should not have a trailing slash. Use /path-to/project instead of /path-to/project/")

    # Ensure input directory exists
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"*** 🚨 Input directory not found: {args.input}")

    # Load and validate metadata
    metadata_file = os.path.join(args.input, "metadata.csv")
    # Check if metadata file exists
    if not os.path.exists(metadata_file):
        raise FileNotFoundError(f"*** 🚨 Metadata file not found in input directory: {metadata_file}")

    # Read metadata
    metadata = pd.read_csv(metadata_file)
    # Ensure 'sample' column exists
    if "sample" not in metadata.columns:
        raise ValueError("*** 🚨 Metadata file must contain a 'sample' column.")

    # Reserved keyword check
    if "sample_id" in metadata.columns:
        raise ValueError("*** 🚨 'sample_id' is a reserved keyword used internally by the pipeline. "
                        "Please remove it from your input file or rename it to avoid conflict.")

    # -------------------------------
    # Species validation
    valid_species = {"hs", "mm"}
    if args.species not in valid_species:
        raise ValueError(f"*** 🚨 Invalid species '{args.species}'. Must be one of {valid_species}.")

    # -------------------------------
    # Doublet detection method checks
    if args.doublet_method == "scrublet" and args.scrublet_cutoff is None:
        raise ValueError("*** 🚨 Argument 'scrublet_cutoff' must be provided when using 'scrublet' for doublet detection.")

    # -------------------------------
    # Batch correction check
    if args.batch_correction is not None: 
        if args.batch_correction.lower() == "harmony":
            if not args.batch_vars:
                raise ValueError(f"*** 🚨 Batch correction selected, but --batch_vars was not provided. "
                                 "Specify column(s) in adata.obs to use for Harmony batch correction. "
                                "(comma-separated if multiple, e.g., 'sample_id,time_point').")

    # -------------------------------
     # Annotation method check
    if args.annotation_method is not None:
        methods = parse_vars(args.annotation_method) # Allow multiple methods
        # Validate each method separately
        valid_methods = {"scimilarity", "celltypist"}
        invalid_methods = set(methods) - valid_methods
        if invalid_methods:
            raise ValueError(f"*** 🚨 Invalid annotation method(s): {', '.join(invalid_methods)}. "
                             f"Choose one or both from ['scimilarity', 'celltypist'].")

        # SCimilarity checks
        if "scimilarity" in methods:
            if args.sci_model_path is None:
                raise ValueError("*** 🚨 SCimilarity annotation requires a model path. Provide --sci_model_path.")        
            if not os.path.exists(args.sci_model_path):
                raise ValueError(f"*** 🚨 The SCimilarity model path '{args.sci_model_path}' does not exist. "
                                 "Provide a valid directory containing the model.")
        # CellTypist checks
        if "celltypist" in methods:
            if args.cty_model_path is None:
                raise ValueError("*** 🚨 CellTypist annotation requires a model path. Provide --cty_model_path.")
            if not os.path.exists(args.cty_model_path):
                raise ValueError(f"*** 🚨 The CellTypist model path '{args.cty_model_path}' does not exist. "
                                 "Provide a valid directory containing the model.")
            if args.cty_model_name is None:
                raise ValueError("*** 🚨 CellTypist annotation requires a model name. Provide --cty_model_name.")
                # Construct full model path and check if the file exists
                model_file = os.path.join(args.cty_model_path, args.cty_model_name + ".pkl")
                if not os.path.exists(model_file):
                    raise FileNotFoundError(f"*** 🚨 CellTypist model file not found: {args.cty_model_name + '.pkl'}. "
                                            "Ensure the model file exists in the specified directory.")
    # -------------------------------
    # TSNE computation check
    valid_options = {"yes", "no"}
    if args.compute_tsne.lower() not in valid_options:
        raise ValueError(f"*** 🚨 Invalid value '{args.compute_tsne}' for --compute_tsne. "
                         f"Expected one of {valid_options}.")

    # -------------------------------

    print(f"*** ✅ All checks passed.")
    return args

    # -------------------------------