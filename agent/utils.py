# utils.py
# -------------------------------

import os
import re
import json
import uuid
import pandas as pd
import gradio as gr
from typing import List, Dict, Any, TypedDict, Literal, Optional
from langchain_openai import AzureChatOpenAI, ChatOpenAI
from langchain_anthropic import ChatAnthropic  
from dotenv import load_dotenv
import google.generativeai as genai
from transformers import pipeline, AutoModelForCausalLM, AutoTokenizer

# -------------------------------

base_path = os.path.dirname(os.path.abspath(__file__))

# -------------------------------

def get_llm_from_env(path_to_env):
    load_dotenv(dotenv_path=os.path.join(path_to_env, ".env"))
    provider = os.getenv("PROVIDER", "").strip()

    if provider == "Azure OpenAI":
        return AzureChatOpenAI(
            deployment_name=os.getenv("AZURE_OPENAI_DEPLOYMENT_NAME"),
            openai_api_version=os.getenv("AZURE_OPENAI_API_VERSION"),
            azure_endpoint=os.getenv("AZURE_OPENAI_ENDPOINT"),
            openai_api_key=os.getenv("AZURE_OPENAI_API_KEY"),
        )
    elif provider == "OpenAI":
        return ChatOpenAI(
            model_name=os.getenv("OPENAI_MODEL", "gpt-4"),
            openai_api_key=os.getenv("OPENAI_API_KEY"),
        )
    elif provider == "Anthropic":
        return ChatAnthropic(api_key=os.getenv("ANTHROPIC_API_KEY"))
    elif provider == "Google":
        # Example using google.generativeai (Gemini)        
        genai.configure(api_key=os.getenv("GOOGLE_API_KEY"))
        model_name = os.getenv("GOOGLE_MODEL", "gemini-pro")
        # You may need to adjust this based on your actual Gemini integration
        return genai.GenerativeModel(model_name)
    elif provider == "Local":
        # Example using Hugging Face Transformers for a local model
        model_path = os.getenv("LOCAL_MODEL_PATH")
        if not model_path:
            raise ValueError("LOCAL_MODEL_PATH must be set for Local provider.")
        tokenizer = AutoTokenizer.from_pretrained(model_path)
        model = AutoModelForCausalLM.from_pretrained(model_path)
        return pipeline("text-generation", model=model, tokenizer=tokenizer)
    else:
        raise ValueError("Invalid or unsupported PROVIDER in .env")

# -------------------------------

def store_metadata_json(metadata: dict, project_dir: Optional[str] = None, filename: Optional[str] = None) -> str:
    try:
        os.makedirs(project_dir, exist_ok=True)
        out_path = os.path.join(project_dir, filename)
        with open(out_path, "w") as f:
            json.dump(metadata, f, indent=2)

        return f"✅ Metadata saved to `{out_path}`"
    except Exception as e:
        return f"❌ Failed to save metadata: {str(e)}"

# -------------------------------

def store_metadata_csv(metadata: dict, project_dir: str, filename: Optional[str] = None) -> str:
    try:
        os.makedirs(project_dir, exist_ok=True)

        # Build DataFrame
        df = pd.DataFrame(metadata)

        # Drop fully empty columns (NaN or blank)
        df.dropna(axis=1, how="all", inplace=True)
        df = df.loc[:, df.apply(lambda col: any(str(cell).strip() != "" for cell in col))]

        out_path = os.path.join(project_dir, filename)
        df.to_csv(out_path, index=False)

        return f"✅ Metadata saved to `{out_path}`"
    except Exception as e:
        return f"❌ Failed to save metadata: {str(e)}"
        return f"❌ Failed to fetch and save metadata: {str(e)}"

# -------------------------------

def csv_filename(filename: Optional[str] = None) -> str:
    # Set filename
    if filename:
        if not filename.endswith(".csv"):
            filename += ".csv"
    else:
        uid = str(uuid.uuid4())[:8]
        filename = f"metadata_{uid}.csv"
    return filename

# -------------------------------

def json_filename(filename: Optional[str] = None) -> str:
    if filename:
        if not filename.endswith(".json"):
            filename += ".json"
    else:
        uid = str(uuid.uuid4())[:8]
        filename = f"metadata_{uid}.json"
    return filename

# -------------------------------

def txt_filename(filename: Optional[str] = None) -> str:
    if filename:
        if not filename.endswith(".txt"):
            filename += ".txt"
    else:
        uid = str(uuid.uuid4())[:8]
        filename = f"metadata_{uid}.txt"
    return filename

# -------------------------------

gr_css = """
#chatbot_aes {
    padding: 20px;
    border-radius: 8px;
    font-size: 8px !important;
    line-height: 1.6;
}

.chat-font {
    padding: 20px;
    border-radius: 8px;
    font-size: 30px !important;
    line-height: 1.6;
}

#logs_terminal_panel .accordion-header {
    font-weight: bold;
    font-size: 12px;
}

#logs_terminal_panel { background-color: #ccece6; }

#logs_browser_panel { background-color: #e6f7ff; }

#logs_history_panel { background-color: #fbeee6; }

#log_viewer_aes textarea {
    background-color: #1435F3 !important;
    color: white !important;
    font-family: monospace;
    font-size: 11px !important;
}

#terminal_aes textarea {
    background-color: black !important;
    color: white !important;
    height: 300px !important;        /* Fixed height */
    overflow-y: auto !important;     /* Enable vertical scrolling */
    resize: none !important;         /* Optional: prevent manual resize */
    font-family: monospace;
    font-size: 14px !important;
    padding: 0.5em;
}

#pdf_upload_aes, 
#pdf_upload_aes * {
    color: black !important;
}

#pdf_upload_aes button {
    background-color: #FF5F1F !important;
    font-weight: bold;

#fixed_top_row .gr-box,
#fixed_top_row .gr-form,
#fixed_top_row .gr-textbox textarea,
#fixed_top_row {
    height: 350px !important; /* Or any fixed height */
    overflow: hidden;         /* Optional: prevents stretching */
    align-items: stretch !important; /* Force columns to fill height */
}
}
"""

# -------------------------------

chatbot_theme = gr.themes.Default(
    primary_hue="blue",
    secondary_hue="slate",
    neutral_hue="gray",
    radius_size="lg",
    text_size="md",  # Options: "sm", "md", "lg", "xl"
    font=["Inter", "sans-serif"]
).set(
    body_background_fill="#f6f9fb",
    body_text_color="#1a237e",
    button_primary_background_fill="#1976d2",
    button_primary_background_fill_hover="#1565c0",
    button_primary_text_color="#fff",
    button_secondary_background_fill="#e3eafc",
    button_secondary_text_color="#1976d2",
    input_background_fill="#fff",
    input_border_color="#90caf9",
    input_border_width="2px",
    block_shadow="0 4px 24px 0 rgba(30, 136, 229, 0.08)",
    loader_color="#1976d2",
    slider_color="#1976d2" ,
)

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

    # -------------------------------