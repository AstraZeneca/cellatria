# utils.py
# -------------------------------

import os
import re
import json
import uuid
import subprocess
import threading
import time
import tempfile
from datetime import datetime
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
    load_dotenv(dotenv_path=os.path.join(path_to_env, ".env"), override=True)
    provider = os.getenv("PROVIDER", "").strip()

    if provider == "Azure OpenAI":
        return AzureChatOpenAI(
            deployment_name=os.getenv("AZURE_OPENAI_DEPLOYMENT_NAME"),
            openai_api_version=os.getenv("AZURE_OPENAI_API_VERSION"),
            azure_endpoint=os.getenv("AZURE_OPENAI_ENDPOINT"),
            openai_api_key=os.getenv("AZURE_OPENAI_API_KEY"),
            temperature = os.getenv("AZURE_OPENAI_TEMPERATURE"),
            top_p = os.getenv("AZURE_OPENAI_TOP_P")   
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

#logs_terminal_panel { background-color: #ccece6; }   /* teal-green */
#logs_browser_panel  { background-color: #e6f7ff; }   /* light blue */
#logs_history_panel  { background-color: #fbeee6; }   /* light peach */
#agent_backend_panel { background-color: #fff9db; }   /* soft yellow */

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

#btn_export button, #btn_llm_meta button {
  background-color: #2d6cdf !important;
  color: #ffffff !important;
  border-color: #2d6cdf !important;
}
#btn_export button:hover, #btn_llm_meta button:hover {
  background-color: #1f56b3 !important;
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
    if args.batch_correction:
        method = args.batch_correction.strip().lower()
        if method not in {"harmony", "scvi"}:
            raise ValueError(
                "*** 🚨 Invalid --batch_correction value. "
                "Choose one of: 'harmony', 'scvi', or omit the flag."
            )

        if not args.batch_vars or not args.batch_vars.strip():
            raise ValueError(
                "*** 🚨 Batch correction selected, but --batch_vars was not provided. "
                "For Harmony, provide one or more columns (comma-separated). "
                "For scVI, provide exactly one column."
            )

        # Parse list (without accessing adata yet)
        btchvrs = [v.strip() for v in args.batch_vars.split(",") if v.strip()]

        if method == "harmony":
            if len(btchvrs) < 1:
                raise ValueError(
                    "*** 🚨 Harmony requires at least one batch variable in --batch_vars "
                    "(e.g., 'sample_id' or 'donor_id,sample_id')."
                )

        elif method == "scvi":
            if len(btchvrs) != 1:
                raise ValueError(
                    "*** 🚨 scVI supports exactly one batch variable in --batch_vars "
                    "(e.g., '--batch_vars donor_id')."
                )

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
# Configuration
LOG_PATH = "/tmp/cellatria_log.txt"

# Logging
def log_status(message: str):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(LOG_PATH, "a") as f:
        # f.write(f"[{timestamp}] {message}\n")
        f.write(f"{message}\n")

def read_log():
    try:
        with open(LOG_PATH, "r") as f:
            return f.read()
    except FileNotFoundError:
        return "📁 No logs yet."

# -------------------------------

# Terminal Session
class TerminalSession:
    def __init__(self):
        self.process = subprocess.Popen(
            ["bash"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1
        )
        self.lock = threading.Lock()
        self.latest_output = ""
        self.history = ""  

        # Start thread to continuously read output
        threading.Thread(target=self._read_output, daemon=True).start()

    def _read_output(self):
        for line in self.process.stdout:
            with self.lock:
                self.latest_output += line

    def run_command(self, cmd):
        if self.process.poll() is not None:
            return "❌ Shell has exited."

        with self.lock:
            self.latest_output = ""

        self.process.stdin.write(cmd + "\n")
        self.process.stdin.flush()
        
        time.sleep(0.3)
        
        with self.lock:
            result = self.latest_output.strip()
        return result

# -------------------------------
# Store transcript outside the session
terminal_transcript = ""

# Strip ANSI codes before displaying
ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')

def clean_ansi(text):
    return ansi_escape.sub('', text)

# -------------------------------
# Initialize the terminal session
terminal = TerminalSession()

# -------------------------------

# Define the function used in Gradio UI
def terminal_interface(cmd):
    global terminal_transcript
    raw_output = terminal.run_command(cmd)
    cleaned_output = clean_ansi(raw_output).strip()

    silent_commands = [
        "cd ",       # change directory
        "export ",   # set environment variable
        "unset ",    # unset environment variable
        "alias ",    # define an alias
        "unalias ",  # remove an alias
        "set ",      # set shell options or positional parameters
        "shopt ",    # shell options (bash-specific)
        "trap ",     # set signal handlers
        "pushd ",    # push directory onto stack
        "popd",      # pop directory from stack
        "exec ",     # replace shell with command
        "true",      # does nothing, returns 0
        "false",     # does nothing, returns 1
        ":",         # no-op command
        "clear",     # clears the screen
        "reset",     # resets the terminal
        "wait",      # waits for background jobs
        "disown",    # removes jobs from shell's job table
        "bg",        # resume job in background
        "fg",        # resume job in foreground
        "jobs",      # list background jobs
        "readonly ", # mark variables as read-only
        "declare ",  # variable declaration (bash)
        "typeset ",  # synonym for declare (ksh/zsh)
        "let ",      # arithmetic evaluation
        "source ",   # source a script
        ".",         # shorthand for source
    ]
    is_silent = any(cmd.strip().startswith(sc) for sc in silent_commands)
    
    if cleaned_output:
        terminal_transcript += f"\n$ {cmd}\n{cleaned_output}\n---"
    elif is_silent:
        terminal_transcript += f"\n$ {cmd}\n---"  # Silent: do not show fake "(no output)"
    else:
        # For commands like `ls` or `cat` returning true empty output
        terminal_transcript += f"\n$ {cmd}\n---"

    return terminal_transcript.strip(), ""  # keep input cleared

# -------------------------------
# File Browser Handlers
def fb_list_subdirs_and_files(path):
    try:
        items = os.listdir(path)
        dirs = sorted([item for item in items if os.path.isdir(os.path.join(path, item))])
        files = sorted([item for item in items if os.path.isfile(os.path.join(path, item))])
        return dirs, files, None
    except Exception as e:
        return [], [], f"❌ Error: {str(e)}"

def fb_get_dropdown_choices(path):
    dirs, _, _ = fb_list_subdirs_and_files(path)
    dirs = [f"📁 {d}" for d in dirs]
    parent = os.path.dirname(path.rstrip("/"))
    if parent and os.path.abspath(parent) != os.path.abspath(path):
        return [".. (Up)"] + dirs
    return dirs

def fb_initial_refresh(path):
    choices = fb_get_dropdown_choices(path)
    dirs, files, error = fb_list_subdirs_and_files(path)
    if error:
        file_display_val = f"<span style='color:red'>{error}</span>"
    else:
        file_display_val = "\n".join(f"📄 {f}" for f in files) or "No files in this directory."
    current_path = f"**Current Path:** `{path}`  \n**Folders:** {len(dirs)} | **Files:** {len(files)}"
    return (
        gr.Dropdown(choices=choices, value=None, interactive=bool(choices)),
        file_display_val,
        path,
        current_path
    )

def fb_navigate_subdir(subdir, base):
    if subdir and subdir.startswith("📁 "):
        subdir = subdir[2:]
    if subdir == ".. (Up)":
        new_path = os.path.dirname(base.rstrip("/"))
        if not new_path:
            new_path = base
    elif subdir:
        new_path = os.path.join(base, subdir)
    else:
        new_path = base
    choices = fb_get_dropdown_choices(new_path)
    dirs, files, error = fb_list_subdirs_and_files(new_path)
    if error:
        file_display_val = f"<span style='color:red'>{error}</span>"
    else:
        file_display_val = "\n".join(f"📄 {f}" for f in files) or "No files in this directory."
    current_path = f"**Current Path:** `{new_path}`  \n**Folders:** {len(dirs)} | **Files:** {len(files)}"
    return (
        gr.Dropdown(choices=choices, value=None, interactive=bool(choices)),
        file_display_val,
        new_path,
        current_path
    )
            
# -------------------------------

# Chat Export Feature
def export_chat_history(state_data):
    uid = uuid.uuid4().hex[:8]
    filename = f"chat_{uid}.json"
    filepath = os.path.join(tempfile.gettempdir(), filename)
    with open(filepath, "w") as f:
        json.dump(state_data, f, indent=2)
    return filepath

# -------------------------------

# Provide export function inside closure
def export_llm_metadata(llm_meta):
    uid = uuid.uuid4().hex[:8]
    filename = f"llm_{uid}.json"
    filepath = os.path.join(tempfile.gettempdir(), filename)
    with open(filepath, "w") as f:
        json.dump(llm_meta, f, indent=2)
    return filepath

# -------------------------------
# Initial message
initial_message = {
    "role": "assistant",
    "content": (
        "👋 Hello! I'm **cellAtria**, your assistant for analyzing single-cell RNA-seq (scRNA-seq) datasets.\n\n"
        "Here's what I can help you with:\n"
        "1. Extract structured metadata from scientific articles (PDF or URL).\n"
        "2. Store and organize metadata in structured project directories.\n"
        "3. Access public databases (currently support GEO) and fetch associated sample metadata.\n"
        "4. Download and organize scRNA-seq datasets.\n"        
        "5. Trigger CellExpress standardized single-cell data processing.\n\n"
        "To see a list of all available actions, type `'help'`.\n"
        "**Let's see how far we can fly together.** 🕊️\n\n"
        "To get started, let's set up your working directory. I can show you the current one or help create a new path. You can change your working directory anytime you wish.\n"
        f"📂 **Current Working Directory:** `{os.getcwd()}`"
    )
}

# -------------------------------

import json, time, uuid, psutil
from datetime import datetime
from langchain_core.messages import HumanMessage
import tiktoken

def count_tokens(text: str, model_name="gpt-4"):
    try:
        enc = tiktoken.encoding_for_model(model_name)
    except KeyError:
        enc = tiktoken.get_encoding("cl100k_base")
    return len(enc.encode(text))

def run_graph_with_metrics(graph, human_messages: list, config: dict):
    import psutil, time, json
    from datetime import datetime
    from langchain_core.messages import HumanMessage

    proc = psutil.Process()
    proc.cpu_percent(interval=None)  # Prime CPU meter
    start_ts = time.perf_counter()

    # Prepare input
    partial_state = {"messages": human_messages}
    full_text = "\n".join([m.content for m in human_messages])
    msg_bytes = len(full_text.encode("utf-8"))
    msg_tokens = count_tokens(full_text)

    steps = []
    try:
        print(f"\n🧾 Processing message block: {full_text[:80]}...")
        for step in graph.stream(partial_state, config=config):
            steps.append(step)
            print(step)  # ✅ Real-time streaming output

        status = "success"
        output_text = "".join(str(s) for s in steps)
        output_bytes = len(output_text.encode("utf-8"))
        output_tokens = count_tokens(output_text)

    except Exception as e:
        print(f"🚨 Error: {e}")
        status = "error"
        output_bytes = 0
        output_tokens = 0

    duration = round(time.perf_counter() - start_ts, 2)
    message_metrics = {
        "input_size_bytes": msg_bytes,
        "input_token_count": msg_tokens,
        "output_size_bytes": output_bytes,
        "output_token_count": output_tokens,
        "duration_sec": duration,
        "status": status
    }

    return steps, message_metrics

# -------------------------------

def get_llm_metadata(llm):
    """
    Return a provider-agnostic dict describing the LLM that was actually constructed.
    Tries to be robust across Azure OpenAI, OpenAI, Anthropic, Google/Gemini, and HF/local.
    """
    meta = {
        "provider": None,
        "model": None,
        "version": None,
        "temperature": None,
        "top_p": None,
        "raw_class": llm.__class__.__name__,
        "raw_module": llm.__class__.__module__,
    }

    # ---- Azure OpenAI / LangChain AzureChatOpenAI ----
    if hasattr(llm, "deployment_name") or hasattr(llm, "azure_endpoint"):
        meta["provider"] = "azure-openai"
        # Azure: deployment_name is effectively the model we care about
        meta["model"] = getattr(llm, "deployment_name", None)
        meta["version"] = getattr(llm, "openai_api_version", None)
        meta["temperature"] = getattr(llm, "temperature", None)
        meta["top_p"] = getattr(llm, "top_p", None)
        return meta

    # ---- OpenAI (non-Azure) ----
    # e.g. langchain_openai.ChatOpenAI
    if hasattr(llm, "model_name"):
        meta["provider"] = "openai"
        meta["model"] = getattr(llm, "model_name", None)
        meta["temperature"] = getattr(llm, "temperature", None)
        meta["top_p"] = getattr(llm, "top_p", None)
        # OpenAI rarely exposes a separate "version" here, so leave None
        return meta

    # ---- Anthropic ----
    # e.g. langchain_anthropic.ChatAnthropic
    if llm.__class__.__module__.startswith("langchain_anthropic") or llm.__class__.__name__.lower().startswith("chatanthropic"):
        meta["provider"] = "anthropic"
        # Anthropic models are usually specified as .model or .model_name
        model = getattr(llm, "model", None) or getattr(llm, "model_name", None)
        meta["model"] = model
        meta["temperature"] = getattr(llm, "temperature", None)
        # Anthropic may not have top_p in the same way
        return meta

    # ---- Google / Gemini (very wrapper-dependent) ----
    # if the user returned a google.generativeai model
    if llm.__class__.__module__.startswith("google") or "genai" in llm.__class__.__module__.lower():
        meta["provider"] = "google-genai"
        # many Gemini wrappers have .model_name or .model
        model = getattr(llm, "model_name", None) or getattr(llm, "model", None)
        meta["model"] = model
        return meta

    # ---- Hugging Face / local pipeline ----
    # most HF pipelines have .model and .tokenizer attrs
    if hasattr(llm, "model") and hasattr(llm, "__call__"):
        meta["provider"] = "local"
        # try to get model id
        try:
            meta["model"] = getattr(llm.model, "name_or_path", None)
        except Exception:
            meta["model"] = None
        return meta

    # ---- fallback ----
    meta["provider"] = "unknown"
    return meta

# -------------------------------