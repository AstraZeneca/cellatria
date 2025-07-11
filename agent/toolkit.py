# toolkit.py
# -------------------------------

from langchain_core.tools import tool
from langchain.agents import Tool
import requests
from bs4 import BeautifulSoup
from bs4.element import Comment
from typing import List, Dict, Any, TypedDict, Literal, Optional, Union
from pydantic import BaseModel, Field, ValidationError
import fitz
import os
import json
import csv
import re
import signal
import GEOparse
import pandas as pd
import requests
import shutil
import psutil
import glob
import uuid
import argparse
import subprocess
from utils import (csv_filename, json_filename, txt_filename, base_path, checks_args)

# -------------------------------
# Extracts structured metadata from a scientific article at a given URL.
# -------------------------------

class ArticleURL(BaseModel):
    url: str

@tool(args_schema=ArticleURL)
def fetch_article_metadata_url(url: str) -> List[str]:
    """
    Fetch visible text from a scientific article. If blocked, ask user to upload PDF or paste text manually.
    After metadata extraction offer data storage in JSON or TXT format.
    """
    def visible_texts(element):
        return element.parent.name not in ["style", "script", "head", "meta", "[document]"] and not isinstance(element, Comment)

    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, "html.parser")
        texts = soup.findAll(string=True)
        visible = filter(visible_texts, texts)
        lines = [t.strip() for t in visible if t.strip()]
        return lines

    except requests.exceptions.HTTPError as http_err:
        if response.status_code == 403:
            return [
                "This site does not allow automated access (HTTP 403).",
                "Please upload the PDF."
            ]
        return [f"❌ HTTP error {response.status_code}: {str(http_err)}"]

    except Exception as e:
        return [f"❌ Failed to fetch article: {str(e)}"]

# -------------------------------
# Holds article metadata for review and edits before saving.
# -------------------------------

# Define schema
class ArticleMetadataSchema(BaseModel):
    Project: str
    Title: str
    Authors: str
    Tissue: str
    Species: str
    Disease: str
    Data_Modality: str
    Data_Availability: str
    Publication_Date: str
    Publisher: str
    Conflicts_of_Interest: str

# Cache for metadata
art_metadata_cache = {}

@tool(args_schema=ArticleMetadataSchema)
def refine_article_metadata(**kwargs) -> str:
    """
    Temporarily registers metadata for review and approval.
    This does not persist the data — it only prepares it for final storage if confirmed.
    """
    art_metadata_cache.clear()
    art_metadata_cache.update(kwargs)  # `kwargs` is already a dict
    return "✅ Metadata prepared for review. Awaiting user confirmation to store."

# -------------------------------
# Saves extracted article metadata to a file (JSON or TXT).
# -------------------------------

class MetadataStorageInput(BaseModel):
    art_metadata_cache: ArticleMetadataSchema
    path: str
    filename: Optional[str] = None
    suffix: Optional[str] = "json"  # default to JSON

@tool(args_schema=MetadataStorageInput)
def store_article_metadata_file(art_metadata_cache: ArticleMetadataSchema, path: str, filename: Optional[str] = None, suffix: str = "json") -> str:
    """
    Store art_metadata_cache to file in the specified format and directory.
    Supported formats: JSON or TXT.
    If `path` is provided, save output.
    If `filename` is not provided, generate a UUID-based filename.
    If `suffix` is 'json', save in JSON format.
    If `suffix` is 'txt', save in TXT format.
    """
    
    if not art_metadata_cache:
        return "⚠️ No art_metadata_cache available in memory to store."

    try:
        os.makedirs(path, exist_ok=True)
        if suffix == "json":
            filename = json_filename(filename)
            out_path = os.path.join(path, filename)
            with open(out_path, "w") as f:
                json.dump(art_metadata_cache.dict(), f, indent=2)

        elif suffix == "txt":
            filename = txt_filename(filename)
            out_path = os.path.join(path, filename)
            with open(out_path, "w") as f:
                for key, val in art_metadata_cache.dict().items():
                    f.write(f"{key}: {val}\n")

        return f"✅ Metadata saved to `{out_path}`"

    except Exception as e:
        return f"❌ Failed to save cached metadata: {str(e)}"

# -------------------------------  
# Extracts structured metadata from an uploaded scientific article PDF
# -------------------------------   

class ArticlePDF(BaseModel):
    path: str

@tool(args_schema=ArticlePDF)
def fetch_article_metadata_pdf(path: str) -> List[str]:
    """
    Fetch visible text from a scientific PDF article. returning paragraphs for LLM to process.
    After metadata extraction offer data storage in JSON or TXT format.
    """
    try:
        doc = fitz.open(path)
        lines = []
        for page in doc:
            text = page.get_text()
            if text.strip():
                lines.extend([line.strip() for line in text.split("\n") if line.strip()])
        doc.close()

        return lines or ["⚠️ PDF parsed, but no readable text found."]
    
    except Exception as e:
        return [f"❌ Failed to extract text from PDF: {str(e)}"]

# -------------------------------   
# Retrieves metadata for a GEO accession (e.g., GSE123456).
# -------------------------------   

class GEOInput(BaseModel):
    gse_id: str

@tool(args_schema=GEOInput)
def fetch_geo_metadata(gse_id: str) -> dict:
    """
    Fetch GEO metadata including title, organism, GSM sample descriptions, and linked BioSample and SRA (SRX) accessions from Relations.
    After metadata extraction offer data storage in JSON or CSV format.
    """

    try:
        geo_url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}"
        response = requests.get(geo_url, timeout=10)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, "html.parser")

        title_cell = soup.find("td", string="Title")
        organism_cell = soup.find("td", string="Organism")
        geo_title = title_cell.find_next_sibling("td").get_text(strip=True) if title_cell else "N/A"
        geo_organism = organism_cell.find_next_sibling("td").get_text(strip=True) if organism_cell else "N/A"

        # Extract GSM sample metadata
        sample_rows = []
        for row in soup.find_all("tr"):
            cols = row.find_all("td")
            if len(cols) >= 2:
                gsm_id = cols[0].get_text(strip=True)
                gsm_desc = cols[1].get_text(strip=True)
                if re.match(r"GSM\d{6,}", gsm_id):
                    gsm_url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm_id}"
                    gsm_response = requests.get(gsm_url, timeout=10)
                    gsm_soup = BeautifulSoup(gsm_response.text, "html.parser")

                    biosample = None
                    sra_ids = []

                    relations_header = gsm_soup.find("strong", string="Relations")
                    if relations_header:
                        relations_table = relations_header.find_parent("table")
                        if relations_table:
                            for row in relations_table.find_all("tr"):
                                cols = row.find_all("td")
                                if len(cols) == 2:
                                    label = cols[0].get_text(strip=True)
                                    link = cols[1].find("a")
                                    value = link.get_text(strip=True) if link else cols[1].get_text(strip=True)
                                    if label == "BioSample":
                                        biosample = value
                                    elif label == "SRA":
                                        sra_ids.append(value)

                    sample_rows.append({
                        "sample": gsm_id,
                        "description": gsm_desc,
                        "biosample_id": biosample if biosample is not None else "Unavailable",
                        "sra_ids": ", ".join(sra_ids) if sra_ids else "Unavailable"
                    })

        metadata = {
                "geo_title": geo_title,
                "geo_organism": geo_organism,
                "geo_link": geo_url,
                "geo_samples": sample_rows
        }

        return metadata

    except Exception as e:
        return {"error": str(e)}

# ------------------------------- 
# Holds GEO metadata for review and edits before saving
# -------------------------------   

class GeoSampleMetadata(BaseModel):
    sample: str
    description: str
    biosample_id: Optional[str]
    sra_ids: Optional[str]  # comma-separated string

class GeoMetadataSchema(BaseModel):
    geo_title: str
    geo_organism: str
    geo_link: str
    geo_samples: List[GeoSampleMetadata]
    
geo_metadata_cache = {}
@tool(args_schema=GeoMetadataSchema)
def refine_geo_metadata(**kwargs) -> str:
    """
    Temporarily registers metadata for review and approval.
    This does not persist the data — it only prepares it for final storage if confirmed.
    """
    geo_metadata_cache.clear()
    geo_metadata_cache.update(kwargs)  # `kwargs` is already a dict
    return "✅ Metadata prepared for review. Awaiting user confirmation to store."

# -------------------------------  
# Saves extracted GEO metadata to a file (JSON or CSV).
# -------------------------------   

class MetadataStorageInput(BaseModel):
    geo_metadata_cache: GeoMetadataSchema
    path: str
    filename: Optional[str] = None
    suffix: Optional[str] = "json"  # default to JSON

@tool(args_schema=MetadataStorageInput)
def store_geo_metadata_file(geo_metadata_cache: GeoMetadataSchema, path: str, filename: Optional[str] = None, suffix: str = "json") -> str:
    """
    Store geo_metadata_cache to file in the specified format and directory.
    Supported formats: JSON or CSV.
    If `path` is provided, save output.
    If `filename` is not provided, generate a UUID-based filename.
    If `suffix` is 'json', save in JSON format.
    If `suffix` is 'txt', save in CSV format.
    """
    
    if not geo_metadata_cache:
        return "⚠️ No geo_metadata_cache available in memory to store."

    try:
        os.makedirs(path, exist_ok=True)
        if suffix == "json":
            filename = json_filename(filename)
            out_path = os.path.join(path, filename)
            with open(out_path, "w") as f:
                json.dump(geo_metadata_cache.dict(), f, indent=2)

        elif suffix == "csv":
            filename = csv_filename(filename)
            out_path = os.path.join(path, filename)            
            with open(out_path, "w") as f:
                samples = [s.dict() for s in geo_metadata_cache.geo_samples]
                writer = csv.DictWriter(f, fieldnames=samples[0].keys())
                writer.writeheader()
                writer.writerows(samples)

        return f"✅ Metadata saved to `{out_path}`"

    except Exception as e:
        return f"❌ Failed to save cached metadata: {str(e)}"

# -------------------------------
# Downloads all supplementary files for a GEO dataset.
# -------------------------------

class DownloadDatasetArgs(BaseModel):
    gse_id: str = Field(..., description="GEO accession ID (e.g., GSE204716)")
    path: str = Field(..., description="Base directory to save dataset files")

@tool(args_schema=DownloadDatasetArgs)
def download_geo(gse_id: str, path: str) -> str:
    """
    Downloads GEO dataset using GEOparse and saves each sample's supplementary files
    in a structured directory under a peroject directory.
    """
    try:
        # Step 1: Download GSE object
        gse = GEOparse.get_GEO(geo=gse_id, destdir=path, how="full")

        # Step 2: Track the original working directory
        original_dir = os.getcwd()

        # Step 3: Download each GSM's supplementary files into its own folder
        for i, (gsm_id, gsm) in enumerate(gse.gsms.items(), start=1):
            sample_dir = os.path.join(path, gsm_id)
            os.makedirs(sample_dir, exist_ok=True)
            os.chdir(sample_dir)

            print(f"{i}. ⬇️  Downloading for {gsm_id}...")
            gsm.download_supplementary_files()

            # Flatten any extracted subdirectories
            for item in os.listdir('.'):
                if os.path.isdir(item):
                    subdir = os.path.join(sample_dir, item)
                    for file in os.listdir(subdir):
                        shutil.move(os.path.join(subdir, file), os.path.join(sample_dir, file))
                    os.rmdir(subdir)

            os.chdir(original_dir)

        return f"✅ Downloaded {len(gse.gsms)} samples to `{path}`"

    except Exception as e:
        return f"❌ Download failed: {str(e)}"

# -------------------------------
# Downloads supplementary files for a list of GSM accessions.
# -------------------------------

class DownloadSamplesArgs(BaseModel):
    gsm_ids: list[str] = Field(..., description="List of GSM IDs, e.g., ['GSM1234567', 'GSM1234568']")
    path: str = Field(..., description="Base directory to save all sample folders")

@tool(args_schema=DownloadSamplesArgs)
def download_gsm(gsm_ids: list[str], path: str) -> str:
    """
    Downloads supplementary files for multiple GSM accessions.
    Each GSM is stored in its own subfolder under <path>/<gsm_id>.
    Removes the .txt descriptor after download.
    """
    
    results = []
    original_dir = os.getcwd()          # remember where we started

    for gsm_id in gsm_ids:
        try:
            sample_dir = os.path.join(path, gsm_id)
            os.makedirs(sample_dir, exist_ok=True)
            os.chdir(sample_dir)        # <-- work inside the GSM folder

            gsm = GEOparse.get_GEO(geo=gsm_id, destdir=sample_dir, how="full")
            gsm.download_supplementary_files()   # now lands inside sample_dir

            # ── flatten everything under sample_dir
            for root, _, files in os.walk(sample_dir):
                for file in files:
                    src = os.path.join(root, file)
                    dst = os.path.join(sample_dir, file)
                    if src != dst:
                        shutil.move(src, dst)

            # remove empty dirs
            for root, dirs, _ in os.walk(sample_dir, topdown=False):
                for d in dirs:
                    dir_path = os.path.join(root, d)
                    if not os.listdir(dir_path):
                        os.rmdir(dir_path)

            # remove descriptor txt
            txt_file = os.path.join(sample_dir, f"{gsm_id}.txt")
            if os.path.exists(txt_file):
                os.remove(txt_file)

            results.append(f"✅ {gsm_id}")
        except Exception as e:
            results.append(f"❌ {gsm_id}: {e}")
        finally:
            os.chdir(original_dir)      # return to parent before next GSM

    return "\n".join(results)

# -------------------------------
# Downloads a file from a direct URL to your project folder.
# -------------------------------

class DownloadFilesToPathArgs(BaseModel):
    urls: List[str] = Field(..., description="List of direct URLs to files to download.")
    path: str = Field(..., description="Target directory to save all downloaded files.")

@tool(args_schema=DownloadFilesToPathArgs)
def download_file(urls: List[str], path: str) -> str:
    """
    Downloads multiple files from given URLs and saves them to a shared local directory.
    """
    os.makedirs(path, exist_ok=True)
    results = []

    for url in urls:
        try:
            filename = os.path.basename(url.split("?")[0])
            full_path = os.path.join(path, filename)

            with requests.get(url, stream=True, timeout=20) as r:
                r.raise_for_status()
                with open(full_path, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)

            results.append(f"✅ `{filename}` saved to `{full_path}`")
        except Exception as e:
            results.append(f"❌ Failed to download `{url}` → {str(e)}")

    return "\n".join(results)

# -------------------------------  
# Creates a new directory at the specified path.
# -------------------------------       

class MakeDirectoryArgs(BaseModel):
    path: str

@tool(args_schema=MakeDirectoryArgs)
def make_directory(path: str) -> str:
    """
    Creates a directory at the specified path if it doesn't already exist.
    """
    try:
        os.makedirs(path, exist_ok=True)
        return f"✅ Created directory: {path}"
    except Exception as e:
        return f"❌ Failed to create directory: {str(e)}"

# -------------------------------   
# Moves a file from one location to another.
# -------------------------------   

class MoveFileArgs(BaseModel):
    src: str
    dest: str

@tool(args_schema=MoveFileArgs)
def move_file(src: str, dest: str) -> str:
    """
    Moves a file from source to destination.
    """
    try:
        shutil.move(src, dest)
        return f"✅ Moved file from {src} to {dest}"
    except Exception as e:
        return f"❌ Failed to move file: {str(e)}"

# -------------------------------  
# Renames a file or moves it to a new name/location.
# -------------------------------  

class RenameFileArgs(BaseModel):
    src: str = Field(..., description="Full path to the existing file to rename.")
    dest: str = Field(..., description="Full path with the new file name.")

@tool(args_schema=RenameFileArgs)
def rename_file(src: str, dest: str) -> str:
    """
    Renames a file from `src` to `dest`.
    """
    try:
        if not os.path.isfile(src):
            return f"❌ Source file does not exist: {src}"
        os.rename(src, dest)
        return f"✅ Renamed file from `{src}` to `{dest}`"
    except Exception as e:
        return f"❌ Failed to rename file: {str(e)}" 

# -------------------------------   
# Shows the contents of a directory in tree format.
# -------------------------------  

class ListDirectoryArgs(BaseModel):
    path: str

def generate_tree(path: str, prefix: str = "") -> List[str]:
    entries = sorted(os.listdir(path))
    lines = []

    for index, entry in enumerate(entries):
        full_path = os.path.join(path, entry)
        connector = "└── " if index == len(entries) - 1 else "├── "
        lines.append(prefix + connector + entry)
        if os.path.isdir(full_path):
            extension = "    " if index == len(entries) - 1 else "│   "
            lines.extend(generate_tree(full_path, prefix + extension))

    return lines

@tool(args_schema=ListDirectoryArgs)
def list_directory(path: str) -> List[str]:
    """
    Returns a tree-style listing of the contents of a directory.
    """
    try:
        if not os.path.exists(path):
            return [f"❌ Path does not exist: {path}"]
        if not os.path.isdir(path):
            return [f"❌ Path is not a directory: {path}"]

        return generate_tree(path)
    except Exception as e:
        return [f"❌ Failed to list directory: {str(e)}"]

# -------------------------------   
# Deletes a file or directory (recursively if needed).
# -------------------------------   

class RemovePathArgs(BaseModel):
    path: str

@tool(args_schema=RemovePathArgs)
def remove_file_or_dir(path: str) -> str:
    """
    Deletes a file or directory at the given path.
    """
    try:
        if os.path.isdir(path):
            shutil.rmtree(path)
        else:
            os.remove(path)
        return f"🗑️ Removed: {path}"
    except Exception as e:
        return f"❌ Failed to remove: {str(e)}"

# -------------------------------  
# Shows the size of a file in megabytes.
# -------------------------------  

class FileSizeArgs(BaseModel):
    path: str

@tool(args_schema=FileSizeArgs)
def get_file_size(path: str) -> str:
    """
    Returns the size of the specified file in MB.
    """
    try:
        size_mb = os.path.getsize(path) / (1024 * 1024)
        return f"{size_mb:.2f} MB"
    except Exception as e:
        return f"❌ Failed to get file size: {str(e)}"

# -------------------------------  
# Creates and saves a CSV file from tabular data.
# -------------------------------  

# Define the schema
class CustomCsvArgs(BaseModel):
    data: List[Dict[str, str]] = Field(..., description="A list of dictionaries where each dictionary is a row in the CSV and keys are column names.")
    path: Optional[str] = None
    filename: Optional[str] = None

@tool(args_schema=CustomCsvArgs)
def create_custom_csv(data: List[Dict[str, str]], path: Optional[str] = None, filename: Optional[str] = None) -> str:
    """
    Creates a customizable metadata file.
    Accepts a list of dictionaries where each dictionary is a row and keys are column names.
    If `filename` is not provided, generate a UUID-based filename.
    Saves the CSV inside the provided project directory.
    """
    try:
        os.makedirs(path, exist_ok=True)

        # Set filename
        filename = csv_filename(filename)
        output_path = os.path.join(path, filename)

        df = pd.DataFrame(data)
        df.to_csv(output_path, index=False)

        return f"✅ Custom sample metadata saved to: {output_path}"
    except Exception as e:
        return f"❌ Failed to save custom sample metadata: {str(e)}"

# -------------------------------
# Creates and saves a JSON file from structured data.
# -------------------------------

class CustomJsonArgs(BaseModel):
    data: List[Dict[str, str]] = Field(..., description="A list of dictionaries to be saved as JSON.")
    path: Optional[str] = None
    filename: Optional[str] = None

@tool(args_schema=CustomJsonArgs)
def create_custom_json(data: List[Dict[str, str]], path: Optional[str] = None, filename: Optional[str] = None) -> str:
    """
    Creates a customizable metadata file in JSON format.
    """
    try:
        os.makedirs(path, exist_ok=True)
        filename = json_filename(filename)
        output_path = os.path.join(path, filename)
        with open(output_path, "w") as f:
            json.dump(data, f, indent=2)
        return f"✅ Custom sample metadata saved to: {output_path}"
    except Exception as e:
        return f"❌ Failed to save custom sample metadata: {str(e)}"

# -------------------------------
# Creates and saves a plain text file from a string.
# -------------------------------

class CustomTxtArgs(BaseModel):
    data: str = Field(..., description="The text content to be saved in the TXT file.")
    path: Optional[str] = None
    filename: Optional[str] = None

@tool(args_schema=CustomTxtArgs)
def create_custom_txt(data: str, path: Optional[str] = None, filename: Optional[str] = None) -> str:
    """
    Creates a customizable metadata file in TXT format.
    """
    try:
        os.makedirs(path, exist_ok=True)
        filename = txt_filename(filename)
        output_path = os.path.join(path, filename)
        with open(output_path, "w") as f:
            f.write(data)
        return f"✅ Custom sample metadata saved to: {output_path}"
    except Exception as e:
        return f"❌ Failed to save custom sample metadata: {str(e)}"

# -------------------------------
# Displays the contents of a CSV, JSON, or TXT file.
# -------------------------------

class FileCache(BaseModel):
    raw: Union[pd.DataFrame, dict, str]
    instructions: Optional[str] = None
    class Config:
        arbitrary_types_allowed = True

class InspectFileArgs(BaseModel):
    path: str = Field(..., description="Directory containing the metadata file.")
    filename: str = Field(..., description="Name of the file to inspect.")
    suffix: Optional[str] = Field(None, description="Optional file extension: 'csv', 'json', or 'txt'.")

@tool(args_schema=InspectFileArgs)
def inspect_file(path: str, filename: str, suffix: Optional[str] = None) -> FileCache:
    """
    Reads a file from disk and returns a FileCache object containing its raw contents.
    Supports CSV, JSON, and TXT formats.
    """
    try:
        full_path = os.path.join(path, filename)
        if not os.path.isfile(full_path):
            raise FileNotFoundError(f"File not found: {full_path}")

        ext = (suffix or os.path.splitext(full_path)[1][1:]).lower()

        if ext == "csv":
            df = pd.read_csv(full_path)
            return FileCache(raw=df)
        elif ext == "json":
            with open(full_path, "r") as f:
                return FileCache(raw=json.load(f))
        elif ext == "txt":
            with open(full_path, "r") as f:
                return FileCache(raw=f.read())
        else:
            raise ValueError(f"Unsupported file extension: {ext}. Use csv, json, or txt.")
    except Exception as e:
        raise RuntimeError(f"❌ Failed to read file: {str(e)}")

# -------------------------------  
# Displays your current working directory.
# -------------------------------   

class gwdArgs(BaseModel):
    pass

@tool(args_schema=gwdArgs)
def get_working_directory() -> str:
    """
    Returns the current working directory.
    """
    return os.getcwd()  

# -------------------------------    
# Sets or changes your working directory.
# -------------------------------     

wd_cache = {"path": None}

class swdArgs(BaseModel):
    path: str = Field(..., description="Absolute path to set as the working directory.")

@tool(args_schema=swdArgs)
def set_working_directory(path: str) -> str:
    """
    Sets the working directory to the specified path.
    Creates the directory if it doesn't exist.
    Caches the path globally for other tools to access.
    """
    try:
        os.makedirs(path, exist_ok=True)
        os.chdir(path)
        wd_cache["path"] = os.getcwd()
        return f"✅ Working directory set to: {wd_cache['path']}"
    except Exception as e:
        return f"❌ Failed to set working directory: {str(e)}"

# -------------------------------
# Standardizes 10X Genomics filenames in a directory.
# -------------------------------

class Fix3gz10xArgs(BaseModel):
    path: str

@tool(args_schema=Fix3gz10xArgs)
def fix_10x_file_format(path: str) -> str:
    """
    Recursively checks and fixes subdirectories under `path` to ensure each contains
    10x Genomics files: `matrix.mtx.gz`, `features.tsv.gz`, and `barcodes.tsv.gz`.
    Files ending in these names will be renamed accordingly.
    """
    expected_files = {
        "matrix.mtx.gz",
        "features.tsv.gz",
        "barcodes.tsv.gz"
    }

    for subdir, dirs, files in os.walk(path):
        renamed = []
        missing = []
        found_suffixes = {f for f in files if any(f.endswith(suffix) for suffix in expected_files)}

        for expected in expected_files:
            matched = [f for f in files if f.endswith(expected)]
            if not matched:
                missing.append(f"- Missing `{expected}`")
            else:
                actual = matched[0]
                if actual != expected:
                    old_path = os.path.join(subdir, actual)
                    new_path = os.path.join(subdir, expected)
                    try:
                        os.rename(old_path, new_path)
                        renamed.append(f"- Renamed `{actual}` → `{expected}`")
                    except Exception as e:
                        renamed.append(f"- ❌ Failed to rename `{actual}`: {e}")

        if not found_suffixes:
            continue

    return f"✅ Validation and renaming complete under `{path}`."

# -------------------------------
# Shows an overview of CellExpress pipeline capabilities.
# -------------------------------

class cellExpressDesc(BaseModel):
    pass

@tool(args_schema=cellExpressDesc)
def get_cellexpress_info() -> str:
    """
    Use this for an overview of the CellExpress pipeline capabilities.
    """
    with open(os.path.abspath(os.path.join(base_path, "..", "cellexpress", "README.md")), "r") as f:
        return f.read()

# -------------------------------
# Sets or updates a CellExpress pipeline argument.
# -------------------------------

class CellExpressCache:
    def __init__(self):
        self._store: Dict[str, Any] = {}

    def set_args(self, new_args: Dict[str, Any]):
        self._store.update(new_args)

    def get_args(self) -> Dict[str, Any]:
        return self._store.copy()

    def clear(self):
        self._store.clear()

    def missing_required_fields(self, required_fields=["input", "project", "species", "tissue", "disease"]):
        return [field for field in required_fields if field not in self._store or self._store[field] in [None, ""]]

cellexpress_cache = CellExpressCache()

class ConfigArgInput(BaseModel):
    args: Dict[str, str] = Field(..., description="Dictionary of argument names and values to set/update.")

@tool(args_schema=ConfigArgInput)
def configure_cellexpress(args: Dict[str, str]) -> str:
    """
    Store or update multiple CellExpress arguments. Validates against full argument schema.
    Accept single argument prompts as well as multiple arguments.
    """
    cellexpress_cache.set_args(args)

    # Try full schema validation
    try:
        all_args = cellexpress_cache.get_args()
        validated = CellExpressArgs(**all_args)
        return "✅ All required arguments set. Ready to run CellExpress."
    except ValidationError as e:
        missing = cellexpress_cache.missing_required_fields()
        return f"Still missing or invalid: {missing}\nValidation errors:\n{e}"

# -------------------------------
# Previews the current CellExpress CLI configuration.
# -------------------------------

class EmptyInput(BaseModel):
    pass

@tool(args_schema=EmptyInput)
def preview_cellexpress_config() -> str:
    """
    Show current configuration of CellExpress arguments after validation.
    """
    args = cellexpress_cache.get_args()
    if not args:
        return "⚠️ No arguments configured yet."

    try:
        validated = CellExpressArgs(**args)
    except ValidationError as e:
        return f"⚠️ Invalid config:\n{e}"

    formatted = "\n".join([f"--{k} {v}" for k, v in validated.dict().items() if v is not None])
    return f"📋 Current CellExpress CLI configuration:\n```bash\n{formatted}\n```"

# -------------------------------
# Clears all cached CellExpress configuration arguments.
# -------------------------------

class EmptyInput(BaseModel):
    pass

@tool(args_schema=EmptyInput)
def reset_cellexpress_config() -> str:
    """
    Clears all cached CellExpress arguments.
    """
    cellexpress_cache.clear()
    return "🧹 Cleared all cached CellExpress arguments."

# -------------------------------
# Validates the current CellExpress configuration.
# -------------------------------

class EmptyInput(BaseModel):
    pass

@tool(args_schema=EmptyInput)
def validate_cellexpress_config() -> str:
    """
    Performs full semantic validation of CellExpress arguments using the main pipeline's internal checks.
    """
    try:
        args = CellExpressArgs(**cellexpress_cache.get_args())
    except ValidationError as e:
        return f"❌ Pydantic schema validation failed:\n{e}"

    try:
        namespace_args = argparse.Namespace(**args.dict())
        checks_args(namespace_args)
        return "✅ All arguments passed semantic validation. Ready to run CellExpress."
    except Exception as e:
        return f"❌ Semantic validation failed:\n{e}"

 # -------------------------------       

class CellExpressArgs(BaseModel):
    # General arguments
    input: str = Field(..., description="Input directory containing the sample folders and metadata.csv.")
    project: str = Field(..., description="Project name.")
    species: str = Field(..., description="Species: 'hs' for human, 'mm' for mouse.")
    tissue: str = Field(..., description="Tissue name.")
    disease: str = Field(..., description="Disease name.")

    # Quality control arguments
    min_umi_per_cell: Optional[int] = Field(750, description="Minimum UMI counts required per cell to pass filtering (default: 750).")
    max_umi_per_cell: Optional[int] = Field(None, description="Maximum UMI counts per cell (default: None).") 
    min_genes_per_cell: Optional[int] = Field(250, description="Minimum genes per cell (default: 250).") 
    max_genes_per_cell: Optional[int] = Field(None, description="Maximum genes per cell (default: None).") 
    min_cell: Optional[int] = Field(3, description="Genes detected in at least this many cells (default: 3).") 
    max_mt_percent: Optional[float] = Field(15.0, description="Max mitochondrial gene percentage (default: 15).") 
    doublet_method: Optional[str] = Field(None, description="Method for doublet identification. Options: 'scrublet' (default: None).") 
    scrublet_cutoff: Optional[float] = Field(0.25, description="Maximum allowed doublet score. Applies only if doublet_method='scrublet'. (default: 0.25).") 

    # Analysis arguments
    norm_target_sum: Optional[float] = Field(1e4, description="Target sum for total counts per cell during normalization (default: 1e4).") 
    n_top_genes: Optional[int] = Field(2000, description="Number of top highly variable genes to retain (default: 2000).") 
    regress_out: Optional[str] = Field("no", description="Regress out total counts and mitochondrial percentage before scaling (default: 'no').")
    scale_max_value: Optional[float] = Field(10, description="Clip (truncate) to this value after scaling (default: 10).") 
    n_pcs: Optional[int] = Field(30, description="Number of principal components to compute (default: 30).") 
    batch_correction: Optional[str] = Field(None, description="Method for batch correction: Options: 'harmony' (default: None).")
    batch_vars: Optional[str] = Field(None, description="Column(s) in adata.obs to use for Harmony batch correction (default: None).")
    n_neighbors: Optional[int] = Field(15, description="Number of neighbors for k-NN graph (default: 15).") 
    resolution: Optional[float] = Field(0.6, description="Resolution parameter for Leiden clustering (default: 0.6).")
    compute_tsne: Optional[str] = Field("no", description="Option to compute t-SNE embedding in addition to UMAP (default: 'no').")

    # Annotation
    annotation_method: Optional[str] = Field(None, description="Annotation method(s) to use: 'scimilarity', 'celltypist', or both (default: 'None').")
    sci_model_path: Optional[str] = Field(None, description="Path to SCimilarity model directory, if used (default: 'None').") 
    cty_model_path: Optional[str] = Field(None, description="Path to directory containing CellTypist `.pkl` models (default: 'None').") 
    cty_model_name: Optional[str] = Field(None, description="Name of the CellTypist model to use for annotation - without '.pkl' extension (default: None).")

    # Differentially expressed analysis
    pval_threshold: Optional[float] = Field(0.05, description="P-value threshold for filtering marker genes (default: 0.05).")
    logfc_threshold: Optional[float] = Field(0.25, description="Log fold-change threshold for filtering marker genes (default: 0.25).")
    dea_method: Optional[str] = Field("wilcoxon", description="Statistical test used for marker genes identification (default: 'wilcoxon').")
    top_n_deg_leidn: Optional[int] = Field(100, description="Number of top marker genes to return per 'leiden' groups (default: 100).")
    top_n_deg_scim: Optional[int] = Field(100, description="Number of top marker genes to return per 'scimilarity' annotated cells (default: 100).")
    top_n_deg_cltpst: Optional[int] = Field(100, description="Number of top marker genes to return per 'celltypist' annotated cells (default: 100).")
    pts_threshold: Optional[float] = Field(0.1, description="Minimum fraction of cells expressing a gene for it to be considered a marker (default: 0.1).")

    # Documentation/metadata
    doc_url: Optional[str] = Field(None, description="URL to the documentation associated with the data (default: None).")
    data_url: Optional[str] = Field(None, description="URL to the publicly available dataset or repository (default: None).")

    # QC-only mode
    only_qc: Optional[str] = Field("no", description="Set to 'yes' to generate an overview of quality control (QC) metrics (default: 'no').")
    fix_gene_names: Optional[str] = Field("no", description="Fix gene names in `adata.var_names` if Ensembl IDs are detected (default: 'no').")

    # Runtime/environment
    limit_threads: Optional[int] = Field(1, description="Apply thread limits to avoid OpenBLAS/OMP crashes (default: 1).")

# -------------------------------
# Launches the CellExpress pipeline with current settings.
# -------------------------------

@tool(args_schema=CellExpressArgs)
def run_cellexpress(**kwargs) -> str:
    """
    Launches CellExpress with the provided arguments. Runs in detached mode with logging and PID tracking.
    """
    try:
        args = CellExpressArgs(**cellexpress_cache.get_args())  # Validate arguments
        args_dict = args.dict()
    except ValidationError as e:
        return f"❌ Cannot launch CellExpress. Argument validation failed:\n{str(e)}"

    input_path = args_dict.get("input")
    if input_path is None:
        raise ValueError("❌ Missing required argument: --input")
    if not os.path.exists(input_path):
        raise ValueError(f"❌ The specified input path does not exist: {input_path}")

    env = os.environ.copy()
    env["PYTHONPATH"] = os.path.abspath(os.path.join(base_path, "..", "cellexpress"))
    base_cmd = ["python", "-m", "main"]

    for key, val in args_dict.items():
        if val is not None:
            if key in {"project", "tissue", "disease"}:
                val = str(val).replace(" ", "_")
            base_cmd.append(f"--{key}")
            base_cmd.append(str(val))

    try:
        # Construct log file path
        uid = uuid.uuid4().hex[:8]
        log_path = os.path.join(input_path, f"cellexpress_{uid}.log")
        err_path = os.path.join(input_path, f"cellexpress_{uid}.err")

        with open(log_path, "w") as log_file, open(err_path, "w") as err_file:
            proc = subprocess.Popen(base_cmd, stdout=log_file, stderr=err_file, env=env)

        return (
            f"🚀 CellExpress job submitted successfully!\n\n"
            f"🔧 PID: `{proc.pid}`\n"
            f"📄 Log file: `{log_path}`\n"
            f"🧾 CLI: `{' '.join(base_cmd)}`"
        )
    except Exception as e:
        return f"❌ Failed to launch CellExpress: {str(e)}"

# -------------------------------
# Checks the status of a CellExpress job by PID.
# -------------------------------

class PIDInput(BaseModel):
    pid: int = Field(..., description="Process ID (PID) of the CellExpress job to check")

@tool(args_schema=PIDInput)
def check_cellexpress_status(pid: int) -> str:
    """
    Checks the status of a CellExpress job by its PID.
    Returns whether the job is running or has completed.
    """
    try:
        proc = psutil.Process(pid)
        if proc.is_running() and proc.status() != psutil.STATUS_ZOMBIE:
            cmdline = " ".join(proc.cmdline())
            return (
                f"✅ Job is still running.\n"
                f"🔧 PID: `{pid}`\n"
                f"🧾 Command: `{cmdline}`"
            )
        else:
            return f"⚠️ PID `{pid}` is no longer running (may have completed or exited)."
    except psutil.NoSuchProcess:
        return f"❌ No process found with PID `{pid}`."
    except Exception as e:
        return f"❌ Failed to check status for PID `{pid}`: {str(e)}"

# -------------------------------
# Stops a running CellExpress job by its PID.
# -------------------------------

class TerminateJobArgs(BaseModel):
    pid: int = Field(..., description="Process ID (PID) of the running CellExpress job to terminate.")

@tool(args_schema=TerminateJobArgs)
def terminate_cellexpress_job(pid: int) -> str:
    """
    Terminates a running CellExpress job given its PID.
    """
    try:
        os.kill(pid, signal.SIGTERM)
        return f"🛑 Sent termination signal (SIGTERM) to process {pid}."
    except ProcessLookupError:
        return f"⚠️ No process found with PID {pid}. It may have already exited."
    except PermissionError:
        return f"❌ Permission denied to terminate process {pid}."
    except Exception as e:
        return f"❌ Failed to terminate process {pid}: {str(e)}"

# -------------------------------
# Shows the latest lines from a CellExpress log file.
# -------------------------------

class ReviewLogInput(BaseModel):
    path: str = Field(..., description="Directory containing the CellExpress log file.")
    filename: str = Field(..., description="Name of the log file (e.g., cellexpress_1234abcd.log).")
    n_lines: int = Field(default=30, description="Number of lines to include from the end of the log.")

@tool(args_schema=ReviewLogInput)
def review_cellexpress_log(path: str, filename: str = None, n_lines: int = 30) -> str:
    """
    Reads the tail of a CellExpress log file and asks the LLM to interpret the current stage.
    """
    full_path = os.path.join(path, filename)

    if not os.path.isfile(full_path):
        return f"❌ File not found: {full_path}"

    try:
        with open(full_path, "r") as f:
            lines = f.readlines()

        tail = lines[-n_lines:] if len(lines) >= n_lines else lines
        log_tail = "".join(tail)

        return (
            f"📄 **Recent CellExpress log output:**\n"
            f"```\n{log_tail}\n```\n\n"
            f"💡 Please analyze this log output and summarize the current stage or task being performed by the pipeline."
        )

    except Exception as e:
        return f"❌ Failed to read log file: {str(e)}"

# -------------------------------

class MathInput(BaseModel):
    x: float
    y: float

@tool(args_schema=MathInput)
def multiply_xy(x: float, y: float) -> float:
    """Multiply two floats together."""
    return x * y

@tool(args_schema=MathInput)
def add_xy(x: float, y: float) -> float:
    """Add two floats together."""
    return x + y

@tool(args_schema=MathInput)
def exponentiate_xy(x: float, y: float) -> float:
    """Exponentiate base x to power y."""
    return x ** y

# -------------------------------

tools = [
    multiply_xy, 
    add_xy, 
    exponentiate_xy, 
    fetch_article_metadata_url,
    fetch_article_metadata_pdf,
    store_article_metadata_file,
    store_geo_metadata_file,
    refine_article_metadata,
    refine_geo_metadata,
    fetch_geo_metadata,
    download_geo,
    download_gsm,
    download_file,
    make_directory,
    move_file,
    rename_file,
    list_directory,
    inspect_file,
    create_custom_csv,
    create_custom_json,
    create_custom_txt,
    remove_file_or_dir,
    get_file_size,
    get_working_directory,
    set_working_directory,
    fix_10x_file_format,
    get_cellexpress_info,
    configure_cellexpress,
    preview_cellexpress_config,
    reset_cellexpress_config,
    validate_cellexpress_config,
    run_cellexpress,
    review_cellexpress_log,
    terminate_cellexpress_job,
    check_cellexpress_status
]

# -------------------------------