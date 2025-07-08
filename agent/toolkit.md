# cellAtria Tool Reference Guide

**Metadata Analysis (7 tools)**
- `fetch_article_metadata_url`: Extracts structured metadata from a scientific article at a given URL.
- `fetch_article_metadata_pdf`: Extracts structured metadata from an uploaded scientific article PDF.
- `store_article_metadata_file`: Saves extracted article metadata to a file (JSON or TXT).
- `refine_article_metadata`: Holds article metadata for review and edits before saving.
- `fetch_geo_metadata`: Retrieves metadata for a GEO accession (e.g., GSE123456).
- `store_geo_metadata_file`: Saves extracted GEO metadata to a file (JSON or CSV).
- `refine_geo_metadata`: Holds GEO metadata for review and edits before saving.

**data Downloads & Organization (12 tools)**
- `download_geo`: Downloads all supplementary files for a GEO dataset.
- `download_gsm`: Downloads supplementary files for a list of GSM accessions.
- `download_file`: Downloads a file from a direct URL to your project folder.
- `make_directory`: Creates a new directory at the specified path.
- `list_directory`: Shows the contents of a directory in tree format.
- `inspect_file`: Displays the contents of a CSV, JSON, or TXT file.
- `get_file_size`: Shows the size of a file in megabytes.
- `move_file`: Moves a file from one location to another.
- `rename_file`: Renames a file or moves it to a new name/location.
- `remove_file_or_dir`: Deletes a file or directory (recursively if needed).
- `get_working_directory`: Displays your current working directory.
- `set_working_directory`: Sets or changes your working directory.

**File Preparation (4 tools)**
- `fix_10x_file_format`: Standardizes 10X Genomics filenames in a directory.
- `create_custom_csv`: Creates and saves a CSV file from tabular data.
- `create_custom_json`: Creates and saves a JSON file from structured data.
- `create_custom_txt`: Creates and saves a plain text file from a string.

**CellEpress Pipeline (9 tools)**
- `get_cellexpress_info`: Shows an overview of CellExpress pipeline capabilities.
- `configure_cellexpress`: Sets or updates a CellExpress pipeline argument.
- `preview_cellexpress_config`: Previews the current CellExpress CLI configuration.
- `reset_cellexpress_config`: Clears all cached CellExpress configuration arguments.
- `validate_cellexpress_config`: Validates the current CellExpress configuration.
- `run_cellexpress`: Launches the CellExpress pipeline with current settings.
- `review_cellexpress_log`: Shows the latest lines from a CellExpress log file.
- `terminate_cellexpress_job`: Stops a running CellExpress job by its PID.
- `check_cellexpress_status`: Checks the status of a CellExpress job by PID.