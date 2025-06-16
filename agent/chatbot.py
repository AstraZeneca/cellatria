#!/usr/bin/env python3
# -------------------------------

import os
import sys
import argparse
from base import create_cellatria

# -------------------------------

HELP_TEXT = """
cellatria - Agentic Triage of Regulated single-cell data Ingestion and Analysis

Usage:
  cellatria [--help]

Environment variables:
  ENV_PATH           Path to your environment directory (default: /data)
  PROVIDER           LLM provider ("Azure OpenAI", "OpenAI", "Anthropic", "Google", "Local")
  ...                (see .env for full list)

Example:
  docker run -it --rm -p 7860:7860 -v /path/to/data:/data -e ENV_PATH=/data cellatria:latest
"""

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("--help", "-h", action="store_true")
parser.add_argument("--env_path", type=str, default=None)
args, unknown = parser.parse_known_args()

if args.help:
    print(HELP_TEXT)
    sys.exit(0)

# -------------------------------

env_path = args.env_path
graph, cellatria = create_cellatria(env_path)
cellatria.launch(server_name="0.0.0.0", server_port=7860, share=False)

# -------------------------------