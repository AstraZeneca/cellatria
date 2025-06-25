#!/usr/bin/env python3
# -------------------------------

import os
import sys
import argparse

# -------------------------------

HELP_TEXT = """
cellAtria - Agentic Triage of Regulated single-cell data Ingestion and Analysis
Version: 1.0.0
"""

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("--version", "-v", action="store_true")
parser.add_argument("--env_path", type=str, default="/mnt/work/projects")
args, unknown = parser.parse_known_args()

if args.version:
    print(HELP_TEXT)
    sys.exit(0)

# -------------------------------

from base import create_cellatria
graph, cellatria = create_cellatria(args.env_path)

# -------------------------------

print("\n")
print("✅ cellAtria successfully initialized.")
print("📍 Copy and paste the link below in your browser to interact with the agent:")
print("👉 http://0.0.0.0:7860\n")

# -------------------------------

cellatria.launch(server_name="0.0.0.0", server_port=7860, share=False, inbrowser=False, quiet=True)

# -------------------------------