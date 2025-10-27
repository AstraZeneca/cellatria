import os
import sys

# -------------------------------

# Import the function to create the cellatria app and graph
from base import create_cellatria

# -------------------------------

# Initialize the agent graph and Gradio Blocks app
env_path = "/mnt/work/projects" # user provides to docker
graph, cellatria = create_cellatria(env_path)

# -------------------------------

# Get Domino environment variables for run ID, project owner, and project name
RUN_ID = os.environ.get("DOMINO_RUN_ID")
PROJECT_OWNER = os.environ.get("DOMINO_PROJECT_OWNER")
PROJECT_NAME = os.environ.get("DOMINO_PROJECT_NAME")

# -------------------------------

# Launch the Gradio app with a custom root path and server settings for Domino
cellatria.launch(root_path="https://domino.astrazeneca.net/{}/{}/r/notebookSession/{}".format(PROJECT_OWNER, PROJECT_NAME, RUN_ID), server_name="0.0.0.0", server_port=8888, share = False)