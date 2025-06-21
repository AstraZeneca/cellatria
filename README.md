<p align="center">
  <a href="#"><img src="http://www.repostatus.org/badges/latest/active.svg" alt="Project Status"/></a>
  <a href="#"><img src="https://img.shields.io/badge/lifecycle-Stable-brightgreen.svg" alt="Lifecycle"/></a>
  <a href="https://langchain-ai.github.io/langgraph/"><img src="https://img.shields.io/badge/docs-latest-brightgreen" alt="Docs"/></a>
    <a href="#"><img src="https://img.shields.io/badge/made%20with-Python-830051?style=flat&logo=python&logoColor=white" alt="Made with Python"/></a>
  <a href="#"><img src="https://img.shields.io/badge/container-Docker-830051?style=flat&logo=docker&logoColor=white" alt="Docker"/></a>
  <a href="#"><img src="https://img.shields.io/badge/platform-GitHub-830051?style=flat&logo=github&logoColor=white" alt="Platform"/></a>
  <a href="https://github.com/langchain-ai/langgraph"><img src="https://img.shields.io/badge/built%20with-LangGraph-830051?style=flat&logo=python&logoColor=white" alt="LangGraph"/></a>
  <a href="#"><img src="https://img.shields.io/badge/agentic-AI%20Agent-830051?style=flat&logo=robotframework&logoColor=white" alt="Agentic"/></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-830051.svg" alt="License"/></a>
  <a href="https://github.com/nourin-nn/cellatria/actions/workflows/docker.yml"><img src="https://github.com/nourin-nn/cellatria/actions/workflows/docker.yml/badge.svg" alt="Docker Build Status"/></a>
</p>

<!-- Version Banner -->
<p align="center" width="100%">
<img width="15%" src="https://img.shields.io/badge/release-v1.0.0-brightgreen.svg?style=for-the-badge" alt="Release v1.0.0"/>
</p>
<p align="center" width="100%">
  <img width="55%" src="cellatria_git_logo.png"> 
</p>

**CellAtria** is a modular, agent-driven platform designed to automate the end-to-end ingestion, curation, and analysis of single-cell RNA sequencing (scRNA-seq) datasets. By integrating large language models (LLMs) with domain-specific bioinformatics toolchains, CellAtria streamlines the full lifecycle of single-cell studies — from literature parsing and metadata extraction to data acquisition and pipeline execution — all accessible through a natural language interface.


---

## 📘 Key Features
<details>
<br>

- Accepts primary research articles as **PDFs** or **URLs**.
- Extracts structured metadata such as sample annotations, organism, tissue type, and GEO accession identifiers.
- Resolves **GSE (study-level)** and **GSM (sample-level)** dependencies across GEO and organizes raw data accordingly.
- Orchestrates full ingestion pipelines and triggers **CellExpress** — an integrated, containerized scRNA-seq analysis framework.
- Empowers users to interact with data and tools via natural language, abstracting away scripting complexity.
- Supports metadata introspection, file transfers, directory traversal, and summarization tools.
- All actions are composed into reusable graph-based tools that operate as callable agent nodes.

> 📍 Additional details on the underlying toolkits and LLM initialization logic can be found in the [system prompts reference](https://github.com/nourin-nn/cellatria/blob/main/agent/system_prompts.md)

</details>

---

## 🚀 Getting Started
<details>

### 1️⃣ Prerequisites

- **Docker**: Install [Docker](https://docs.docker.com/get-docker/) and ensure the Docker daemon is running.
- **Data Directory**: Prepare a working directory to store your datasets and outputs.
- **Environment Configuration**: Provide a `.env` file with credentials and runtime configuration (see [Environment File Setup](#env_setup)).

---

### 2️⃣ Docker Images

Pull the latest cellAtria Docker image from **GitHub Container Registry (GHCR)**:

```bash
docker pull ghcr.io/nourin-nn/cellatria:v1.0.0
```

---

### 3️⃣  Launching CellAtria
Start the agent with the following command (replace paths with your actual directories)::

```bash
docker run --platform=linux/amd64 -it --rm \
  -p 7860:7860 \
  -v /path/to/your/project/directory:/data \
  -v /path/to/your/env/directory:/envdir \
  cellatria:v1.0.0 cellatria \
  --env_path /envdir
```

Command Breakdown:

- `-p 7860:7860`: Exposes the Gradio UI on port 7860.
- `-v /path/to/your/project/directory:/data`: Mounts your project directory into the container.
- `-v /path/to/your/env/directory:/envdir`: Mounts your `.env` directory for configuration.
- `cellatria:v1.0.0 cellatria`: Specifies the Docker image and the entrypoint command to launch the app inside the container.
- `--env_path /envdir`: Tells cellAtria where to find the `.env` file for provider setup.

> 📍 Once launched, the agent will initialize and provide a local or proxied URL for interaction.  Simply open the link printed in your terminal to begin using cellAtria through your browser.

</details>

---

<a name="env_setup"></a>
## ⚙️ Environment File Setup

<details>

### Quick Start

CellAtria requires a `.env` file to configure access to your chosen LLM provider and local runtime paths.

> 📍 Download the template [`.env`](./path/to/env_template.env) and replace it with the actual path to your `.env` in the repository.

### Supported LLM Backends

- `azure`: Azure OpenAI (enterprise-grade access to GPT models)
- `openai`: Standard OpenAI API (e.g., GPT-4, GPT-3.5)
- `anthropic`: Claude models via the Anthropic API
- `google`: Gemini models via Google Cloud / Vertex AI
- `local`: Offline models (e.g., Llama.cpp, Ollama, Hugging Face)

> 📍 Set the `PROVIDER` variable in your `.env` file to one of the supported values above. Only one provider can be active at a time.

### Instructions

1. Copy the `.env` template into your environment directory (e.g., `/envdir/.env`).
2. Set `PROVIDER=your_choice` in the file.
3. Fill in the required fields for your selected provider.

> 📍 You only need to configure the block for the provider you're using. The rest can remain commented.

</details>

---

## 📬 Contact

<details>
<br>

For help and questions please contact the [cellatria's maintenance team](mailto:ni.nouri@gmail.com).

</details>

---