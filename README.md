[![Project Status](http://www.repostatus.org/badges/latest/active.svg)]()
[![Lifecycle](https://img.shields.io/badge/lifecycle-Stable-brightgreen.svg)]()
[![Made with Python](https://img.shields.io/badge/made%20with-Python-3776AB?style=flat&logo=python&logoColor=white)]()
[![Docker](https://img.shields.io/badge/container-Docker-2496ED?style=flat&logo=docker&logoColor=white)]()
[![Platform](https://img.shields.io/badge/platform-GitHub-black?style=flat&logo=github&logoColor=white)]()
[![Docs](https://img.shields.io/badge/docs-latest-blue)](https://langchain-ai.github.io/langgraph/)
[![LangGraph](https://img.shields.io/badge/built%20with-LangGraph-6A5ACD?style=flat&logo=python&logoColor=white)](https://github.com/langchain-ai/langgraph)
[![Agentic](https://img.shields.io/badge/agentic-AI%20Agent-FFB300?style=flat&logo=robotframework&logoColor=white)]()
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

<!-- Version Banner -->
<img src="https://img.shields.io/badge/release-v1.0.1-blue.svg?style=for-the-badge" alt="Release v1.0.0"/>

<p align="center" width="100%">
  <img width="70%" src="cellatria_git_logo.png"> 
</p>

**CellAtria** is an intelligent, agent-driven platform for automating the end-to-end ingestion, preparation, and analysis of single-cell RNA-seq (scRNA-seq) datasets. It integrates large language models (LLMs) with domain-specific computational toolchains to streamline single-cell workflows - from literature parsing and metadata extraction to data acquisition and pipeline execution - via a natural language interface.

---

## 📌 Key Features
<details>

### 🔍 Literature-Guided Initialization
- Accepts a **URL** or **PDF** of a primary article.
- Extracts structured metadata (e.g., sample annotations, GEO accessions) from manuscripts.
- Enables **zero-shot** dataset discovery and processing through document parsing.

### 📂 Metadata-Aware Data Acquisition
- Supports **GSE-level (study-wide)** and **GSM-level (sample-specific)** retrieval from public repositories.
- Automatically resolves GEO relationships and organizes files in compliant directory schemas.
- Prepares `metadata.csv` aligned with pipeline requirements.

### 🧬 Integrated Analysis with CellExpress
- Executes the co-developed **CellExpress** pipeline—an end-to-end, containerized scRNA-seq analysis framework.
- Performs normalization, HVG selection, batch correction, clustering, marker gene detection, and cell type annotation.
- Supports tissue-agnostic (SCimilarity) and tissue-specific (CellTypist) annotations.

### 🧠 Conversational Workflow Orchestration
- Facilitates **multi-turn, context-aware reasoning** to walk users through dataset processing.
- Automatically executes validated commands, logging each step for reproducibility.
- Offers fine-grained control over tool behavior through structured dialogue.

### 🔧 Tool-Driven Modular Architecture
- Tools are embedded as graph nodes and accessed via natural language.
- Supports metadata inspection, file downloads, directory traversal, report summarization, and more.
- Agent actions are both **traceable** and **auditable**, supporting regulatory workflows.

</details>

---

## 📘 Getting Started
<details>
<summary>Expand for step-by-step instructions</summary>

### 1️⃣ Prerequisites

- **Docker**: Install [Docker](https://docs.docker.com/get-docker/) and ensure it is running.
- **Data Directory**: Prepare a directory with your input data and a valid `.env` configuration file (see [Configuration](#configuration)).

---

### 2️⃣ Launch CellAtria via Docker

Run the following command in your terminal (replace `/path/to/your/project/directory` and `/path/to/your/env/directory` with your actual directories):

```bash
docker run --platform=linux/amd64 -it --rm \
  -p 7860:7860 \
  -v /path/to/your/project/directory:/data \
  -v /path/to/your/env/directory:/envdir \
  cellatria:v1.0.0 cellatria --env_path /envdir
```

**Command Breakdown:**
- `-p 7860:7860` maps the app port to your host.
- `-v /path/to/your/project/directory:/data` mounts your project directory as `/data` in the container.
- `-v /path/to/your/env/directory:/envdir` mounts your environment directory (with `.env`) as `/envdir`.
- `cellatria --env_path /envdir` launches the agent with your environment directory.

</details>

---

## ⚙️ Environment File Setup
<details>

Before running CellAtria, you must provide a `.env` file containing your configuration and API keys.  
This file tells CellAtria which LLM provider to use and how to connect to it.

- **Download the `.env` template:**  
  [CellAtria .env Template](./path/to/your/env_template.env)  
  *(Replace with the actual path or link to your template file in the repository)*

### Compatible LLM Providers

CellAtria supports seamless integration with the following large language model (LLM) providers:

- **Azure OpenAI**  
  Use enterprise-grade Azure OpenAI endpoints for secure, scalable access to GPT models.
- **OpenAI**  
  Connect directly to OpenAI’s public API for models like GPT-4 and GPT-3.5.
- **Anthropic**  
  Leverage Claude models via the Anthropic API.
- **Google Gemini / Vertex AI**  
  Access Google’s Gemini models through the Google Cloud API.
- **Local**  
  Run local models (e.g., Llama.cpp, Ollama, Hugging Face) for private, offline inference.

> **Note:**  
> Only one provider can be active at a time. Set the `PROVIDER` variable in your `.env` to your desired backend.

---

### Instructions

1. **Download or copy** the `.env` template from the link above into your environment directory (e.g., `/envdir/.env`).
2. **Set** the `PROVIDER` variable to match your desired LLM backend (see list above).
3. **Fill in** the required fields for your chosen provider (see comments in the template).
4. **Keep your API keys secure**—do not share your `.env` file publicly.

> **Tip:**  
> You only need to fill in the section for your selected provider.

For more details on each provider’s configuration, see the [Configuration Guide](#configuration-guide).

</details>

---