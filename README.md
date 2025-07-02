<p align="center">
  <a href="#"><img src="https://img.shields.io/badge/made%20with-Python-830051?style=flat&logo=python&logoColor=white" alt="Made with Python"/></a>
  <a href="#"><img src="https://img.shields.io/badge/container-Docker-830051?style=flat&logo=docker&logoColor=white" alt="Docker"/></a>
  <a href="#"><img src="https://img.shields.io/badge/platform-GitHub-830051?style=flat&logo=github&logoColor=white" alt="Platform"/></a>
  <a href="https://github.com/langchain-ai/langgraph"><img src="https://img.shields.io/badge/built%20with-LangGraph-830051?style=flat&logo=python&logoColor=white" alt="LangGraph"/></a>
  <a href="#"><img src="https://img.shields.io/badge/agentic-AI%20Agent-830051?style=flat&logo=robotframework&logoColor=white" alt="Agentic"/></a>
  <br>
  <a href="#"><img src="http://www.repostatus.org/badges/latest/active.svg" alt="Project Status"/></a>
  <a href="#"><img src="https://img.shields.io/badge/lifecycle-Stable-brightgreen.svg" alt="Lifecycle"/></a>
  <a href="#"><img src="https://img.shields.io/badge/docs-latest-brightgreen?style=flat" alt="Docs"/></a>
  <a href="#"><img src="https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat" alt="Contributions welcome"/></a>
</p>

<!-- Version Banner -->
<p align="center" width="100%">
  <a href="https://github.com/nourin-nn/cellatria/actions/workflows/docker.yml"><img src="https://github.com/nourin-nn/cellatria/actions/workflows/docker.yml/badge.svg" alt="Docker Build Status"/></a>
  <br>
  <img width="15%" src="https://img.shields.io/badge/release-v1.0.0-4444AA.svg?style=for-the-badge" alt="Release v1.0.0"/>
</p>
<p align="center" width="100%">
  <img width="100%" src="cellatria_git_logo.png"> 
</p>

---

## ✨ Introduction
<details>
<br>

**cellAtria** is a modular, agent-driven platform designed to automate the end-to-end metadata curation, data ingestion, and single-cell RNA sequencing (scRNA-seq) analysis. By integrating large language models (LLMs) with domain-specific bioinformatics toolchains, cellAtria streamlines the full lifecycle of single-cell studies — from literature parsing and metadata extraction to data acquisition and pipeline execution — all accessible through a natural language interface.

> **cellAtria** empowers users to interact with data and tools via natural language, abstracting away scripting complexity.

</details>

---

## 💡 Key Features
<details>
<br>

- Accepts primary research articles as **PDFs** or **URLs**.
- Extracts structured metadata such as sample annotations, organism, tissue type, and GEO (Gene Expression Omnibus) accession identifiers.
- Resolves **GSE (study-level)** and **GSM (sample-level)** dependencies across GEO and organizes data accordingly.
- Orchestrates full ingestion pipelines and triggers [**CellExpress**](#cellexpress) — an integrated, containerized scRNA-seq analysis framework.
- Supports metadata introspection, file transfers, directory traversal, and summarization tools.
- All actions are composed into reusable graph-based tools that operate as callable agent nodes.

> Additional details on the underlying toolkits and LLM initialization logic can be found in the [system prompts reference](https://github.com/nourin-nn/cellatria/blob/main/agent/system_prompts.md)

</details>

---

## 🚀  Getting Started
<details>

### (1) Prerequisites

- **Docker**: Install [Docker](https://docs.docker.com/get-docker/) and ensure the Docker daemon is running.
- **Environment Configuration**: Provide a `.env` file with credentials and parameters (see [LLM Configuration](#env_setup)).

---

### (2) Docker Images

The **cellAtria** GitHub repository is equipped with a GitHub Actions workflow that automatically builds and publishes a Docker image to the [GitHub Container Registry](https://github.com/nourin-nn/cellatria/pkgs/container/cellatria) upon each successful release or update.

Pull the latest **cellAtria** Docker image using:

```bash
# Run this command in your terminal
docker pull ghcr.io/nourin-nn/cellatria:v1.0.0
```

> This container includes all necessary dependencies to launch the **cellAtria** agent.

---

### (3)  Launching Agent
Start the agent with the following command (replace paths with your actual directories)::

```bash
# Run this command in your terminal
docker run -it --rm \
  -p 7860:7860 \
  -v /path/to/your/project/directory:/data \
  -v /path/to/your/env/directory:/envdir \
  ghcr.io/nourin-nn/cellatria:v1.0.0 cellatria \
  --env_path /envdir
```

Command Breakdown:

- `-p 7860:7860`: Exposes the agent user interface (UI) on port 7860.
- `-v /path/to/your/project/directory:/data`: Mounts your project directory into the container.
- `-v /path/to/your/env/directory:/envdir`: Mounts your `.env` directory for configuration (see [LLM Configuration](#env_setup)).
- `ghcr.io/nourin-nn/cellatria:v1.0.0 cellatria`: Specifies the Docker image and the entrypoint command to launch the app inside the container.
- `--env_path /envdir`: Tells agent where to find the `.env` file for provider setup.

> macOS users with Apple Silicon (M1/M2): You may encounter a warning due to platform mismatch. To ensure compatibility, add `--platform=linux/amd64` when running the container. 

> Once launched, the agent will initialize and provide a local URL for interaction. Simply open the link printed in your terminal to begin using cellAtria through your browser.

---

**Mounting a Working Directory:**

When running the container, any host directory you want to access must be explicitly mounted using Docker’s `-v` (volume) flag. The container only has access to the paths you specify during mounting.

For example, the following command:

```bash
-v /absolute/path/on/host:/data
```

makes the contents of `/absolute/path/on/host` on your host machine available inside the container at `/data`.

> If you set a working directory inside the container (e.g., `my_project`), ensure that all command-line flags referencing it use the container's path, i.e., `/data/my_project`.

> Accessing paths outside the mounted host directory will not be possible from within the container.

</details>

---

<a name="env_setup"></a>
## 🛠️ LLM Configuration

<details>

### Quick Start

cellAtria requires a `.env` file to configure access to your chosen LLM provider.

> You can download the template [`.env`](https://github.com/nourin-nn/cellatria/blob/main/.env), fill in the necessary credentials and parameters. Ensure the directory containing the `.env` file is mounted into the container.

### Supported LLM Backends

- `azure`: Azure OpenAI (enterprise-grade access to GPT models)
- `openai`: Standard OpenAI API (e.g., GPT-4, GPT-3.5)
- `anthropic`: Claude models via the Anthropic API
- `google`: Gemini models via Google Cloud / Vertex AI
- `local`: Offline models (e.g., Llama.cpp, Ollama, Hugging Face)

> Set the `PROVIDER` variable in your `.env` file to one of the supported values above. Only one provider can be active at a time. You only need to configure the block for the provider you're using. The rest can remain commented.

</details>

---

<a name="cellexpress"></a>
## 🚂 CellExpress Engine
<details>
<br>

**CellExpress** is a companion pipeline embedded within the **cellAtria** framework. It delivers a reproducible and automated workflow for processing single-cell RNA-seq datasets (scRNA-seq) — from raw count matrices to comprehensive cell type annotations and report generation.

- Designed to lower bioinformatics barriers, **CellExpress** implements a comprehensive set of state-of-the-art, Scanpy-based processing stages, including quality control (performed globally or per sample), data transformation (including normalization, highly variable gene selection, and scaling), dimensionality reduction (UMAP and t-SNE), graph-based clustering, and marker gene identification. Additional tools are integrated to support advanced analysis tasks, including doublet detection, batch correction, and automated cell type annotation using both tissue-agnostic and tissue-specific models. All analytical steps are executed sequentially under centralized control, with parameters fully configurable via a comprehensive input schema. 

> For full details, usage instructions, and configuration options, refer to the [CellExpress README](https://github.com/nourin-nn/cellatria/blob/main/cellexpress/README.md).

</details>

---

## 🛠️ Computing Environment

<details>
<br>

The `Dockerfile` defines the dedicated computing environment for executing **cellatria** and the co-developed **CellExpress** pipelie in a consistent and reproducible manner. 
It includes all required Python and R dependencies, along with support for HTML reporting and visualization. 
Built on an Ubuntu-based system, the environment also provides essential system-level packages to support end-to-end 
pipeline execution. 

</details>

---

## 🧠 Usage Intuition
<details>
<br>

While **cellAtria** supports flexible, user-driven interactions, its functionality is governed by an underlying **execution narrative** — a structured flow of modular actions that define how tasks are interpreted, routed, and executed. Users may invoke any module independently; however, for optimal results and seamless orchestration, we recommend following the intended workflow trajectory below.

**cellAtria's internal logic integrates:**

1. **Document Parsing** — Extracts structured metadata from narrative-formatted scientific documents (article URL or PDF).  
2. **Accession Resolution** — Identifies relevant GEO (Gene Expression Omnibus) accession IDs from parsed metadata.  
3. **Dataset Retrieval** — Downloads datasets directly from public repositories.  
4. **File & Data Organization** — Structures downloaded content into a consistent directory schema.  
5. **Pipeline Configuration** — Prepares **CellExpress** arguments and environmental parameters for execution.  
6. **CellExpress Execution** — Launches the standardized single-cell analysis pipeline.  

> This modular, agent-guided framework allows users to begin at any point while preserving logical consistency across steps.

</details>

---

## 📬 Contact

<details>
<br>

We welcome community contributions! If you'd like to help improve **cellAtria**, please feel free to 
suggest enhancements.


| Role         | Name               | Contact                                     |
|--------------|--------------------|---------------------------------------------|
| Maintainer   | Nima Nouri         | [ni.nouri@gmail.com](mailto:ni.nouri@gmail.com) | 

</details>

---