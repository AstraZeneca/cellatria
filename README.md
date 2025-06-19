# CellAtria: Agentic Triage of Regulated single-cell data Ingestion and Analysis

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

- **Docker**: Ensure [Docker](https://docs.docker.com/get-docker/) is installed and running on your system.
- **Data Directory**: Prepare a directory containing your input data and a valid `.env` configuration file (see [Configuration](#configuration) for details).

---

### 2️⃣ Launch CellAtria via Docker

Run the following command in your terminal, replacing `/Users/kjlg643/work/projects` with the path to your project directory:

```bash
docker run --platform=linux/amd64 -it --rm \
  -p 7860:7860 \
  -v $(pwd):/data \
  -v /Users/kjlg643/work/projects:/envdir \
  cellatria:v1.0.0 cellatria --env_path /envdir

</details>

---