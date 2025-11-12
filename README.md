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
  <img src="https://img.shields.io/badge/Maturity%20Level-ML--0-brightgreen" alt="Maturity level-0"/>
  <br>
  <a href="https://github.com/AstraZeneca/cellatria/actions/workflows/docker.yml"><img src="https://github.com/AstraZeneca/cellatria/actions/workflows/docker.yml/badge.svg?branch=main" alt="cellatria image"/></a>
  <br>
</p>

<!-- Version Banner -->
<p align="center" width="100%">  
  <img width="15%" src="https://img.shields.io/badge/release-v1.0.0-4444AA.svg?style=for-the-badge" alt="Release v1.0.0"/>
</p>

<!-- CellAtria Infographic -->
<p align="center" width="100%">
  <img width="100%" src="cellatria_git_fig0.png"> 
</p>

---

## ✨ Introduction
<details>
<br>


**CellAtria** is an agentic AI system that enables **full-lifecycle, document-to-analysis automation** in single-cell research. It integrates natural language interaction with a robust, graph-based, multi-actor execution framework. The system orchestrates diverse tasks, ranging from literature parsing and metadata extraction to dataset retrieval and downstream scRNA-seq analysis via the co-developed [**CellExpress**](#cellexpress) pipeline.

> Through its comprehensive interface, **CellAtria** empowers users to engage with a language model augmented by task-specific tools. This eliminates the need for manual command-line operations, accelerating data onboarding and the reuse of public single-cell resources.

<p align="center" width="100%">
  <img width="55%" src="cellatria_git_fig1.png"> 
</p>
<p align="center">
  <small><em><strong>Language model-mediated orchestration of toolchains</strong>. Upon receiving a user prompt, the CellAtria interface transfers the request to the LLM agent, which interprets intent and autonomously invokes relevant tools. Outputs are returned through the interface, completing a full cycle of context-aware execution.</em></small>
</p>


</details>

---

## 💡 Key Features
<details>
<br>

- **Flexible Input**: Accepts primary research articles as **PDFs** or **URLs** for seamless integration.
- **Automated Metadata Extraction**: Extracts structured metadata, including sample annotations, organism, tissue type, and GEO (Gene Expression Omnibus) accession identifiers.
- **Intelligent Data Retrieval**: Resolves and organizes GEO datasets by accessing both **GSE (study-level)** and **GSM (sample-level)** records, ensuring structured and comprehensive data retrieval.
- **Integrated Analysis Pipeline**: Orchestrates full pipeline configuration and launches [**CellExpress**](#cellexpress), a containerized framework for standardized scRNA-seq analysis, ensuring reproducible results.
- **Enhanced User Control**: Enables metadata editing, secure file transfers, and direct file system management within the agent session.
- **Modular & Reusable Architecture**: Composes all core actions into reusable, graph-based tools that serve as callable agent nodes, fostering extensibility.

> Additional details on the underlying toolkits can be found in the [toolkit reference](https://github.com/AstraZeneca/cellatria/blob/main/agent/toolkit.md)

</details>

---

## 🚀  Getting Started
<details>
<br>

Users can initiate CellAtria agent in just three simple steps:

1. Open a Terminal,
2. Pull the preconfigured Docker image (which includes all necessary dependencies), and
3. Run the container, specifying the path to the file that contains credentials for your preferred LLM provider.

This minimal setup ensures consistent execution across environments—whether on a laptop, workstation, or HPC—without requiring users to manage dependencies manually. Detailed instructions follow below.

---

### (1) Prerequisites

- **Docker**: Install [Docker](https://docs.docker.com/get-docker/) and ensure the Docker daemon is running.
- **Environment Configuration**: Provide a `.env` file with credentials and parameters (see [LLM Configuration](#env_setup) section below).

---

### (2) Docker Images

The **CellAtria** repository includes a GitHub Actions workflow that builds and publishes a preconfigured Docker image to the [GitHub Container Registry](https://github.com/AstraZeneca/cellatria/pkgs/container/cellatria).

Pull the latest **CellAtria** Docker image using:

```bash
docker pull ghcr.io/astrazeneca/cellatria:v1.0.0
```

> This image contains all dependencies needed to run  the **CellAtria** agent in a consistent environment.

If you prefer to build the image locally instead of pulling from GHCR:

**Clone the repository**
```bash
git clone https://github.com/AstraZeneca/cellatria.git
cd cellatria
```

**Build the Docker image**
```bash
docker build -t ghcr.io/astrazeneca/cellatria:v1.0.0 .
```

**Verify the image**
```bash
docker images | grep cellatria
```

---

### (3)  Launching Agent
Start the agent with the following command (replace paths with your actual directories):

```bash
# Run this command in your terminal
docker run -it --rm \
  -p 7860:7860 \
  -v /path/to/your/project/directory:/data \
  -v /path/to/your/env/directory:/envdir \
  ghcr.io/astrazeneca/cellatria:v1.0.0 cellatria \
  --env_path /envdir
```

Command Breakdown:

- `-p 7860:7860`: Exposes the agent user interface (UI) on port 7860.
- `-v /path/to/your/project/directory:/data`: Mounts your project directory into the container.
- `-v /path/to/your/env/directory:/envdir`: Mounts your `.env` directory for configuration (see [LLM Configuration](#env_setup) section below).
- `ghcr.io/astrazeneca/cellatria:v1.0.0 cellatria`: Specifies the Docker image and the entrypoint command to launch the app inside the container.
- `--env_path /envdir`: Tells agent where to find the `.env` file for provider setup.

> macOS users with Apple Silicon (M1/M2): You may encounter a warning due to platform mismatch. To ensure compatibility, add `--platform=linux/amd64` when running the container (i.e., `docker run --platform=linux/amd64 -it --rm`). 

Once launched, the agent will initialize and provide a local URL for interaction. Simply open the link printed in your terminal to begin using CellAtria through your browser.

---

**Mounting a Working Directory:**

When running the container, any host directory you want the container to access must be explicitly mounted using Docker’s `-v` (volume) flag. The container can only see and interact with the directories you specify at runtime.

For example, the following command:

```bash
-v /absolute/path/on/host:/data
```

makes the contents of `/absolute/path/on/host` on your host machine available inside the container at `/data`.

> If you set a working directory inside the container (e.g., `my_project`), make sure to reference it using the container’s path — for instance: `/data/my_project`. Attempting to access files or directories outside the mounted path from within the container will fail, as they are not visible to the container’s filesystem.

</details>

---

<a name="env_setup"></a>
## 🛠️ LLM Configuration

<details>

### Quick Start

CellAtria requires a `.env` file to configure access to your chosen LLM provider. You can download the template [`.env`](https://github.com/AstraZeneca/cellatria/blob/main/.env), fill in the necessary credentials and parameters. Ensure the directory containing the `.env` file is mounted into the container.

### Supported LLM Backends

- `azure`: Azure OpenAI (enterprise-grade access to GPT models)
- `openai`: Standard OpenAI API (e.g., GPT-4, GPT-3.5)
- `anthropic`: Claude models via the Anthropic API
- `google`: Gemini models via Google Cloud / Vertex AI

> Set the `PROVIDER` variable in your `.env` file to one of the supported values above. Only one provider can be active at a time. 

> You only need to configure the block for the provider you're using. The rest can remain commented.

</details>

---

<a name="cellexpress"></a>
## 🚂 CellExpress Engine
<details>
<br>

**CellExpress** is a companion pipeline embedded within the **CellAtria** framework. It delivers a reproducible and automated workflow for processing single-cell RNA-seq datasets (scRNA-seq) - from raw count matrices to comprehensive cell type annotations and report generation.

> Designed to lower bioinformatics barriers, **CellExpress** implements a comprehensive set of state-of-the-art, Scanpy-based processing stages, including quality control (performed globally or per sample), data transformation (including normalization, highly variable gene selection, and scaling), dimensionality reduction (UMAP and t-SNE), graph-based clustering, and marker gene identification. Additional tools are integrated to support advanced analysis tasks, including doublet detection, batch correction, and automated cell type annotation using both tissue-agnostic and tissue-specific models. All analytical steps are executed sequentially under centralized control, with parameters fully configurable via a comprehensive input schema. 

---

### Run CellExpress in Standalone Mode

**CellExpress** is a fully standalone pipeline for comprehensive scRNA-seq data analysis. It can be orchestrated either through an agentic system - as incorporated into the **CellAtria** framework - or via direct command-line execution.

To execute the CellExpress pipeline directly using Docker, use the following command:

```bash
# Run this command in your terminal
docker run -it --rm \
  -v /path/to/your/local/data:/data \
  ghcr.io/astrazeneca/cellatria:v1.0.0 cellexpress \
    --input /data \
    --project your_project_name \
    --species `species` \
    --tissue `tissue` \
    --disease `disease` \
    [--additional `options`...]
```

Command Breakdown:

- `-v /path/to/your/local/data:/data`: Mounts your project directory into the container.
- `ghcr.io/astrazeneca/cellatria:v1.0.0 cellexpress`: Specifies the Docker image and the entrypoint command to launch **CellExpress** inside the container.
-  [--additional `options`...]: arguments to configure pipeline.

> macOS users with Apple Silicon (M1/M2): You may encounter a warning due to platform mismatch. To ensure compatibility, add `--platform=linux/amd64` when running the container (i.e., `docker run --platform=linux/amd64 -it --rm`). 

> For full details, usage instructions, and configuration options, refer to the [CellExpress README](https://github.com/AstraZeneca/cellatria/blob/main/cellexpress/README.md).

> See [CellExpress in Action](https://github.com/AstraZeneca/cellatria/blob/main/docs/cellespress_usage.md) for a Docker-based execution walkthrough.

</details>

---

## 🛠️ Computing Environment

<details>
<br>

The `Dockerfile` defines the dedicated computing environment for executing **CellAtria** and the co-developed **CellExpress** pipelie in a consistent and reproducible manner. 
It includes all required Python and R dependencies, along with support for HTML reporting and visualization. 
Built on an Ubuntu-based system, the environment also provides essential system-level packages to support end-to-end 
pipeline execution. 

</details>

---

## 🧠 Usage Intuition
<details>
<br>

While **CellAtria** supports flexible, user-driven interactions, its functionality is governed by an underlying **execution narrative** — a structured flow of modular actions that define how tasks are interpreted, routed, and executed. Users may invoke any module independently; however, for optimal results and seamless orchestration, we recommend following the intended workflow trajectory below.

**CellAtria's internal logic integrates the following key stages:**

1.  **Document Parsing** - Extracts structured metadata from narrative-formatted scientific documents (article URL or PDF).
2.  **Accession Resolution** - Identifies relevant GEO (Gene Expression Omnibus) accession IDs from parsed metadata.
3.  **Dataset Retrieval** - Downloads datasets directly from public repositories.
4.  **File & Data Organization** - Structures downloaded content into a consistent directory schema for analysis.
5.  **Pipeline Configuration** - Prepares **CellExpress** arguments and environmental parameters for execution.
6.  **CellExpress Execution** - Launches the standardized single-cell analysis pipeline in a detached mode. 

> This modular, agent-guided framework allows users to begin at any point while preserving logical consistency across steps.

> See the [CellAtria Walkthrough](https://github.com/AstraZeneca/cellatria/blob/main/docs/cellatria_walkthrough_git.pdf)
 for a step-by-step visual tutorial.

</details>

---

## 📖 Related Publication
<details>
<br>

If you use this repository, please cite:

> Nima Nouri, et al. (2025). An Agentic AI Framework for Ingestion and Standardization of Single-Cell RNA-seq Data Analysis. *bioRxiv*. https://doi.org/10.1101/2025.07.31.667880

```
@article{nouri2025agentic,
  title={An Agentic AI Framework for Ingestion and Standardization of Single-Cell RNA-seq Data Analysis},
  author={Nouri, Nima and Artzi, Ronen and Savova, Virginia},
  journal={bioRxiv},
  year={2025},
  publisher={Cold Spring Harbor Laboratory}
}
```

</details>

---

## ⚖️ Ethical Use, Data Governance, and Compliance
<details>
<br>

As automated frameworks increasingly interface with regulated data sources, ensuring ethical and compliant data handling becomes essential. When extracting metadata from scientific literature, CellAtria’s schema captures repository identifiers, publisher information, and author-reported conflicts of interest, providing legal and provenance context for each dataset. For proprietary or sensitive data, the Docker-based architecture of CellAtria enables reproducible execution within secure, auditable, and GxP-compliant environments, a critical requirement in clinical research. Nonetheless, data governance policies supersede automation, and it remains the user’s responsibility to verify licensing terms and consent scope. While agentic systems can augment user capabilities, they are not a substitute for ethical oversight, which is preserved through a maintained human-in-the-loop paradigm.

</details>

---

## 📬 Contact
<details>
<br>

| Role         | Name               | Contact                                     |
|--------------|--------------------|---------------------------------------------|
| Author/Maintainer   | Nima Nouri         | [nima.nouri@astrazeneca.com](mailto:nima.nouri@astrazeneca.com) | 

</details>

---