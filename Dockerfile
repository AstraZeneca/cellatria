# ==============================================================================
# CellAtria Docker Image
# ==============================================================================
# This Dockerfile creates a comprehensive environment for CellAtria, an agentic
# triage system for regulated single-cell data ingestion and analysis.
# ==============================================================================

# Get base image 
FROM ubuntu:22.04
# Set non-interactive frontend
ARG DEBIAN_FRONTEND=noninteractive
# -----------------------------------
# OCI-compliant image metadata
LABEL org.opencontainers.image.title="cellAtria"
LABEL org.opencontainers.image.version="v1.0.0"
LABEL org.opencontainers.image.description="CellAtria: Agentic Triage of Regulated single-cell data Ingestion and Analysis."
LABEL org.opencontainers.image.authors="Nima Nouri <nima.nouri@astrazeneca.com>"
LABEL org.opencontainers.image.source="https://github.com/azu-oncology-rd/cellatria"
# -----------------------------------
# System package installation
RUN apt-get update --fix-missing && \ 
		apt-get install -y --no-install-recommends --fix-missing \
		pkgconf \
		build-essential \
        python3-pip \
        python3-dev \
        python3-setuptools \
        python-is-python3 \
        apt-utils \
        ed \
        less \
        locales \
        vim-tiny \
        wget \
        tree \
        ca-certificates \
        apt-transport-https \
        sshpass \
        openssh-client \
        curl \
        cmake \
        git \
        bzip2 \
        vim \
        nano \
        acl \
        software-properties-common \
        fonts-liberation \
        gcc \
        g++ \
        gfortran \
        libssl-dev \
        libcurl4-openssl-dev \
        libhdf5-dev \
        libmpich-dev \
        libblosc-dev \
        liblzma-dev \
        libzstd-dev \
        libopenblas-dev \        
        libxml2-dev \
        libnlopt-dev \
        libicu-dev \
        libjpeg62 \
        libgeos-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        libglpk-dev \
        libcairo2-dev \
        libffi-dev \
        libgit2-dev \
        libblas-dev \
        liblapack-dev \
        gnupg2 \
        lsb-release \
        dirmngr \
        gdebi-core \
        psmisc \
        jq \
        ninja-build \
        clang \ 
        clang-tidy \
        && \
        apt-get clean && \
        rm -rf /var/lib/apt/lists/*
# -----------------------------------
# Install Pandoc for R Markdown document conversion and report generation
RUN apt-get update && apt-get install -y pandoc
# -----------------------------------
# Python Package Installation
RUN pip install --no-cache-dir --upgrade setuptools wheel
RUN pip install --no-cache-dir numpy pandas scipy
RUN pip install --no-cache-dir plotly matplotlib seaborn anndata scanpy tqdm scikit-learn h5py networkx scrublet nose annoy 'zarr<3'
RUN pip install --no-cache-dir torch torchvision torchaudio torchsummary torchopt entmax
RUN pip install --no-cache-dir harmonypy igraph leidenalg celltypist scimilarity
RUN pip install --no-cache-dir fa2_modified
RUN pip install --no-cache-dir --upgrade openai langchain langchainhub langchain-community langchain-openai langchain-core langchain-anthropic
RUN pip install --no-cache-dir --upgrade transformers accelerate gradio langgraph PyMuPDF GEOparse beautifulsoup4 google-generativeai
# -----------------------------------
# R Installation and Configuration
# Update indices
RUN apt update -qq
# Install two helper packages
RUN apt install apt-transport-https software-properties-common dirmngr
RUN curl -sSL \
'http://keyserver.ubuntu.com/pks/lookup?op=get&search=0xE298A3A825C0D65DFD57CBB651716619E084DAB9' \
| apt-key add -
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
# add the R repo from CRAN -- adjust 'jammy'
RUN add-apt-repository 'deb http://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/'
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y \
    r-base \
    r-base-core \
    r-base-dev \
    r-recommended
# -----------------------------------
# Install R packages
RUN Rscript -e "options(repos='http://cran.rstudio.com/'); install.packages('devtools', clean=TRUE)"
RUN Rscript -e "options(repos='http://cran.rstudio.com/'); install.packages(c('progress','Rcpp','Rcpp11','RcppAnnoy'), clean=TRUE)"
ARG R_DEPS="c('ggplot2', 'dplyr', 'gtools', 'grid', 'gridtext', \
                'jsonlite', 'kableExtra', 'DT', 'scales','RColorBrewer', \
                'plotly', 'visNetwork', 'ggrepel', 'gtools', 'viridis', \
                'gridExtra', 'tidyr', 'DescTools')"	
RUN Rscript -e "options(repos='http://cran.rstudio.com/'); install.packages(${R_DEPS}, clean=TRUE)"
# -----------------------------------
# Copy the CellAtria application files into the container
# and set up the working directory structure
# -----------------------------------
# Copy all files into Docker
RUN mkdir -p /opt/cellatria
WORKDIR /opt/cellatria
COPY . /opt/cellatria/
# -----------------------------------
# Make cellatria CLI callable via `cellatria`
RUN chmod +x /opt/cellatria/agent/chatbot.py
RUN ln /opt/cellatria/agent/chatbot.py /usr/local/bin/cellatria
# -----------------------------------
# Make cellexpress CLI callable via `cellexpress`
RUN chmod +x /opt/cellatria/cellexpress/main.py
RUN ln -s /opt/cellatria/cellexpress/main.py /usr/local/bin/cellexpress
# -----------------------------------
# The VOLUME instruction and the -v option to docker run 
VOLUME /data
WORKDIR /data
# -----------------------------------
# Expose the port used by Gradio
EXPOSE 7860
# -----------------------------------
# Configure Python paths and data locations
ENV PYTHONPATH=/opt/cellatria/agent
ENV ENV_PATH=/data
# -----------------------------------
# Default command launches the CellAtria chatbot interface
CMD ["/usr/local/bin/cellatria"]
# -----------------------------------