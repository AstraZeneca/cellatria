# getting bas image
FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive
# ===================================
# installing python and all required packages
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
# ===================================
# Checking python installation
RUN echo "Python3 version:" && python3 --version && \
    echo "Pip3 version:" && pip3 --version && \
    echo "Python3 path: $(which python3)" && \
    echo "Pip3 path: $(which pip3)"
RUN echo "Python version:" && python --version && \
    echo "Pip version:" && pip --version && \
    echo "Python path: $(which python)" && \
    echo "Pip path: $(which pip)"
# =====================================================
RUN pip install --no-cache-dir --upgrade setuptools wheel
RUN pip install --no-cache-dir numpy pandas scipy
RUN pip install --no-cache-dir plotly matplotlib seaborn anndata scanpy tqdm scikit-learn h5py networkx scrublet nose annoy 'zarr<3'
RUN pip install --no-cache-dir torch torchvision torchaudio torchsummary torchopt entmax
RUN pip install --no-cache-dir harmonypy igraph leidenalg celltypist scimilarity
RUN pip install --no-cache-dir fa2_modified
RUN pip install --no-cache-dir --upgrade openai langchain langchainhub langchain-community langchain-openai langchain-core langchain-anthropic
RUN pip install --no-cache-dir --upgrade transformers accelerate gradio langgraph PyMuPDF GEOparse beautifulsoup4 google-generativeai
# =====================================================
### Verify installation
RUN python -c "import sys, torch; \
    print(f'Python installed successfully! Version: {sys.version}'); \
    print(f'Torch installed successfully! Version: {torch.__version__}')"
RUN python3 -c "missing=[]; packages=['setuptools','wheel','numpy','pandas','scipy','plotly','matplotlib','seaborn','anndata','scanpy','tqdm','sklearn','h5py','networkx','scrublet','nose','annoy','zarr','torch','torchvision','torchaudio','torchsummary','torchopt','entmax','harmonypy','igraph','leidenalg','celltypist','scimilarity','fa2_modified','openai','langchain','langchainhub','langchain_community','langchain_openai','transformers','accelerate','langchain_core','gradio','langgraph','fitz','GEOparse', 'bs4']; import importlib; [missing.append(pkg) for pkg in packages if (lambda: (importlib.import_module(pkg), False)[1] if hasattr(importlib, 'import_module') else True)()]; print('✅ All Python packages installed successfully' if not missing else f'❌ Missing packages: {missing}')"
# =====================================================  
# installing R and all required libraries
# update indices
RUN apt update -qq
# install two helper packages we need
RUN apt install apt-transport-https software-properties-common dirmngr
RUN curl -sSL \
'http://keyserver.ubuntu.com/pks/lookup?op=get&search=0xE298A3A825C0D65DFD57CBB651716619E084DAB9' \
| apt-key add -
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
# add the R 4.0 repo from CRAN -- adjust 'jammy'
RUN add-apt-repository 'deb http://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/'
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y \
    r-base \
    r-base-core \
	r-base-dev \
	r-recommended
# ===================================
### Install repos from CRAN
RUN Rscript -e "options(repos='http://cran.rstudio.com/'); install.packages('devtools', clean=TRUE)"
RUN Rscript -e "options(repos='http://cran.rstudio.com/'); install.packages(c('progress','Rcpp','Rcpp11','RcppAnnoy'), clean=TRUE)"
ARG R_DEPS="c('ggplot2', 'dplyr', 'gtools', 'grid', 'gridtext', \
                'jsonlite', 'kableExtra', 'DT', 'scales','RColorBrewer', \
                'plotly', 'visNetwork', 'ggrepel', 'gtools', 'viridis', \
                'gridExtra', 'tidyr')"	
RUN Rscript -e "options(repos='http://cran.rstudio.com/'); install.packages(${R_DEPS}, clean=TRUE)"
# =====================================================
# Clean up to reduce image size #
RUN apt clean && rm -rf /var/lib/apt/lists/*
# =====================================================
RUN mkdir -p /opt/cellatria
WORKDIR /opt/cellatria
COPY . /opt/cellatria/
RUN chmod +x /opt/cellatria/agent/chatbot.py
RUN ln /opt/cellatria/agent/chatbot.py /usr/local/bin/cellatria
# ===================================
# The VOLUME instruction and the -v option to docker run tell docker to store 
VOLUME /data
WORKDIR /data
# ===================================
# Expose the port used by Gradio
EXPOSE 7860
# ===================================
# Set maintainer
LABEL maintainer="Nima Nouri <ni.nouri@gmail.com>" \
      description=""
LABEL version="1.0.0"
# ===================================
ENV PYTHONPATH=/opt/cellatria/agent
ENV ENV_PATH=/data
CMD ["/usr/local/bin/cellatria"]
# ===================================
# To run your agent:
# docker run -it --rm -p 7860:7860 -v /path/to/data:/data cellatria:latest