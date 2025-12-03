# Dockerfile for Mechanism-First Causal Graphs Pipeline
# Ensures reproducibility for Nature-tier publication

FROM continuumio/miniconda3:23.5.2-0

LABEL maintainer="your.email@institution.edu"
LABEL description="Mechanism-First Causal Graphs for Noncoding GWAS"
LABEL version="1.0.0"

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    git \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    wget \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    tabix \
    bcftools \
    && rm -rf /var/lib/apt/lists/*

# Create working directory
WORKDIR /app

# Copy environment file
COPY environment.yml /app/environment.yml

# Create conda environment
RUN conda env create -f environment.yml && \
    conda clean -afy

# Activate environment
SHELL ["conda", "run", "-n", "mechanism-gwas", "/bin/bash", "-c"]

# Install R packages (SuSiE, COLOC) within conda environment
RUN R -e "install.packages(c('susieR', 'coloc', 'data.table', 'dplyr', 'tidyr'), repos='https://cloud.r-project.org/')"

# Copy source code
COPY src/ /app/src/
COPY scripts/ /app/scripts/
COPY workflow/ /app/workflow/
COPY config/ /app/config/

# Install package in development mode
COPY setup.py pyproject.toml /app/
RUN pip install -e .

# Create data and results directories
RUN mkdir -p /app/data /app/results /app/logs

# Set entrypoint
ENTRYPOINT ["conda", "run", "-n", "mechanism-gwas"]
CMD ["snakemake", "--help"]
