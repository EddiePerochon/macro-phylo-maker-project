FROM rocker/r-ver:4.4.1

ENV DOCKER_API_VERSION=1.43
ENV DEBIAN_FRONTEND=noninteractive
ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=true

# System dependencies for R packages, plotting, devtools, ChronoSTA Python support,
# and Docker CLI access from inside the container.
RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    git \
    curl \
    wget \
    ca-certificates \
    build-essential \
    gfortran \
    make \
    cmake \
    pkg-config \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libcairo2-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff5-dev \
    libwebp-dev \
    libuv1-dev \
    libglpk-dev \
    libgmp3-dev \
    libmpfr-dev \
    libgsl-dev \
    python3 \
    python3-pip \
    python3-venv \
    docker.io \
    && rm -rf /var/lib/apt/lists/*

# R packages used across the MacroPhyloMaker workflow.
# The README already lists devtools, ape, phytools, here, readr and additional
# workflow packages such as phangorn, stringr, progress, igraph, and MonoPhy.
RUN Rscript -e 'install.packages(c( \
    "devtools", \
    "ape", \
    "phytools", \
    "phangorn", \
    "stringr", \
    "progress", \
    "igraph", \
    "MonoPhy", \
    "here", \
    "readr", \
    "data.table", \
    "dplyr", \
    "tidyr", \
    "tibble", \
    "ggplot2", \
    "remotes" \
  ), repos = "https://cloud.r-project.org", Ncpus = parallel::detectCores())'

# Python packages needed for the ChronoSTA-enabled step.
# Install them into an isolated virtual environment to avoid system-Python issues.
RUN python3 -m venv /opt/chronosta-venv && \
    /opt/chronosta-venv/bin/python -m pip install --upgrade pip setuptools wheel && \
    /opt/chronosta-venv/bin/python -m pip install --no-cache-dir \
      biopython \
      pandas \
      numpy \
      scipy \
      matplotlib

ENV PATH="/opt/chronosta-venv/bin:${PATH}"
ENV CHRONOSTA_PYTHON="/opt/chronosta-venv/bin/python"

WORKDIR /work

CMD ["R"]
