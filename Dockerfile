# ============================================================
# NanoHIV-DR pipeline container
# ============================================================
#
# Target architecture:
#     linux/amd64
#
# Contains:
#
#   Oxford Nanopore
#     - Dorado
#     - minimap2
#     - samtools
#     - seqtk
#     - NanoQ
#     - Medaka
#
#   Host removal
#     - SanitizeMe
#
#   HIV analysis
#     - CodFreq
#     - BCFtools
#
#   Reporting
#     - Python
#     - Bash
#     - NanoHIV-DR REPORT scripts
#
# ============================================================

FROM mambaorg/micromamba:latest

# ============================================================
# Base configuration
# ============================================================

USER root

ENV DEBIAN_FRONTEND=noninteractive
ENV MAMBA_ROOT_PREFIX=/opt/conda

ENV PATH="/opt/conda/bin:/opt/conda/envs/codfreq/bin:/opt/dorado/bin:/opt/codfreq/bin:/opt/REPORT:/opt/REPORT/bin:${PATH}"

# ============================================================
# System packages
# ============================================================

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        bash \
        wget \
        curl \
        git \
        python3 \
        python3-pip \
        python3-dev \
        unzip \
        tar \
        gzip \
        bzip2 \
        xz-utils \
        ca-certificates \
        coreutils \
        procps \
        grep \
        sed \
        gawk \
        findutils \
        build-essential \
        gcc \
        g++ \
        make \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev \
        libcurl4-openssl-dev \
        pkg-config \
    && \
    rm -rf /var/lib/apt/lists/*

# ============================================================
# Bioinformatics software
# ============================================================
#
# These programs are installed in the base micromamba
# environment and therefore are available to Nextflow
# processes using the container.
#
# ============================================================

RUN micromamba install -y \
    -n base \
    -c conda-forge \
    -c bioconda \
        minimap2 \
        samtools \
        seqtk \
        nanoq \
        medaka \
        ivar \
        bcftools \
        pysam \
        sanitizeme \
    && \
    micromamba clean --all --yes
# ============================================================
# Apple Silicon compatibility
# ============================================================
#
# pyabpoa bundled with Medaka triggers SIGILL when the
# linux/amd64 container is run under emulation on Apple Silicon.
# NanoHIV-DR uses medaka_consensus, which functions without
# pyabpoa.
#
# ============================================================

RUN python3 -m pip uninstall -y pyabpoa

# ============================================================
# CodFreq environment
# ============================================================
#
# CodFreq is kept in its own Python environment because it
# has specific Python/build dependencies.
#
# ============================================================

RUN micromamba create -y \
    -n codfreq \
    -c conda-forge \
        python=3.11 \
        pip \
        setuptools \
        wheel \
        cython=0.29.35 \
        gcc \
        gxx \
        make \
        pkg-config \
    && \
    micromamba clean --all --yes

# ============================================================
# Download CodFreq source
# ============================================================

RUN git clone \
    --depth 1 \
    https://github.com/hivdb/codfreq.git \
    /tmp/codfreq

# ============================================================
# CodFreq Python dependencies
# ============================================================
#
# These versions follow the dependency versions used by the
# CodFreq project.
#
# ============================================================

RUN micromamba run -n codfreq \
    pip install \
    --no-cache-dir \
        "pysam==0.21.0" \
        "cutadapt==4.4" \
        "orjson==3.9.1" \
        "click==8.1.3" \
        "dnaio==0.10.0" \
        "more-itertools==9.1.0" \
        "pafpy==0.2.0" \
        "pygments==2.15.1" \
        "pyyaml==6.0.1" \
        "tqdm==4.65.0" \
        "types-setuptools==67.8.0.0" \
        "typing-extensions==4.7.1" \
        "urllib3==2.0.4" \
        "xopen==1.7.0"

# ============================================================
# post-align
# ============================================================

RUN micromamba run -n codfreq \
    pip install \
    --no-cache-dir \
    "https://github.com/hivdb/post-align/archive/cb28b83e2d9639960533f805f7dc2612ea63ddc6.zip"

# ============================================================
# Install CodFreq
# ============================================================
#
# Do NOT assume a particular executable name such as
# sam2codfreq.
#
# The package installs its own command-line scripts.
#
# ============================================================

RUN cd /tmp/codfreq && \
    micromamba run -n codfreq \
    pip install \
    --no-cache-dir \
    --ignore-installed \
    .

# ============================================================
# CodFreq executable discovery
# ============================================================
#
# Show exactly which CodFreq commands were installed.
# This avoids the previous build failure caused by assuming
# that sam2codfreq exists.
#
# ============================================================

RUN echo "============================================================" && \
    echo "CodFreq installed executables" && \
    echo "============================================================" && \
    ls -lah /opt/conda/envs/codfreq/bin/ && \
    echo "" && \
    find /opt/conda/envs/codfreq/bin \
        -maxdepth 1 \
        -type f \
        -perm /111 \
        -printf '%f\n' \
        | sort

# ============================================================
# Verify CodFreq Python package
# ============================================================

RUN micromamba run -n codfreq \
    python -c \
    "import codfreq; print('CodFreq Python package:', codfreq.__file__)"

# ============================================================
# Verify CodFreq command-line programs
# ============================================================
#
# fastq2codfreq is the command expected for the FASTQ branch.
#
# We test for it without hard-coding sam2codfreq.
#
# ============================================================

RUN if [ -x /opt/conda/envs/codfreq/bin/fastq2codfreq ]; then \
        echo "fastq2codfreq found"; \
        /opt/conda/envs/codfreq/bin/fastq2codfreq --help >/dev/null; \
    else \
        echo "ERROR: fastq2codfreq was not installed"; \
        echo "Installed CodFreq commands:"; \
        find /opt/conda/envs/codfreq/bin \
            -maxdepth 1 \
            -type f \
            -perm /111 \
            -printf '%f\n' \
            | sort; \
        exit 1; \
    fi

# ============================================================
# Make CodFreq executable available in PATH
# ============================================================

RUN mkdir -p /opt/codfreq/bin && \
    ln -sf \
        /opt/conda/envs/codfreq/bin/fastq2codfreq \
        /opt/codfreq/bin/fastq2codfreq

# ============================================================
# Remove CodFreq source
# ============================================================

RUN rm -rf /tmp/codfreq

# ============================================================
# Additional CodFreq runtime tools
# ============================================================
#
# The current CodFreq workflow uses alignment/processing
# utilities including minimap2, samtools and fastp.
#
# minimap2 and samtools are already available in the base
# environment.
#
# Install fastp into the base environment.
#
# ============================================================

RUN micromamba install -y \
    -n base \
    -c conda-forge \
    -c bioconda \
        fastp \
    && \
    micromamba clean --all --yes

# ============================================================
# Dorado (downloaded at build time; pinned for reproducibility)
# ============================================================
#
# The official Linux x64 Dorado tarball bundles the CUDA
# libraries and also runs on CPU (there is no separate CPU-only
# build). Set DORADO_VERSION to the version you validated
# against; check your local version with:  dorado --version
#
# Override at build time, e.g.:
#   docker buildx build --build-arg DORADO_VERSION=2.0.0 (or another version) ...
#
ARG DORADO_VERSION=2.1.0

RUN set -eux; \
    url="https://cdn.oxfordnanoportal.com/software/analysis/dorado-${DORADO_VERSION}-linux-x64.tar.gz"; \
    wget -q "$url" -O /tmp/dorado.tar.gz; \
    mkdir -p /opt/dorado; \
    tar -xzf /tmp/dorado.tar.gz -C /opt/dorado --strip-components=1; \
    rm -f /tmp/dorado.tar.gz; \
    chmod +x /opt/dorado/bin/dorado

# ============================================================
# NanoHIV-DR reporting software
# ============================================================

COPY REPORT /opt/REPORT

RUN chmod +x /opt/REPORT/preprocessing.sh && \
    if [ -d /opt/REPORT/bin ]; then \
        find /opt/REPORT/bin \
            -type f \
            -name "*.py" \
            -exec chmod +x {} \; ; \
    fi

# ============================================================
# REPORT Python environment
# ============================================================
#
# The legacy SierraPy client pins dependencies that are not
# compatible with Python 3.12. Keep reporting isolated in a
# dedicated Python 3.10 environment.
#
# ============================================================

RUN micromamba create -y \
    -n report \
    -c conda-forge \
    python=3.10 \
    pip \
    && \
    micromamba clean --all --yes

RUN if [ -f /opt/REPORT/requirements.txt ]; then \
        micromamba run -n report python -m pip install \
            --no-cache-dir \
            -r /opt/REPORT/requirements.txt; \
    fi

# ============================================================
# Final container checks
# ============================================================

RUN echo "============================================================" && \
    echo " NanoHIV-DR container verification" && \
    echo "============================================================" && \
    echo "" && \
    echo "System Python:" && \
    python3 --version && \
    echo "" && \
    echo "CodFreq Python:" && \
    /opt/conda/envs/codfreq/bin/python --version && \
    echo "" && \
    echo "CodFreq package:" && \
    /opt/conda/envs/codfreq/bin/python \
        -c "import codfreq; print(codfreq.__file__)" && \
    echo "" && \
    echo "fastq2codfreq:" && \
    /opt/conda/envs/codfreq/bin/fastq2codfreq --help >/dev/null && \
    echo "fastq2codfreq OK" && \
    echo "" && \
    echo "Dorado:" && \
    /opt/dorado/bin/dorado --version && \
    echo "" && \
    echo "Minimap2:" && \
    micromamba run -n base minimap2 --version && \
    echo "" && \
    echo "Samtools:" && \
    micromamba run -n base samtools --version | head -n 1 && \
    echo "" && \
    echo "Seqtk:" && \
    micromamba run -n base seqtk 2>&1 || true && \
    echo "" && \
    echo "NanoQ:" && \
    micromamba run -n base nanoq --version && \
    echo "" && \
    echo "Medaka:" && \
    micromamba run -n base medaka --version && \
    micromamba run -n base python3 -c "import medaka; from medaka import prediction; print('Medaka prediction import OK')" && \
    echo "" && \
    echo "BCFtools:" && \
    micromamba run -n base bcftools --version | head -n 1 && \
    echo "" && \
    echo "SanitizeMe:" && \
    micromamba run -n base \
        SanitizeMe_CLI.py -h >/dev/null && \
    echo "SanitizeMe_CLI.py OK" && \
    echo "" && \
    echo "fastp:" && \
    micromamba run -n base fastp --version && \
    echo "" && \
    echo "REPORT:" && \
    test -x /opt/REPORT/preprocessing.sh && \
    echo "/opt/REPORT/preprocessing.sh OK" && \
    echo "" && \
    echo "REPORT Python:" && \
    micromamba run -n report python --version && \
    echo "" && \
    echo "SierraPy:" && \
    micromamba run -n report sierrapy --help >/dev/null && \
    echo "SierraPy OK" && \
    echo "============================================================"

# ============================================================
# Working directory
# ============================================================

WORKDIR /data

# ============================================================
# Default command
# ============================================================

CMD ["bash"]
