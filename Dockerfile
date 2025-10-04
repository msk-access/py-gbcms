# syntax=docker/dockerfile:1
FROM python:3.11-slim

# Install OS-level build dependencies commonly needed for cyvcf2, pysam, numpy, etc.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libsqlite3-dev \
    libgdbm-dev \
    libreadline-dev \
    libffi-dev \
    libxml2-dev \
    libxslt1-dev \
    gfortran \
    libopenblas-dev \
    liblapack-dev \
    git \
    autoconf \
    && rm -rf /var/lib/apt/lists/*

# Install uv (optional helper) — you can remove if you prefer pip directly
RUN pip install --no-cache-dir uv

WORKDIR /app

# Copy only the files needed for installation first to leverage Docker cache
COPY pyproject.toml pyproject.lock* README.md LICENSE* /app/
COPY src/ /app/src/

# Create a virtualenv, activate it and install the package with extras
RUN uv venv .venv && \
    /bin/bash -lc "source .venv/bin/activate && uv pip install --no-cache-dir '.[all]'" 

# Ensure the venv bin is first in PATH
ENV PATH="/app/.venv/bin:${PATH}"

# Working directory for running
WORKDIR /data

# Entrypoint/command defaults
ENTRYPOINT ["gbcms"]
CMD ["--help"]

LABEL maintainer="MSK-ACCESS <access@mskcc.org>"
LABEL description="Python implementation of GetBaseCountsMultiSample (gbcms) for calculating base counts in BAM files"
LABEL org.opencontainers.image.source="https://github.com/msk-access/py-gbcms"
LABEL org.opencontainers.image.documentation="https://github.com/msk-access/py-gbcms/blob/main/README.md"
LABEL org.opencontainers.image.licenses="AGPL-3.0"
LABEL org.opencontainers.image.base.image="python:3.11-slim"