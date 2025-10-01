# Alternative Dockerfile using Ubuntu base
# This is larger but may be more familiar to some users

# Multi-stage build for py-gbcms
FROM ubuntu:22.04 as builder

# Prevent interactive prompts during build
ENV DEBIAN_FRONTEND=noninteractive

# Install Python and build dependencies
RUN apt-get update && apt-get install -y \
    python3.11 \
    python3.11-dev \
    python3-pip \
    build-essential \
    gcc \
    g++ \
    make \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libhts-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Create symlinks for python/pip
RUN ln -s /usr/bin/python3.11 /usr/bin/python && \
    ln -s /usr/bin/pip3 /usr/bin/pip

# Install uv
RUN pip install --no-cache-dir uv

# Set working directory
WORKDIR /build

# Copy project files
COPY pyproject.toml README.md LICENSE ./
COPY src/ ./src/

# Install the package with all optional dependencies
# This includes cyvcf2 (fast VCF parsing) and Ray (distributed computing)
RUN uv pip install --system --no-cache ".[all]"

# Final stage
FROM ubuntu:22.04

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install Python and runtime dependencies
RUN apt-get update && apt-get install -y \
    python3.11 \
    python3-pip \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    libcurl4 \
    libssl3 \
    libhts3 \
    samtools \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Create symlinks
RUN ln -s /usr/bin/python3.11 /usr/bin/python

# Copy installed packages from builder
COPY --from=builder /usr/local/lib/python3.11/dist-packages /usr/local/lib/python3.11/dist-packages
COPY --from=builder /usr/local/bin/gbcms /usr/local/bin/gbcms
COPY --from=builder /usr/local/bin/getbasecounts /usr/local/bin/getbasecounts

# Create working directory
WORKDIR /data

# Set entrypoint
ENTRYPOINT ["gbcms"]

# Default command (show help)
CMD ["--help"]

# Verify installation
RUN gbcms version

# Metadata
LABEL maintainer="MSK-ACCESS <access@mskcc.org>"
LABEL description="Python implementation of GetBaseCountsMultiSample (gbcms) for calculating base counts in BAM files"
LABEL version="2.0.0"
LABEL org.opencontainers.image.source="https://github.com/msk-access/getbasecounts"
LABEL org.opencontainers.image.documentation="https://github.com/msk-access/getbasecounts/blob/main/README.md"
LABEL org.opencontainers.image.licenses="Apache-2.0"
LABEL base.image="ubuntu:22.04"
