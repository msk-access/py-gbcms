# Dockerfile for py-gbcms
# 
# Build: docker build -t py-gbcms .
# Run:   docker run py-gbcms --help

# Stage 1: Builder (with Rust compilation)
FROM python:3.11-bookworm AS builder

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

WORKDIR /app

# Install system dependencies for Rust compilation
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    pkg-config \
    git \
    cmake \
    clang \
    libclang-dev \
    llvm-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

ENV LIBCLANG_PATH="/usr/lib/llvm-14/lib"
ENV BINDGEN_EXTRA_CLANG_ARGS="-I/usr/lib/llvm-14/lib/clang/14.0.6/include"

# Install Rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Install maturin
RUN pip install --no-cache-dir maturin

# Copy source
COPY . /app

# Build unified wheel with maturin (includes both Python and Rust)
RUN maturin build --release --out /app/dist

# Stage 2: Runtime (slim image)
FROM python:3.11-slim-bookworm

WORKDIR /app

# Install runtime dependencies for htslib
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4 \
    libssl3 \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    && rm -rf /var/lib/apt/lists/*

# Copy and install the unified wheel
COPY --from=builder /app/dist/*.whl /app/dist/
RUN pip install --no-cache-dir /app/dist/*.whl

# Verify installation
RUN python -c "from gbcms import _rs; import gbcms; print(f'gbcms {gbcms.__version__} ready')"

ENTRYPOINT ["gbcms"]
CMD ["--help"]