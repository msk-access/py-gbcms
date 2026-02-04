# Dockerfile for py-gbcms
# 
# For development/local builds (includes Rust compilation):
#   docker build -t py-gbcms .
#
# For CI builds (uses pre-built wheel - faster):
#   docker build --build-arg WHEEL_PATH=dist/*.whl -t py-gbcms .

ARG WHEEL_PATH=""

# Stage 1: Builder (only used if no pre-built wheel provided)
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
    libssl-dev

ENV LIBCLANG_PATH="/usr/lib/llvm-14/lib"
ENV BINDGEN_EXTRA_CLANG_ARGS="-I/usr/lib/llvm-14/lib/clang/14.0.6/include"

# Install Rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Install build tools
RUN pip install --no-cache-dir maturin build

# Copy source
COPY . /app

# Build gbcms_rs wheel
WORKDIR /app/rust
RUN maturin build --release --out /app/dist

# Stage 2: Runtime
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

# Copy wheels from builder or use pre-built
COPY --from=builder /app/dist /app/dist
COPY . /app

# Install Rust extension wheel first, then Python package
RUN pip install --no-cache-dir /app/dist/*.whl
RUN pip install --no-cache-dir .

# Verify installation
RUN python -c "import gbcms; import gbcms_rs; print(f'gbcms {gbcms.__version__} ready')"

ENTRYPOINT ["gbcms"]
CMD ["--help"]