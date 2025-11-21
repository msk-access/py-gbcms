# Installation

gbcms requires a Python environment and the Rust toolchain to build the core engine.

## Prerequisites

*   **Python**: 3.10 or higher
*   **Rust**: Latest stable version (install via [rustup](https://rustup.rs/))

## Installing from Source

1.  **Clone the repository**:
    ```bash
    git clone https://github.com/msk-access/py-gbcms.git
    cd py-gbcms
    ```

2.  **Create a virtual environment** (recommended):
    ```bash
    python -m venv .venv
    source .venv/bin/activate
    ```

3.  **Install dependencies and build**:
    This command will install Python dependencies and compile the Rust extension (`gbcms_rs`).
    ```bash
    pip install .
    ```

    *Note: If you are developing, use `pip install -e .` for editable mode.*

## Verifying Installation

Run the help command to verify that the CLI is accessible:

```bash
gbcms --help
```

You should see the help message for the `gbcms` command.

## Using Docker

Docker is the recommended way to run gbcms to ensure a consistent environment.

### 1. Build the Image

```bash
docker build -t gbcms:latest .
```

### 2. Run with Docker

When running with Docker, you need to mount your data directories so the container can access them.

```bash
docker run --rm \
    -v /path/to/data:/data \
    gbcms:latest \
    run \
    --variants /data/variants.vcf \
    --fasta /data/reference.fa \
    --bam /data/sample.bam \
    --output-dir /data/output
```

**Note**: All paths passed to `gbcms` must be paths *inside* the container (e.g., `/data/...`).

