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
python -m gbcms.cli --help
```

You should see the help message for the `gbcms` command.
