# Installation

## Requirements

- Python 3.9 or later
- Rust toolchain (`cargo`, `rustc` 1.80+) when compiling from source
- Python dependencies: `numpy`, `pyyaml`, `msgpack`

## From source

Installing from source requires a working Rust compiler (`rustc` and `cargo`):

```bash
git clone https://github.com/ProteinDF/ProteinDF_bridge.git
cd ProteinDF_bridge
python -m venv .venv
source .venv/bin/activate
pip install .
```

For development (editable mode):

```bash
pip install -e .
```

## Pre-built wheels

Pre-compiled wheels for Linux (x86_64, aarch64) and macOS (Apple Silicon) will be provided via GitHub Releases (planned in PR#50).
