# Installation

## Requirements

- Python 3.9 or later
- Rust toolchain (`cargo`, `rustc`, latest stable) when compiling from source
- Python dependencies: `numpy`, `pyyaml`, `msgpack`

## Pre-built wheels (Recommended)

Pre-built wheels are available on GitHub Releases for Linux (x86_64, aarch64 with manylinux2014) and macOS (Apple Silicon). No Rust toolchain is required when installing from a pre-built wheel.

To let `pip` automatically select the compatible wheel for your environment:

```bash
pip install proteindf_bridge --find-links https://github.com/ProteinDF/ProteinDF_bridge/releases/expanded_assets/<tag>
```

> **Note**: Be sure to point to `expanded_assets/<tag>` rather than `tag/<tag>`. The standard release tag page on GitHub renders release assets asynchronously, preventing `pip` from finding them, whereas `expanded_assets` serves the direct anchor links that `pip` can scrape.

Alternatively, install a specific wheel directly by URL:

```bash
pip install https://github.com/ProteinDF/ProteinDF_bridge/releases/download/<tag>/proteindf_bridge-<version>-<platform_tag>.whl
```

## From source

Installing from source requires a working Rust compiler (`rustc` and `cargo`, latest stable):

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
