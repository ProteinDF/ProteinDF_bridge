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

> **Note**: `expanded_assets/<tag>` is an undocumented (non-public) GitHub web endpoint that returns the HTML fragment of release assets, which `pip` can scrape (unlike `releases/tag/<tag>` which loads assets asynchronously). Since this is not an officially documented feature of GitHub, installing via the direct wheel URL below is also supported and provided as a reliable alternative.

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
