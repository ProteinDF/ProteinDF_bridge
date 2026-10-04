# ProteinDF_bridge -- bridge scripts between ProteinDF and molecular data formats

ProteinDF_bridge is a Python library for reading, manipulating, and converting various molecular structure formats (PDB, mmCIF, MOL2, GRO, Amber PRMTOP, etc.) for use with ProteinDF.

## Requirements

- Python 3.9 or later
- Rust toolchain (`cargo`, `rustc`, latest stable) when installing from source
- Python packages:
  - `numpy`
  - `pyyaml`
  - `msgpack`

## Installation

### Pre-built Wheels (Recommended)

Pre-built wheels are available on GitHub Releases for Linux (x86_64, aarch64 with manylinux2014) and macOS (Apple Silicon). No Rust toolchain is required when installing from a wheel.

To let `pip` automatically select the matching wheel for your platform:

```bash
pip install proteindf_bridge --find-links https://github.com/ProteinDF/ProteinDF_bridge/releases/expanded_assets/<tag>
```

> **Note**: Use the `expanded_assets/<tag>` endpoint instead of `releases/tag/<tag>`. GitHub Releases loads assets dynamically on tag pages which prevents `pip` from scraping link targets, whereas `expanded_assets` provides the direct HTML link list that `pip` can parse.

Alternatively, you can specify the direct URL to the appropriate wheel file:

```bash
pip install https://github.com/ProteinDF/ProteinDF_bridge/releases/download/<tag>/proteindf_bridge-<version>-<platform_tag>.whl
```

### From Source

Installing from source compiles the bundled Rust extension (`proteindf_bridge.rs`) and requires a working Rust toolchain (`rustc` and `cargo`, latest stable).

```bash
git clone https://github.com/ProteinDF/ProteinDF_bridge.git
cd ProteinDF_bridge
pip install .
```

For development (editable install):

```bash
pip install -e .
```

## Quick Start

```python
from proteindf_bridge import BioPdb, AtomGroup
import proteindf_bridge.rs as rs

# Load structure from PDB
pdb = BioPdb()
pdb.load("protein.pdb")

# Get AtomGroup representation
ag = pdb.get_atomgroup()
print(f"Total atoms: {ag.get_number_of_all_atoms()}")
```

## Running Tests

```bash
# Python test suite
python -m unittest discover -s tests

# Rust crate test suite
cargo test --manifest-path rust/Cargo.toml
```

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for the branching model (GitFlow) and release process.

## License

ProteinDF_bridge is licensed under the GNU General Public License v3.0 (GPLv3).
