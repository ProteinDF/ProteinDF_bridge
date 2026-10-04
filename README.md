# ProteinDF_bridge -- bridge scripts between ProteinDF and molecular data formats

ProteinDF_bridge is a Python library for reading, manipulating, and converting various molecular structure formats (PDB, mmCIF, MOL2, GRO, Amber PRMTOP, etc.) for use with ProteinDF.

## Requirements

- Python 3.9 or later
- Rust toolchain (`cargo`, `rustc` 1.80+) when installing from source
- Python packages:
  - `numpy`
  - `pyyaml`
  - `msgpack`

## Installation

### From Source

Installing from source compiles the bundled Rust extension (`proteindf_bridge.rs`) and requires a working Rust toolchain. Pre-built wheels for common platforms will be provided on GitHub Releases (PR#50).

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
