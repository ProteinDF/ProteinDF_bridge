# ProteinDF_bridge -- bridge scripts between ProteinDF and molecular data formats

ProteinDF_bridge is a Python library for reading, manipulating, and converting various molecular structure formats (PDB, mmCIF, MOL2, GRO, Amber PRMTOP, etc.) for use with ProteinDF.

## Requirements

- Python 3.8 or later
- Python packages:
  - `numpy`
  - `pyyaml`
  - `msgpack`

## Installation

### Clone repository

```bash
git clone https://github.com/ProteinDF/ProteinDF_bridge.git
cd ProteinDF_bridge
```

### Install using pip

```bash
pip install .
```

For development (editable install):

```bash
pip install -e .
```

## Quick Start

```python
from proteindf_bridge import BioPdb, AtomGroup

# Load structure from PDB
pdb = BioPdb()
pdb.load("protein.pdb")

# Get AtomGroup representation
ag = pdb.get_atomgroup()
print(f"Total atoms: {ag.get_number_of_all_atoms()}")
```

## Running Tests

```bash
pytest
```

## License

ProteinDF_bridge is licensed under the GNU General Public License v3.0 (GPLv3).
