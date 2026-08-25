# ProteinDF_bridge

`ProteinDF_bridge` is a Python library and set of command-line tools that
bridge the [ProteinDF](https://proteindf.github.io/) quantum chemistry
engine with common molecular structure formats and other data/packages.

It centers on an in-memory molecular representation (`AtomGroup`, `Atom`,
`Bond`) that can be serialized to/from a compact binary format (msgpack,
referred to as "brd") or YAML, and converted to and from external formats
such as PDB, Amber `prmtop`, GROMACS `.gro`, MOL2, mmCIF, and XYZ.

```{toctree}
:maxdepth: 2
:caption: Contents

installation
usage
api/modules
```

## Indices and tables

- {ref}`genindex`
- {ref}`modindex`
- {ref}`search`
