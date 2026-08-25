# Command-line tools

`ProteinDF_bridge` installs a set of small, single-purpose CLI scripts under
`scripts/`. Each one is a thin wrapper around `proteindf_bridge` that reads
one format and writes another, or performs one structure-editing operation
on a bridge (`.brd`) file. Run any of them with `-h` for full options.

## Format conversion

| Script | Description |
| --- | --- |
| `pdb2brd.py` | translate from PDB to bridge file |
| `brd2pdb.py` | transform bridge file to PDB file |
| `gro2brd.py` | translate from gro to bridge file |
| `brd2gro.py` | output gro format from bridge file |
| `xyz2brd.py` | transform XYZ file to bridge file |
| `brd2xyz.py` | transform bridge file to XYZ file |
| `mmcif2txt.py` | parse mmCIF file |
| `mmcif2mol2.py` | parse mmCIF file to mol2 file |
| `read_amber_prmtop.py` | parse Amber prmtop |
| `mpac2yml.py` / `yml2mpac.py` | convert between MsgPack and YAML bridge representations |
| `mpac2txt.py` | display file formatted by MsgPack using YAML |
| `brd2txt.py` | print molecular bridge file |
| `db2txt.py` | print DB (sqlite3) file |

## Structure editing

| Script | Description |
| --- | --- |
| `brd-select.py` | bridge file selector |
| `brd-select-path.py` | bridge file selector (by path) |
| `brd-divide.py` | bridge file divider |
| `brd-divide-mainchain.py` | bridge file divider (main chain) |
| `brd-restructure.py` | restructure brd file by reference file |
| `brd-renumber-resid.py` | renumber resid in bridge file |
| `brd-setup-bond.py` | setup bonds |
| `brd-show-bonds.py` | show bonds |
| `brd-show-res.py` | print residues in the bridge file |
| `remove_wat.py` | remove water molecules in bridge file |
| `reorder.py` | reorder protein |
| `neutralize.py` | neutralize protein |
| `crystallize.py` | crystallize molecules |
| `superposer.py` | superpose |

## Analysis

| Script | Description |
| --- | --- |
| `brd-box.py` | calc box size |
| `brd-density.py` | calc density |
| `brd-formula.py` | print molecular formula |

## Development

| Script | Description |
| --- | --- |
| `doctest_runner.py` | run doctests across `proteindf_bridge` |
| `module_inspect.py` | inspect module contents |
