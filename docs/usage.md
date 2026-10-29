# Command-line tools

`ProteinDF_bridge` installs a set of small, single-purpose CLI scripts as
entry-point console scripts (placed in your Python environment's `bin/` directory,
accessible via `PATH`). The underlying implementations reside in
`proteindf_bridge.cli`. Each command is a thin wrapper around `proteindf_bridge`
that reads one format and writes another, or performs one structure-editing
operation on a bridge (`.brd`) file. Run any of them with `-h` for full options.

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
| `brd-ramachandran.py` | calculate Ramachandran angles (phi/psi) |

## Development

| Script | Description |
| --- | --- |
| `doctest_runner.py` | run doctests across `proteindf_bridge` |
| `module_inspect.py` | inspect module contents |

## Python: checking missing atoms

You can inspect an `AtomGroup` hierarchy for missing atoms against standard
Chemical Component Dictionary (CCD) residue definitions using
`AtomGroup.check_missing_atoms()`. This inspection is non-destructive and does
not modify the tree. By default, only heavy atoms are verified; pass
`include_hydrogens=True` to inspect hydrogen atoms as well.

### Basic usage

```python
import proteindf_bridge as pb

pdb = pb.Pdb("proteindf_bridge/data/1hls.pdb")
ag = pdb.get_atomgroup()

# Inspect heavy atoms
report = ag.check_missing_atoms()
print(report.is_complete)    # True
print(report.total_missing)  # 0
print(str(report))           # MissingAtomsReport: complete (0 missing atoms)

# Inspect individual incomplete residues if any
for item in report.incomplete:
    print(f"{item.path} ({item.comp_id}): missing {item.missing}")
```

### Components without CCD templates

If the structure contains residues or ligands not found in the CCD template
database, their check status is set to `"no_template"`. These components are
summarized in `report.no_template_components` (a mapping from component ID to
count). To provide custom component definitions, load the built-in database,
register your definitions with `add_from_file()`, and pass `db` to
`check_missing_atoms()`:

```python
import proteindf_bridge.rs as rs

db = rs.CcdTemplateDb.builtin()
db.add_from_file("custom_components.cif")
report = ag.check_missing_atoms(db=db)
```

### Optional and unexpected atoms

Certain atoms are treated as optional and do not cause a residue to be marked
incomplete if absent:
- C-terminal `OXT` for amino acids
- 5'-terminal phosphate atoms (`P`, `OP1`, `OP2`, `OP3`) for nucleic acids
- Thiol hydrogen `HG` on cysteine residues participating in disulfide bonds (in hydrogen mode)

Atoms present in the residue that are not defined in the template are reported
in `ResidueAtomCheck.unexpected`. If all required template atoms are present, the
residue status remains `"complete"` even when unexpected atoms exist.

### Hydrogen inspection caveats

When checking hydrogens (`include_hydrogens=True`), false positives may occur
due to physiological or non-standard protonation states (such as neutral vs.
charged forms of HIS, ASP, or GLU). In such cases, hydrogen atoms expected by the
standard neutral template may be absent from the structure.

## Python: completing missing atoms

You can reconstruct and add missing **heavy atoms** to incomplete residues in an
`AtomGroup` hierarchy using `AtomGroup.complete_missing_atoms()`. Missing heavy
atoms are generated from ideal coordinates in the Chemical Component Dictionary
(CCD) residue template and inserted into the residue. Existing atoms are never
moved or removed. Hydrogen atoms are not reconstructed by this method, and
entirely missing residues (where no anchor atoms exist) are not handled.

### Basic usage

The recommended workflow is to inspect for missing atoms with
`check_missing_atoms()`, reconstruct missing heavy atoms with
`complete_missing_atoms()`, and finally add missing hydrogens with
`add_missing_hydrogens()`:

```python
import proteindf_bridge as pb

pdb = pb.Pdb("proteindf_bridge/data/1hls.pdb")
ag = pdb.get_atomgroup()

# Simulate missing heavy atoms
ag["model_1"]["A"]["1"].del_atom("4_O")
ag["model_1"]["A"]["3"].del_atom("35_CG2")

# 1. Inspect missing atoms before completion
check_before = ag.check_missing_atoms()
print(check_before.is_complete)    # False
print(check_before.total_missing)  # 2

# 2. Reconstruct missing heavy atoms
report = ag.complete_missing_atoms()
print(report.total_added)          # 2
print(len(report.completed))       # 2
print(len(report.skipped))         # 0
print(str(report))
# CompletionReport: 2 heavy atom(s) added across 2 completed residue(s), 0 skipped residue(s)

# 3. Add missing hydrogens after heavy atoms are complete
ag.add_missing_hydrogens()

# 4. Verify that the structure is now complete
check_after = ag.check_missing_atoms()
print(check_after.is_complete)     # True
```

### Completion statuses

Each inspected residue in `report.residues` is assigned one of the following
statuses:
- `"completed"`: Missing heavy atoms were successfully reconstructed and added.
- `"unchanged"`: The residue already contains all required heavy atoms.
- `"no_template"`: The component ID was not found in the CCD template database.
- `"insufficient_anchors"`: Fewer than 3 matching heavy atoms are available for
  superposition, or existing anchor atoms are collinear.
- `"poor_fit"`: Anchor superposition RMSD exceeded `max_rmsd` (default: 1.0 Å),
  indicating possible atom misnaming or incorrect element assignment.

Residues that could not be completed are collected in `report.skipped`.

### Limitations and accuracy

Reconstruction superimposes ideal coordinates from CCD templates onto existing
anchor atoms. When distal side-chain atoms are missing, their torsional angles
cannot be determined from the remaining anchors and may differ from the true
conformation. This method does not perform protein structure prediction; its goal
is to produce a complete set of atoms that satisfies structural inspection.

### Order of operations

Always run `complete_missing_atoms()` before adding hydrogen atoms
(`add_missing_hydrogens()`). Reconstructing heavy atoms first ensures that
subsequent hydrogen placement can use the complete heavy-atom framework to position
all attached hydrogens correctly.

## Python: adding missing hydrogens

Hydrogen atoms can be added to an `AtomGroup` using `add_missing_hydrogens()` based
on ideal coordinates from CCD residue templates:

```python
# Add missing hydrogens using automatic histidine tautomer detection (default)
report = ag.add_missing_hydrogens()

# Or specify the histidine protonation mode explicitly
report = ag.add_missing_hydrogens(his="hip")
```

### Histidine protonation handling

The `his` parameter controls how histidine residues are protonated:

- `"auto"` (default): Automatically determines the protonation state (HID, HIE, HIP,
  or metal-coordinated) based on metal coordination, residue name, existing hydrogens,
  and salt bridges (in priority order).
- `"hip"`: Traditional full protonation behavior where both `HD1` and `HE2` are added
  to every histidine residue.

Any value other than `"auto"` or `"hip"` raises a `ValueError`.

