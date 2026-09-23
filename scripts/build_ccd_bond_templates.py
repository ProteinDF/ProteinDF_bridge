#!/usr/bin/env python3
# SPDX-FileCopyrightText: The ProteinDF development team
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Extracts standard amino acid, nucleic acid, and water bond+geometry templates from
wwPDB CCD and packages them into a MessagePack binary database for proteindf-bridge.

Each atom entry includes its element symbol and idealized 3D coordinates
(falling back to model coordinates when idealized ones are absent), used both for
bond-order resolution (`AtomGroup::apply_ccd_bond_templates`) and for CCD-reference
hydrogen addition (`RUST_PORT_SPEC.md` §3.16).

Target components:
- 20 standard amino acids: ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE,
                           LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL
- 8 standard nucleic acids: DA, DC, DG, DT, A, C, G, U
- Water: HOH
"""

import os
import re
import sys
import urllib.request
from pathlib import Path

try:
    import msgpack
except ImportError:
    print("msgpack not found, please run with .venv/bin/python", file=sys.stderr)
    sys.exit(1)

COMPONENTS = [
    # 20 standard amino acids
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
    # 8 standard nucleic acids (DNA & RNA)
    "DA", "DC", "DG", "DT",
    "A", "C", "G", "U",
    # Water
    "HOH",
]

RCSB_URL_TEMPLATE = "https://files.rcsb.org/ligands/view/{comp_id}.cif"

def fetch_ccd_cif(comp_id: str) -> str:
    url = RCSB_URL_TEMPLATE.format(comp_id=comp_id)
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "ProteinDF-Bridge-CCD-Builder/1.0"}
    )
    with urllib.request.urlopen(req) as resp:
        return resp.read().decode("utf-8")

def parse_ccd_cif(cif_text: str, comp_id: str) -> dict:
    lines = cif_text.splitlines()
    atoms = []
    bonds = []

    # Map mmCIF value_order to integer bond order (consistent with format/mmcif.rs)
    order_map = {
        "SING": 1,
        "DOUB": 2,
        "TRIP": 3,
        "QUAD": 4,
        "AROM": 1,
    }

    i = 0
    n = len(lines)
    while i < n:
        line = lines[i].strip()
        if line == "loop_":
            i += 1
            headers = []
            while i < n and lines[i].strip().startswith("_"):
                headers.append(lines[i].strip())
                i += 1
            
            # Check which loop block this is
            is_atom_loop = any(h.startswith("_chem_comp_atom.") for h in headers)
            is_bond_loop = any(h.startswith("_chem_comp_bond.") for h in headers)

            rows = []
            while i < n:
                row_line = lines[i].strip()
                if not row_line or row_line.startswith("#"):
                    i += 1
                    continue
                if row_line == "loop_" or row_line.startswith("_") or row_line.startswith("data_"):
                    break
                # Parse whitespace-separated tokens (handling quoted tokens)
                tokens = []
                for m in re.finditer(r'("[^"]*"|\'[^\']*\'|\S+)', row_line):
                    tok = m.group(1)
                    if (tok.startswith('"') and tok.endswith('"')) or (tok.startswith("'") and tok.endswith("'")):
                        tok = tok[1:-1]
                    tokens.append(tok)
                if len(tokens) == len(headers):
                    rows.append(dict(zip(headers, tokens)))
                i += 1

            if is_atom_loop:
                atom_id_col = "_chem_comp_atom.atom_id"
                type_symbol_col = "_chem_comp_atom.type_symbol"
                ideal_cols = (
                    "_chem_comp_atom.pdbx_model_Cartn_x_ideal",
                    "_chem_comp_atom.pdbx_model_Cartn_y_ideal",
                    "_chem_comp_atom.pdbx_model_Cartn_z_ideal",
                )
                model_cols = (
                    "_chem_comp_atom.model_Cartn_x",
                    "_chem_comp_atom.model_Cartn_y",
                    "_chem_comp_atom.model_Cartn_z",
                )

                def parse_xyz(r, cols):
                    values = [r.get(c) for c in cols]
                    if any(v is None or v in ("?", ".") for v in values):
                        return None
                    try:
                        return [float(v) for v in values]
                    except ValueError:
                        return None

                for r in rows:
                    if atom_id_col not in r:
                        continue
                    element = r.get(type_symbol_col, "X")
                    if element == "D":
                        element = "H"
                    ideal_xyz = parse_xyz(r, ideal_cols) or parse_xyz(r, model_cols)
                    atoms.append({
                        "name": r[atom_id_col],
                        "element": element,
                        "ideal_xyz": ideal_xyz,
                    })

            elif is_bond_loop:
                a1_col = "_chem_comp_bond.atom_id_1"
                a2_col = "_chem_comp_bond.atom_id_2"
                order_col = "_chem_comp_bond.value_order"
                for r in rows:
                    if a1_col in r and a2_col in r:
                        a1 = r[a1_col]
                        a2 = r[a2_col]
                        vo = r.get(order_col, "SING")
                        bo = order_map.get(vo, 0)
                        bonds.append([a1, a2, bo])
        else:
            i += 1

    return {
        "comp_id": comp_id,
        "atoms": atoms,
        "bonds": bonds,
    }

def main():
    templates = {}
    print(f"Fetching and parsing {len(COMPONENTS)} CCD components from RCSB...")

    for comp_id in COMPONENTS:
        print(f"  Fetching {comp_id}...", end="", flush=True)
        cif_text = fetch_ccd_cif(comp_id)
        parsed = parse_ccd_cif(cif_text, comp_id)
        templates[comp_id] = parsed
        print(f" OK: {len(parsed['atoms'])} atoms, {len(parsed['bonds'])} bonds")

    output_dir = Path(__file__).resolve().parent.parent / "rust" / "crates" / "proteindf-bridge" / "src" / "data"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "ccd_bond_templates.msgpack"

    print(f"Writing MessagePack database to {output_path}...")
    with open(output_path, "wb") as fp:
        msgpack.pack(templates, fp)

    size = output_path.stat().st_size
    print(f"Done! Output size: {size} bytes ({size / 1024:.1f} KB)")

if __name__ == "__main__":
    main()
