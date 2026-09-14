// Copyright (C) 2015 The ProteinDF development team.
// see also AUTHORS and README if provided.
//
// This file is a part of the ProteinDF software package.
//
// The ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

use pdf_bridge::atom::Atom;
use pdf_bridge::atom_group::AtomGroup;
use pdf_bridge::ion_pair::IonPair;
use pdf_bridge::position::Position;

#[test]
fn test_ion_pair_detection() {
    // Ported from tests/test_ionpair.py
    let mut protein = AtomGroup::with_name("protein");
    let mut chain_a = AtomGroup::with_name("A");

    // GLU (anionic residue: CD, OE1, OE2)
    let mut glu = AtomGroup::with_name("1");
    glu.name = "GLU".to_string();
    glu.set_atom(
        "CD",
        Atom::new_with_pos("C", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    glu.set_atom(
        "OE1",
        Atom::new_with_pos("O", Position::new(0.0, 1.0, 0.0)).unwrap(),
    );
    glu.set_atom(
        "OE2",
        Atom::new_with_pos("O", Position::new(1.0, 0.0, 0.0)).unwrap(),
    );
    chain_a.set_group("1", glu);

    // LYS (cationic residue: NZ) - distance < 4.0 Å
    let mut lys = AtomGroup::with_name("2");
    lys.name = "LYS".to_string();
    lys.set_atom(
        "NZ",
        Atom::new_with_pos("N", Position::new(0.5, 2.0, 0.0)).unwrap(),
    );
    chain_a.set_group("2", lys);

    protein.set_group("A", chain_a);

    let ip = IonPair::new(&protein);
    let pairs = ip.get_ion_pairs();

    assert_eq!(pairs.len(), 1);
    assert_eq!(pairs[0].anion_type, "GLU");
    assert_eq!(pairs[0].cation_type, "LYS");
    assert_eq!(pairs[0].anion_path, "/A/1/");
    assert_eq!(pairs[0].cation_path, "/A/2/");
}

#[test]
fn test_asp_and_arg_detection() {
    let mut protein = AtomGroup::new();
    let mut chain_a = AtomGroup::new();

    // ASP at origin: CG(0,0,0), OD1(0,1,0), OD2(1,0,0) -> center at (1/3, 1/3, 0)
    let mut asp = AtomGroup::new();
    asp.name = "ASP".to_string();
    asp.set_atom(
        "CG",
        Atom::new_with_pos("C", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    asp.set_atom(
        "OD1",
        Atom::new_with_pos("O", Position::new(0.0, 1.0, 0.0)).unwrap(),
    );
    asp.set_atom(
        "OD2",
        Atom::new_with_pos("O", Position::new(1.0, 0.0, 0.0)).unwrap(),
    );
    chain_a.set_group("1", asp);

    // ARG near ASP: CZ(1,1,0), NH1(1,2,0), NH2(2,1,0)
    let mut arg = AtomGroup::new();
    arg.name = "ARG".to_string();
    arg.set_atom(
        "CZ",
        Atom::new_with_pos("C", Position::new(1.0, 1.0, 0.0)).unwrap(),
    );
    arg.set_atom(
        "NH1",
        Atom::new_with_pos("N", Position::new(1.0, 2.0, 0.0)).unwrap(),
    );
    arg.set_atom(
        "NH2",
        Atom::new_with_pos("N", Position::new(2.0, 1.0, 0.0)).unwrap(),
    );
    chain_a.set_group("2", arg);

    protein.set_group("A", chain_a);

    let pairs = IonPair::find_ion_pairs(&protein);
    // ARG has 3 sites: ARG (center), ARG1 (NH1), ARG2 (NH2)
    // All 3 should be within 4.0 Å of ASP center (1/3, 1/3, 0)
    assert_eq!(pairs.len(), 3);
    let cation_types: Vec<&str> = pairs.iter().map(|p| p.cation_type.as_str()).collect();
    assert!(cation_types.contains(&"ARG"));
    assert!(cation_types.contains(&"ARG1"));
    assert!(cation_types.contains(&"ARG2"));
    assert!(pairs.iter().all(|p| p.anion_type == "ASP"));
}

#[test]
fn test_ctm_and_ntm_detection() {
    let mut protein = AtomGroup::new();
    let mut chain_a = AtomGroup::new();

    // C-terminal residue ALA with OXT at (0,0,0)
    let mut cterm = AtomGroup::new();
    cterm.name = "ALA".to_string();
    cterm.set_atom(
        "C",
        Atom::new_with_pos("C", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    cterm.set_atom(
        "O",
        Atom::new_with_pos("O", Position::new(0.0, 1.0, 0.0)).unwrap(),
    );
    cterm.set_atom(
        "OXT",
        Atom::new_with_pos("O", Position::new(1.0, 0.0, 0.0)).unwrap(),
    );
    chain_a.set_group("1", cterm);

    // N-terminal residue GLY with H3 and N at (0, 2, 0)
    let mut nterm = AtomGroup::new();
    nterm.name = "GLY".to_string();
    nterm.set_atom(
        "N",
        Atom::new_with_pos("N", Position::new(0.0, 2.0, 0.0)).unwrap(),
    );
    nterm.set_atom(
        "H3",
        Atom::new_with_pos("H", Position::new(0.0, 2.5, 0.0)).unwrap(),
    );
    chain_a.set_group("2", nterm);

    protein.set_group("A", chain_a);

    let pairs = IonPair::find_ion_pairs(&protein);
    assert_eq!(pairs.len(), 1);
    assert_eq!(pairs[0].anion_type, "CTM");
    assert_eq!(pairs[0].cation_type, "NTM");
}

#[test]
fn test_distance_threshold_boundary() {
    // Test boundary: 3.99 Å (should match) vs 4.01 Å (should not match)
    let mut protein = AtomGroup::new();
    let mut chain_a = AtomGroup::new();

    // LYS at origin
    let mut lys = AtomGroup::new();
    lys.name = "LYS".to_string();
    lys.set_atom(
        "NZ",
        Atom::new_with_pos("N", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    chain_a.set_group("1", lys);

    // GLU1 center at (3.99, 0, 0) -> distance 3.99 < 4.0
    let mut glu1 = AtomGroup::new();
    glu1.name = "GLU".to_string();
    glu1.set_atom(
        "CD",
        Atom::new_with_pos("C", Position::new(3.99, 0.0, 0.0)).unwrap(),
    );
    glu1.set_atom(
        "OE1",
        Atom::new_with_pos("O", Position::new(3.99, 0.0, 0.0)).unwrap(),
    );
    glu1.set_atom(
        "OE2",
        Atom::new_with_pos("O", Position::new(3.99, 0.0, 0.0)).unwrap(),
    );
    chain_a.set_group("2", glu1);

    // GLU2 center at (4.01, 0, 0) -> distance 4.01 > 4.0
    let mut glu2 = AtomGroup::new();
    glu2.name = "GLU".to_string();
    glu2.set_atom(
        "CD",
        Atom::new_with_pos("C", Position::new(4.01, 0.0, 0.0)).unwrap(),
    );
    glu2.set_atom(
        "OE1",
        Atom::new_with_pos("O", Position::new(4.01, 0.0, 0.0)).unwrap(),
    );
    glu2.set_atom(
        "OE2",
        Atom::new_with_pos("O", Position::new(4.01, 0.0, 0.0)).unwrap(),
    );
    chain_a.set_group("3", glu2);

    protein.set_group("A", chain_a);

    let pairs = IonPair::find_ion_pairs(&protein);
    assert_eq!(pairs.len(), 1);
    assert_eq!(pairs[0].anion_path, "/A/2/");
    assert_eq!(pairs[0].cation_path, "/A/1/");
}

#[test]
fn test_missing_atom_safe_skip() {
    // Unlike Python which raises KeyError on incomplete residues,
    // Rust safely skips incomplete residues without crashing the whole analysis.
    let mut protein = AtomGroup::new();
    let mut chain_a = AtomGroup::new();

    // 1. Incomplete GLU: missing OE2 (only has CD and OE1)
    let mut incomplete_glu = AtomGroup::new();
    incomplete_glu.name = "GLU".to_string();
    incomplete_glu.set_atom(
        "CD",
        Atom::new_with_pos("C", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    incomplete_glu.set_atom(
        "OE1",
        Atom::new_with_pos("O", Position::new(1.0, 0.0, 0.0)).unwrap(),
    );
    // Note: OE2 is deliberately missing
    chain_a.set_group("1", incomplete_glu);

    // 2. Complete ASP: CG, OD1, OD2
    let mut asp = AtomGroup::new();
    asp.name = "ASP".to_string();
    asp.set_atom(
        "CG",
        Atom::new_with_pos("C", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    asp.set_atom(
        "OD1",
        Atom::new_with_pos("O", Position::new(1.0, 0.0, 0.0)).unwrap(),
    );
    asp.set_atom(
        "OD2",
        Atom::new_with_pos("O", Position::new(0.0, 1.0, 0.0)).unwrap(),
    );
    chain_a.set_group("2", asp);

    // 3. Complete LYS: NZ
    let mut lys = AtomGroup::new();
    lys.name = "LYS".to_string();
    lys.set_atom(
        "NZ",
        Atom::new_with_pos("N", Position::new(0.5, 0.5, 0.0)).unwrap(),
    );
    chain_a.set_group("3", lys);

    protein.set_group("A", chain_a);

    // Should not crash, and should detect the valid ASP-LYS pair while skipping incomplete GLU
    let pairs = IonPair::find_ion_pairs(&protein);
    assert_eq!(pairs.len(), 1);
    assert_eq!(pairs[0].anion_path, "/A/2/");
    assert_eq!(pairs[0].anion_type, "ASP");
    assert_eq!(pairs[0].cation_path, "/A/3/");
    assert_eq!(pairs[0].cation_type, "LYS");
}
