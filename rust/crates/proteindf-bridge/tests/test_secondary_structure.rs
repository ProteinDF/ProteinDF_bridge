// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::path::PathBuf;

use proteindf_bridge::atom::Atom;
use proteindf_bridge::atom_group::AtomGroup;
use proteindf_bridge::format::pdb::Pdb;
use proteindf_bridge::position::Position;
use proteindf_bridge::secondary_structure::{calc_secondary_structure, SsCode};

fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

fn make_atom(name: &str, symbol: &str, pos: Position) -> Atom {
    let mut a = Atom::new_with_pos(symbol, pos).unwrap();
    a.name = name.to_string();
    a
}

#[test]
fn test_secondary_structure_1hls_chain_a() {
    let path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&path, None).expect("failed to load 1hls.pdb");
    let ag = pdb.get_atomgroup(None, None).expect("failed to get ag");

    let model_1 = ag.get_group("model_1").expect("model_1 not found");
    let chain_a = model_1.get_group("A").expect("chain A not found");

    let ss_list = calc_secondary_structure(chain_a);
    assert_eq!(ss_list.len(), 21);

    // Expected benchmark values from docs/rust-port-handoff.md:
    let expected = [
        ("1", "GLY", SsCode::Loop),
        ("2", "ILE", SsCode::Loop),
        ("3", "VAL", SsCode::Helix),
        ("4", "GLU", SsCode::Helix),
        ("5", "GLN", SsCode::Helix),
        ("6", "CYS", SsCode::Helix),
        ("7", "CYS", SsCode::Loop),
        ("8", "THR", SsCode::Loop),
        ("9", "SER", SsCode::Loop),
        ("10", "ILE", SsCode::Loop),
        ("11", "CYS", SsCode::Loop),
        ("12", "SER", SsCode::Loop),
        ("13", "LEU", SsCode::Loop),
        ("14", "TYR", SsCode::Loop),
        ("15", "GLN", SsCode::Loop),
        ("16", "LEU", SsCode::Loop),
        ("17", "GLU", SsCode::Helix),
        ("18", "ASN", SsCode::Helix),
        ("19", "TYR", SsCode::Helix),
        ("20", "CYS", SsCode::Loop),
        ("21", "ASN", SsCode::Loop),
    ];

    for (i, (exp_key, exp_name, exp_code)) in expected.iter().enumerate() {
        assert_eq!(ss_list[i].residue_key, *exp_key);
        assert_eq!(ss_list[i].residue_name, *exp_name);
        assert_eq!(
            ss_list[i].code, *exp_code,
            "Mismatch at residue {}: expected {}, got {}",
            exp_key, exp_code, ss_list[i].code
        );
    }
}

#[test]
fn test_secondary_structure_1hls_chain_b() {
    let path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&path, None).expect("failed to load 1hls.pdb");
    let ag = pdb.get_atomgroup(None, None).expect("failed to get ag");

    let model_1 = ag.get_group("model_1").expect("model_1 not found");
    let chain_b = model_1.get_group("B").expect("chain B not found");

    let ss_list = calc_secondary_structure(chain_b);
    assert_eq!(ss_list.len(), 30);

    // Expected benchmark values from docs/rust-port-handoff.md:
    let expected = [
        ("1", "PHE", SsCode::Loop),
        ("2", "VAL", SsCode::Loop),
        ("3", "ASN", SsCode::Loop),
        ("4", "GLN", SsCode::Loop),
        ("5", "HIS", SsCode::Loop),
        ("6", "LEU", SsCode::Loop),
        ("7", "CYS", SsCode::Loop),
        ("8", "GLY", SsCode::Loop),
        ("9", "SER", SsCode::Helix),
        ("10", "HIS", SsCode::Helix),
        ("11", "LEU", SsCode::Helix),
        ("12", "VAL", SsCode::Helix),
        ("13", "GLU", SsCode::Helix),
        ("14", "ALA", SsCode::Helix),
        ("15", "LEU", SsCode::Helix),
        ("16", "HIS", SsCode::Helix),
        ("17", "LEU", SsCode::Helix),
        ("18", "VAL", SsCode::Helix),
        ("19", "CYS", SsCode::Helix),
        ("20", "GLY", SsCode::Loop),
        ("21", "GLU", SsCode::Loop),
        ("22", "ARG", SsCode::Loop),
        ("23", "GLY", SsCode::Loop),
        ("24", "PHE", SsCode::Loop),
        ("25", "PHE", SsCode::Loop),
        ("26", "TYR", SsCode::Loop),
        ("27", "THR", SsCode::Loop),
        ("28", "PRO", SsCode::Loop),
        ("29", "LYS", SsCode::Loop),
        ("30", "THR", SsCode::Loop),
    ];

    for (i, (exp_key, exp_name, exp_code)) in expected.iter().enumerate() {
        assert_eq!(ss_list[i].residue_key, *exp_key);
        assert_eq!(ss_list[i].residue_name, *exp_name);
        assert_eq!(
            ss_list[i].code, *exp_code,
            "Mismatch at residue {}: expected {}, got {}",
            exp_key, exp_code, ss_list[i].code
        );
    }
}

#[test]
fn test_secondary_structure_antiparallel_beta_strand_synthetic() {
    // Synthetic 2-strand antiparallel beta-sheet configuration
    // Strand 1: residues 1, 2, 3
    // Strand 2: residues 11, 12, 13
    // Residues 2 and 12 form mutual antiparallel hydrogen bonds:
    // HBond(12 -> 2) and HBond(2 -> 12).
    let mut chain = AtomGroup::with_name("A");

    // Residue 1: provides C_prev for residue 2
    let mut res1 = AtomGroup::with_name("ALA");
    res1.set_atom("N", make_atom("N", "N", Position::new(2.0, -1.5, 0.0)));
    res1.set_atom("CA", make_atom("CA", "C", Position::new(2.0, -1.0, 0.0)));
    res1.set_atom("C", make_atom("C", "C", Position::new(2.0, -0.6, 0.0)));
    res1.set_atom("O", make_atom("O", "O", Position::new(1.0, -0.6, 0.0)));
    chain.set_group("1", res1);

    // Residue 2: Acceptor (C=O at (0,0)->(0, 1.23)) and Donor (N at (3, 0), H at (3, 1.01))
    let mut res2 = AtomGroup::with_name("VAL");
    res2.set_atom("N", make_atom("N", "N", Position::new(3.0, 0.0, 0.0)));
    res2.set_atom("CA", make_atom("CA", "C", Position::new(4.0, -0.6, 0.0)));
    res2.set_atom("C", make_atom("C", "C", Position::new(0.0, 0.0, 0.0)));
    res2.set_atom("O", make_atom("O", "O", Position::new(0.0, 1.23, 0.0)));
    chain.set_group("2", res2);

    // Residue 3: tail of strand 1
    let mut res3 = AtomGroup::with_name("LEU");
    res3.set_atom("N", make_atom("N", "N", Position::new(-1.0, 0.0, 0.0)));
    res3.set_atom("CA", make_atom("CA", "C", Position::new(-2.0, 0.0, 0.0)));
    res3.set_atom("C", make_atom("C", "C", Position::new(-3.0, 0.0, 0.0)));
    res3.set_atom("O", make_atom("O", "O", Position::new(-3.0, 1.0, 0.0)));
    chain.set_group("3", res3);

    // Residue 11: provides C_prev for residue 12
    let mut res11 = AtomGroup::with_name("ILE");
    res11.set_atom("N", make_atom("N", "N", Position::new(-1.0, 4.5, 0.0)));
    res11.set_atom("CA", make_atom("CA", "C", Position::new(-1.0, 4.0, 0.0)));
    res11.set_atom("C", make_atom("C", "C", Position::new(-1.0, 3.5, 0.0)));
    res11.set_atom("O", make_atom("O", "O", Position::new(-2.0, 3.5, 0.0)));
    chain.set_group("11", res11);

    // Residue 12: Donor (N at (0, 2.9), H at (0, 1.89)) and Acceptor (C=O at (3, 2.9)->(3, 1.67))
    let mut res12 = AtomGroup::with_name("TYR");
    res12.set_atom("N", make_atom("N", "N", Position::new(0.0, 2.9, 0.0)));
    res12.set_atom("CA", make_atom("CA", "C", Position::new(1.0, 3.5, 0.0)));
    res12.set_atom("C", make_atom("C", "C", Position::new(3.0, 2.9, 0.0)));
    res12.set_atom("O", make_atom("O", "O", Position::new(3.0, 1.67, 0.0)));
    chain.set_group("12", res12);

    // Residue 13: tail of strand 2
    let mut res13 = AtomGroup::with_name("PHE");
    res13.set_atom("N", make_atom("N", "N", Position::new(4.0, 2.9, 0.0)));
    res13.set_atom("CA", make_atom("CA", "C", Position::new(5.0, 2.9, 0.0)));
    res13.set_atom("C", make_atom("C", "C", Position::new(6.0, 2.9, 0.0)));
    res13.set_atom("O", make_atom("O", "O", Position::new(6.0, 4.0, 0.0)));
    chain.set_group("13", res13);

    let ss_list = calc_secondary_structure(&chain);
    let get_code = |key: &str| ss_list.iter().find(|s| s.residue_key == key).unwrap().code;

    assert_eq!(
        get_code("2"),
        SsCode::Strand,
        "Residue 2 must be detected as Strand (E)"
    );
    assert_eq!(
        get_code("12"),
        SsCode::Strand,
        "Residue 12 must be detected as Strand (E)"
    );
}

#[test]
fn test_unsorted_chain_keys() {
    let path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&path, None).expect("failed to load 1hls.pdb");
    let ag = pdb.get_atomgroup(None, None).expect("failed to get ag");

    let model_1 = ag.get_group("model_1").expect("model_1 not found");
    let chain_a = model_1.get_group("A").expect("chain A not found");

    let ss_sorted = calc_secondary_structure(chain_a);

    // Build scrambled chain with deliberately reversed keys
    let mut chain_scrambled = AtomGroup::with_name("A");
    let mut keys = chain_a.get_group_list();
    keys.reverse();
    for k in keys {
        let res = chain_a.get_group(&k).unwrap();
        chain_scrambled.set_group(&k, res.clone());
    }

    let ss_scrambled = calc_secondary_structure(&chain_scrambled);
    assert_eq!(ss_sorted, ss_scrambled);
}
