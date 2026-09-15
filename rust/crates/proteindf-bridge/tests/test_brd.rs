// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::path::PathBuf;

use proteindf_bridge::atom::Atom;
use proteindf_bridge::atom_group::{AtomGroup, BondRecord};
use proteindf_bridge::brd::{load_atomgroup, load_brd_yui, save_atomgroup, save_brd_yui};
use proteindf_bridge::position::Position;

fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

#[test]
fn test_load_conformer_fixtures() {
    let conformers = [
        "ACE_ALA_NME_trans1.brd",
        "ACE_ALA_NME_trans2.brd",
        "ACE_ALA_NME_cis1.brd",
        "ACE_ALA_NME_cis2.brd",
    ];

    for name in conformers {
        let path = test_data_dir().join(name);
        let ag = load_atomgroup(&path).unwrap_or_else(|e| panic!("failed to load {name}: {e}"));

        assert_eq!(
            ag.get_number_of_groups(),
            3,
            "failed group count for {name}"
        );
        let groups = ag.get_group_list();
        assert_eq!(groups, vec!["1", "2", "3"]);

        let atom_list = ag.get_atom_list();
        assert_eq!(atom_list.len(), 22, "failed atom count for {name}");

        // Check the first three atoms for trans1/cis conformers
        let paths: Vec<String> = atom_list.iter().map(|a| a.path.clone()).collect();
        assert_eq!(
            &paths[..3],
            &["/1/1_CH3", "/1/2_H1", "/1/3_H2"],
            "atom paths mismatch for {name}"
        );
    }

    // Check ACE_ALA_NME.brd (different atom naming order)
    let path = test_data_dir().join("ACE_ALA_NME.brd");
    let ag = load_atomgroup(&path).expect("failed to load ACE_ALA_NME.brd");
    assert_eq!(ag.get_number_of_groups(), 3);
    assert_eq!(ag.name, "_");
    let atom_list = ag.get_atom_list();
    assert_eq!(atom_list.len(), 22);
    let paths: Vec<String> = atom_list.iter().map(|a| a.path.clone()).collect();
    assert_eq!(
        &paths[..3],
        &["/1/1_HH31", "/1/2_CH3", "/1/3_HH32"],
        "atom paths mismatch for ACE_ALA_NME.brd"
    );
}

#[test]
fn test_load_nml_fixtures() {
    for name in ["NML.brd", "NML_trans.brd"] {
        let path = test_data_dir().join(name);
        let ag = load_atomgroup(&path).unwrap_or_else(|e| panic!("failed to load {name}: {e}"));

        assert_eq!(ag.get_number_of_groups(), 0);
        let atom_list = ag.get_atom_list();
        assert_eq!(atom_list.len(), 12, "failed atom count for {name}");

        let paths: Vec<String> = atom_list.iter().map(|a| a.path.clone()).collect();
        assert_eq!(paths[0], "/0");
        assert_eq!(paths[1], "/1");
        assert_eq!(paths[2], "/2");
    }
}

#[test]
fn test_atom_raw_data_roundtrip() {
    let mut atom = Atom::new();
    atom.set_atomic_number(6);
    atom.name = "CA".to_string();
    atom.charge = 0.25;
    atom.xyz = Position::new(1.0, 2.0, 3.0);
    atom.force = Position::new(0.01, -0.02, 0.03);

    let raw = atom.get_raw_data();
    let reconstructed = Atom::from_raw_data(&raw).expect("failed to reconstruct atom");

    assert_eq!(reconstructed.atomic_number(), 6);
    assert_eq!(reconstructed.name, "CA");
    assert!((reconstructed.charge - 0.25).abs() < 1e-6);
    assert_eq!(reconstructed.xyz, Position::new(1.0, 2.0, 3.0));
    assert_eq!(reconstructed.force, Position::new(0.01, -0.02, 0.03));
}

#[test]
fn test_atomgroup_raw_data_roundtrip() {
    let mut root = AtomGroup::with_name("protein");

    let mut res1 = AtomGroup::with_name("ALA");
    let mut a1 = Atom::new();
    a1.set_atomic_number(6);
    a1.name = "CA".to_string();
    a1.xyz = Position::new(1.0, 2.0, 3.0);
    res1.set_atom("CA", a1);

    let mut a2 = Atom::new();
    a2.set_atomic_number(7);
    a2.name = "N".to_string();
    a2.xyz = Position::new(0.0, 1.0, 2.0);
    res1.set_atom("N", a2);

    res1.set_bonds(vec![BondRecord {
        atom1_path: "/CA".to_string(),
        atom2_path: "/N".to_string(),
        order: 1,
    }]);

    root.set_group("1", res1);

    // Save and load via MessagePack file
    let tmp_dir = std::env::temp_dir();
    let tmp_file = tmp_dir.join("test_roundtrip.brd");
    save_atomgroup(&root, &tmp_file).expect("failed to save atomgroup");

    let loaded = load_atomgroup(&tmp_file).expect("failed to load atomgroup");
    assert_eq!(loaded.name, "protein");
    assert_eq!(loaded.get_number_of_groups(), 1);

    let sub = loaded.get_group("1").expect("missing group 1");
    assert_eq!(sub.name, "ALA");
    assert_eq!(sub.get_number_of_atoms(), 2);

    let ca = sub.get_atom("CA").expect("missing CA");
    assert_eq!(ca.atomic_number(), 6);
    assert_eq!(ca.xyz, Position::new(1.0, 2.0, 3.0));

    assert_eq!(sub.bonds().len(), 1);
    assert_eq!(sub.bonds()[0].order, 1);

    let _ = std::fs::remove_file(tmp_file);
}

#[test]
fn test_yui_format_uncompressed_and_zstd() {
    let path = test_data_dir().join("ACE_ALA_NME_trans1.brd");
    let original = load_atomgroup(&path).expect("failed to load fixture");

    let tmp_dir = std::env::temp_dir();

    // 1. Test uncompressed YUI format
    let yui_raw = tmp_dir.join("test_trans1_uncompressed.brd");
    save_brd_yui(&original, &yui_raw, false).expect("failed to save uncompressed YUI brd");

    let loaded_raw = load_brd_yui(&yui_raw).expect("failed to load uncompressed YUI brd");
    assert_eq!(loaded_raw.get_number_of_groups(), 3);
    assert_eq!(loaded_raw.get_atom_list().len(), 22);

    // 2. Test Zstd-compressed YUI format
    let yui_zstd = tmp_dir.join("test_trans1_zstd.brd");
    save_brd_yui(&original, &yui_zstd, true).expect("failed to save zstd YUI brd");

    let loaded_zstd = load_brd_yui(&yui_zstd).expect("failed to load zstd YUI brd");
    assert_eq!(loaded_zstd.get_number_of_groups(), 3);
    assert_eq!(loaded_zstd.get_atom_list().len(), 22);

    // Coordinate comparison between original and loaded_zstd
    let orig_atoms = original.clone().get_atom_list();
    let zstd_atoms = loaded_zstd.get_atom_list();
    for (a1, a2) in orig_atoms.iter().zip(zstd_atoms.iter()) {
        assert_eq!(a1.path, a2.path);
        assert_eq!(a1.atomic_number(), a2.atomic_number());
        assert!((a1.xyz.distance_from(&a2.xyz)).abs() < 1e-8);
    }

    // 3. Test invalid magic header rejection
    let bad_header_file = tmp_dir.join("test_bad_magic.brd");
    std::fs::write(&bad_header_file, b"BAD\0\x01\x00data").unwrap();
    let err = load_brd_yui(&bad_header_file);
    assert!(err.is_err());

    let _ = std::fs::remove_file(yui_raw);
    let _ = std::fs::remove_file(yui_zstd);
    let _ = std::fs::remove_file(bad_header_file);
}
