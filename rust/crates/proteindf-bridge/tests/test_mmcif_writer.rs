// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

//! Tests for PDBx/mmCIF structure writer (`RUST_PORT_SPEC.md` §3.17, PR#42).

use std::collections::HashMap;
use std::path::PathBuf;

use proteindf_bridge::atom::Atom;
use proteindf_bridge::atom_group::AtomGroup;
use proteindf_bridge::ccd_templates::CcdTemplateDb;
use proteindf_bridge::format::mmcif_writer::{
    is_standard_residue, parse_model_key, parse_residue_key, quote_cif_value,
};
use proteindf_bridge::format::{MmcifWriteOptions, Pdb, SimpleMmcif};
use proteindf_bridge::position::Position;

fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

fn count_all_hydrogens(ag: &AtomGroup) -> usize {
    let mut count = 0;
    for (_m_key, m) in ag.groups() {
        for (_c_key, c) in m.groups() {
            for (_r_key, r) in c.groups() {
                for (_a_key, a) in r.atoms() {
                    if a.atomic_number() == 1 {
                        count += 1;
                    }
                }
            }
        }
    }
    count
}

/// Verifies that two AtomGroups match in models, chains, residues, atom names, elements, charges, and coordinates.
fn assert_atomgroups_match_roundtrip(orig: &AtomGroup, reloaded: &AtomGroup) {
    assert_eq!(
        orig.get_number_of_groups(),
        reloaded.get_number_of_groups(),
        "model count mismatch"
    );
    for (model_key, orig_model) in orig.groups() {
        let rel_model = reloaded
            .get_group(model_key)
            .unwrap_or_else(|| panic!("missing model {model_key} in reloaded"));
        assert_eq!(
            orig_model.get_number_of_groups(),
            rel_model.get_number_of_groups(),
            "chain count mismatch in {model_key}"
        );

        for (chain_key, orig_chain) in orig_model.groups() {
            let rel_chain = rel_model
                .get_group(chain_key)
                .unwrap_or_else(|| panic!("missing chain {chain_key} in {model_key}"));
            assert_eq!(
                orig_chain.get_number_of_groups(),
                rel_chain.get_number_of_groups(),
                "residue count mismatch in {model_key}/{chain_key}"
            );

            for (res_key, orig_res) in orig_chain.groups() {
                let rel_res = rel_chain.get_group(res_key).unwrap_or_else(|| {
                    panic!("missing residue {res_key} in {model_key}/{chain_key}")
                });
                assert_eq!(
                    orig_res.name, rel_res.name,
                    "residue name mismatch for {model_key}/{chain_key}/{res_key}"
                );
                assert_eq!(
                    orig_res.get_number_of_atoms(),
                    rel_res.get_number_of_atoms(),
                    "atom count mismatch in {model_key}/{chain_key}/{res_key}"
                );

                let mut orig_atoms: HashMap<&str, &Atom> = HashMap::new();
                for (_k, a) in orig_res.atoms() {
                    orig_atoms.insert(&a.name, a);
                }
                let mut rel_atoms: HashMap<&str, &Atom> = HashMap::new();
                for (_k, a) in rel_res.atoms() {
                    rel_atoms.insert(&a.name, a);
                }

                assert_eq!(
                    orig_atoms.len(),
                    rel_atoms.len(),
                    "atom count mismatch in {model_key}/{chain_key}/{res_key}"
                );

                for (name, orig_atom) in orig_atoms {
                    let rel_atom = rel_atoms.get(name).unwrap_or_else(|| {
                        panic!("missing atom {name} in {model_key}/{chain_key}/{res_key}")
                    });
                    assert_eq!(
                        orig_atom.atomic_number(),
                        rel_atom.atomic_number(),
                        "atomic number mismatch for {model_key}/{chain_key}/{res_key}/{name}"
                    );
                    assert!(
                        (orig_atom.charge - rel_atom.charge).abs() < 1e-4,
                        "charge mismatch for {model_key}/{chain_key}/{res_key}/{name}: orig {} vs rel {}",
                        orig_atom.charge,
                        rel_atom.charge
                    );
                    assert!(
                        (orig_atom.xyz.x - rel_atom.xyz.x).abs() <= 1e-3,
                        "coord x mismatch for {model_key}/{chain_key}/{res_key}/{name}: orig {} vs rel {}",
                        orig_atom.xyz.x,
                        rel_atom.xyz.x
                    );
                    assert!(
                        (orig_atom.xyz.y - rel_atom.xyz.y).abs() <= 1e-3,
                        "coord y mismatch for {model_key}/{chain_key}/{res_key}/{name}: orig {} vs rel {}",
                        orig_atom.xyz.y,
                        rel_atom.xyz.y
                    );
                    assert!(
                        (orig_atom.xyz.z - rel_atom.xyz.z).abs() <= 1e-3,
                        "coord z mismatch for {model_key}/{chain_key}/{res_key}/{name}: orig {} vs rel {}",
                        orig_atom.xyz.z,
                        rel_atom.xyz.z
                    );
                }
            }
        }
    }
}

// ========================================================================
// 1. Real Data Roundtrip Tests
// ========================================================================

#[test]
fn test_roundtrip_real_1hls() {
    let path = test_data_dir().join("1HLS.cif");
    let cif = SimpleMmcif::from_file(&path).expect("failed to load 1HLS.cif");
    let orig = cif
        .get_structure_atomgroup(None, None)
        .expect("failed to get structure atomgroup");

    // 1HLS has 20 NMR models
    assert_eq!(orig.get_number_of_groups(), 20);

    let mut buf = Vec::new();
    let opts = MmcifWriteOptions {
        data_block_name: "1HLS".to_string(),
        charge_to_b_factor: false,
    };
    SimpleMmcif::write_structure(&orig, &mut buf, &opts).expect("failed to write 1HLS mmCIF");

    let cif_str = String::from_utf8(buf).expect("invalid UTF-8");
    let reloaded_cif = SimpleMmcif::from_str(&cif_str).expect("failed to reload 1HLS mmCIF");
    let reloaded = reloaded_cif
        .get_structure_atomgroup(None, None)
        .expect("failed to get reloaded atomgroup");

    assert_atomgroups_match_roundtrip(&orig, &reloaded);
}

#[test]
fn test_save_structure_real_1hls() {
    let path = test_data_dir().join("1HLS.cif");
    let cif = SimpleMmcif::from_file(&path).expect("failed to load 1HLS.cif");
    let orig = cif
        .get_structure_atomgroup(None, None)
        .expect("failed to get structure atomgroup");

    let out_path = std::env::temp_dir().join("1hls_roundtrip_test.cif");
    let opts = MmcifWriteOptions {
        data_block_name: "1HLS".to_string(),
        charge_to_b_factor: false,
    };
    SimpleMmcif::save_structure(&orig, &out_path, &opts).expect("failed to save 1HLS to file");

    let reloaded_cif = SimpleMmcif::from_file(&out_path).expect("failed to reload saved file");
    let reloaded = reloaded_cif
        .get_structure_atomgroup(None, None)
        .expect("failed to get reloaded atomgroup");

    assert_atomgroups_match_roundtrip(&orig, &reloaded);
    let _ = std::fs::remove_file(&out_path);
}

#[test]
fn test_roundtrip_real_2fb4() {
    let path = test_data_dir().join("2FB4.cif");
    let cif = SimpleMmcif::from_file(&path).expect("failed to load 2FB4.cif");
    let orig = cif
        .get_structure_atomgroup(None, None)
        .expect("failed to get structure atomgroup");

    let mut buf = Vec::new();
    let opts = MmcifWriteOptions {
        data_block_name: "2FB4".to_string(),
        charge_to_b_factor: false,
    };
    SimpleMmcif::write_structure(&orig, &mut buf, &opts).expect("failed to write 2FB4 mmCIF");

    let cif_str = String::from_utf8(buf).expect("invalid UTF-8");
    let reloaded_cif = SimpleMmcif::from_str(&cif_str).expect("failed to reload 2FB4 mmCIF");
    let reloaded = reloaded_cif
        .get_structure_atomgroup(None, None)
        .expect("failed to get reloaded atomgroup");

    assert_atomgroups_match_roundtrip(&orig, &reloaded);
}

#[test]
fn test_roundtrip_real_2mgo() {
    let path = test_data_dir().join("2MGO.cif");
    let cif = SimpleMmcif::from_file(&path).expect("failed to load 2MGO.cif");
    let orig = cif
        .get_structure_atomgroup(None, None)
        .expect("failed to get structure atomgroup");

    // 2MGO has 20 NMR models
    assert_eq!(orig.get_number_of_groups(), 20);

    let mut buf = Vec::new();
    let opts = MmcifWriteOptions {
        data_block_name: "2MGO".to_string(),
        charge_to_b_factor: false,
    };
    SimpleMmcif::write_structure(&orig, &mut buf, &opts).expect("failed to write 2MGO mmCIF");

    let cif_str = String::from_utf8(buf).expect("invalid UTF-8");
    let reloaded_cif = SimpleMmcif::from_str(&cif_str).expect("failed to reload 2MGO mmCIF");
    let reloaded = reloaded_cif
        .get_structure_atomgroup(None, None)
        .expect("failed to get reloaded atomgroup");

    assert_atomgroups_match_roundtrip(&orig, &reloaded);
}

#[test]
fn test_roundtrip_real_3i3z() {
    let path = test_data_dir().join("3I3Z.cif");
    let cif = SimpleMmcif::from_file(&path).expect("failed to load 3I3Z.cif");
    let orig = cif
        .get_structure_atomgroup(None, None)
        .expect("failed to get structure atomgroup");

    let mut buf = Vec::new();
    let opts = MmcifWriteOptions {
        data_block_name: "3I3Z".to_string(),
        charge_to_b_factor: false,
    };
    SimpleMmcif::write_structure(&orig, &mut buf, &opts).expect("failed to write 3I3Z mmCIF");

    let cif_str = String::from_utf8(buf).expect("invalid UTF-8");
    let reloaded_cif = SimpleMmcif::from_str(&cif_str).expect("failed to reload 3I3Z mmCIF");
    let reloaded = reloaded_cif
        .get_structure_atomgroup(None, None)
        .expect("failed to get reloaded atomgroup");

    assert_atomgroups_match_roundtrip(&orig, &reloaded);
}

// ========================================================================
// 2. Hydrogenation Pipeline Roundtrip Test
// ========================================================================

#[test]
fn test_hydrogenation_pipeline_roundtrip() {
    let pdb_path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&pdb_path, None).expect("failed to read 1hls.pdb");
    let mut structure = pdb
        .get_atomgroup(None, None)
        .expect("failed to convert to AtomGroup");

    structure.setup().expect("setup failed");

    let db = CcdTemplateDb::global();
    let report = structure
        .add_missing_hydrogens(db)
        .expect("hydrogenation failed on 1hls");
    assert!(report.total_added_hydrogens > 0);

    let initial_h_count = count_all_hydrogens(&structure);
    assert!(initial_h_count > 0);

    // Save and reload via mmCIF
    let mut buf = Vec::new();
    let opts = MmcifWriteOptions::default();
    SimpleMmcif::write_structure(&structure, &mut buf, &opts).expect("failed to write mmCIF");

    let cif_str = String::from_utf8(buf).expect("invalid UTF-8");
    let reloaded_cif = SimpleMmcif::from_str(&cif_str).expect("failed to reload mmCIF");
    let reloaded = reloaded_cif
        .get_structure_atomgroup(None, None)
        .expect("failed to get structure atomgroup");

    let reloaded_h_count = count_all_hydrogens(&reloaded);
    assert_eq!(initial_h_count, reloaded_h_count, "hydrogen count mismatch");

    assert_atomgroups_match_roundtrip(&structure, &reloaded);
}

// ========================================================================
// 3. Synthetic Data Tests
// ========================================================================

#[test]
fn test_synthetic_negative_residue_and_insertion_code() {
    let mut root = AtomGroup::new();
    let mut model = AtomGroup::with_name("model_1");
    let mut chain = AtomGroup::with_name("A");

    // Residue with negative number and insertion code: -1A
    let mut res1 = AtomGroup::with_name("ALA");
    let mut atom1 = Atom::new();
    atom1.name = "CA".to_string();
    atom1.set_atomic_number(6);
    atom1.xyz = Position::new(1.0, 2.0, 3.0);
    res1.set_atom("1_CA", atom1);
    chain.set_group("-1A", res1);

    // Residue with positive number and insertion code: 52B
    let mut res2 = AtomGroup::with_name("GLY");
    let mut atom2 = Atom::new();
    atom2.name = "CA".to_string();
    atom2.set_atomic_number(6);
    atom2.xyz = Position::new(4.0, 5.0, 6.0);
    res2.set_atom("2_CA", atom2);
    chain.set_group("52B", res2);

    model.set_group("A", chain);
    root.set_group("model_1", model);

    let mut buf = Vec::new();
    SimpleMmcif::write_structure(&root, &mut buf, &MmcifWriteOptions::default()).unwrap();
    let cif_str = String::from_utf8(buf).unwrap();

    let reloaded = SimpleMmcif::from_str(&cif_str)
        .unwrap()
        .get_structure_atomgroup(None, None)
        .unwrap();

    let rel_chain = reloaded
        .get_group("model_1")
        .unwrap()
        .get_group("A")
        .unwrap();
    assert!(rel_chain.has_group("-1A"));
    assert!(rel_chain.has_group("52B"));
    assert_atomgroups_match_roundtrip(&root, &reloaded);
}

#[test]
fn test_synthetic_quoting_and_roundtrip() {
    let mut root = AtomGroup::new();
    let mut model = AtomGroup::with_name("model_1");
    let mut chain = AtomGroup::with_name("A");

    // Residue with quote in atom name (like nucleic acid O5')
    let mut res = AtomGroup::with_name("DA");

    let mut a1 = Atom::new();
    a1.name = "O5'".to_string();
    a1.set_atomic_number(8);
    a1.xyz = Position::new(1.0, 0.0, 0.0);
    res.set_atom("1_O5'", a1);

    // Reserved keyword-like atom name
    let mut a2 = Atom::new();
    a2.name = "data_foo".to_string();
    a2.set_atomic_number(6);
    a2.xyz = Position::new(2.0, 0.0, 0.0);
    res.set_atom("2_data_foo", a2);

    // Atom name starting with underscore
    let mut a3 = Atom::new();
    a3.name = "_X".to_string();
    a3.set_atomic_number(6);
    a3.xyz = Position::new(3.0, 0.0, 0.0);
    res.set_atom("3__X", a3);

    chain.set_group("1", res);
    model.set_group("A", chain);
    root.set_group("model_1", model);

    let mut buf = Vec::new();
    SimpleMmcif::write_structure(&root, &mut buf, &MmcifWriteOptions::default()).unwrap();
    let cif_str = String::from_utf8(buf).unwrap();

    // Check that O5' is quoted with double quotes
    assert!(cif_str.contains("\"O5'\""));
    // Check that data_foo is quoted with single quotes
    assert!(cif_str.contains("'data_foo'"));
    // Check that _X is quoted with single quotes
    assert!(cif_str.contains("'_X'"));

    let reloaded = SimpleMmcif::from_str(&cif_str)
        .unwrap()
        .get_structure_atomgroup(None, None)
        .unwrap();

    assert_atomgroups_match_roundtrip(&root, &reloaded);
}

#[test]
fn test_synthetic_mixed_quotes_error() {
    let mut root = AtomGroup::new();
    let mut model = AtomGroup::with_name("model_1");
    let mut chain = AtomGroup::with_name("A");
    let mut res = AtomGroup::with_name("ALA");

    let mut atom = Atom::new();
    atom.name = "A'\"B".to_string(); // contains both ' and "
    atom.set_atomic_number(6);
    atom.xyz = Position::new(0.0, 0.0, 0.0);
    res.set_atom("1_A'\"B", atom);

    chain.set_group("1", res);
    model.set_group("A", chain);
    root.set_group("model_1", model);

    let mut buf = Vec::new();
    let result = SimpleMmcif::write_structure(&root, &mut buf, &MmcifWriteOptions::default());
    assert!(result.is_err(), "expected error for mixed quotes");
    let err_msg = result.unwrap_err().to_string();
    assert!(err_msg.contains("both single and double quotes"));
}

#[test]
fn test_synthetic_schema_violation_error() {
    let mut root = AtomGroup::new();
    let mut model = AtomGroup::with_name("model_1");
    let mut chain = AtomGroup::with_name("A");

    // Atom placed directly under chain instead of residue (depth violation)
    let mut direct_atom = Atom::new();
    direct_atom.name = "CA".to_string();
    direct_atom.set_atomic_number(6);
    chain.set_atom("1_CA", direct_atom);

    model.set_group("A", chain);
    root.set_group("model_1", model);

    let mut buf = Vec::new();
    let result = SimpleMmcif::write_structure(&root, &mut buf, &MmcifWriteOptions::default());
    assert!(result.is_err(), "expected schema violation error");
    let err_msg = result.unwrap_err().to_string();
    assert!(err_msg.contains("schema validation failed"));
}

#[test]
fn test_synthetic_unparseable_residue_key_error() {
    let mut root = AtomGroup::new();
    let mut model = AtomGroup::with_name("model_1");
    let mut chain = AtomGroup::with_name("A");
    let mut res = AtomGroup::with_name("ALA");

    let mut atom = Atom::new();
    atom.name = "CA".to_string();
    atom.set_atomic_number(6);
    res.set_atom("1_CA", atom);

    // Invalid residue key: letters only, cannot extract integer
    chain.set_group("INVALID", res);
    model.set_group("A", chain);
    root.set_group("model_1", model);

    let mut buf = Vec::new();
    let result = SimpleMmcif::write_structure(&root, &mut buf, &MmcifWriteOptions::default());
    assert!(
        result.is_err(),
        "expected error for unparseable residue key"
    );
    let err_msg = result.unwrap_err().to_string();
    assert!(err_msg.contains("does not start with a valid integer"));
}

#[test]
fn test_synthetic_unparseable_model_key_error() {
    let mut root = AtomGroup::new();
    let mut model = AtomGroup::with_name("invalid_model");
    let mut chain = AtomGroup::with_name("A");
    let mut res = AtomGroup::with_name("ALA");

    let mut atom = Atom::new();
    atom.name = "CA".to_string();
    atom.set_atomic_number(6);
    res.set_atom("1_CA", atom);

    chain.set_group("1", res);
    model.set_group("A", chain);
    root.set_group("invalid_model", model);

    let mut buf = Vec::new();
    let result = SimpleMmcif::write_structure(&root, &mut buf, &MmcifWriteOptions::default());
    assert!(result.is_err(), "expected error for unparseable model key");
    let err_msg = result.unwrap_err().to_string();
    assert!(err_msg.contains("does not start with 'model_'"));
}

#[test]
fn test_synthetic_non_integer_charge_handling() {
    let mut root = AtomGroup::new();
    let mut model = AtomGroup::with_name("model_1");
    let mut chain = AtomGroup::with_name("A");
    let mut res = AtomGroup::with_name("ALA");

    let mut atom = Atom::new();
    atom.name = "CA".to_string();
    atom.set_atomic_number(6);
    atom.charge = 0.35; // non-integer partial charge
    atom.xyz = Position::new(1.0, 2.0, 3.0);
    res.set_atom("1_CA", atom);

    chain.set_group("1", res);
    model.set_group("A", chain);
    root.set_group("model_1", model);

    // 1. charge_to_b_factor = false -> error
    let mut buf = Vec::new();
    let opts_err = MmcifWriteOptions {
        charge_to_b_factor: false,
        ..Default::default()
    };
    let result = SimpleMmcif::write_structure(&root, &mut buf, &opts_err);
    assert!(
        result.is_err(),
        "expected error for non-integer charge when charge_to_b_factor is false"
    );
    let err_msg = result.unwrap_err().to_string();
    assert!(err_msg.contains("non-integer formal charge"));

    // 2. charge_to_b_factor = true -> success, written to B_iso_or_equiv, formal charge is '?'
    let mut buf_ok = Vec::new();
    let opts_ok = MmcifWriteOptions {
        charge_to_b_factor: true,
        ..Default::default()
    };
    SimpleMmcif::write_structure(&root, &mut buf_ok, &opts_ok)
        .expect("write should succeed with charge_to_b_factor = true");

    let cif_str = String::from_utf8(buf_ok).unwrap();
    // B_iso_or_equiv should be 0.35, formal charge should be ?
    assert!(cif_str.contains(" 0.35 ? 1 ALA A CA 1"));
}

#[test]
fn test_synthetic_group_pdb_classification() {
    let mut root = AtomGroup::new();
    let mut model = AtomGroup::with_name("model_1");
    let mut chain = AtomGroup::with_name("A");

    // Standard amino acid -> ATOM
    let mut res_ala = AtomGroup::with_name("ALA");
    let mut a_ala = Atom::new();
    a_ala.name = "CA".to_string();
    a_ala.set_atomic_number(6);
    res_ala.set_atom("1_CA", a_ala);
    chain.set_group("1", res_ala);

    // Standard nucleic acid -> ATOM
    let mut res_da = AtomGroup::with_name("DA");
    let mut a_da = Atom::new();
    a_da.name = "C1'".to_string();
    a_da.set_atomic_number(6);
    res_da.set_atom("2_C1'", a_da);
    chain.set_group("2", res_da);

    // Non-standard component / ligand (LIG) -> HETATM
    let mut res_lig = AtomGroup::with_name("LIG");
    let mut a_lig = Atom::new();
    a_lig.name = "C1".to_string();
    a_lig.set_atomic_number(6);
    res_lig.set_atom("3_C1", a_lig);
    chain.set_group("3", res_lig);

    // Water (HOH) -> HETATM
    let mut res_hoh = AtomGroup::with_name("HOH");
    let mut a_hoh = Atom::new();
    a_hoh.name = "O".to_string();
    a_hoh.set_atomic_number(8);
    res_hoh.set_atom("4_O", a_hoh);
    chain.set_group("4", res_hoh);

    model.set_group("A", chain);
    root.set_group("model_1", model);

    let mut buf = Vec::new();
    SimpleMmcif::write_structure(&root, &mut buf, &MmcifWriteOptions::default()).unwrap();
    let cif_str = String::from_utf8(buf).unwrap();

    let lines: Vec<&str> = cif_str
        .lines()
        .filter(|l| l.starts_with("ATOM") || l.starts_with("HETATM"))
        .collect();
    assert_eq!(lines.len(), 4);
    assert!(lines[0].starts_with("ATOM")); // ALA
    assert!(lines[1].starts_with("ATOM")); // DA
    assert!(lines[2].starts_with("HETATM")); // LIG
    assert!(lines[3].starts_with("HETATM")); // HOH

    // Also check label_seq_id: ATOM has sequence number, HETATM has '.'
    assert!(lines[0].contains(" 1 ? ")); // label_seq_id = 1 for ALA
    assert!(lines[1].contains(" 2 ? ")); // label_seq_id = 2 for DA
    assert!(lines[2].contains(" . ? ")); // label_seq_id = . for LIG
    assert!(lines[3].contains(" . ? ")); // label_seq_id = . for HOH
}

// ========================================================================
// 4. Large Structure (>100,000 atoms) and Benchmark
// ========================================================================

fn build_large_synthetic_protein(total_residues: usize, atoms_per_residue: usize) -> AtomGroup {
    let mut root = AtomGroup::new();
    let mut model = AtomGroup::with_name("model_1");
    let mut chain = AtomGroup::with_name("A");

    let mut serial: usize = 0;
    for r in 1..=total_residues {
        let res_key = r.to_string();
        let mut residue = AtomGroup::with_name("ALA");
        for a in 0..atoms_per_residue {
            serial += 1;
            let mut atom = Atom::new();
            atom.name = format!("C{a}");
            atom.set_atomic_number(6);
            atom.xyz = Position::new(r as f64, a as f64, 0.0);
            residue.set_atom(&format!("{serial}_{}", atom.name), atom);
        }
        chain.set_group(&res_key, residue);
    }
    model.set_group("A", chain);
    root.set_group("model_1", model);
    root
}

#[test]
fn test_large_structure_over_100k_atoms() {
    // 10,001 residues (exceeding PDB 9,999 residue limit) * 10 atoms = 100,010 atoms (exceeding 99,999 atom limit)
    let ag = build_large_synthetic_protein(10_001, 10);
    assert_eq!(ag.get_number_of_all_atoms(), 100_010);

    let mut buf = Vec::new();
    SimpleMmcif::write_structure(&ag, &mut buf, &MmcifWriteOptions::default())
        .expect("writing 100k atoms should succeed");

    let cif_str = String::from_utf8(buf).expect("valid UTF-8");
    let reloaded_cif =
        SimpleMmcif::from_str(&cif_str).expect("reloading 100k atoms should succeed");
    let reloaded = reloaded_cif
        .get_structure_atomgroup(None, None)
        .expect("extracting atomgroup should succeed");

    assert_eq!(reloaded.get_number_of_all_atoms(), 100_010);

    // Verify coordinates of the last atom
    let last_orig = ag
        .get_group("model_1")
        .unwrap()
        .get_group("A")
        .unwrap()
        .get_group("10001")
        .unwrap()
        .get_atom("C9")
        .unwrap();
    let last_rel = reloaded
        .get_group("model_1")
        .unwrap()
        .get_group("A")
        .unwrap()
        .get_group("10001")
        .unwrap()
        .get_atom("C9")
        .unwrap();

    assert_eq!(last_orig.atomic_number(), last_rel.atomic_number());
    assert!((last_orig.xyz.x - last_rel.xyz.x).abs() < 1e-3);
    assert!((last_orig.xyz.y - last_rel.xyz.y).abs() < 1e-3);
    assert!((last_orig.xyz.z - last_rel.xyz.z).abs() < 1e-3);
}

#[test]
#[ignore = "benchmark for 1,000,000 atoms (run with -- --ignored --nocapture)"]
fn benchmark_write_1_million_atoms() {
    // 100,000 residues * 10 atoms = 1,000,000 atoms
    eprintln!("Building 1,000,000 atom synthetic structure...");
    let build_start = std::time::Instant::now();
    let ag = build_large_synthetic_protein(100_000, 10);
    assert_eq!(ag.get_number_of_all_atoms(), 1_000_000);
    eprintln!("Built 1,000,000 atoms in {:?}", build_start.elapsed());

    // Benchmark writing to a null / buffered stream
    let write_start = std::time::Instant::now();
    let mut sink = std::io::BufWriter::with_capacity(1024 * 1024, std::io::sink());
    SimpleMmcif::write_structure(&ag, &mut sink, &MmcifWriteOptions::default())
        .expect("write failed");
    use std::io::Write;
    sink.flush().expect("flush failed");
    let elapsed = write_start.elapsed();
    eprintln!("Time to write 1,000,000 atoms: {:?}", elapsed);
    println!("BENCHMARK_1M_ATOMS_TIME: {:?}", elapsed);
}

// ========================================================================
// 5. Unit Tests for Helpers
// ========================================================================

#[test]
fn test_parse_residue_key_cases() {
    assert_eq!(parse_residue_key("1").unwrap(), (1, None));
    assert_eq!(parse_residue_key("-1").unwrap(), (-1, None));
    assert_eq!(parse_residue_key("52A").unwrap(), (52, Some("A")));
    assert_eq!(parse_residue_key("-10B").unwrap(), (-10, Some("B")));
    assert_eq!(parse_residue_key("+3C").unwrap(), (3, Some("C")));

    assert!(parse_residue_key("").is_err());
    assert!(parse_residue_key("A").is_err());
    assert!(parse_residue_key("-A").is_err());
    assert!(parse_residue_key("+").is_err());
}

#[test]
fn test_parse_model_key_cases() {
    assert_eq!(parse_model_key("model_1").unwrap(), 1);
    assert_eq!(parse_model_key("model_20").unwrap(), 20);

    assert!(parse_model_key("model_0").is_err());
    assert!(parse_model_key("model_-1").is_err());
    assert!(parse_model_key("model_abc").is_err());
    assert!(parse_model_key("1").is_err());
    assert!(parse_model_key("").is_err());
}

#[test]
fn test_quote_cif_value_cases() {
    assert_eq!(quote_cif_value("ALA").unwrap(), "ALA");
    assert_eq!(quote_cif_value("O5'").unwrap(), "\"O5'\"");
    assert_eq!(quote_cif_value("A B").unwrap(), "'A B'");
    assert_eq!(quote_cif_value("").unwrap(), "''");
    assert_eq!(quote_cif_value("_tag").unwrap(), "'_tag'");
    assert_eq!(quote_cif_value("data_foo").unwrap(), "'data_foo'");
    assert_eq!(quote_cif_value("loop_").unwrap(), "'loop_'");
    assert_eq!(quote_cif_value(".").unwrap(), "'.'");
    assert_eq!(quote_cif_value("?").unwrap(), "'?'");

    // Both single and double quotes -> error
    assert!(quote_cif_value("a'b\"c").is_err());
}

#[test]
fn test_is_standard_residue() {
    assert!(is_standard_residue("ALA"));
    assert!(is_standard_residue("GLY"));
    assert!(is_standard_residue("DA"));
    assert!(is_standard_residue("U"));

    assert!(!is_standard_residue("HOH"));
    assert!(!is_standard_residue("LIG"));
    assert!(!is_standard_residue("UNK"));
}
