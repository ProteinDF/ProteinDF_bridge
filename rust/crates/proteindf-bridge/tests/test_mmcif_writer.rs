// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

//! Tests for PDBx/mmCIF structure writer (`RUST_PORT_SPEC.md` §3.17, PR#42).

use std::collections::HashMap;
use std::path::PathBuf;

use proteindf_bridge::atom::Atom;
use proteindf_bridge::atom_group::AtomGroup;
use proteindf_bridge::ccd_templates::CcdTemplateDb;
use proteindf_bridge::format::mmcif_writer::{
    determine_cif_quote_kind, is_standard_residue, parse_model_key, parse_residue_key,
    quote_cif_value, validate_data_block_name, CifQuoteKind,
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
                // Verify uniqueness of atom names within the original residue
                assert_eq!(
                    orig_atoms.len(),
                    orig_res.get_number_of_atoms(),
                    "duplicate atom name found in original residue {model_key}/{chain_key}/{res_key}"
                );

                let mut rel_atoms: HashMap<&str, &Atom> = HashMap::new();
                for (_k, a) in rel_res.atoms() {
                    rel_atoms.insert(&a.name, a);
                }
                // Verify uniqueness of atom names within the reloaded residue
                assert_eq!(
                    rel_atoms.len(),
                    rel_res.get_number_of_atoms(),
                    "duplicate atom name found in reloaded residue {model_key}/{chain_key}/{res_key}"
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
fn test_synthetic_empty_chain_id_roundtrip() {
    let mut root = AtomGroup::new();
    let mut model = AtomGroup::with_name("model_1");
    // Empty chain ID is represented by group key "_"
    let mut chain = AtomGroup::with_name("_");

    let mut res = AtomGroup::with_name("GLY");
    let mut atom = Atom::new();
    atom.name = "CA".to_string();
    atom.set_atomic_number(6);
    atom.xyz = Position::new(1.0, 2.0, 3.0);
    res.set_atom("1_CA", atom);
    chain.set_group("1", res);

    model.set_group("_", chain);
    root.set_group("model_1", model);

    let mut buf = Vec::new();
    SimpleMmcif::write_structure(&root, &mut buf, &MmcifWriteOptions::default()).unwrap();
    let cif_str = String::from_utf8(buf).unwrap();

    // Verify that chain ID in _atom_site.label_asym_id and auth_asym_id is written as '.'
    assert!(
        cif_str.contains(" GLY . ? 1 ? "),
        "expected '.' for empty chain label_asym_id"
    );
    assert!(
        cif_str.contains(" GLY . CA 1"),
        "expected '.' for empty chain auth_asym_id"
    );

    let reloaded = SimpleMmcif::from_str(&cif_str)
        .unwrap()
        .get_structure_atomgroup(None, None)
        .unwrap();

    let rel_model = reloaded.get_group("model_1").unwrap();
    assert!(
        rel_model.has_group("_"),
        "reloaded model should have chain key '_'"
    );
    assert_atomgroups_match_roundtrip(&root, &reloaded);
}

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
    assert!(buf.is_empty(), "buffer must remain empty on error");
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
    assert!(buf.is_empty(), "buffer must remain empty on error");
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
    assert!(buf.is_empty(), "buffer must remain empty on error");
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
    assert!(buf.is_empty(), "buffer must remain empty on error");
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
    atom.charge = -0.8341; // 4-decimal partial charge (RESP-like)
    atom.xyz = Position::new(1.0, 2.0, 3.0);
    res.set_atom("1_CA", atom);

    chain.set_group("1", res);
    model.set_group("A", chain);
    root.set_group("model_1", model);

    // 1. charge_to_b_factor = false -> error and buffer is empty
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
    assert!(buf.is_empty(), "buffer must remain empty on error");

    // 2. charge_to_b_factor = true -> success, written to B_iso_or_equiv with 4 decimals
    let mut buf_ok = Vec::new();
    let opts_ok = MmcifWriteOptions {
        charge_to_b_factor: true,
        ..Default::default()
    };
    SimpleMmcif::write_structure(&root, &mut buf_ok, &opts_ok)
        .expect("write should succeed with charge_to_b_factor = true");

    let cif_str = String::from_utf8(buf_ok).unwrap();
    // B_iso_or_equiv should preserve 4 decimals: -0.8341
    assert!(
        cif_str.contains(" -0.8341 ? 1 ALA A CA 1"),
        "expected -0.8341 in B_iso_or_equiv; got: {cif_str}"
    );
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
// 4. Atomic Error Handling Tests (Review Item 1)
// ========================================================================

#[test]
fn test_atomic_write_error_leaves_buffer_empty() {
    let mut root = AtomGroup::new();
    let mut model = AtomGroup::with_name("model_1");
    let mut chain = AtomGroup::with_name("A");
    let mut res = AtomGroup::with_name("ALA");

    // Atom 1: valid
    let mut a1 = Atom::new();
    a1.name = "N".to_string();
    a1.set_atomic_number(7);
    a1.charge = 0.0;
    res.set_atom("1_N", a1);

    // Atom 2: invalid partial charge (triggers mid-residue error)
    let mut a2 = Atom::new();
    a2.name = "CA".to_string();
    a2.set_atomic_number(6);
    a2.charge = 0.25;
    res.set_atom("2_CA", a2);

    chain.set_group("1", res);
    model.set_group("A", chain);
    root.set_group("model_1", model);

    let mut buf = Vec::new();
    let res = SimpleMmcif::write_structure(&root, &mut buf, &MmcifWriteOptions::default());
    assert!(res.is_err(), "expected error for partial charge");
    assert!(
        buf.is_empty(),
        "write_structure must write nothing to buffer when pre-validation fails"
    );
}

#[test]
fn test_atomic_save_error_leaves_destination_untouched() {
    let mut root = AtomGroup::new();
    let mut model = AtomGroup::with_name("model_1");
    let mut chain = AtomGroup::with_name("A");
    let mut res = AtomGroup::with_name("ALA");

    // Valid atom 1, invalid atom 2
    let mut a1 = Atom::new();
    a1.name = "N".to_string();
    a1.set_atomic_number(7);
    res.set_atom("1_N", a1);

    let mut a2 = Atom::new();
    a2.name = "CA".to_string();
    a2.set_atomic_number(6);
    a2.charge = 0.5;
    res.set_atom("2_CA", a2);

    chain.set_group("1", res);
    model.set_group("A", chain);
    root.set_group("model_1", model);

    let dest_path = std::env::temp_dir().join(format!(
        "atomic_test_{}.cif",
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));

    // Case 1: Destination does not exist -> must not be created
    let res = SimpleMmcif::save_structure(&root, &dest_path, &MmcifWriteOptions::default());
    assert!(res.is_err());
    assert!(
        !dest_path.exists(),
        "destination file must not exist after failed save_structure"
    );

    // Case 2: Destination already exists with prior content -> must remain unmodified
    std::fs::write(&dest_path, "PREVIOUS CONTENT").unwrap();
    let res = SimpleMmcif::save_structure(&root, &dest_path, &MmcifWriteOptions::default());
    assert!(res.is_err());
    let content = std::fs::read_to_string(&dest_path).unwrap();
    assert_eq!(
        content, "PREVIOUS CONTENT",
        "existing file content must not be modified after failed save_structure"
    );

    let _ = std::fs::remove_file(&dest_path);
}

// ========================================================================
// 5. Large Structure (>100,000 atoms) and Benchmark
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
// 6. Unit Tests for Helpers (Review Items 4 & 5)
// ========================================================================

#[test]
fn test_data_block_name_validation() {
    assert_eq!(validate_data_block_name("structure").unwrap(), "structure");
    assert_eq!(
        validate_data_block_name("data_structure").unwrap(),
        "structure"
    );
    assert_eq!(validate_data_block_name("1HLS").unwrap(), "1HLS");

    // Empty
    assert!(validate_data_block_name("").is_err());
    assert!(validate_data_block_name("data_").is_err());

    // Whitespace / Control
    assert!(validate_data_block_name("a b").is_err());
    assert!(validate_data_block_name("data_a\nb").is_err());

    // Invalid characters
    assert!(validate_data_block_name("data_foo#bar").is_err());
    assert!(validate_data_block_name("data_'foo'").is_err());
    assert!(validate_data_block_name("data_\"foo\"").is_err());
    assert!(validate_data_block_name("data_foo;bar").is_err());
}

#[test]
fn test_determine_cif_quote_kind_cases() {
    assert_eq!(determine_cif_quote_kind("ALA").unwrap(), CifQuoteKind::None);
    assert_eq!(
        determine_cif_quote_kind("O5'").unwrap(),
        CifQuoteKind::Double
    );
    assert_eq!(
        determine_cif_quote_kind("A B").unwrap(),
        CifQuoteKind::Single
    );
    assert_eq!(determine_cif_quote_kind("").unwrap(), CifQuoteKind::Single);
    assert_eq!(
        determine_cif_quote_kind("_tag").unwrap(),
        CifQuoteKind::Single
    );
    assert_eq!(
        determine_cif_quote_kind("data_foo").unwrap(),
        CifQuoteKind::Single
    );
    assert_eq!(
        determine_cif_quote_kind("loop_").unwrap(),
        CifQuoteKind::Single
    );
    assert_eq!(determine_cif_quote_kind(".").unwrap(), CifQuoteKind::Single);
    assert_eq!(determine_cif_quote_kind("?").unwrap(), CifQuoteKind::Single);

    // Both quotes -> Error
    assert!(determine_cif_quote_kind("a'b\"c").is_err());
}

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

// ============================================================================
// PR#43: _struct_conn writing and reader covale support tests
// ============================================================================

/// Helper to normalize an atom path to `/{model}/{chain}/{res}/{atom_name}` without serial ID prefix.
fn normalize_atom_path(ag: &AtomGroup, path: &str) -> String {
    let parts = AtomGroup::divide_path(path);
    if parts.len() < 4 {
        return path.to_string();
    }
    let model_key = &parts[0];
    let chain_key = &parts[1];
    let res_key = &parts[2];
    let atom_key = &parts[3];

    let atom_name = ag
        .get_group(model_key)
        .and_then(|m| m.get_group(chain_key))
        .and_then(|c| c.get_group(res_key))
        .and_then(|r| r.get_atom(atom_key))
        .map(|a| a.name.as_str())
        .unwrap_or(atom_key);

    format!("/{model_key}/{chain_key}/{res_key}/{atom_name}")
}

/// Helper to extract normalized inter-residue bonds `(min_path, max_path, order)` from an AtomGroup.
fn get_inter_residue_bonds(ag: &AtomGroup) -> Vec<(String, String, usize)> {
    let mut bonds = Vec::new();
    let all = ag.get_bond_list_ref();
    for b in all {
        let parts1 = AtomGroup::divide_path(&b.atom1_path);
        let parts2 = AtomGroup::divide_path(&b.atom2_path);
        if parts1.len() >= 4 && parts2.len() >= 4 {
            // Check if they belong to different residues (different chain or different residue key)
            if parts1[1] != parts2[1] || parts1[2] != parts2[2] {
                let norm1 = normalize_atom_path(ag, &b.atom1_path);
                let norm2 = normalize_atom_path(ag, &b.atom2_path);
                let pair = if norm1 <= norm2 {
                    (norm1, norm2, b.order)
                } else {
                    (norm2, norm1, b.order)
                };
                bonds.push(pair);
            }
        }
    }
    bonds.sort();
    bonds.dedup();
    bonds
}

/// PR#43 Criterion 1: Real data reading test with covale connections.
///
/// Fixture entry: 1WCT (Conotoxin derivative with O-glycosylation and modified residues)
/// Acquisition date: 2026-10-03
/// Source URL: https://files.rcsb.org/download/1WCT.cif
/// Baseline values counted directly from 1WCT.cif via grep:
///   grep -c -E "^covale[0-9]+" tests/data/1WCT.cif => 8 covale records
///   grep -c -E "^disulf[0-9]+" tests/data/1WCT.cif => 2 disulf records
/// Expected inter-residue bond pairs:
///   covale1: /model_1/A/1/C    - /model_1/A/2/N    (CGU 1 to CYS 2, order 1)
///   covale2: /model_1/A/3/C    - /model_1/A/4/N    (CYS 3 to CGU 4, order 1)
///   covale3: /model_1/A/4/C    - /model_1/A/5/N    (CGU 4 to ASP 5, order 1)
///   covale4: /model_1/A/6/C    - /model_1/A/7/N    (GLY 6 to BTR 7, order 1)
///   covale5: /model_1/A/7/C    - /model_1/A/8/N    (BTR 7 to CYS 8, order 1)
///   covale6: /model_1/A/10/OG1 - /model_1/B/1/C1   (THR 10 to NGA 1 O-Glycosylation, order 1)
///   covale7: /model_1/A/12/C   - /model_1/A/13/N   (ALA 12 to HYP 13, order 1)
///   covale8: /model_1/B/1/O3   - /model_1/B/2/C1   (NGA 1 to GAL 2 glycan bond, order 1)
///   disulf1: /model_1/A/2/SG   - /model_1/A/8/SG   (CYS 2 to CYS 8, order 1)
///   disulf2: /model_1/A/3/SG   - /model_1/A/9/SG   (CYS 3 to CYS 9, order 1)
#[test]
fn test_real_covale_loading_1wct() {
    let cif_path = test_data_dir().join("1WCT.cif");
    let cif = SimpleMmcif::from_file(&cif_path).expect("failed to load 1WCT.cif");

    // 1. Verify StructConnRecord parsing from table
    let records = cif.get_struct_conn_records();
    let covale_records: Vec<_> = records
        .iter()
        .filter(|r| r.conn_type_id == "covale")
        .collect();
    let disulf_records: Vec<_> = records
        .iter()
        .filter(|r| r.conn_type_id == "disulf")
        .collect();

    assert_eq!(
        covale_records.len(),
        8,
        "covale record count mismatch against baseline"
    );
    assert_eq!(
        disulf_records.len(),
        2,
        "disulf record count mismatch against baseline"
    );

    // Verify pdbx_value_order parsing
    let covale6 = covale_records
        .iter()
        .find(|r| r.id == "covale6")
        .expect("covale6 not found");
    assert_eq!(covale6.pdbx_value_order.as_deref(), Some("sing"));
    assert_eq!(covale6.bond_order(), 1);

    let covale8 = covale_records
        .iter()
        .find(|r| r.id == "covale8")
        .expect("covale8 not found");
    assert_eq!(covale8.pdbx_value_order.as_deref(), Some("sing"));
    assert_eq!(covale8.bond_order(), 1);

    // 2. Verify AtomGroup bond linking
    let ag = cif
        .get_structure_atomgroup(Some(1), None)
        .expect("failed to get AtomGroup for 1WCT");
    let inter_bonds = get_inter_residue_bonds(&ag);
    assert_eq!(
        inter_bonds.len(),
        10,
        "total inter-residue bond count (8 covale + 2 disulf) mismatch"
    );

    // Verify O-glycosylation bond exists (THR 10 OG1 - NGA 1 C1)
    let has_glyco = inter_bonds.iter().any(|(p1, p2, order)| {
        *order == 1
            && ((p1 == "/model_1/A/10/OG1" && p2 == "/model_1/B/1/C1")
                || (p2 == "/model_1/A/10/OG1" && p1 == "/model_1/B/1/C1"))
    });
    assert!(
        has_glyco,
        "O-glycosylation bond (THR 10 OG1 - NGA 1 C1) not found in AtomGroup"
    );

    // Verify glycan-glycan bond exists (NGA 1 O3 - GAL 2 C1)
    let has_glycan_link = inter_bonds.iter().any(|(p1, p2, order)| {
        *order == 1
            && ((p1 == "/model_1/B/1/O3" && p2 == "/model_1/B/2/C1")
                || (p2 == "/model_1/B/1/O3" && p1 == "/model_1/B/2/C1"))
    });
    assert!(
        has_glycan_link,
        "glycan-glycan bond (NGA 1 O3 - GAL 2 C1) not found in AtomGroup"
    );

    // Verify disulfide bonds exist
    let has_ss1 = inter_bonds.iter().any(|(p1, p2, _)| {
        (p1 == "/model_1/A/2/SG" && p2 == "/model_1/A/8/SG")
            || (p2 == "/model_1/A/2/SG" && p1 == "/model_1/A/8/SG")
    });
    let has_ss2 = inter_bonds.iter().any(|(p1, p2, _)| {
        (p1 == "/model_1/A/3/SG" && p2 == "/model_1/A/9/SG")
            || (p2 == "/model_1/A/3/SG" && p1 == "/model_1/A/9/SG")
    });
    assert!(has_ss1, "disulfide bond CYS 2 - CYS 8 not found");
    assert!(has_ss2, "disulfide bond CYS 3 - CYS 9 not found");
}

/// PR#43 Criterion 2: Coexistence of file-derived inter-residue bonds with ag.setup().
///
/// Verifies for 1WCT.cif, 1HLS.cif, and 2FB4.cif that:
/// - File-derived inter-residue bonds are strictly preserved without duplication.
/// - CCD templates add intra-residue chemical bonds, increasing total bond count.
#[test]
fn test_setup_coexistence_with_file_bonds() {
    let test_cases = ["1WCT.cif", "1HLS.cif", "2FB4.cif"];

    for filename in test_cases {
        let cif_path = test_data_dir().join(filename);
        let cif = SimpleMmcif::from_file(&cif_path)
            .unwrap_or_else(|e| panic!("failed to load {filename}: {e}"));
        let mut ag = cif
            .get_structure_atomgroup(Some(1), None)
            .unwrap_or_else(|e| panic!("failed to get structure atomgroup for {filename}: {e}"));

        let initial_inter_bonds = get_inter_residue_bonds(&ag);
        let initial_total_bonds = ag.get_bond_list_ref().len();

        // Perform setup()
        ag.setup()
            .unwrap_or_else(|e| panic!("ag.setup() failed for {filename}: {e}"));

        let post_inter_bonds = get_inter_residue_bonds(&ag);
        let post_total_bonds = ag.get_bond_list_ref().len();

        // 1. Total bond count must increase (CCD templates add intra-residue bonds)
        assert!(
            post_total_bonds > initial_total_bonds,
            "{filename}: setup() should add intra-residue bonds (initial={initial_total_bonds}, post={post_total_bonds})"
        );

        // 2. All initial file-derived inter-residue bonds must still be present with identical order
        for init_bond in &initial_inter_bonds {
            assert!(
                post_inter_bonds.contains(init_bond),
                "{filename}: file-derived bond {:?} was lost after setup()",
                init_bond
            );
        }

        // 3. No duplicate bonds between any atom pair
        let all_post_bonds = ag.get_bond_list_ref();
        let mut seen_pairs = std::collections::HashSet::new();
        for b in &all_post_bonds {
            let key = if b.atom1_path <= b.atom2_path {
                (b.atom1_path.clone(), b.atom2_path.clone())
            } else {
                (b.atom2_path.clone(), b.atom1_path.clone())
            };
            assert!(
                seen_pairs.insert(key.clone()),
                "{filename}: duplicate bond found for atom pair {:?}",
                key
            );
        }
    }
}

/// PR#43 Criterion 3: Roundtrip test for inter-residue bonds in 1WCT.cif (covale + disulf).
#[test]
fn test_roundtrip_inter_residue_bonds_1wct() {
    let cif_path = test_data_dir().join("1WCT.cif");
    let cif = SimpleMmcif::from_file(&cif_path).expect("failed to load 1WCT.cif");
    let ag = cif
        .get_structure_atomgroup(Some(1), None)
        .expect("failed to get AtomGroup");

    let orig_inter_bonds = get_inter_residue_bonds(&ag);
    assert_eq!(orig_inter_bonds.len(), 10);

    // Write to memory buffer
    let mut buf = Vec::new();
    let opts = MmcifWriteOptions::default();
    SimpleMmcif::write_structure(&ag, &mut buf, &opts).expect("failed to write 1WCT mmCIF");

    let text = String::from_utf8(buf).expect("invalid utf-8 output");
    assert!(
        text.contains("_struct_conn.id"),
        "_struct_conn loop should be present"
    );
    assert!(text.contains("disulf1"));
    assert!(text.contains("covale1"));

    // Reload mmCIF from text
    let reloaded_cif = SimpleMmcif::from_str(&text).expect("failed to reload written mmCIF");
    let reloaded_ag = reloaded_cif
        .get_structure_atomgroup(Some(1), None)
        .expect("failed to get reloaded AtomGroup");

    let reloaded_inter_bonds = get_inter_residue_bonds(&reloaded_ag);

    assert_eq!(
        orig_inter_bonds.len(),
        reloaded_inter_bonds.len(),
        "roundtrip inter-residue bond count mismatch"
    );

    for (p1, p2, order) in &orig_inter_bonds {
        let found = reloaded_inter_bonds.iter().any(|(rp1, rp2, rorder)| {
            *rorder == *order && ((rp1 == p1 && rp2 == p2) || (rp1 == p2 && rp2 == p1))
        });
        assert!(
            found,
            "original bond ({p1}, {p2}, order {order}) missing in reloaded structure"
        );
    }
}

/// PR#43 Criterion 3: Roundtrip test for disulfide bonds in 1HLS.cif (NMR 20 models, disulf).
#[test]
fn test_roundtrip_disulf_1hls() {
    let cif_path = test_data_dir().join("1HLS.cif");
    let cif = SimpleMmcif::from_file(&cif_path).expect("failed to load 1HLS.cif");
    let ag = cif
        .get_structure_atomgroup(None, None)
        .expect("failed to get AtomGroup");

    let orig_inter_bonds = get_inter_residue_bonds(&ag);
    assert_eq!(
        orig_inter_bonds.len(),
        60,
        "1HLS has 20 models with 3 disulfide bonds each = 60 bonds"
    );

    let mut buf = Vec::new();
    let opts = MmcifWriteOptions::default();
    SimpleMmcif::write_structure(&ag, &mut buf, &opts).expect("failed to write 1HLS mmCIF");

    let text = String::from_utf8(buf).expect("invalid utf-8 output");
    assert!(text.contains("_struct_conn.id"));
    assert!(text.contains("disulf1"));

    let reloaded_cif = SimpleMmcif::from_str(&text).expect("failed to reload 1HLS mmCIF");
    let reloaded_ag = reloaded_cif
        .get_structure_atomgroup(None, None)
        .expect("failed to get reloaded AtomGroup");

    let reloaded_inter_bonds = get_inter_residue_bonds(&reloaded_ag);
    assert_eq!(orig_inter_bonds, reloaded_inter_bonds);
}

/// PR#43 Criterion 3: Roundtrip test for disulfide bonds in 2FB4.cif (insertion codes, disulf).
#[test]
fn test_roundtrip_disulf_2fb4() {
    let cif_path = test_data_dir().join("2FB4.cif");
    let cif = SimpleMmcif::from_file(&cif_path).expect("failed to load 2FB4.cif");
    let ag = cif
        .get_structure_atomgroup(None, None)
        .expect("failed to get AtomGroup");

    let orig_inter_bonds = get_inter_residue_bonds(&ag);
    assert_eq!(orig_inter_bonds.len(), 6, "2FB4 has 6 disulfide bonds");

    let mut buf = Vec::new();
    let opts = MmcifWriteOptions::default();
    SimpleMmcif::write_structure(&ag, &mut buf, &opts).expect("failed to write 2FB4 mmCIF");

    let text = String::from_utf8(buf).expect("invalid utf-8 output");
    assert!(text.contains("_struct_conn.id"));
    assert!(text.contains("disulf1"));

    let reloaded_cif = SimpleMmcif::from_str(&text).expect("failed to reload 2FB4 mmCIF");
    let reloaded_ag = reloaded_cif
        .get_structure_atomgroup(None, None)
        .expect("failed to get reloaded AtomGroup");

    let reloaded_inter_bonds = get_inter_residue_bonds(&reloaded_ag);
    assert_eq!(orig_inter_bonds, reloaded_inter_bonds);
}

/// PR#43 Criterion 4: Synthetic data test:
/// - Standard peptide bonds (ALA C - GLY N) are excluded from _struct_conn.
/// - Standard nucleic backbone bonds (DA O3' - DT P) are excluded from _struct_conn.
/// - Intra-residue bonds (CA - CB) are excluded from _struct_conn.
/// - Double bond (order: 2) between different residues is exported as `doub` and preserved on roundtrip.
/// - Non-standard residue covalent bond is exported as `covale`.
#[test]
fn test_synthetic_struct_conn_filtering_and_order() {
    let mut root = AtomGroup::new();
    let mut model = AtomGroup::new();
    let mut chain_a = AtomGroup::new();

    // Residue 1: ALA (standard amino acid)
    let mut r1 = AtomGroup::new();
    r1.name = "ALA".to_string();
    let mut a1_ca = Atom::new();
    a1_ca.name = "CA".to_string();
    a1_ca.xyz = Position::new(0.0, 0.0, 0.0);
    let mut a1_c = Atom::new();
    a1_c.name = "C".to_string();
    a1_c.xyz = Position::new(1.0, 0.0, 0.0);
    r1.set_atom("1_CA", a1_ca.clone());
    r1.set_atom("1_C", a1_c.clone());
    // Intra-residue bond CA-C
    r1.add_bond(&a1_ca, &a1_c, 1);

    // Residue 2: GLY (standard amino acid)
    let mut r2 = AtomGroup::new();
    r2.name = "GLY".to_string();
    let mut a2_n = Atom::new();
    a2_n.name = "N".to_string();
    a2_n.xyz = Position::new(2.0, 0.0, 0.0);
    let mut a2_ca = Atom::new();
    a2_ca.name = "CA".to_string();
    a2_ca.xyz = Position::new(3.0, 0.0, 0.0);
    r2.set_atom("2_N", a2_n.clone());
    r2.set_atom("2_CA", a2_ca.clone());

    // Residue 3: DA (standard DNA nucleotide)
    let mut r3 = AtomGroup::new();
    r3.name = "DA".to_string();
    let mut a3_o3 = Atom::new();
    a3_o3.name = "O3'".to_string();
    a3_o3.xyz = Position::new(4.0, 0.0, 0.0);
    r3.set_atom("3_O3'", a3_o3.clone());

    // Residue 4: DT (standard DNA nucleotide)
    let mut r4 = AtomGroup::new();
    r4.name = "DT".to_string();
    let mut a4_p = Atom::new();
    a4_p.name = "P".to_string();
    a4_p.xyz = Position::new(5.0, 0.0, 0.0);
    r4.set_atom("4_P", a4_p.clone());

    // Residue 5: LIG (non-standard ligand) with a double bond to GLY CA
    let mut r5 = AtomGroup::new();
    r5.name = "LIG".to_string();
    let mut a5_x = Atom::new();
    a5_x.name = "X1".to_string();
    a5_x.xyz = Position::new(6.0, 0.0, 0.0);
    r5.set_atom("5_X1", a5_x.clone());

    chain_a.set_group("1", r1);
    chain_a.set_group("2", r2);
    chain_a.set_group("3", r3);
    chain_a.set_group("4", r4);
    chain_a.set_group("5", r5);

    model.set_group("A", chain_a);
    root.set_group("model_1", model);

    // Add inter-residue bonds using atoms resolved from the hierarchy (so their paths are populated)
    let a1_c = root.get_atom_by_path("/model_1/A/1/1_C").unwrap().clone();
    let a2_n = root.get_atom_by_path("/model_1/A/2/2_N").unwrap().clone();
    let a3_o3 = root.get_atom_by_path("/model_1/A/3/3_O3'").unwrap().clone();
    let a4_p = root.get_atom_by_path("/model_1/A/4/4_P").unwrap().clone();
    let a2_ca = root.get_atom_by_path("/model_1/A/2/2_CA").unwrap().clone();
    let a5_x = root.get_atom_by_path("/model_1/A/5/5_X1").unwrap().clone();

    let model_mut = root.get_group_mut("model_1").unwrap();

    // 1. Peptide bond: ALA 1 C - GLY 2 N (should be excluded from _struct_conn)
    model_mut.add_bond(&a1_c, &a2_n, 1);

    // 2. Nucleic backbone bond: DA 3 O3' - DT 4 P (should be excluded from _struct_conn)
    model_mut.add_bond(&a3_o3, &a4_p, 1);

    // 3. Inter-residue double bond: GLY 2 CA - LIG 5 X1 (order: 2, should be exported as covale with doub)
    model_mut.add_bond(&a2_ca, &a5_x, 2);

    let mut buf = Vec::new();
    let opts = MmcifWriteOptions::default();
    SimpleMmcif::write_structure(&root, &mut buf, &opts)
        .expect("failed to write synthetic structure");

    let text = String::from_utf8(buf).expect("invalid utf-8");

    // Extract the _struct_conn section from the written CIF text
    let struct_conn_section = if let Some((_, conn_part)) = text.split_once("_struct_conn.id") {
        conn_part
    } else {
        panic!("_struct_conn section not found in written output:\n{text}");
    };

    // Peptide bond (ALA C - GLY N) and nucleic backbone (O3' - P) must NOT be in _struct_conn
    assert!(
        !struct_conn_section.contains("ALA"),
        "peptide bond should be excluded from _struct_conn"
    );
    assert!(
        !struct_conn_section.contains("O3'"),
        "nucleic backbone bond should be excluded from _struct_conn"
    );

    // The double bond must be present as covale1 and doub
    assert!(
        struct_conn_section.contains("covale1 covale doub"),
        "expected 'covale1 covale doub' in _struct_conn section:\n{struct_conn_section}"
    );
    assert!(struct_conn_section.contains("LIG"));

    // Reload and verify bond order is preserved as 2 (doub)
    let reloaded_cif = SimpleMmcif::from_str(&text).expect("failed to reload synthetic mmCIF");
    let reloaded_ag = reloaded_cif
        .get_structure_atomgroup(Some(1), None)
        .expect("failed to get reloaded synthetic AtomGroup");

    let inter_bonds = get_inter_residue_bonds(&reloaded_ag);
    assert_eq!(
        inter_bonds.len(),
        1,
        "only the LIG-GLY bond should exist as an inter-residue bond"
    );
    assert_eq!(inter_bonds[0].2, 2, "bond order 2 (doub) must be preserved");
}

/// Verifies AtomGroup::get_bond_list_ref produces identical records to get_bond_list.
#[test]
fn test_get_bond_list_ref_matches_get_bond_list() {
    let cif_path = test_data_dir().join("1WCT.cif");
    let cif = SimpleMmcif::from_file(&cif_path).expect("failed to load 1WCT.cif");
    let mut ag = cif
        .get_structure_atomgroup(Some(1), None)
        .expect("failed to get AtomGroup");

    let ref_bonds = ag.get_bond_list_ref();
    let mut_bonds = ag.get_bond_list();

    assert_eq!(ref_bonds, mut_bonds);
}

/// PR#42 / PR#43 invariant: If validation fails before writing (e.g., due to invalid characters in struct_conn),
/// nothing should be written to the writer and it remains empty.
#[test]
fn test_struct_conn_validation_failure_leaves_buffer_empty() {
    let mut root = AtomGroup::new();
    root.name = "root".to_string();
    let mut model = AtomGroup::new();
    model.name = "model_1".to_string();
    let mut chain = AtomGroup::new();
    chain.name = "A".to_string();
    let mut r1 = AtomGroup::new();
    r1.name = "RES1".to_string();
    let mut r2 = AtomGroup::new();
    r2.name = "RES2".to_string();

    let mut a1 = Atom::new();
    a1.name = "CA".to_string();
    a1.set_atomic_number(6);
    a1.xyz = Position::new(0.0, 0.0, 0.0);

    // Introduce mixed quotes (' and ") in the atom name which cannot be safely quoted in mmCIF
    let mut a2 = Atom::new();
    a2.name = "N'\"BAD".to_string();
    a2.set_atomic_number(7);
    a2.xyz = Position::new(1.0, 0.0, 0.0);

    r1.set_atom("1_CA", a1);
    r2.set_atom("2_N'\"BAD", a2);
    chain.set_group("1", r1);
    chain.set_group("2", r2);
    model.set_group("A", chain);
    root.set_group("model_1", model);

    let a1_ref = root.get_atom_by_path("/model_1/A/1/1_CA").unwrap().clone();
    let a2_ref = root
        .get_atom_by_path("/model_1/A/2/2_N'\"BAD")
        .unwrap()
        .clone();
    root.get_group_mut("model_1")
        .unwrap()
        .add_bond(&a1_ref, &a2_ref, 1);

    let mut buf = Vec::new();
    let result = SimpleMmcif::write_structure(&root, &mut buf, &MmcifWriteOptions::default());
    assert!(
        result.is_err(),
        "validation should fail due to mixed quotes in atom name for struct_conn"
    );
    assert!(
        buf.is_empty(),
        "buffer must remain completely empty on validation failure"
    );
}
