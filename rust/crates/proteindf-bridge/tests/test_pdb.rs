// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::path::PathBuf;

use proteindf_bridge::atom_group::AtomGroup;
use proteindf_bridge::format::pdb::Pdb;

fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

#[test]
fn test_2mgo_hierarchy_and_counts() {
    let path = test_data_dir().join("2MGO.pdb");
    let pdb = Pdb::from_file(&path, None).expect("failed to load 2MGO.pdb");

    let mut ag = pdb
        .get_atomgroup(None, None)
        .expect("failed to get AtomGroup");

    // 20 models
    assert_eq!(ag.get_number_of_groups(), 20);

    // Check model_1 structure
    let m1 = ag.get_group("model_1").expect("model_1 not found");
    assert_eq!(m1.get_number_of_groups(), 1); // chain A

    let chain_a = m1.get_group("A").expect("chain A not found");
    assert_eq!(chain_a.get_number_of_groups(), 9); // 9 residues

    // Check each residue atom count and name (matching Python tests)
    let res_checks = [
        ("1", "CYS", 12),
        ("2", "TYR", 21),
        ("3", "ILE", 19),
        ("4", "GLN", 17),
        ("5", "ASN", 14),
        ("6", "CYS", 10),
        ("7", "PRO", 14),
        ("8", "LEU", 19),
        ("9", "GLY", 8),
    ];

    for (seq, expected_name, expected_atom_count) in res_checks {
        let res = chain_a
            .get_group(seq)
            .unwrap_or_else(|| panic!("residue {seq} not found"));
        assert_eq!(res.name, expected_name, "residue {seq} name mismatch");
        assert_eq!(
            res.get_number_of_atoms(),
            expected_atom_count,
            "residue {seq} atom count mismatch"
        );
    }

    // Total atom count across the entire ensemble (20 models * 134 atoms)
    assert_eq!(ag.get_atom_list().len(), 2680);

    // Verify SSBOND disulfide bond linking:
    // Each of the 20 models should have 1 disulfide bond routed within chain A (CYS1 SG <-> CYS6 SG).
    let bond_list = ag.get_bond_list();
    assert_eq!(
        bond_list.len(),
        20,
        "expected 20 SSBOND linkages across 20 models"
    );

    for i in 1..=20 {
        let expected_a1 = format!("/model_{i}/A/1/6_SG");
        let expected_a2 = format!("/model_{i}/A/6/89_SG");
        let found = bond_list.iter().any(|b| {
            (b.atom1_path == expected_a1 && b.atom2_path == expected_a2)
                || (b.atom1_path == expected_a2 && b.atom2_path == expected_a1)
        });
        assert!(
            found,
            "SSBOND for model_{i} between CYS1 SG and CYS6 SG not found in bond list"
        );
    }
}

#[test]
fn test_1hls_real_pdb() {
    let path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&path, None).expect("failed to load 1hls.pdb");

    let ag = pdb
        .get_atomgroup(None, None)
        .expect("failed to get AtomGroup");

    // 1 model: model_1
    assert_eq!(ag.get_number_of_groups(), 1);
    let m1 = ag.get_group("model_1").expect("model_1 not found");

    // 2 chains: A and B
    assert_eq!(m1.get_number_of_groups(), 2);
    let chain_a = m1.get_group("A").expect("chain A not found");
    let chain_b = m1.get_group("B").expect("chain B not found");

    // Residue counts
    assert_eq!(chain_a.get_number_of_groups(), 21);
    assert_eq!(chain_b.get_number_of_groups(), 30);

    // Atom counts per chain
    assert_eq!(chain_a.get_atom_list().len(), 312);
    assert_eq!(chain_b.get_atom_list().len(), 470);

    // Total atoms
    let atoms = ag.get_atom_list();
    assert_eq!(atoms.len(), 782);

    // Verify first atom properties (matches Python: N N 0.0 (-2.414, 8.071, 6.020))
    let first = &atoms[0];
    assert_eq!(first.name, "N");
    assert_eq!(first.symbol().unwrap(), "N");
    assert_eq!(first.charge, 0.0);
    assert!((first.xyz.x - (-2.414)).abs() < 1e-5);
    assert!((first.xyz.y - 8.071).abs() < 1e-5);
    assert!((first.xyz.z - 6.020).abs() < 1e-5);
}

#[test]
fn test_3i3zh_real_pdb() {
    let path = test_data_dir().join("3i3zH.pdb");
    let pdb = Pdb::from_file(&path, None).expect("failed to load 3i3zH.pdb");

    let ag = pdb
        .get_atomgroup(None, None)
        .expect("failed to get AtomGroup");

    // 1 model: model_1
    assert_eq!(ag.get_number_of_groups(), 1);
    let m1 = ag.get_group("model_1").expect("model_1 not found");

    // 2 chains: A and B
    assert_eq!(m1.get_number_of_groups(), 2);
    let chain_a = m1.get_group("A").expect("chain A not found");
    let chain_b = m1.get_group("B").expect("chain B not found");

    // Residue counts
    assert_eq!(chain_a.get_number_of_groups(), 21);
    assert_eq!(chain_b.get_number_of_groups(), 30);

    // Atom counts per chain
    assert_eq!(chain_a.get_atom_list().len(), 312);
    assert_eq!(chain_b.get_atom_list().len(), 468);

    // Total atoms
    let atoms = ag.get_atom_list();
    assert_eq!(atoms.len(), 780);

    // First atom: H1 H 0.0 (-7.107, -26.611, 8.235)
    let first = &atoms[0];
    assert_eq!(first.name, "H1");
    assert_eq!(first.symbol().unwrap(), "H");
    assert_eq!(first.charge, 0.0);
    assert!((first.xyz.x - (-7.107)).abs() < 1e-5);
    assert!((first.xyz.y - (-26.611)).abs() < 1e-5);
    assert!((first.xyz.z - 8.235).abs() < 1e-5);
}

#[test]
fn test_set_by_atomgroup_roundtrip() {
    let path = test_data_dir().join("1hls.pdb");
    let pdb_orig = Pdb::from_file(&path, None).expect("failed to load 1hls.pdb");
    let ag_orig = pdb_orig
        .get_atomgroup(None, None)
        .expect("get_atomgroup failed");

    let mut pdb_rebuilt = Pdb::new(None);
    pdb_rebuilt
        .set_by_atomgroup(&ag_orig, false)
        .expect("set_by_atomgroup failed");

    let text = pdb_rebuilt.get_text();
    assert!(!text.is_empty());
    assert!(text.contains("MODEL        1"));
    assert!(text.contains("ATOM  "));
    assert!(text.contains("TER   "));

    // Parse the generated text back and verify counts
    let pdb_parsed = Pdb::from_str(&text, None).expect("failed to parse generated text");
    let ag_parsed = pdb_parsed
        .get_atomgroup(None, None)
        .expect("get_atomgroup failed");

    assert_eq!(
        ag_parsed.get_atom_list().len(),
        ag_orig.get_atom_list().len()
    );
    let m_orig = ag_orig.get_group("model_1").unwrap();
    let m_parsed = ag_parsed.get_group("model_1").unwrap();
    assert_eq!(
        m_parsed.get_number_of_groups(),
        m_orig.get_number_of_groups()
    );
}

#[test]
fn test_amber_mode_modpdb() {
    let pdb = Pdb::new(Some("AMBER"));
    let mut ag = AtomGroup::new();

    // Create model -> chain -> residue HIS with HD1, HD2, HE1, HE2 -> HIP
    let mut model = AtomGroup::new();
    let mut chain = AtomGroup::new();
    let mut res_his = AtomGroup::new();
    res_his.name = "HIS".to_string();

    for name in &["HD1", "HD2", "HE1", "HE2"] {
        let mut atm = proteindf_bridge::atom::Atom::new();
        atm.name = (*name).to_string();
        res_his.set_atom(name, atm);
    }

    chain.set_group("1", res_his);
    model.set_group("A", chain);
    ag.set_group("model_1", model);

    let mod_ag = pdb.get_modpdb_atomgroup(&ag);
    let his_res = mod_ag
        .get_group("model_1")
        .unwrap()
        .get_group("A")
        .unwrap()
        .get_group("1")
        .unwrap();
    assert_eq!(his_res.name, "HIP");
}

#[test]
fn test_renumber() {
    let path = test_data_dir().join("1hls.pdb");
    let mut pdb = Pdb::from_file(&path, None).expect("failed to load");
    pdb.renumber();

    if let Some(records) = pdb.data().get(&1) {
        for (i, rec) in records.iter().enumerate() {
            assert_eq!(rec.serial, i + 1);
        }
    }
}

#[test]
fn test_error_handling() {
    let mut pdb = Pdb::new(None);
    let res = pdb.load("non_existent_file.pdb");
    assert!(res.is_err());

    let invalid_content = "ATOM  abcde";
    let mut pdb_err = Pdb::new(None);
    let res2 = pdb_err.parse_str(invalid_content);
    assert!(res2.is_err());
}

#[test]
fn test_occupancy_temp_factor_invalid_error() {
    // Valid 80-column line for reference:
    // "ATOM      1  N   CYS A   1       4.874   2.855   0.366  1.00  0.00           N  "
    // Column 55-60 (0-indexed 54..60): occupancy "  1.00"
    // Column 61-66 (0-indexed 60..66): temp_factor "  0.00"

    // 1. Invalid occupancy (non-empty, non-numeric) should propagate error
    let invalid_occ =
        "ATOM      1  N   CYS A   1       4.874   2.855   0.366  XXXX  0.00           N  ";
    let mut pdb = Pdb::new(None);
    let err = pdb.parse_str(invalid_occ).unwrap_err();
    assert!(
        err.to_string().contains("ATOM occupancy"),
        "error should mention ATOM occupancy: {err}"
    );

    // 2. Invalid temp_factor (non-empty, non-numeric) should propagate error
    let invalid_temp =
        "ATOM      1  N   CYS A   1       4.874   2.855   0.366  1.00  YYYY           N  ";
    let mut pdb = Pdb::new(None);
    let err = pdb.parse_str(invalid_temp).unwrap_err();
    assert!(
        err.to_string().contains("ATOM temp_factor"),
        "error should mention ATOM temp_factor: {err}"
    );

    // 3. Blank occupancy and temp_factor should use default values (1.0, 0.0) without error
    let blank_occ_temp =
        "ATOM      1  N   CYS A   1       4.874   2.855   0.366                      N  ";
    let mut pdb = Pdb::new(None);
    pdb.parse_str(blank_occ_temp)
        .expect("blank occupancy and temp_factor should be valid");
    let records = pdb.data().get(&1).unwrap();
    assert_eq!(records[0].occupancy, 1.0);
    assert_eq!(records[0].temp_factor, 0.0);
}

#[test]
fn test_2mgo_real_pdb_conect() {
    let path = test_data_dir().join("2MGO.pdb");
    let pdb = Pdb::from_file(&path, None).expect("failed to load 2MGO.pdb");

    // Check CONECT record parsing from real PDB fixture
    assert_eq!(pdb.conects().len(), 1);
    assert_eq!(pdb.conects(), &[(6, 89)]);

    let mut ag = pdb
        .get_atomgroup(None, None)
        .expect("failed to get AtomGroup");

    // Verify deduplication between SSBOND and CONECT:
    // Exactly 1 bond per model (20 total across 20 models), not duplicated to 40.
    let bond_list = ag.get_bond_list();
    assert_eq!(
        bond_list.len(),
        20,
        "expected 20 bonds across 20 models (no duplicate between SSBOND and CONECT)"
    );
}

#[test]
fn test_conect_multiple_partners_synthetic() {
    // Synthetic PDB with central carbon (serial 1) bonded to 4 hydrogens (serials 2, 3, 4, 5)
    // in a single CONECT line: CONECT    1    2    3    4    5
    // along with reverse CONECT lines to test deduplication.
    let pdb_content = "\
ATOM      1  C   MOL A   1       0.000   0.000   0.000  1.00  0.00           C  
ATOM      2  H1  MOL A   1       1.000   0.000   0.000  1.00  0.00           H  
ATOM      3  H2  MOL A   1       0.000   1.000   0.000  1.00  0.00           H  
ATOM      4  H3  MOL A   1       0.000   0.000   1.000  1.00  0.00           H  
ATOM      5  H4  MOL A   1      -1.000   0.000   0.000  1.00  0.00           H  
CONECT    1    2    3    4    5
CONECT    2    1
CONECT    3    1
CONECT    4    1
CONECT    5    1
";
    let pdb = Pdb::from_str(pdb_content, None).expect("failed to parse pdb string");

    // 4 unique bonds: (1, 2), (1, 3), (1, 4), (1, 5)
    assert_eq!(pdb.conects().len(), 4);
    assert_eq!(pdb.conects(), &[(1, 2), (1, 3), (1, 4), (1, 5)]);

    let mut ag = pdb
        .get_atomgroup(None, None)
        .expect("failed to get AtomGroup");

    let bond_list = ag.get_bond_list();
    assert_eq!(bond_list.len(), 4);

    for record in &bond_list {
        let (a1, a2) = ag.resolve_bond(record).expect("bond must resolve");
        match (a1.name.as_str(), a2.name.as_str()) {
            ("C", "H1") | ("H1", "C") => assert_eq!(record.order, 1),
            ("C", "H2") | ("H2", "C") => assert_eq!(record.order, 1),
            ("C", "H3") | ("H3", "C") => assert_eq!(record.order, 1),
            ("C", "H4") | ("H4", "C") => assert_eq!(record.order, 1),
            other => panic!("unexpected bond pair: {:?}", other),
        }
    }
}

#[test]
fn test_conect_invalid_error() {
    let invalid_conect = "CONECT   XX    1\n";
    let mut pdb = Pdb::new(None);
    let err = pdb.parse_str(invalid_conect).unwrap_err();
    assert!(err.to_string().contains("CONECT serial"));

    let invalid_partner = "CONECT    1   YY\n";
    let mut pdb2 = Pdb::new(None);
    let err2 = pdb2.parse_str(invalid_partner).unwrap_err();
    assert!(err2.to_string().contains("CONECT partner serial"));
}
