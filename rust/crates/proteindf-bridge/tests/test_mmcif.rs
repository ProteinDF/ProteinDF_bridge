// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::path::PathBuf;

use proteindf_bridge::format::pdb::Pdb;
use proteindf_bridge::format::SimpleMmcif;

#[test]
fn test_load_mmcif_basic() {
    let mmcif_content = r#"data_test_entry
_entry.id test_entry
_cell.length_a 10.0
_cell.length_b 20.0
_cell.length_c 30.0
"#;

    let mut cif = SimpleMmcif::new();
    cif.load_from_str(mmcif_content).unwrap();

    assert!(cif.data().contains_key("data_test_entry"));
    let block = cif.get_data_block("data_test_entry").unwrap();
    assert_eq!(
        block.key_values.get("_entry.id").map(|s| s.as_str()),
        Some("test_entry")
    );
    assert_eq!(
        block.key_values.get("_cell.length_a").map(|s| s.as_str()),
        Some("10.0")
    );
    assert_eq!(
        block.key_values.get("_cell.length_b").map(|s| s.as_str()),
        Some("20.0")
    );
    assert_eq!(
        block.key_values.get("_cell.length_c").map(|s| s.as_str()),
        Some("30.0")
    );
}

#[test]
fn test_ala_ccd_atomgroup() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/ALA.cif");
    let cif = SimpleMmcif::from_file(&path).expect("failed to load ALA.cif");

    let mol_names = cif.get_molecule_names();
    assert_eq!(mol_names, vec!["data_ALA".to_string()]);

    let ag = cif
        .get_atomgroup("data_ALA")
        .expect("failed to get ALA atomgroup");

    // Check molecule/residue name
    assert_eq!(ag.name, "ALA");

    // Check atom count: 13 atoms
    assert_eq!(ag.get_number_of_atoms(), 13);

    // Check expected coordinates for key atoms (from pdbx_model_Cartn_*_ideal)
    let n = ag.get_atom("N").expect("N atom missing");
    assert_eq!(n.name, "N");
    assert_eq!(n.symbol().unwrap(), "N");
    assert!((n.xyz.x - (-0.966)).abs() < 1e-4);
    assert!((n.xyz.y - 0.493).abs() < 1e-4);
    assert!((n.xyz.z - 1.500).abs() < 1e-4);

    let ca = ag.get_atom("CA").expect("CA atom missing");
    assert_eq!(ca.name, "CA");
    assert_eq!(ca.symbol().unwrap(), "C");
    assert!((ca.xyz.x - 0.257).abs() < 1e-4);
    assert!((ca.xyz.y - 0.418).abs() < 1e-4);
    assert!((ca.xyz.z - 0.692).abs() < 1e-4);

    let c = ag.get_atom("C").expect("C atom missing");
    assert_eq!(c.name, "C");
    assert_eq!(c.symbol().unwrap(), "C");
    assert!((c.xyz.x - (-0.094)).abs() < 1e-4);
    assert!((c.xyz.y - 0.017).abs() < 1e-4);
    assert!((c.xyz.z - (-0.716)).abs() < 1e-4);

    let o = ag.get_atom("O").expect("O atom missing");
    assert_eq!(o.name, "O");
    assert_eq!(o.symbol().unwrap(), "O");
    assert!((o.xyz.x - (-1.056)).abs() < 1e-4);
    assert!((o.xyz.y - (-0.682)).abs() < 1e-4);
    assert!((o.xyz.z - (-0.923)).abs() < 1e-4);

    let cb = ag.get_atom("CB").expect("CB atom missing");
    assert_eq!(cb.symbol().unwrap(), "C");
    assert!((cb.xyz.x - 1.204).abs() < 1e-4);
    assert!((cb.xyz.y - (-0.620)).abs() < 1e-4);
    assert!((cb.xyz.z - 1.296).abs() < 1e-4);

    let oxt = ag.get_atom("OXT").expect("OXT atom missing");
    assert_eq!(oxt.symbol().unwrap(), "O");
    assert!((oxt.xyz.x - 0.661).abs() < 1e-4);
    assert!((oxt.xyz.y - 0.439).abs() < 1e-4);
    assert!((oxt.xyz.z - (-1.742)).abs() < 1e-4);

    // Check bonds: 12 bonds total
    let bonds = ag.bonds();
    assert_eq!(bonds.len(), 12);

    // Verify bond orders: C=O is bond order 2, all other 11 bonds are order 1
    let mut double_bonds = 0;
    let mut single_bonds = 0;
    for bond in bonds {
        let p1 = bond.atom1_path.trim_start_matches('/');
        let p2 = bond.atom2_path.trim_start_matches('/');
        if bond.order == 2 {
            double_bonds += 1;
            // Should be between C and O
            assert!(
                (p1 == "C" && p2 == "O") || (p1 == "O" && p2 == "C"),
                "Unexpected double bond between {} and {}",
                p1,
                p2
            );
        } else if bond.order == 1 {
            single_bonds += 1;
        }
    }
    assert_eq!(double_bonds, 1);
    assert_eq!(single_bonds, 11);
}

#[test]
fn test_deuterium_and_coordinate_fallback() {
    let mmcif_content = r#"data_D_test
_chem_comp.id D_test
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.model_Cartn_x
_chem_comp_atom.model_Cartn_y
_chem_comp_atom.model_Cartn_z
D_test D1 D 1.5 2.5 3.5
D_test C1 C 0.0 0.0 0.0
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
D_test D1 C1 SING
"#;

    let mut cif = SimpleMmcif::new();
    cif.load_from_str(mmcif_content).unwrap();
    let ag = cif.get_atomgroup("data_D_test").unwrap();

    assert_eq!(ag.name, "D_test");
    assert_eq!(ag.get_number_of_atoms(), 2);

    let d1 = ag.get_atom("D1").unwrap();
    // Deuterium "D" should be mapped to Hydrogen "H" (symbol H, atomic number 1)
    assert_eq!(d1.symbol().unwrap(), "H");
    assert_eq!(d1.atomic_number(), 1);
    // Should fallback to model_Cartn_* coordinates (1.5, 2.5, 3.5)
    assert!((d1.xyz.x - 1.5).abs() < 1e-6);
    assert!((d1.xyz.y - 2.5).abs() < 1e-6);
    assert!((d1.xyz.z - 3.5).abs() < 1e-6);

    let bonds = ag.bonds();
    assert_eq!(bonds.len(), 1);
    assert_eq!(bonds[0].order, 1);
}

#[test]
fn test_error_handling() {
    let cif = SimpleMmcif::new();
    let result = cif.get_atomgroup("non_existent");
    assert!(result.is_err());
}

#[test]
fn test_1hls_cif_matches_pdb() {
    let data_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let cif_path = data_dir.join("1HLS.cif");
    let pdb_path = data_dir.join("1hls.pdb");

    let cif = SimpleMmcif::from_file(&cif_path).expect("failed to load 1HLS.cif");
    let pdb = Pdb::from_file(&pdb_path, None).expect("failed to load 1hls.pdb");

    // Select model 1 from both
    let cif_ag = cif
        .get_structure_atomgroup(Some(1), None)
        .expect("failed to get cif atomgroup");
    let pdb_ag = pdb
        .get_atomgroup(Some(1), None)
        .expect("failed to get pdb atomgroup");

    // Check model 1 exists in both
    let cif_m1 = cif_ag.get_group("model_1").expect("cif model_1 not found");
    let pdb_m1 = pdb_ag.get_group("model_1").expect("pdb model_1 not found");

    // Chains: both should have chain A and B
    assert_eq!(cif_m1.get_group_list(), pdb_m1.get_group_list());

    let cif_chain_a = cif_m1.get_group("A").unwrap();
    let pdb_chain_a = pdb_m1.get_group("A").unwrap();
    let cif_chain_b = cif_m1.get_group("B").unwrap();
    let pdb_chain_b = pdb_m1.get_group("B").unwrap();

    // Residue counts: A has 21 residues, B has 30 residues
    assert_eq!(cif_chain_a.get_number_of_groups(), 21);
    assert_eq!(pdb_chain_a.get_number_of_groups(), 21);
    assert_eq!(cif_chain_b.get_number_of_groups(), 30);
    assert_eq!(pdb_chain_b.get_number_of_groups(), 30);

    // Atom counts per chain: A has 312 atoms, B has 470 atoms
    assert_eq!(cif_chain_a.get_atom_list().len(), 312);
    assert_eq!(pdb_chain_a.get_atom_list().len(), 312);
    assert_eq!(cif_chain_b.get_atom_list().len(), 470);
    assert_eq!(pdb_chain_b.get_atom_list().len(), 470);

    // Total atoms: 782
    let cif_atoms = cif_ag.get_atom_list();
    let pdb_atoms = pdb_ag.get_atom_list();
    assert_eq!(cif_atoms.len(), 782);
    assert_eq!(pdb_atoms.len(), 782);

    // Verify first atom properties
    let cif_first = &cif_atoms[0];
    let pdb_first = &pdb_atoms[0];
    assert_eq!(cif_first.name, "N");
    assert_eq!(cif_first.symbol().unwrap(), "N");
    assert_eq!(cif_first.charge, 0.0);
    assert_eq!(pdb_first.name, "N");
    assert_eq!(pdb_first.symbol().unwrap(), "N");
    assert_eq!(pdb_first.charge, 0.0);

    // Coordinate agreement across all atoms (tolerance 1e-3)
    for (i, (c_atom, p_atom)) in cif_atoms.iter().zip(pdb_atoms.iter()).enumerate() {
        assert_eq!(
            c_atom.name, p_atom.name,
            "atom {i} name mismatch: cif={} pdb={}",
            c_atom.name, p_atom.name
        );
        assert_eq!(
            c_atom.atomic_number(),
            p_atom.atomic_number(),
            "atom {i} atomic number mismatch"
        );
        assert!(
            (c_atom.xyz.x - p_atom.xyz.x).abs() < 1e-3,
            "atom {i} x mismatch: cif={} pdb={}",
            c_atom.xyz.x,
            p_atom.xyz.x
        );
        assert!(
            (c_atom.xyz.y - p_atom.xyz.y).abs() < 1e-3,
            "atom {i} y mismatch: cif={} pdb={}",
            c_atom.xyz.y,
            p_atom.xyz.y
        );
        assert!(
            (c_atom.xyz.z - p_atom.xyz.z).abs() < 1e-3,
            "atom {i} z mismatch: cif={} pdb={}",
            c_atom.xyz.z,
            p_atom.xyz.z
        );
    }

    // Also test that cif.get_atomgroup("data_1HLS") dispatches to structure parser (loads all 20 models)
    let dispatched = cif.get_atomgroup("data_1HLS").expect("dispatch failed");
    assert_eq!(dispatched.get_number_of_groups(), 20);
    assert_eq!(dispatched.get_atom_list().len(), 15640);
    let disp_m1 = dispatched.get_group("model_1").expect("model_1 not found");
    assert_eq!(disp_m1.get_atom_list().len(), 782);
}

#[test]
fn test_2mgo_cif_matches_pdb() {
    let data_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let cif_path = data_dir.join("2MGO.cif");
    let pdb_path = data_dir.join("2MGO.pdb");

    let cif = SimpleMmcif::from_file(&cif_path).expect("failed to load 2MGO.cif");
    let pdb = Pdb::from_file(&pdb_path, None).expect("failed to load 2MGO.pdb");

    // Compare all 20 models
    let mut cif_ag = cif
        .get_structure_atomgroup(None, None)
        .expect("failed to get cif atomgroup");
    let mut pdb_ag = pdb
        .get_atomgroup(None, None)
        .expect("failed to get pdb atomgroup");

    assert_eq!(cif_ag.get_number_of_groups(), 20);
    assert_eq!(pdb_ag.get_number_of_groups(), 20);

    // Check model 1 structure
    let cif_m1 = cif_ag.get_group("model_1").expect("cif model_1 not found");
    let pdb_m1 = pdb_ag.get_group("model_1").expect("pdb model_1 not found");

    assert_eq!(cif_m1.get_number_of_groups(), 1); // chain A
    assert_eq!(pdb_m1.get_number_of_groups(), 1);

    let cif_chain_a = cif_m1.get_group("A").unwrap();
    let pdb_chain_a = pdb_m1.get_group("A").unwrap();
    assert_eq!(cif_chain_a.get_number_of_groups(), 9); // 9 residues
    assert_eq!(pdb_chain_a.get_number_of_groups(), 9);

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
        let c_res = cif_chain_a
            .get_group(seq)
            .unwrap_or_else(|| panic!("cif residue {seq} not found"));
        let p_res = pdb_chain_a
            .get_group(seq)
            .unwrap_or_else(|| panic!("pdb residue {seq} not found"));

        assert_eq!(c_res.name, expected_name);
        assert_eq!(p_res.name, expected_name);
        assert_eq!(c_res.get_number_of_atoms(), expected_atom_count);
        assert_eq!(p_res.get_number_of_atoms(), expected_atom_count);
    }

    // Check total atoms across all 20 models (20 * 134 = 2680)
    assert_eq!(cif_ag.get_atom_list().len(), 2680);
    assert_eq!(pdb_ag.get_atom_list().len(), 2680);

    // Verify SSBOND disulfide bond linking from _struct_conn matches pdb ssbonds
    let cif_bonds = cif_ag.get_bond_list();
    let pdb_bonds = pdb_ag.get_bond_list();
    assert_eq!(cif_bonds.len(), 20);
    assert_eq!(pdb_bonds.len(), 20);

    // In model 1, atom serials match PDB exactly (6 and 89)
    let m1_found = cif_bonds.iter().any(|b| {
        (b.atom1_path == "/model_1/A/1/6_SG" && b.atom2_path == "/model_1/A/6/89_SG")
            || (b.atom1_path == "/model_1/A/6/89_SG" && b.atom2_path == "/model_1/A/1/6_SG")
    });
    assert!(m1_found, "model_1 SSBOND exact path mismatch");

    // Across all 20 models, each model has a disulfide bond between CYS1 SG and CYS6 SG
    // (Note: In mmCIF, atom id is globally sequential across models, avoiding PDB's 99,999 limit)
    for i in 1..=20 {
        let prefix1 = format!("/model_{i}/A/1/");
        let prefix2 = format!("/model_{i}/A/6/");
        let found = cif_bonds.iter().any(|b| {
            let matches_forward = b.atom1_path.starts_with(&prefix1)
                && b.atom1_path.ends_with("_SG")
                && b.atom2_path.starts_with(&prefix2)
                && b.atom2_path.ends_with("_SG");
            let matches_backward = b.atom2_path.starts_with(&prefix1)
                && b.atom2_path.ends_with("_SG")
                && b.atom1_path.starts_with(&prefix2)
                && b.atom1_path.ends_with("_SG");
            matches_forward || matches_backward
        });
        assert!(
            found,
            "cif SSBOND for model_{i} between CYS1 SG and CYS6 SG not found"
        );
    }
}

#[test]
fn test_3i3z_cif_hierarchy_and_altloc() {
    let data_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let cif_path = data_dir.join("3I3Z.cif");
    let cif = SimpleMmcif::from_file(&cif_path).expect("failed to load 3I3Z.cif");

    let block = cif.get_data_block("data_3I3Z").unwrap();
    let all_records = block.get_atom_site_records().expect("records error");
    assert_eq!(
        all_records.len(),
        480,
        "total raw atom_site records in 3I3Z.cif"
    );

    // Count altlocs in raw records
    let count_alt_a = all_records.iter().filter(|r| r.label_alt_id == "A").count();
    let count_alt_b = all_records.iter().filter(|r| r.label_alt_id == "B").count();
    assert_eq!(count_alt_a, 20, "expected 20 altloc A records");
    assert_eq!(count_alt_b, 20, "expected 20 altloc B records");

    // Parse with default altloc ("A" / blank retained, "B" filtered out)
    let ag = cif
        .get_structure_atomgroup(None, None)
        .expect("failed to parse 3I3Z.cif with default altloc");

    assert_eq!(ag.get_number_of_groups(), 1); // 1 model: model_1
    let m1 = ag.get_group("model_1").expect("model_1 missing");

    // 2 chains: A and B (auth_asym_id merges water molecules into A and B; C/D do not exist)
    assert_eq!(m1.get_number_of_groups(), 2);
    let chain_a = m1.get_group("A").expect("chain A missing");
    let chain_b = m1.get_group("B").expect("chain B missing");

    // Residue counts with auth_seq_id:
    // Chain A: 21 polymer residues + 23 water residues = 44 residues
    // Chain B: 30 polymer residues + 34 water residues = 64 residues
    assert_eq!(chain_a.get_number_of_groups(), 44, "chain A residue count");
    assert_eq!(chain_b.get_number_of_groups(), 64, "chain B residue count");

    // Atom counts with default altloc (20 altloc B atoms excluded: 18 from B polymer, 2 from B water)
    // Chain A: 186 atoms (163 polymer + 23 water, no altloc B)
    // Chain B: 274 atoms (258 - 18 = 240 polymer + 36 - 2 = 34 water)
    assert_eq!(chain_a.get_atom_list().len(), 186, "chain A atom count");
    assert_eq!(
        chain_b.get_atom_list().len(),
        274,
        "chain B atom count with altloc A"
    );

    // Total atoms: 480 - 20 = 460
    assert_eq!(ag.get_atom_list().len(), 460);

    // Verify insertion code is captured in AtomSiteRecord
    // (All insertion codes in 3I3Z are absent/empty as documented in known gaps)
    for rec in &all_records {
        assert!(
            rec.pdbx_pdb_ins_code.is_empty(),
            "expected empty insertion code in 3I3Z.cif"
        );
    }
}

#[test]
fn test_mmcif_invalid_coordinate_error() {
    let invalid_cif = r#"data_invalid
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 N N . GLY A not_a_number 8.071 6.020
"#;

    let cif = SimpleMmcif::from_str(invalid_cif).expect("failed to tokenize/parse mmcif structure");
    let result = cif.get_structure_atomgroup(None, None);
    assert!(
        result.is_err(),
        "expected error on invalid coordinate string 'not_a_number'"
    );

    let err_str = result.err().unwrap().to_string();
    assert!(
        err_str.contains("Cartn_x") && err_str.contains("invalid float"),
        "error message should mention Cartn_x and invalid float: {err_str}"
    );
}

#[test]
fn test_ccd_bond_order_arom_and_quad() {
    let mmcif_content = r#"data_TEST_BONDS
_chem_comp.id TEST_BONDS
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_model_Cartn_x_ideal
_chem_comp_atom.pdbx_model_Cartn_y_ideal
_chem_comp_atom.pdbx_model_Cartn_z_ideal
TEST_BONDS C1 C 0.0 0.0 0.0
TEST_BONDS C2 C 1.4 0.0 0.0
TEST_BONDS Re1 RE 0.0 2.0 0.0
TEST_BONDS Re2 RE 0.0 4.2 0.0
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
TEST_BONDS C1 C2 AROM
TEST_BONDS Re1 Re2 QUAD
"#;

    let mut cif = SimpleMmcif::new();
    cif.load_from_str(mmcif_content).unwrap();
    let ag = cif.get_atomgroup("data_TEST_BONDS").unwrap();

    let bonds = ag.bonds();
    assert_eq!(bonds.len(), 2);

    let mut arom_found = false;
    let mut quad_found = false;

    for bond in bonds {
        let p1 = bond.atom1_path.trim_start_matches('/');
        let p2 = bond.atom2_path.trim_start_matches('/');
        if (p1 == "C1" && p2 == "C2") || (p1 == "C2" && p2 == "C1") {
            assert_eq!(bond.order, 1, "AROM bond order must be 1");
            arom_found = true;
        } else if (p1 == "Re1" && p2 == "Re2") || (p1 == "Re2" && p2 == "Re1") {
            assert_eq!(bond.order, 4, "QUAD bond order must be 4");
            quad_found = true;
        }
    }

    assert!(arom_found, "AROM bond between C1 and C2 should be found");
    assert!(quad_found, "QUAD bond between Re1 and Re2 should be found");
}

#[test]
fn test_ccd_missing_coordinate_error() {
    let mmcif_content = r#"data_MISSING_COORD
_chem_comp.id MISSING_COORD
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_model_Cartn_x_ideal
_chem_comp_atom.pdbx_model_Cartn_y_ideal
_chem_comp_atom.pdbx_model_Cartn_z_ideal
MISSING_COORD C1 C 0.0 0.0 0.0
MISSING_COORD C2 C ? ? ?
"#;

    let mut cif = SimpleMmcif::new();
    cif.load_from_str(mmcif_content).unwrap();
    let result = cif.get_atomgroup("data_MISSING_COORD");

    assert!(
        result.is_err(),
        "expected error when coordinates are missing for atom C2"
    );

    let err_str = result.err().unwrap().to_string();
    assert!(
        err_str.contains("C2") && err_str.contains("Missing or unparseable coordinates"),
        "error message should mention atom 'C2' and missing coordinates: {err_str}"
    );
}

#[test]
fn test_ccd_multiple_data_blocks() {
    let mmcif_content = r#"data_ALA
_chem_comp.id ALA
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_model_Cartn_x_ideal
_chem_comp_atom.pdbx_model_Cartn_y_ideal
_chem_comp_atom.pdbx_model_Cartn_z_ideal
ALA N N -0.966 0.493 1.500
ALA CA C 0.257 0.418 0.692
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
ALA N CA SING

data_BNZ
_chem_comp.id BNZ
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_model_Cartn_x_ideal
_chem_comp_atom.pdbx_model_Cartn_y_ideal
_chem_comp_atom.pdbx_model_Cartn_z_ideal
BNZ C1 C 0.000 1.396 0.000
BNZ C2 C 1.209 0.698 0.000
BNZ C3 C 1.209 -0.698 0.000
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
BNZ C1 C2 AROM
BNZ C2 C3 AROM
"#;

    let mut cif = SimpleMmcif::new();
    cif.load_from_str(mmcif_content).unwrap();

    let mut names = cif.get_molecule_names();
    names.sort();
    assert_eq!(names, vec!["data_ALA".to_string(), "data_BNZ".to_string()]);

    // Check ALA block
    let ag_ala = cif.get_atomgroup("data_ALA").unwrap();
    assert_eq!(ag_ala.name, "ALA");
    assert_eq!(ag_ala.get_number_of_atoms(), 2);
    assert_eq!(ag_ala.bonds().len(), 1);
    assert_eq!(ag_ala.bonds()[0].order, 1);

    // Check BNZ block
    let ag_bnz = cif.get_atomgroup("data_BNZ").unwrap();
    assert_eq!(ag_bnz.name, "BNZ");
    assert_eq!(ag_bnz.get_number_of_atoms(), 3);
    assert_eq!(ag_bnz.bonds().len(), 2);
    assert_eq!(ag_bnz.bonds()[0].order, 1);
    assert_eq!(ag_bnz.bonds()[1].order, 1);
}

#[test]
fn test_mmcif_insertion_code_residues() {
    let mmcif_content = r#"data_INS_CODE_TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_asym_id
_atom_site.auth_comp_id
_atom_site.auth_seq_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . ALA A 52 ? 0.000 0.000 0.000 A ALA 52 N 1
ATOM 2 C CA . ALA A 52 ? 1.000 0.000 0.000 A ALA 52 CA 1
ATOM 3 N N . GLY A 52 A 2.000 0.000 0.000 A GLY 52 N 1
ATOM 4 C CA . GLY A 52 A 3.000 0.000 0.000 A GLY 52 CA 1
ATOM 5 N N . SER A 52 B 4.000 0.000 0.000 A SER 52 N 1
ATOM 6 C CA . SER A 52 B 5.000 0.000 0.000 A SER 52 CA 1
"#;
    let mut cif = SimpleMmcif::new();
    cif.load_from_str(mmcif_content).unwrap();
    let ag = cif.get_atomgroup("INS_CODE_TEST").unwrap();

    let chain_a = ag
        .get_group("model_1")
        .expect("model_1 must exist")
        .get_group("A")
        .expect("chain A must exist");

    // All three residues must exist as separate groups
    assert!(chain_a.has_group("52"), "residue 52 must exist");
    assert!(chain_a.has_group("52A"), "residue 52A must exist");
    assert!(chain_a.has_group("52B"), "residue 52B must exist");

    let res_52 = chain_a.get_group("52").unwrap();
    let res_52a = chain_a.get_group("52A").unwrap();
    let res_52b = chain_a.get_group("52B").unwrap();

    assert_eq!(res_52.name, "ALA");
    assert_eq!(res_52a.name, "GLY");
    assert_eq!(res_52b.name, "SER");

    assert_eq!(res_52.get_number_of_atoms(), 2);
    assert_eq!(res_52a.get_number_of_atoms(), 2);
    assert_eq!(res_52b.get_number_of_atoms(), 2);
}

#[test]
fn test_mmcif_insertion_code_struct_conn() {
    let mmcif_content = r#"data_CONN_TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_asym_id
_atom_site.auth_comp_id
_atom_site.auth_seq_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . CYS A 52 A 0.000 0.000 0.000 A CYS 52 N 1
ATOM 2 S SG . CYS A 52 A 1.000 0.000 0.000 A CYS 52 SG 1
ATOM 3 N N . CYS A 100 ? 5.000 0.000 0.000 A CYS 100 N 1
ATOM 4 S SG . CYS A 100 ? 3.040 0.000 0.000 A CYS 100 SG 1

loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.pdbx_ptnr1_PDB_ins_code
_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_comp_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.pdbx_ptnr2_PDB_ins_code
_struct_conn.ptnr2_auth_asym_id
_struct_conn.ptnr2_auth_comp_id
_struct_conn.ptnr2_auth_seq_id
disulf1 disulf A CYS 52 SG A A CYS 52 A CYS 100 SG ? A CYS 100
"#;
    let mut cif = SimpleMmcif::new();
    cif.load_from_str(mmcif_content).unwrap();

    let mut ag = cif.get_atomgroup("CONN_TEST").unwrap();
    let bonds = ag.get_bond_list();
    assert_eq!(
        bonds.len(),
        1,
        "struct_conn disulfide bond must be established for insertion code residue"
    );
    let (a1, a2) = ag.resolve_bond(&bonds[0]).expect("bond must resolve");
    assert!(
        (a1.path.contains("/52A/") && a2.path.contains("/100/"))
            || (a1.path.contains("/100/") && a2.path.contains("/52A/")),
        "Bond must connect CYS 52A SG and CYS 100 SG, got {} and {}",
        a1.path,
        a2.path
    );
}

#[test]
fn test_2fb4_real_mmcif_insertion_codes_and_struct_conn() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/2FB4.cif");
    let mut cif = SimpleMmcif::new();
    cif.load(&path).expect("failed to load 2FB4.cif");

    let mut ag = cif
        .get_atomgroup("2FB4")
        .expect("failed to get 2FB4 atomgroup");

    let chain_h = ag
        .get_group("model_1")
        .expect("model_1 must exist")
        .get_group("H")
        .expect("chain H must exist");

    // 1. Verify all 5 residues with seq 101 exist as separate groups
    for &(res_key, expected_name, expected_atom_count) in &[
        ("101", "GLY", 4),
        ("101A", "HIS", 10),
        ("101B", "GLY", 4),
        ("101C", "PHE", 11),
        ("101D", "CYS", 6),
    ] {
        assert!(
            chain_h.has_group(res_key),
            "chain H must have residue group {res_key}"
        );
        let res = chain_h.get_group(res_key).unwrap();
        assert_eq!(
            res.name, expected_name,
            "residue {res_key} must be {expected_name}"
        );
        assert_eq!(
            res.get_number_of_atoms(),
            expected_atom_count,
            "residue {res_key} must have {expected_atom_count} atoms"
        );
    }

    // 2. Verify all 5 residues with seq 104 exist as separate groups
    for &(res_key, expected_name, expected_atom_count) in &[
        ("104", "ALA", 5),
        ("104A", "SER", 6),
        ("104B", "CYS", 6),
        ("104C", "PHE", 11),
        ("104D", "GLY", 4),
    ] {
        assert!(
            chain_h.has_group(res_key),
            "chain H must have residue group {res_key}"
        );
        let res = chain_h.get_group(res_key).unwrap();
        assert_eq!(
            res.name, expected_name,
            "residue {res_key} must be {expected_name}"
        );
        assert_eq!(
            res.get_number_of_atoms(),
            expected_atom_count,
            "residue {res_key} must have {expected_atom_count} atoms"
        );
    }

    // 3. Verify disulfide bond between insertion-code residues CYS 101D and CYS 104B is resolved from _struct_conn
    let bonds = ag.get_bond_list();
    let ssbond_101d_104b = bonds.iter().find(|b| {
        (b.atom1_path.contains("/H/101D/") && b.atom2_path.contains("/H/104B/"))
            || (b.atom1_path.contains("/H/104B/") && b.atom2_path.contains("/H/101D/"))
    });
    assert!(
        ssbond_101d_104b.is_some(),
        "disulfide bond between CYS 101D and CYS 104B must be established in 2FB4 mmCIF"
    );
    assert_eq!(ssbond_101d_104b.unwrap().order, 1);
}
