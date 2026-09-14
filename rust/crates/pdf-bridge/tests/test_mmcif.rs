// Copyright (C) 2019 The ProteinDF development team.
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

use std::path::PathBuf;

use pdf_bridge::format::pdb::Pdb;
use pdf_bridge::format::SimpleMmcif;

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
    let all_records = block.get_atom_site_records();
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

    // 4 chains: A, B, C, D
    assert_eq!(m1.get_number_of_groups(), 4);
    let chain_a = m1.get_group("A").expect("chain A missing");
    let chain_b = m1.get_group("B").expect("chain B missing");
    let chain_c = m1.get_group("C").expect("chain C missing");
    let chain_d = m1.get_group("D").expect("chain D missing");

    // Residue counts
    assert_eq!(chain_a.get_number_of_groups(), 21, "chain A residue count");
    assert_eq!(chain_b.get_number_of_groups(), 30, "chain B residue count");
    assert_eq!(
        chain_c.get_number_of_groups(),
        1,
        "chain C (water HETATM) residue count"
    );
    assert_eq!(
        chain_d.get_number_of_groups(),
        1,
        "chain D (water HETATM) residue count"
    );

    // Atom counts with default altloc (20 altloc B atoms excluded: 18 from B, 2 from D)
    // Chain A: 163 atoms
    // Chain B: 258 - 18 = 240 atoms
    // Chain C: 23 atoms
    // Chain D: 36 - 2 = 34 atoms
    assert_eq!(chain_a.get_atom_list().len(), 163, "chain A atom count");
    assert_eq!(
        chain_b.get_atom_list().len(),
        240,
        "chain B atom count with altloc A"
    );
    assert_eq!(chain_c.get_atom_list().len(), 23, "chain C atom count");
    assert_eq!(
        chain_d.get_atom_list().len(),
        34,
        "chain D atom count with altloc A"
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
