// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::path::PathBuf;

use proteindf_bridge::atom::Atom;
use proteindf_bridge::atom_group::AtomGroup;
use proteindf_bridge::bond::Bond;
use proteindf_bridge::ccd_templates::CcdTemplateDb;
use proteindf_bridge::format::Pdb;

#[test]
fn test_ccd_template_db_embedded() {
    let db = CcdTemplateDb::global();
    assert_eq!(db.len(), 29);
    assert!(!db.is_empty());

    // Check ALA
    let ala = db.lookup("ALA").expect("ALA template should exist");
    assert_eq!(ala.comp_id, "ALA");
    assert!(ala.atoms.contains(&"N".to_string()));
    assert!(ala.atoms.contains(&"CA".to_string()));
    assert!(ala.atoms.contains(&"C".to_string()));
    assert!(ala.atoms.contains(&"O".to_string()));
    assert!(ala.atoms.contains(&"CB".to_string()));
    // C=O should have order 2
    let co_bond = ala
        .bonds
        .iter()
        .find(|(a1, a2, _)| (a1 == "C" && a2 == "O") || (a1 == "O" && a2 == "C"));
    assert_eq!(co_bond, Some(&("C".to_string(), "O".to_string(), 2)));

    // Check ARG: Guanidino group CZ=NH2 double bond
    let arg = db.lookup("ARG").expect("ARG template should exist");
    let cz_nh2 = arg
        .bonds
        .iter()
        .find(|(a1, a2, _)| (a1 == "CZ" && a2 == "NH2") || (a1 == "NH2" && a2 == "CZ"));
    assert_eq!(cz_nh2, Some(&("CZ".to_string(), "NH2".to_string(), 2)));

    // Check HOH: Water
    let hoh = db.lookup("HOH").expect("HOH template should exist");
    assert_eq!(hoh.atoms.len(), 3);
    assert_eq!(hoh.bonds.len(), 2);
    for (_, _, order) in &hoh.bonds {
        assert_eq!(*order, 1);
    }
}

#[test]
fn test_apply_ccd_bond_templates_1hls_real_pdb() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/1hls.pdb");
    let mut pdb = Pdb::new(None);
    pdb.load(&path).expect("failed to load 1hls.pdb");

    let mut ag = pdb
        .get_atomgroup(None, None)
        .expect("failed to get 1hls atomgroup");

    // Verify baseline: Bond::setup() alone without templates assigns bond order 1 to all detected bonds
    let mut ag_baseline = ag.clone();
    let mut bond = Bond::new();
    bond.setup(&mut ag_baseline).expect("Bond::setup failed");
    let baseline_bonds = ag_baseline.get_bond_list();
    assert!(!baseline_bonds.is_empty());
    // In pure VDW heuristic, every bond has order 1
    assert!(
        baseline_bonds.iter().all(|b| b.order == 1),
        "Bond::setup alone should only produce order 1 bonds"
    );

    // Now apply CCD bond templates on ag
    let db = CcdTemplateDb::global();
    ag.apply_ccd_bond_templates(db);

    let bonds_after_ccd = ag.get_bond_list();
    assert!(!bonds_after_ccd.is_empty());

    // Verify double bonds in standard residues:
    // 1. Chain A residue 4 is GLU (1hls.pdb chain A has GLU at res 4)
    //    Should have backbone C=O (order 2) and sidechain CD=OE1 (order 2)
    let mut glu_double_bonds = 0;
    for b in &bonds_after_ccd {
        if b.atom1_path.contains("/A/4/") && b.atom2_path.contains("/A/4/") && b.order == 2 {
            glu_double_bonds += 1;
        }
    }
    assert!(
        glu_double_bonds >= 2,
        "GLU residue 4 should have at least 2 double bonds (C=O and CD=OE1/OE2), found {}",
        glu_double_bonds
    );

    // 2. Chain B residue 22 is ARG
    //    Should have backbone C=O (order 2) and guanidino CZ=NH2 (order 2)
    let mut arg_double_bonds = 0;
    for b in &bonds_after_ccd {
        if b.atom1_path.contains("/B/22/") && b.atom2_path.contains("/B/22/") && b.order == 2 {
            arg_double_bonds += 1;
        }
    }
    assert!(
        arg_double_bonds >= 2,
        "ARG residue 22 should have at least 2 double bonds (C=O and CZ=NH2), found {}",
        arg_double_bonds
    );

    // 3. Fallback synergy test: Run Bond::setup afterwards for inter-residue peptide bonds
    //    and verify that existing CCD template bonds (order 2) are not overwritten to order 1
    let mut bond2 = Bond::new();
    bond2
        .setup(&mut ag)
        .expect("Bond::setup after CCD templates failed");
    let final_bonds = ag.get_bond_list();
    let co_bond = final_bonds.iter().find(|b| {
        let is_a1_c = b.atom1_path.contains("/A/4/")
            && (b.atom1_path.ends_with("_C") || b.atom1_path.ends_with("/C"));
        let is_a2_o = b.atom2_path.contains("/A/4/")
            && (b.atom2_path.ends_with("_O") || b.atom2_path.ends_with("/O"));
        let is_a1_o = b.atom1_path.contains("/A/4/")
            && (b.atom1_path.ends_with("_O") || b.atom1_path.ends_with("/O"));
        let is_a2_c = b.atom2_path.contains("/A/4/")
            && (b.atom2_path.ends_with("_C") || b.atom2_path.ends_with("/C"));
        (is_a1_c && is_a2_o) || (is_a1_o && is_a2_c)
    });
    assert!(co_bond.is_some(), "C=O bond should be present");
    assert_eq!(
        co_bond.unwrap().order,
        2,
        "C=O bond order should remain 2 after Bond::setup"
    );
}

#[test]
fn test_apply_ccd_bond_templates_does_not_overwrite_existing_bonds() {
    let mut ag = AtomGroup::new();
    ag.set_path("/model_1/A/1/".to_string());
    ag.name = "ALA".to_string();

    let mut atom_c = Atom::new();
    atom_c.name = "C".to_string();
    let mut atom_o = Atom::new();
    atom_o.name = "O".to_string();

    ag.set_atom("C", atom_c);
    ag.set_atom("O", atom_o);

    // Pre-register C-O bond as order 1 (e.g. from an explicit file source like CONECT)
    let c_ref = ag.get_atom("C").unwrap().clone();
    let o_ref = ag.get_atom("O").unwrap().clone();
    ag.add_bond(&c_ref, &o_ref, 1);

    assert_eq!(ag.bonds().len(), 1);
    assert_eq!(ag.bonds()[0].order, 1);

    // Apply CCD templates (where ALA defines C=O with order 2)
    let db = CcdTemplateDb::global();
    ag.apply_ccd_bond_templates(db);

    // The existing bond must NOT be overwritten; it should stay order 1, and no duplicate added
    assert_eq!(ag.bonds().len(), 1, "no duplicate bond should be added");
    assert_eq!(
        ag.bonds()[0].order,
        1,
        "existing file-derived bond order 1 must not be overwritten"
    );
}

#[test]
fn test_apply_ccd_bond_templates_unknown_component_safe_skip() {
    let mut ag = AtomGroup::new();
    ag.set_path("/model_1/A/1/".to_string());
    ag.name = "XYZ_UNKNOWN".to_string();

    let mut atom_x1 = Atom::new();
    atom_x1.name = "X1".to_string();
    let mut atom_x2 = Atom::new();
    atom_x2.name = "X2".to_string();

    ag.set_atom("X1", atom_x1);
    ag.set_atom("X2", atom_x2);

    let db = CcdTemplateDb::global();
    // Applying templates to an unknown component should not error or panic
    ag.apply_ccd_bond_templates(db);

    // No bonds added since component is unknown
    assert_eq!(ag.bonds().len(), 0);
    assert_eq!(ag.get_number_of_atoms(), 2);
}
