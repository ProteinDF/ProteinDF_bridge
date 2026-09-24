// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::path::PathBuf;

use proteindf_bridge::atom::Atom;
use proteindf_bridge::atom_group::AtomGroup;
use proteindf_bridge::bond::Bond;
use proteindf_bridge::ccd_templates::{CcdAtom, CcdBondTemplate, CcdTemplateDb};
use proteindf_bridge::format::amber_prmtop::AmberPrmtop;
use proteindf_bridge::format::gro::SimpleGro;
use proteindf_bridge::format::mmcif::SimpleMmcif;
use proteindf_bridge::format::mol2::SimpleMol2;
use proteindf_bridge::format::Pdb;
use proteindf_bridge::position::Position;

/// Builds a minimal synthetic `CcdAtom` for tests that only care about atom names/bonds
/// (e.g. merge-precedence tests), not geometry.
fn synthetic_atom(name: &str) -> CcdAtom {
    CcdAtom {
        name: name.to_string(),
        element: "C".to_string(),
        ideal_xyz: None,
    }
}

#[test]
fn test_ccd_template_db_embedded() {
    let db = CcdTemplateDb::global();
    assert_eq!(db.len(), 29);
    assert!(!db.is_empty());

    // Check ALA
    let ala = db.lookup("ALA").expect("ALA template should exist");
    assert_eq!(ala.comp_id, "ALA");
    assert!(ala.get_atom("N").is_some());
    assert!(ala.get_atom("CA").is_some());
    assert!(ala.get_atom("C").is_some());
    assert!(ala.get_atom("O").is_some());
    assert!(ala.get_atom("CB").is_some());
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

/// Reference idealized coordinates for ALA, independently read off the raw CCD CIF
/// (`https://files.rcsb.org/ligands/view/ALA.cif`, `_chem_comp_atom` loop,
/// `pdbx_model_Cartn_{x,y,z}_ideal` columns) rather than derived from the code under test.
#[test]
fn test_ccd_template_ala_hydrogen_ideal_coordinates() {
    let db = CcdTemplateDb::global();
    let ala = db.lookup("ALA").expect("ALA template should exist");

    let cases: &[(&str, &str, (f64, f64, f64))] = &[
        ("N", "N", (-0.966, 0.493, 1.500)),
        ("CA", "C", (0.257, 0.418, 0.692)),
        ("C", "C", (-0.094, 0.017, -0.716)),
        ("O", "O", (-1.056, -0.682, -0.923)),
        ("CB", "C", (1.204, -0.620, 1.296)),
        ("H", "H", (-1.383, -0.425, 1.482)),
        ("H2", "H", (-0.676, 0.661, 2.452)),
        ("HA", "H", (0.746, 1.392, 0.682)),
        ("HB1", "H", (1.459, -0.330, 2.316)),
        ("HB2", "H", (0.715, -1.594, 1.307)),
        ("HB3", "H", (2.113, -0.676, 0.697)),
        ("OXT", "O", (0.661, 0.439, -1.742)),
        ("HXT", "H", (0.435, 0.182, -2.647)),
    ];

    for (name, element, (x, y, z)) in cases {
        let atom = ala
            .get_atom(name)
            .unwrap_or_else(|| panic!("ALA atom '{name}' should exist"));
        assert_eq!(&atom.element, element, "element mismatch for atom {name}");
        let (ax, ay, az) = atom
            .ideal_xyz
            .unwrap_or_else(|| panic!("ALA atom '{name}' should have ideal_xyz"));
        assert!(
            (ax - x).abs() < 1e-6 && (ay - y).abs() < 1e-6 && (az - z).abs() < 1e-6,
            "ideal_xyz mismatch for atom {name}: got ({ax}, {ay}, {az}), expected ({x}, {y}, {z})"
        );
    }

    // Hydrogen classification
    assert!(ala.get_atom("HA").unwrap().is_hydrogen());
    assert!(!ala.get_atom("CA").unwrap().is_hydrogen());
}

/// Every embedded template (all 29 standard components) should expose element symbols
/// and idealized coordinates for its hydrogen atoms, not just bare names.
#[test]
fn test_ccd_template_all_embedded_hydrogens_have_ideal_coordinates() {
    let db = CcdTemplateDb::global();
    let comp_ids = [
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
        "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "DA", "DC", "DG", "DT", "A", "C", "G",
        "U", "HOH",
    ];
    assert_eq!(comp_ids.len(), 29);

    for comp_id in comp_ids {
        let template = db
            .lookup(comp_id)
            .unwrap_or_else(|| panic!("template '{comp_id}' should exist"));
        let hydrogens: Vec<_> = template.atoms.iter().filter(|a| a.is_hydrogen()).collect();
        assert!(
            !hydrogens.is_empty(),
            "template '{comp_id}' should have at least one hydrogen atom"
        );
        for h in &hydrogens {
            assert_eq!(h.element, "H", "hydrogen atom '{}' in {comp_id}", h.name);
            assert!(
                h.ideal_xyz.is_some(),
                "hydrogen atom '{}' in {comp_id} should have ideal_xyz",
                h.name
            );
        }
        // Heavy atoms should also be present with resolvable coordinates.
        let heavy_without_xyz: Vec<_> = template
            .atoms
            .iter()
            .filter(|a| !a.is_hydrogen() && a.ideal_xyz.is_none())
            .map(|a| a.name.clone())
            .collect();
        assert!(
            heavy_without_xyz.is_empty(),
            "template '{comp_id}' has heavy atoms without ideal_xyz: {heavy_without_xyz:?}"
        );
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
    assert!(
        ag.get_bond_list().is_empty(),
        "Raw PDB loader output should have no bonds when file lacks CONECT/SSBOND"
    );

    // Verify baseline: Bond::setup_heuristic() alone without templates assigns bond order 1 to all detected bonds
    let mut ag_baseline = ag.clone();
    let mut bond = Bond::new();
    bond.setup_heuristic(&mut ag_baseline)
        .expect("Bond::setup_heuristic failed");
    let baseline_bonds = ag_baseline.get_bond_list();
    assert!(!baseline_bonds.is_empty());
    // In pure VDW/covalent heuristic, every bond has order 1
    assert!(
        baseline_bonds.iter().all(|b| b.order == 1),
        "Bond::setup_heuristic alone should only produce order 1 bonds"
    );

    // Now apply CCD bond templates on a clean ag_ccd
    let mut ag_ccd = ag.clone();
    let db = CcdTemplateDb::global();
    ag_ccd.apply_ccd_bond_templates(db);

    let bonds_after_ccd = ag_ccd.get_bond_list();
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

    // 3. Fallback synergy test: Run Bond::setup_heuristic afterwards for inter-residue peptide bonds
    //    and verify that existing CCD template bonds (order 2) are not overwritten to order 1
    let mut bond2 = Bond::new();
    bond2
        .setup_heuristic(&mut ag_ccd)
        .expect("Bond::setup_heuristic after CCD templates failed");
    let final_bonds = ag_ccd.get_bond_list();
    let co_bonds: Vec<_> = final_bonds
        .iter()
        .filter(|b| {
            let is_a1_c = b.atom1_path.contains("/A/4/")
                && (b.atom1_path.ends_with("_C") || b.atom1_path.ends_with("/C"));
            let is_a2_o = b.atom2_path.contains("/A/4/")
                && (b.atom2_path.ends_with("_O") || b.atom2_path.ends_with("/O"));
            let is_a1_o = b.atom1_path.contains("/A/4/")
                && (b.atom1_path.ends_with("_O") || b.atom1_path.ends_with("/O"));
            let is_a2_c = b.atom2_path.contains("/A/4/")
                && (b.atom2_path.ends_with("_C") || b.atom2_path.ends_with("/C"));
            (is_a1_c && is_a2_o) || (is_a1_o && is_a2_c)
        })
        .collect();
    assert_eq!(
        co_bonds.len(),
        1,
        "GLU 4 C=O bond should appear exactly once (no duplicate bond record)"
    );
    assert_eq!(
        co_bonds[0].order, 2,
        "C=O bond order should remain 2 after Bond::setup_heuristic"
    );

    // Verify that NO duplicate bond records exist across all final bonds
    let mut seen_pairs = std::collections::HashSet::new();
    for b in &final_bonds {
        let key = if b.atom1_path <= b.atom2_path {
            (&b.atom1_path, &b.atom2_path)
        } else {
            (&b.atom2_path, &b.atom1_path)
        };
        assert!(
            seen_pairs.insert(key),
            "Duplicate bond found between {} and {}",
            b.atom1_path,
            b.atom2_path
        );
    }
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

#[test]
fn test_from_mmcif_block_ala_cif_matches_embedded() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/ALA.cif");
    let mut cif = SimpleMmcif::new();
    cif.load(&path).expect("failed to load ALA.cif");

    let block = cif.get_data_block("ALA").expect("data block ALA not found");
    let parsed_template =
        CcdBondTemplate::from_mmcif_block(block, "ALA").expect("from_mmcif_block failed for ALA");

    assert_eq!(parsed_template.comp_id, "ALA");

    let embedded_db = CcdTemplateDb::global();
    let embedded_ala = embedded_db
        .lookup("ALA")
        .expect("ALA not found in embedded DB");

    // Both should contain 13 atoms
    assert_eq!(parsed_template.atoms.len(), 13);
    assert_eq!(parsed_template.atoms, embedded_ala.atoms);

    // Both should contain 12 bonds with identical topology and bond orders
    assert_eq!(parsed_template.bonds.len(), 12);
    assert_eq!(parsed_template.bonds, embedded_ala.bonds);

    // Verify C=O is a double bond
    let co_bond = parsed_template
        .bonds
        .iter()
        .find(|(a1, a2, _)| (a1 == "C" && a2 == "O") || (a1 == "O" && a2 == "C"));
    assert_eq!(co_bond, Some(&("C".to_string(), "O".to_string(), 2)));
}

#[test]
fn test_from_mmcif_block_synthetic_custom_ligand() {
    // Synthetic CCD CIF data for a non-standard custom ligand "LIG"
    let cif_text = r#"
data_LIG
#
_chem_comp.id                                    LIG
_chem_comp.name                                  "CUSTOM LIGAND"
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
LIG C1 C
LIG O1 O
LIG C2 C
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
LIG C1 O1 DOUB
LIG C1 C2 SING
"#;

    let mut cif = SimpleMmcif::new();
    cif.load_from_str(cif_text)
        .expect("failed to parse synthetic ligand CIF");

    let block = cif.get_data_block("LIG").expect("block LIG not found");
    let template =
        CcdBondTemplate::from_mmcif_block(block, "LIG").expect("from_mmcif_block failed");

    assert_eq!(template.comp_id, "LIG");
    assert_eq!(template.atoms.len(), 3);
    // No coordinate columns in this minimal fixture: ideal_xyz must fall back to None,
    // not error (a template usable for bond-order resolution even without geometry).
    for (name, element) in [("C1", "C"), ("O1", "O"), ("C2", "C")] {
        let atom = template.get_atom(name).unwrap();
        assert_eq!(atom.element, element);
        assert!(atom.ideal_xyz.is_none());
    }
    assert_eq!(
        template.bonds,
        vec![
            ("C1".to_string(), "O1".to_string(), 2),
            ("C1".to_string(), "C2".to_string(), 1),
        ]
    );

    // Insert into a freshly created CcdTemplateDb
    let mut custom_db = CcdTemplateDb::new();
    assert!(custom_db.is_empty());
    custom_db.insert(template);
    assert_eq!(custom_db.len(), 1);

    // Apply to AtomGroup
    let mut ag = AtomGroup::new();
    ag.set_path("/model_1/A/100/".to_string());
    ag.name = "LIG".to_string();

    let mut a_c1 = Atom::new();
    a_c1.name = "C1".to_string();
    let mut a_o1 = Atom::new();
    a_o1.name = "O1".to_string();
    let mut a_c2 = Atom::new();
    a_c2.name = "C2".to_string();

    ag.set_atom("C1", a_c1);
    ag.set_atom("O1", a_o1);
    ag.set_atom("C2", a_c2);

    ag.apply_ccd_bond_templates(&custom_db);

    let bonds = ag.get_bond_list();
    assert_eq!(bonds.len(), 2);

    let c1_o1 = bonds.iter().find(|b| {
        (b.atom1_path.ends_with("/C1") && b.atom2_path.ends_with("/O1"))
            || (b.atom1_path.ends_with("/O1") && b.atom2_path.ends_with("/C1"))
    });
    assert!(c1_o1.is_some());
    assert_eq!(
        c1_o1.unwrap().order,
        2,
        "C1=O1 bond order should be 2 from custom template"
    );

    let c1_c2 = bonds.iter().find(|b| {
        (b.atom1_path.ends_with("/C1") && b.atom2_path.ends_with("/C2"))
            || (b.atom1_path.ends_with("/C2") && b.atom2_path.ends_with("/C1"))
    });
    assert!(c1_c2.is_some());
    assert_eq!(
        c1_c2.unwrap().order,
        1,
        "C1-C2 bond order should be 1 from custom template"
    );
}

/// A user-supplied CCD file (via `from_mmcif_block`) with explicit hydrogens and both
/// idealized and model coordinate columns must: (1) normalize deuterium ("D") to "H",
/// (2) prefer idealized over model coordinates when both are present, and (3) fall back
/// to model coordinates when idealized ones are unresolvable ("?").
#[test]
fn test_from_mmcif_block_ideal_priority_and_deuterium_normalization() {
    let cif_text = r#"
data_LG2
#
_chem_comp.id                                    LG2
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.model_Cartn_x
_chem_comp_atom.model_Cartn_y
_chem_comp_atom.model_Cartn_z
_chem_comp_atom.pdbx_model_Cartn_x_ideal
_chem_comp_atom.pdbx_model_Cartn_y_ideal
_chem_comp_atom.pdbx_model_Cartn_z_ideal
LG2 C1 C 0.100 0.200 0.300 1.100 1.200 1.300
LG2 D1 D 0.400 0.500 0.600 1.400 1.500 1.600
LG2 C2 C 0.700 0.800 0.900 ?     ?     ?
LG2 C3 C 2.100 2.200 2.300 3.100 3.200 ?
"#;

    let mut cif = SimpleMmcif::new();
    cif.load_from_str(cif_text)
        .expect("failed to parse synthetic LG2 CIF");
    let block = cif.get_data_block("LG2").expect("block LG2 not found");
    let template =
        CcdBondTemplate::from_mmcif_block(block, "LG2").expect("from_mmcif_block failed");

    // Idealized coordinates take priority over model coordinates when both are present.
    let c1 = template.get_atom("C1").unwrap();
    assert_eq!(c1.element, "C");
    assert_eq!(c1.ideal_xyz, Some((1.100, 1.200, 1.300)));

    // Deuterium ("D") is normalized to hydrogen ("H"), and its idealized coordinates
    // are used just like any other atom.
    let d1 = template.get_atom("D1").unwrap();
    assert_eq!(d1.element, "H");
    assert!(d1.is_hydrogen());
    assert_eq!(d1.ideal_xyz, Some((1.400, 1.500, 1.600)));

    // Unresolvable ("?") idealized coordinates fall back to model coordinates.
    let c2 = template.get_atom("C2").unwrap();
    assert_eq!(c2.element, "C");
    assert_eq!(c2.ideal_xyz, Some((0.700, 0.800, 0.900)));

    // A partially-unresolvable idealized triple (only the z axis is "?") must fall back
    // to the *whole* model triple, not silently mix idealized x/y with model z into a
    // geometrically meaningless point that belongs to neither conformer.
    let c3 = template.get_atom("C3").unwrap();
    assert_eq!(c3.element, "C");
    assert_eq!(c3.ideal_xyz, Some((2.100, 2.200, 2.300)));
}

#[test]
fn test_from_mmcif_block_duplicate_atom_id_rows() {
    // An exact duplicate row (same element and coordinates) for an already-seen atom
    // name is tolerated silently...
    let exact_dup_cif = r#"
data_DUP1
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
DUP1 C1 C
DUP1 C1 C
"#;
    let mut cif = SimpleMmcif::new();
    cif.load_from_str(exact_dup_cif)
        .expect("failed to parse DUP1 CIF");
    let block = cif.get_data_block("DUP1").expect("block DUP1 not found");
    let template =
        CcdBondTemplate::from_mmcif_block(block, "DUP1").expect("exact duplicate should be OK");
    assert_eq!(template.atoms.len(), 1);

    // ...but a *conflicting* duplicate row (different element for the same atom name) is
    // rejected with an error instead of silently keeping whichever row appeared first.
    let conflicting_dup_cif = r#"
data_DUP2
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
DUP2 C1 C
DUP2 C1 N
"#;
    let mut cif2 = SimpleMmcif::new();
    cif2.load_from_str(conflicting_dup_cif)
        .expect("failed to parse DUP2 CIF");
    let block2 = cif2.get_data_block("DUP2").expect("block DUP2 not found");
    let result = CcdBondTemplate::from_mmcif_block(block2, "DUP2");
    assert!(
        result.is_err(),
        "conflicting duplicate atom_id rows should be rejected, not silently resolved"
    );
}

#[test]
fn test_from_mmcif_block_rejects_atom_site() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/1HLS.cif");
    let mut cif = SimpleMmcif::new();
    cif.load(&path).expect("failed to load 1HLS.cif");

    let block = cif.get_data_block("1HLS").expect("block 1HLS not found");
    assert!(block.has_atom_site());

    let result = CcdBondTemplate::from_mmcif_block(block, "1HLS");
    assert!(
        result.is_err(),
        "from_mmcif_block should reject data blocks with _atom_site"
    );
}

#[test]
fn test_ccd_template_db_merge_precedence() {
    let mut db1 = CcdTemplateDb::new();
    db1.insert(CcdBondTemplate {
        comp_id: "XYZ".to_string(),
        atoms: vec![synthetic_atom("A1"), synthetic_atom("A2")],
        bonds: vec![("A1".to_string(), "A2".to_string(), 1)],
    });

    let mut db2 = CcdTemplateDb::new();
    db2.insert(CcdBondTemplate {
        comp_id: "XYZ".to_string(),
        atoms: vec![synthetic_atom("A1"), synthetic_atom("A2")],
        bonds: vec![("A1".to_string(), "A2".to_string(), 2)], // Updated bond order 2
    });
    db2.insert(CcdBondTemplate {
        comp_id: "NEW".to_string(),
        atoms: vec![synthetic_atom("B1")],
        bonds: vec![],
    });

    // Merge db2 into db1 -> db2 entries should overwrite db1 ("last-write-wins")
    db1.merge(&db2);

    assert_eq!(db1.len(), 2);
    let xyz = db1.lookup("XYZ").expect("XYZ must exist in db1");
    assert_eq!(
        xyz.bonds[0].2, 2,
        "merged entry from db2 should overwrite existing entry with bond order 2"
    );
    assert!(db1.lookup("NEW").is_some());
}

#[test]
fn test_bond_setup_after_ccd_does_not_duplicate_bonds() {
    // Construct an isolated ALA residue with realistic coordinates
    let mut ag = AtomGroup::new();
    ag.set_path("/model_1/A/1/".to_string());
    ag.name = "ALA".to_string();

    let mut n = Atom::new();
    n.name = "N".to_string();
    n.set_atomic_number(7);
    n.xyz = Position::new(0.0, 0.0, 0.0);

    let mut ca = Atom::new();
    ca.name = "CA".to_string();
    ca.set_atomic_number(6);
    ca.xyz = Position::new(1.46, 0.0, 0.0);

    let mut c = Atom::new();
    c.name = "C".to_string();
    c.set_atomic_number(6);
    c.xyz = Position::new(2.0, 1.4, 0.0);

    let mut o = Atom::new();
    o.name = "O".to_string();
    o.set_atomic_number(8);
    o.xyz = Position::new(1.3, 2.4, 0.0);

    let mut cb = Atom::new();
    cb.name = "CB".to_string();
    cb.set_atomic_number(6);
    cb.xyz = Position::new(2.0, -0.7, 1.2);

    ag.set_atom("N", n);
    ag.set_atom("CA", ca);
    ag.set_atom("C", c);
    ag.set_atom("O", o);
    ag.set_atom("CB", cb);

    // 1. Apply CCD templates
    let db = CcdTemplateDb::global();
    ag.apply_ccd_bond_templates(db);

    let initial_bonds = ag.get_bond_list();
    assert_eq!(
        initial_bonds.len(),
        4,
        "ALA has 4 heavy atom bonds in CCD (N-CA, CA-C, C-O, CA-CB)"
    );
    let co_initial = initial_bonds
        .iter()
        .find(|b| {
            (b.atom1_path.ends_with("/C") && b.atom2_path.ends_with("/O"))
                || (b.atom1_path.ends_with("/O") && b.atom2_path.ends_with("/C"))
        })
        .expect("C=O bond should exist");
    assert_eq!(
        co_initial.order, 2,
        "C=O bond should have order 2 in CCD template"
    );

    // 2. Call Bond::setup_heuristic()
    let mut bond = Bond::new();
    bond.setup_heuristic(&mut ag)
        .expect("Bond::setup_heuristic failed");

    let after_bonds = ag.get_bond_list();
    assert_eq!(
        after_bonds.len(),
        initial_bonds.len(),
        "Bond::setup_heuristic must not duplicate already registered bonds"
    );

    let co_bonds: Vec<_> = after_bonds
        .iter()
        .filter(|b| {
            (b.atom1_path.ends_with("/C") && b.atom2_path.ends_with("/O"))
                || (b.atom1_path.ends_with("/O") && b.atom2_path.ends_with("/C"))
        })
        .collect();
    assert_eq!(
        co_bonds.len(),
        1,
        "C=O bond must remain unique (no duplicate bond record)"
    );
    assert_eq!(
        co_bonds[0].order, 2,
        "C=O bond order must remain 2 after Bond::setup_heuristic"
    );

    // Check all bonds are unique
    let mut seen = std::collections::HashSet::new();
    for b in &after_bonds {
        let key = if b.atom1_path <= b.atom2_path {
            (&b.atom1_path, &b.atom2_path)
        } else {
            (&b.atom2_path, &b.atom1_path)
        };
        assert!(seen.insert(key), "Duplicate bond found: {:?}", key);
    }
}

#[test]
fn test_atomgroup_setup_1hls_real_pdb() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/1hls.pdb");
    let mut pdb = Pdb::new(None);
    pdb.load(&path).expect("failed to load 1hls.pdb");
    let mut ag = pdb
        .get_atomgroup(None, None)
        .expect("failed to get 1hls atomgroup");

    // Verify raw output has no bonds, then call smart default AtomGroup::setup()
    assert!(
        ag.get_bond_list().is_empty(),
        "Raw PDB loader output should have no bonds before setup()"
    );
    ag.setup().expect("AtomGroup::setup failed");

    let final_bonds = ag.get_bond_list();

    // 1. Verify CCD template bond orders (e.g. GLU 4 C=O has order 2, unique)
    let glu4_co_bonds: Vec<_> = final_bonds
        .iter()
        .filter(|b| {
            let is_a1_c = b.atom1_path.contains("/A/4/")
                && (b.atom1_path.ends_with("_C") || b.atom1_path.ends_with("/C"));
            let is_a2_o = b.atom2_path.contains("/A/4/")
                && (b.atom2_path.ends_with("_O") || b.atom2_path.ends_with("/O"));
            let is_a1_o = b.atom1_path.contains("/A/4/")
                && (b.atom1_path.ends_with("_O") || b.atom1_path.ends_with("/O"));
            let is_a2_c = b.atom2_path.contains("/A/4/")
                && (b.atom2_path.ends_with("_C") || b.atom2_path.ends_with("/C"));
            (is_a1_c && is_a2_o) || (is_a1_o && is_a2_c)
        })
        .collect();
    assert_eq!(
        glu4_co_bonds.len(),
        1,
        "GLU 4 C=O bond should appear exactly once"
    );
    assert_eq!(
        glu4_co_bonds[0].order, 2,
        "GLU 4 C=O bond should have order 2 from CCD template"
    );

    // 2. Verify ARG 22 in Chain B has multiple double bonds (C=O and CZ=NH2)
    let arg_double_bonds = final_bonds
        .iter()
        .filter(|b| {
            (b.atom1_path.contains("/B/22/") || b.atom2_path.contains("/B/22/")) && b.order == 2
        })
        .count();
    assert!(
        arg_double_bonds >= 2,
        "ARG 22 should have at least 2 double bonds, found {}",
        arg_double_bonds
    );

    // 3. Verify inter-residue peptide bond is established via heuristic fallback
    let peptide_bond = final_bonds.iter().find(|b| {
        let is_res4_c = b.atom1_path.contains("/A/4/")
            && (b.atom1_path.ends_with("_C") || b.atom1_path.ends_with("/C"));
        let is_res5_n = b.atom2_path.contains("/A/5/")
            && (b.atom2_path.ends_with("_N") || b.atom2_path.ends_with("/N"));
        let is_res5_n_rev = b.atom1_path.contains("/A/5/")
            && (b.atom1_path.ends_with("_N") || b.atom1_path.ends_with("/N"));
        let is_res4_c_rev = b.atom2_path.contains("/A/4/")
            && (b.atom2_path.ends_with("_C") || b.atom2_path.ends_with("/C"));
        (is_res4_c && is_res5_n) || (is_res5_n_rev && is_res4_c_rev)
    });
    assert!(
        peptide_bond.is_some(),
        "Inter-residue peptide bond between GLU 4 and PRO 5 should be detected by heuristic"
    );
    assert_eq!(
        peptide_bond.unwrap().order,
        1,
        "Peptide bond should have order 1"
    );

    // 4. Verify that NO duplicate bond records exist across all final bonds
    let mut seen_pairs = std::collections::HashSet::new();
    for b in &final_bonds {
        let key = if b.atom1_path <= b.atom2_path {
            (&b.atom1_path, &b.atom2_path)
        } else {
            (&b.atom2_path, &b.atom1_path)
        };
        assert!(
            seen_pairs.insert(key),
            "Duplicate bond found between {} and {}",
            b.atom1_path,
            b.atom2_path
        );
    }
}

#[test]
fn test_resolve_bonds_preserves_existing_file_bonds() {
    let mut ag = AtomGroup::new();
    ag.set_path("/model_1/A/1/".to_string());
    ag.name = "ALA".to_string();

    let mut n = Atom::new();
    n.name = "N".to_string();
    n.set_atomic_number(7);
    n.xyz = Position::new(0.0, 0.0, 0.0);

    let mut ca = Atom::new();
    ca.name = "CA".to_string();
    ca.set_atomic_number(6);
    ca.xyz = Position::new(1.46, 0.0, 0.0);

    let mut c = Atom::new();
    c.name = "C".to_string();
    c.set_atomic_number(6);
    c.xyz = Position::new(2.0, 1.4, 0.0);

    let mut o = Atom::new();
    o.name = "O".to_string();
    o.set_atomic_number(8);
    o.xyz = Position::new(1.3, 2.4, 0.0);

    let mut cb = Atom::new();
    cb.name = "CB".to_string();
    cb.set_atomic_number(6);
    cb.xyz = Position::new(2.0, -0.7, 1.2);

    ag.set_atom("N", n);
    ag.set_atom("CA", ca);
    ag.set_atom("C", c);
    ag.set_atom("O", o);
    ag.set_atom("CB", cb);

    // Pre-register C-O bond as order 1 (e.g. from an explicit file source like CONECT)
    let c_ref = ag.get_atom("C").unwrap().clone();
    let o_ref = ag.get_atom("O").unwrap().clone();
    ag.add_bond(&c_ref, &o_ref, 1);

    assert_eq!(ag.bonds().len(), 1);
    assert_eq!(ag.bonds()[0].order, 1);

    // Call smart default AtomGroup::setup()
    ag.setup().expect("AtomGroup::setup failed");

    // Verify:
    // 1. Total heavy atom bonds for ALA is 4 (N-CA, CA-C, C-O, CA-CB)
    let bonds = ag.get_bond_list();
    assert_eq!(bonds.len(), 4, "Total bonds should be 4 without duplicates");

    // 2. Pre-registered C-O bond order 1 was preserved (not overwritten by CCD's order 2)
    let co_bonds: Vec<_> = bonds
        .iter()
        .filter(|b| {
            (b.atom1_path.ends_with("/C") && b.atom2_path.ends_with("/O"))
                || (b.atom1_path.ends_with("/O") && b.atom2_path.ends_with("/C"))
        })
        .collect();
    assert_eq!(co_bonds.len(), 1, "C-O bond must be unique");
    assert_eq!(
        co_bonds[0].order, 1,
        "Pre-existing C-O bond order 1 must be strictly preserved"
    );

    // 3. No duplicates
    let mut seen = std::collections::HashSet::new();
    for b in &bonds {
        let key = if b.atom1_path <= b.atom2_path {
            (&b.atom1_path, &b.atom2_path)
        } else {
            (&b.atom2_path, &b.atom1_path)
        };
        assert!(seen.insert(key), "Duplicate bond found: {:?}", key);
    }
}

#[test]
fn test_atomgroup_setup_with_db_custom() {
    let mut ag = AtomGroup::new();
    ag.set_path("/model_1/A/1/".to_string());
    ag.name = "XYZ".to_string();

    let mut a1 = Atom::new();
    a1.name = "A1".to_string();
    a1.set_atomic_number(6);
    a1.xyz = Position::new(0.0, 0.0, 0.0);

    let mut a2 = Atom::new();
    a2.name = "A2".to_string();
    a2.set_atomic_number(8);
    a2.xyz = Position::new(1.23, 0.0, 0.0);

    ag.set_atom("A1", a1);
    ag.set_atom("A2", a2);

    // Custom database with XYZ ligand where A1=A2 is order 2
    let mut custom_db = CcdTemplateDb::new();
    custom_db.insert(CcdBondTemplate {
        comp_id: "XYZ".to_string(),
        atoms: vec![synthetic_atom("A1"), synthetic_atom("A2")],
        bonds: vec![("A1".to_string(), "A2".to_string(), 2)],
    });

    ag.setup_with_db(&custom_db).expect("setup_with_db failed");
    let bonds = ag.get_bond_list();
    assert_eq!(bonds.len(), 1);
    assert_eq!(bonds[0].order, 2);
}

#[test]
fn test_loaders_without_bonds_return_empty_bonds_until_setup() {
    let data_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");

    // 1. PDB: 1hls.pdb (has no CONECT records)
    let pdb_path = data_dir.join("1hls.pdb");
    let mut pdb = Pdb::new(None);
    pdb.load(&pdb_path).expect("failed to load 1hls.pdb");
    let mut ag_pdb = pdb
        .get_atomgroup(None, None)
        .expect("failed to get PDB atomgroup");
    assert!(
        ag_pdb.get_bond_list().is_empty(),
        "PDB loader must return empty bonds when file has no explicit bonds"
    );
    ag_pdb.setup().expect("ag_pdb.setup failed");
    assert!(
        !ag_pdb.get_bond_list().is_empty(),
        "PDB atomgroup should resolve bonds after explicit setup()"
    );

    // 2. mmCIF: macromolecule structure without _struct_conn
    const MMCIF_NO_BONDS: &str = "\
data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_asym_id
_atom_site.auth_comp_id
_atom_site.auth_seq_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C C1 . ETH A 1 0.000 0.000 0.000 A ETH 1 C1 1
ATOM 2 C C2 . ETH A 1 1.540 0.000 0.000 A ETH 1 C2 1
";
    let mut mmcif = SimpleMmcif::new();
    mmcif
        .load_from_str(MMCIF_NO_BONDS)
        .expect("parse mmcif without bonds failed");
    let mut ag_cif = mmcif
        .get_atomgroup("test")
        .expect("failed to get mmCIF atomgroup");
    assert!(
        ag_cif.get_bond_list().is_empty(),
        "mmCIF loader must return empty bonds when file has no explicit bonds"
    );
    ag_cif.setup().expect("ag_cif.setup failed");
    assert_eq!(
        ag_cif.get_bond_list().len(),
        1,
        "mmCIF atomgroup should resolve bonds after explicit setup()"
    );

    // 3. PRMTOP: string with C1 and H1 spaced by 1.09 Å and no BONDS section
    const PRMTOP_NO_BONDS: &str = "\
%VERSION  VERSION_STAMP = V0001.000  DATE = 08/25/26  12:00:00
%FLAG ATOM_NAME
%FORMAT(20a4)
C1  H1  
%FLAG CHARGE
%FORMAT(5E16.8)
 0.00000000E+00 0.00000000E+00
%FLAG ATOMIC_NUMBER
%FORMAT(10I8)
       6       1
";
    const INPCRD_TWO_ATOMS: &str = "\
default_name
    2
   0.0000000   0.0000000   0.0000000   1.0900000   0.0000000   0.0000000
";
    let amber = AmberPrmtop::from_strings(PRMTOP_NO_BONDS, INPCRD_TWO_ATOMS).unwrap();
    let mut ag_amber = amber.get_atomgroup().unwrap();
    assert!(
        ag_amber.get_bond_list().is_empty(),
        "PRMTOP loader must return empty bonds when BONDS section is absent"
    );
    ag_amber.setup().expect("ag_amber.setup failed");
    assert_eq!(
        ag_amber.get_bond_list().len(),
        1,
        "PRMTOP atomgroup should resolve covalent bonds after explicit setup()"
    );

    // 4. GRO: sample.gro (GRO format has no bond records)
    let gro_path = data_dir.join("sample.gro");
    let gro = SimpleGro::from_file(&gro_path).expect("failed to load sample.gro");
    let mut ag_gro = gro.get_atomgroup().expect("failed to get GRO atomgroup");
    assert!(
        ag_gro.get_bond_list().is_empty(),
        "GRO loader must return empty bonds"
    );
    ag_gro.setup().expect("ag_gro.setup failed");
    assert!(
        !ag_gro.get_bond_list().is_empty(),
        "GRO atomgroup should resolve bonds after explicit setup()"
    );

    // 5. MOL2: Mol2 without @<TRIPOS>BOND section
    let mol2_no_bonds_str = "\
@<TRIPOS>MOLECULE
ethane_fragment
2 0 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
      1 C1          0.0000    0.0000    0.0000 C.3       1 ETH       0.0000
      2 C2          1.5400    0.0000    0.0000 C.3       1 ETH       0.0000
";
    let mut mol2_no_bonds = SimpleMol2::new();
    mol2_no_bonds
        .parse_str(mol2_no_bonds_str)
        .expect("parse mol2 without bonds failed");
    let mut ag_mol2 = mol2_no_bonds.get_atomgroup().clone();
    assert!(
        ag_mol2.get_bond_list().is_empty(),
        "MOL2 loader must return empty bonds when @<TRIPOS>BOND is absent"
    );
    ag_mol2.setup().expect("ag_mol2.setup failed");
    let mol2_bonds = ag_mol2.get_bond_list();
    assert_eq!(
        mol2_bonds.len(),
        1,
        "MOL2 atomgroup should resolve bonds after explicit setup()"
    );
}

#[test]
fn test_loaders_preserve_explicit_file_bonds() {
    let data_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");

    // MOL2: sample.mol2 has 8 explicit bonds defined in @<TRIPOS>BOND
    let mol2_path = data_dir.join("sample.mol2");
    let mol2 = SimpleMol2::from_file(&mol2_path).expect("failed to load sample.mol2");
    let mut ag_mol2 = mol2.get_atomgroup().clone();
    let bonds = ag_mol2.get_bond_list();
    assert_eq!(
        bonds.len(),
        8,
        "MOL2 loader must preserve exactly the 8 explicit file-derived bonds without duplication"
    );
    assert!(
        bonds.iter().all(|b| b.order == 1),
        "All explicit bonds should retain order 1"
    );

    // mmCIF: 1HLS.cif has 3 explicit disulfide bonds in _struct_conn per model across 20 models (total 60 bonds)
    let mmcif_path = data_dir.join("1HLS.cif");
    let mut mmcif = SimpleMmcif::new();
    mmcif.load(&mmcif_path).expect("failed to load 1HLS.cif");
    let mut ag_cif = mmcif
        .get_atomgroup("1HLS")
        .expect("failed to get mmCIF atomgroup");
    let cif_bonds = ag_cif.get_bond_list();
    assert_eq!(
        cif_bonds.len(),
        60,
        "mmCIF loader must preserve exactly the 60 explicit _struct_conn disulfide bonds across 20 models without implicit setup"
    );
    let mut model_1 = ag_cif.get_group("model_1").unwrap().clone();
    assert_eq!(
        model_1.get_bond_list().len(),
        3,
        "model_1 must have exactly 3 disulfide bonds"
    );
}
