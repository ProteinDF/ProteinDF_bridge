// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::path::PathBuf;

use proteindf_bridge::atom::Atom;
use proteindf_bridge::atom_group::AtomGroup;
use proteindf_bridge::bond::{Bond, COVALENT_BOND_TOLERANCE};
use proteindf_bridge::error::BridgeError;
use proteindf_bridge::format::pdb::Pdb;
use proteindf_bridge::periodic_table::PeriodicTable;
use proteindf_bridge::position::Position;

fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

#[test]
fn test_covalent_radii_periodic_table() {
    // Verify Cordero et al. (2008) covalent radii values for common biological elements
    assert!((PeriodicTable::covalent_radius("H").unwrap() - 0.31).abs() < 1e-6);
    assert!((PeriodicTable::covalent_radius("C").unwrap() - 0.76).abs() < 1e-6);
    assert!((PeriodicTable::covalent_radius("N").unwrap() - 0.71).abs() < 1e-6);
    assert!((PeriodicTable::covalent_radius("O").unwrap() - 0.66).abs() < 1e-6);
    assert!((PeriodicTable::covalent_radius("P").unwrap() - 1.07).abs() < 1e-6);
    assert!((PeriodicTable::covalent_radius("S").unwrap() - 1.05).abs() < 1e-6);
    assert!((PeriodicTable::covalent_radius("Fe").unwrap() - 1.32).abs() < 1e-6);
    assert!((PeriodicTable::covalent_radius("Cm").unwrap() - 1.69).abs() < 1e-6);

    // Verify tolerance constant is 0.45 Å (OpenBabel standard)
    assert!((COVALENT_BOND_TOLERANCE - 0.45).abs() < 1e-6);

    // Test error cases for index 0 and >96
    assert!(matches!(
        PeriodicTable::covalent_radius(0),
        Err(BridgeError::CovalentRadiusNotFound(0))
    ));
    assert!(matches!(
        PeriodicTable::covalent_radius(97),
        Err(BridgeError::CovalentRadiusNotFound(97))
    ));
}

#[test]
fn test_non_covalent_contact_excluded() {
    let mut ag = AtomGroup::with_name("synthetic_pairs");

    // Case 1: Two Carbon atoms at 2.5 Å distance
    // - Legacy VDW cutoff: 1.70 + 1.70 + 0.40 = 3.80 Å (FALSE POSITIVE BOND)
    // - Modern Covalent cutoff: 0.76 + 0.76 + 0.45 = 1.97 Å (CORRECTLY NOT BONDED)
    let c1 = Atom::new_with_pos("C", Position::new(0.0, 0.0, 0.0)).unwrap();
    let c2 = Atom::new_with_pos("C", Position::new(2.5, 0.0, 0.0)).unwrap();
    ag.set_atom("C1", c1);
    ag.set_atom("C2", c2);

    // Case 2: Nitrogen and Oxygen at typical hydrogen bond distance (2.8 Å)
    // - Legacy VDW cutoff: 1.55 + 1.52 + 0.40 = 3.47 Å (FALSE POSITIVE BOND)
    // - Modern Covalent cutoff: 0.71 + 0.66 + 0.45 = 1.82 Å (CORRECTLY NOT BONDED)
    let n1 = Atom::new_with_pos("N", Position::new(0.0, 10.0, 0.0)).unwrap();
    let o1 = Atom::new_with_pos("O", Position::new(2.8, 10.0, 0.0)).unwrap();
    ag.set_atom("N1", n1);
    ag.set_atom("O1", o1);

    // Case 3: Genuine covalent C-C bond (1.54 Å)
    let c3 = Atom::new_with_pos("C", Position::new(0.0, 20.0, 0.0)).unwrap();
    let c4 = Atom::new_with_pos("C", Position::new(1.54, 20.0, 0.0)).unwrap();
    ag.set_atom("C3", c3);
    ag.set_atom("C4", c4);

    let mut bond = Bond::new();
    bond.setup_heuristic(&mut ag)
        .expect("Bond::setup_heuristic failed");

    let bonds = ag.get_bond_list();

    // Only C3-C4 should be bonded! C1-C2 and N1-O1 must NOT be bonded.
    assert_eq!(
        bonds.len(),
        1,
        "Expected exactly 1 genuine covalent bond, but found {}",
        bonds.len()
    );

    let b = &bonds[0];
    let is_c3_c4 = (b.atom1_path.ends_with("C3") && b.atom2_path.ends_with("C4"))
        || (b.atom1_path.ends_with("C4") && b.atom2_path.ends_with("C3"));
    assert!(is_c3_c4, "The only detected bond should be C3-C4");
}

#[test]
fn test_known_covalent_bonds_1hls() {
    let path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&path, None).expect("failed to load 1hls.pdb");

    let mut ag = pdb
        .get_atomgroup(None, None)
        .expect("failed to get AtomGroup");

    let mut bond = Bond::new();
    bond.setup_heuristic(&mut ag)
        .expect("Bond::setup_heuristic failed on 1hls.pdb");

    let bonds = ag.get_bond_list();
    assert!(!bonds.is_empty(), "Bond list should not be empty");

    // 1. Verify inter-residue peptide bonds (C of residue i to N of residue i+1)
    // Chain A, residue 1 (GIVEQ...) to residue 2 (ILE):
    let peptide_a1_a2 = bonds.iter().find(|b| {
        let is_a1_c = b.atom1_path.contains("/A/1/")
            && (b.atom1_path.ends_with("_C") || b.atom1_path.ends_with("/C"));
        let is_a2_n = b.atom2_path.contains("/A/2/")
            && (b.atom2_path.ends_with("_N") || b.atom2_path.ends_with("/N"));
        let is_a1_n = b.atom1_path.contains("/A/2/")
            && (b.atom1_path.ends_with("_N") || b.atom1_path.ends_with("/N"));
        let is_a2_c = b.atom2_path.contains("/A/1/")
            && (b.atom2_path.ends_with("_C") || b.atom2_path.ends_with("/C"));
        (is_a1_c && is_a2_n) || (is_a1_n && is_a2_c)
    });
    assert!(
        peptide_a1_a2.is_some(),
        "Peptide bond between Chain A residue 1 C and residue 2 N should be detected"
    );

    // Chain B, residue 1 (FVNQH...) to residue 2 (VAL):
    let peptide_b1_b2 = bonds.iter().find(|b| {
        let is_b1_c = b.atom1_path.contains("/B/1/")
            && (b.atom1_path.ends_with("_C") || b.atom1_path.ends_with("/C"));
        let is_b2_n = b.atom2_path.contains("/B/2/")
            && (b.atom2_path.ends_with("_N") || b.atom2_path.ends_with("/N"));
        let is_b1_n = b.atom1_path.contains("/B/2/")
            && (b.atom1_path.ends_with("_N") || b.atom1_path.ends_with("/N"));
        let is_b2_c = b.atom2_path.contains("/B/1/")
            && (b.atom2_path.ends_with("_C") || b.atom2_path.ends_with("/C"));
        (is_b1_c && is_b2_n) || (is_b1_n && is_b2_c)
    });
    assert!(
        peptide_b1_b2.is_some(),
        "Peptide bond between Chain B residue 1 C and residue 2 N should be detected"
    );

    // 2. Verify intra-residue covalent bonds (e.g. Chain A residue 4 GLU: CA-CB, CB-CG, CD-OE1)
    let glu4_ca_cb = bonds.iter().find(|b| {
        let is_ca = b.atom1_path.contains("/A/4/")
            && (b.atom1_path.ends_with("_CA") || b.atom1_path.ends_with("/CA"));
        let is_cb = b.atom2_path.contains("/A/4/")
            && (b.atom2_path.ends_with("_CB") || b.atom2_path.ends_with("/CB"));
        let is_ca_2 = b.atom2_path.contains("/A/4/")
            && (b.atom2_path.ends_with("_CA") || b.atom2_path.ends_with("/CA"));
        let is_cb_2 = b.atom1_path.contains("/A/4/")
            && (b.atom1_path.ends_with("_CB") || b.atom1_path.ends_with("/CB"));
        (is_ca && is_cb) || (is_ca_2 && is_cb_2)
    });
    assert!(
        glu4_ca_cb.is_some(),
        "CA-CB bond in GLU residue 4 should be detected"
    );

    // 3. Verify disulfide bond (S-S) detection:
    // 1hls has known intra-chain disulfide bond A6 (CYS) - A11 (CYS)
    let ssbond_a6_a11 = bonds.iter().find(|b| {
        let is_a6 = b.atom1_path.contains("/A/6/")
            && (b.atom1_path.ends_with("_SG") || b.atom1_path.ends_with("/SG"));
        let is_a11 = b.atom2_path.contains("/A/11/")
            && (b.atom2_path.ends_with("_SG") || b.atom2_path.ends_with("/SG"));
        let is_a6_2 = b.atom2_path.contains("/A/6/")
            && (b.atom2_path.ends_with("_SG") || b.atom2_path.ends_with("/SG"));
        let is_a11_2 = b.atom1_path.contains("/A/11/")
            && (b.atom1_path.ends_with("_SG") || b.atom1_path.ends_with("/SG"));
        (is_a6 && is_a11) || (is_a6_2 && is_a11_2)
    });
    assert!(
        ssbond_a6_a11.is_some(),
        "Disulfide bond between Chain A CYS6 SG and CYS11 SG should be detected"
    );
}
