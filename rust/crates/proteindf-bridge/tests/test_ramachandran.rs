// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::path::PathBuf;

use proteindf_bridge::atom::Atom;
use proteindf_bridge::atom_group::AtomGroup;
use proteindf_bridge::format::pdb::Pdb;
use proteindf_bridge::position::{dihedral_angle, Position};
use proteindf_bridge::ramachandran::calc_phi_psi;

fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

fn make_backbone_atom(name: &str, pos: Position) -> Atom {
    let symbol = match name {
        "N" => "N",
        "CA" | "C" => "C",
        other => panic!("unexpected backbone atom name: {}", other),
    };
    let mut a = Atom::new_with_pos(symbol, pos).unwrap();
    a.name = name.to_string();
    a
}

#[test]
fn test_dihedral_angle_geometric_sanity() {
    let p1 = Position::new(0.0, 1.0, 0.0);
    let p2 = Position::new(0.0, 0.0, 0.0);
    let p3 = Position::new(1.0, 0.0, 0.0);

    // Cis conformation: dihedral angle is 0.0 degrees
    let p4_cis = Position::new(1.0, 1.0, 0.0);
    let angle_cis = dihedral_angle(&p1, &p2, &p3, &p4_cis);
    assert!(
        angle_cis.abs() < 1e-10,
        "cis angle expected ~0, got {}",
        angle_cis
    );

    // Trans conformation: dihedral angle is 180.0 (or -180.0) degrees
    let p4_trans = Position::new(1.0, -1.0, 0.0);
    let angle_trans = dihedral_angle(&p1, &p2, &p3, &p4_trans);
    assert!(
        (angle_trans.abs() - 180.0).abs() < 1e-10,
        "trans angle expected ~180, got {}",
        angle_trans
    );

    // +90 degrees: rotated towards -Z in right-handed system looking along +X (p2->p3)
    let p4_pos90 = Position::new(1.0, 0.0, -1.0);
    let angle_pos90 = dihedral_angle(&p1, &p2, &p3, &p4_pos90);
    assert!(
        (angle_pos90 - 90.0).abs() < 1e-10,
        "+90 angle expected 90, got {}",
        angle_pos90
    );

    // -90 degrees: rotated towards +Z in right-handed system looking along +X (p2->p3)
    let p4_neg90 = Position::new(1.0, 0.0, 1.0);
    let angle_neg90 = dihedral_angle(&p1, &p2, &p3, &p4_neg90);
    assert!(
        (angle_neg90 - (-90.0)).abs() < 1e-10,
        "-90 angle expected -90, got {}",
        angle_neg90
    );

    // Degenerate collinear case: does not crash, returns 0.0
    let p1_col = Position::new(0.0, 0.0, 0.0);
    let p2_col = Position::new(1.0, 0.0, 0.0);
    let p3_col = Position::new(2.0, 0.0, 0.0);
    let p4_col = Position::new(3.0, 0.0, 0.0);
    let angle_col = dihedral_angle(&p1_col, &p2_col, &p3_col, &p4_col);
    assert!(angle_col.abs() < 1e-10);
}

#[test]
fn test_ramachandran_real_fixture_1hls() {
    let path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&path, None).expect("failed to load 1hls.pdb");
    let ag = pdb.get_atomgroup(None, None).expect("failed to get ag");

    let model_1 = ag.get_group("model_1").expect("model_1 not found");
    let chain_a = model_1.get_group("A").expect("chain A not found");

    let angles = calc_phi_psi(chain_a);
    assert_eq!(angles.len(), 21);

    // First residue: phi is None, psi is Some
    let first = &angles[0];
    assert_eq!(first.residue_key, "1");
    assert_eq!(first.residue_name, "GLY");
    assert!(first.phi.is_none());
    assert!(first.psi.is_some());

    // Last residue: phi is Some, psi is None
    let last = &angles[20];
    assert_eq!(last.residue_key, "21");
    assert_eq!(last.residue_name, "ASN");
    assert!(last.phi.is_some());
    assert!(last.psi.is_none());

    // Verification against independently computed ground truth in docs/rust-port-handoff.md:
    // Residue 4: phi=70.5993, psi=3.3945
    let res4 = angles.iter().find(|a| a.residue_key == "4").unwrap();
    assert_eq!(res4.residue_name, "GLU");
    let phi4 = res4.phi.unwrap();
    let psi4 = res4.psi.unwrap();
    assert!(
        (phi4 - 70.5993).abs() < 1e-3,
        "Residue 4 phi mismatch: expected 70.5993, got {}",
        phi4
    );
    assert!(
        (psi4 - 3.3945).abs() < 1e-3,
        "Residue 4 psi mismatch: expected 3.3945, got {}",
        psi4
    );

    // Residue 5: phi=121.6311, psi=26.5618
    let res5 = angles.iter().find(|a| a.residue_key == "5").unwrap();
    assert_eq!(res5.residue_name, "GLN");
    let phi5 = res5.phi.unwrap();
    let psi5 = res5.psi.unwrap();
    assert!(
        (phi5 - 121.6311).abs() < 1e-3,
        "Residue 5 phi mismatch: expected 121.6311, got {}",
        phi5
    );
    assert!(
        (psi5 - 26.5618).abs() < 1e-3,
        "Residue 5 psi mismatch: expected 26.5618, got {}",
        psi5
    );

    // Residue 10: phi=84.5507, psi=-99.6314
    let res10 = angles.iter().find(|a| a.residue_key == "10").unwrap();
    assert_eq!(res10.residue_name, "ILE");
    let phi10 = res10.phi.unwrap();
    let psi10 = res10.psi.unwrap();
    assert!(
        (phi10 - 84.5507).abs() < 1e-3,
        "Residue 10 phi mismatch: expected 84.5507, got {}",
        phi10
    );
    assert!(
        (psi10 - (-99.6314)).abs() < 1e-3,
        "Residue 10 psi mismatch: expected -99.6314, got {}",
        psi10
    );
}

#[test]
fn test_ramachandran_missing_atoms_safe_skip() {
    let mut chain = AtomGroup::with_name("A");

    // Residue 1: complete backbone
    let mut res1 = AtomGroup::with_name("1");
    res1.name = "ALA".to_string();
    res1.set_atom("N", make_backbone_atom("N", Position::new(0.0, 0.0, 0.0)));
    res1.set_atom("CA", make_backbone_atom("CA", Position::new(1.0, 0.0, 0.0)));
    res1.set_atom("C", make_backbone_atom("C", Position::new(1.5, 1.0, 0.0)));
    chain.set_group("1", res1);

    // Residue 2: missing C atom (incomplete backbone)
    let mut res2 = AtomGroup::with_name("2");
    res2.name = "GLY".to_string();
    res2.set_atom("N", make_backbone_atom("N", Position::new(2.5, 1.0, 0.0)));
    res2.set_atom("CA", make_backbone_atom("CA", Position::new(3.0, 2.0, 0.0)));
    // C is missing!
    chain.set_group("2", res2);

    // Residue 3: complete backbone
    let mut res3 = AtomGroup::with_name("3");
    res3.name = "VAL".to_string();
    res3.set_atom("N", make_backbone_atom("N", Position::new(4.5, 1.0, 0.0)));
    res3.set_atom("CA", make_backbone_atom("CA", Position::new(5.5, 1.0, 0.0)));
    res3.set_atom("C", make_backbone_atom("C", Position::new(6.0, 2.0, 0.0)));
    chain.set_group("3", res3);

    // Residue 4: complete backbone
    let mut res4 = AtomGroup::with_name("4");
    res4.name = "LEU".to_string();
    res4.set_atom("N", make_backbone_atom("N", Position::new(7.0, 2.0, 0.0)));
    res4.set_atom("CA", make_backbone_atom("CA", Position::new(8.0, 2.0, 0.0)));
    res4.set_atom("C", make_backbone_atom("C", Position::new(8.5, 3.0, 0.0)));
    chain.set_group("4", res4);

    let angles = calc_phi_psi(&chain);

    // Residue 2 was skipped because of missing C
    assert_eq!(angles.len(), 3);
    assert_eq!(angles[0].residue_key, "1");
    assert_eq!(angles[1].residue_key, "3");
    assert_eq!(angles[2].residue_key, "4");

    // Residue 1: phi is None; psi is None because contiguous residue 2 is missing C
    assert!(angles[0].phi.is_none());
    assert!(angles[0].psi.is_none());

    // Residue 3: phi is None because preceding residue 2 is missing; psi is Some because residue 4 is contiguous
    assert!(angles[1].phi.is_none());
    assert!(angles[1].psi.is_some());

    // Residue 4: phi is Some because preceding residue 3 is contiguous; psi is None because it is the last residue
    assert!(angles[2].phi.is_some());
    assert!(angles[2].psi.is_none());
}
