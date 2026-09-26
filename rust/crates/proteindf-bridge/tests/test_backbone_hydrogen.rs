// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::path::PathBuf;

use proteindf_bridge::atom::Atom;
use proteindf_bridge::atom_group::AtomGroup;
use proteindf_bridge::backbone_hydrogen::{
    add_backbone_hydrogens_to_residue, add_backbone_hydrogens_to_residue_in_place,
    build_backbone_amide_hydrogen, build_nterm_hydrogens, STANDARD_AMIDE_NH_BOND_LENGTH,
};
use proteindf_bridge::format::Pdb;
use proteindf_bridge::position::Position;

fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

/// Helper: calculates Euclidean distance between two positions.
fn distance(p1: &Position, p2: &Position) -> f64 {
    let diff = *p1 - *p2;
    diff.length()
}

/// Helper: calculates angle between vector BA and BC in radians.
fn angle_rad(a: &Position, b: &Position, c: &Position) -> f64 {
    let v1 = *a - *b;
    let v2 = *c - *b;
    let len1 = v1.length();
    let len2 = v2.length();
    if len1 < 1e-12 || len2 < 1e-12 {
        return 0.0;
    }
    let cos = (v1.dot(&v2) / (len1 * len2)).clamp(-1.0, 1.0);
    cos.acos()
}

/// Helper: creates a copy of a residue group containing only heavy atoms (strips hydrogens).
fn strip_hydrogens(res: &AtomGroup) -> AtomGroup {
    let mut res_heavy = AtomGroup::with_name(&res.name);
    for (k, atom) in res.atoms() {
        if atom.atomic_number() != 1 {
            res_heavy.set_atom(k, atom.clone());
        }
    }
    res_heavy
}

/// Helper: retrieves a residue group by chain and residue key from a loaded structure.
fn get_residue<'a>(protein: &'a AtomGroup, chain_id: &str, res_key: &str) -> &'a AtomGroup {
    let model = protein
        .get_group("model_1")
        .or_else(|| protein.get_group("1"))
        .expect("model group must exist");
    let chain = model
        .get_group(chain_id)
        .unwrap_or_else(|| panic!("chain {chain_id} must exist"));
    chain
        .get_group(res_key)
        .unwrap_or_else(|| panic!("residue {res_key} in chain {chain_id} must exist"))
}

// 1. Geometric validation: Hand-derived planar geometry test
// Tests that build_backbone_amide_hydrogen strictly obeys:
// - Exact bond length STANDARD_AMIDE_NH_BOND_LENGTH (1.01 A)
// - Lies in the plane of (C_prev, N, CA)
// - Exactly bisects the angle between (C_prev -> N) and (CA -> N)
#[test]
fn test_backbone_amide_geometry_hand_derived() {
    // Construct planar coordinates in the XY plane:
    // N at origin (0, 0, 0)
    // C_prev at (-1.33, 0.0, 0.0) -> vector from N to C_prev is along -X
    // CA at (0.7, 1.2, 0.0) -> in XY plane
    let n = Position::new(0.0, 0.0, 0.0);
    let c_prev = Position::new(-1.33, 0.0, 0.0);
    let ca = Position::new(0.7, 1.2, 0.0);

    let h_opt = build_backbone_amide_hydrogen(&c_prev, &n, &ca);
    assert!(
        h_opt.is_some(),
        "Hydrogen position should be successfully computed"
    );
    let h = h_opt.unwrap();

    // 1. Bond length N-H must match standard literature value exactly (1.01 A)
    let bond_len = distance(&n, &h);
    assert!(
        (bond_len - STANDARD_AMIDE_NH_BOND_LENGTH).abs() < 1e-6,
        "Bond length N-H was {bond_len}, expected {STANDARD_AMIDE_NH_BOND_LENGTH}"
    );

    // 2. Planarity: H must lie in the XY plane (z ~ 0.0)
    assert!(
        h.z.abs() < 1e-6,
        "H.z was {}, expected 0.0 (coplanar with C_prev, N, CA)",
        h.z
    );

    // 3. Angle bisector: angle(C_prev, N, H) == angle(CA, N, H)
    let ang_cn_h = angle_rad(&c_prev, &n, &h);
    let ang_can_h = angle_rad(&ca, &n, &h);
    assert!(
        (ang_cn_h - ang_can_h).abs() < 1e-6,
        "N-H must bisect the C-N-CA angle: ang(C-N-H)={:.4} rad, ang(CA-N-H)={:.4} rad",
        ang_cn_h,
        ang_can_h
    );
}

// 2. Real PDB validation against NMR structure with observed hydrogens (1hls.pdb)
// In 1hls.pdb (NMR structure of insulin), experimental hydrogens are present.
// We verify that calculated backbone amide H directions qualitatively match the
// experimental coordinates for internal residues.
#[test]
fn test_backbone_amide_hydrogen_real_fixture_1hls() {
    let pdb_path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&pdb_path, None).expect("failed to load 1hls.pdb");
    let protein = pdb.get_atomgroup(None, None).expect("failed to get ag");

    // Test Chain A residues 2..21 (insulin A chain)
    // Predecessor of A/2 is A/1
    for resid in 2..=21 {
        let prev_key = (resid - 1).to_string();
        let curr_key = resid.to_string();

        let prev_res = get_residue(&protein, "A", &prev_key);
        let curr_res = get_residue(&protein, "A", &curr_key);

        // PRO residues have no backbone amide H
        if curr_res.name.eq_ignore_ascii_case("PRO") {
            continue;
        }

        let c_prev = &prev_res.get_atom("C").expect("prev C exists").xyz;
        let n_curr = &curr_res.get_atom("N").expect("curr N exists").xyz;
        let ca_curr = &curr_res.get_atom("CA").expect("curr CA exists").xyz;

        let calc_h =
            build_backbone_amide_hydrogen(c_prev, n_curr, ca_curr).expect("calc_h should succeed");

        // If experimental H is present in 1hls, check distance and angle consistency
        if let Some(exp_h) = curr_res.get_atom("H") {
            let dist = distance(&calc_h, &exp_h.xyz);
            // Experimental NMR model coordinates fluctuate, but calculated H should be close (~0.1 - 0.4 A)
            assert!(
                dist < 0.5,
                "Residue A/{} calculated H deviates too far from experimental NMR H: dist = {:.3} A",
                resid,
                dist
            );

            // Vector direction comparison
            let exp_vec = exp_h.xyz - *n_curr;
            let calc_vec = calc_h - *n_curr;
            let cos_sim = exp_vec.dot(&calc_vec) / (exp_vec.length() * calc_vec.length());
            assert!(
                cos_sim > 0.90,
                "Residue A/{} N-H orientation mismatch: cos_sim = {:.3}",
                resid,
                cos_sim
            );
        }
    }
}

// 3. Proline exclusion test:
// In an internal chain, adding backbone hydrogen to a PRO residue must do nothing (return 0 added).
#[test]
fn test_backbone_hydrogen_proline_skipped() {
    let pdb_path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&pdb_path, None).expect("failed to load 1hls.pdb");
    let protein = pdb.get_atomgroup(None, None).expect("failed to get ag");

    // In 1hls, Chain B residue 28 is PRO
    let prev_res = get_residue(&protein, "B", "27"); // THR
    let pro_res = get_residue(&protein, "B", "28"); // PRO

    let mut pro_copy = pro_res.clone();
    // Strip any existing H
    pro_copy.remove_atom("H");

    let report = add_backbone_hydrogens_to_residue_in_place(&mut pro_copy, Some(prev_res))
        .expect("PRO backbone hydrogenation should succeed");

    assert_eq!(
        report.added_hydrogens, 0,
        "Proline must NOT receive a backbone amide hydrogen"
    );
    assert!(!pro_copy.has_atom("H"));
}

// 4. N-terminal hydrogenation test:
// When prev_residue is None, the residue is treated as an N-terminus and receives
// ammonium hydrogens (H1, H2, H3).
#[test]
fn test_backbone_hydrogen_n_terminus() {
    let pdb_path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&pdb_path, None).expect("failed to load 1hls.pdb");
    let protein = pdb.get_atomgroup(None, None).expect("failed to get ag");

    // In 1hls, Chain A residue 1 is GLY (N-terminal)
    let gly1 = get_residue(&protein, "A", "1");

    let mut gly_stripped = strip_hydrogens(gly1);

    assert!(!gly_stripped.has_atom("H1"));
    assert!(!gly_stripped.has_atom("H2"));
    assert!(!gly_stripped.has_atom("H3"));

    let report = add_backbone_hydrogens_to_residue_in_place(&mut gly_stripped, None)
        .expect("N-terminal hydrogenation should succeed");

    assert_eq!(
        report.added_hydrogens, 3,
        "N-terminal GLY should receive 3 hydrogens"
    );
    assert!(gly_stripped.has_atom("H1"));
    assert!(gly_stripped.has_atom("H2"));
    assert!(gly_stripped.has_atom("H3"));

    // Check bond lengths N-H
    let n = gly_stripped.get_atom("N").unwrap();
    for h_name in ["H1", "H2", "H3"] {
        let h = gly_stripped.get_atom(h_name).unwrap();
        let bond_len = distance(&n.xyz, &h.xyz);
        assert!(
            (bond_len - 1.0).abs() < 1e-4,
            "N-{h_name} bond length was {bond_len}, expected 1.0 A"
        );
    }
}

// 5. N-terminal PRO test:
// When an N-terminal residue is PRO, it receives 2 hydrogens (H1, H2) corresponding to secondary amine.
#[test]
fn test_backbone_hydrogen_n_terminal_proline() {
    let mut pro = AtomGroup::with_name("PRO");
    pro.set_atom(
        "N",
        Atom::new_with_pos("N", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    pro.set_atom(
        "CA",
        Atom::new_with_pos("C", Position::new(1.4, 0.0, 0.0)).unwrap(),
    );

    let report = add_backbone_hydrogens_to_residue_in_place(&mut pro, None)
        .expect("N-terminal PRO hydrogenation should succeed");

    assert_eq!(
        report.added_hydrogens, 2,
        "N-terminal PRO should receive 2 hydrogens"
    );
    assert!(pro.has_atom("H1"));
    assert!(pro.has_atom("H2"));
    assert!(!pro.has_atom("H3"));
}

// 6. Degenerate and missing coordinate error handling:
// Verifies proper BridgeError propagation when required atoms are missing or coincident.
#[test]
fn test_backbone_hydrogen_error_handling() {
    let mut curr = AtomGroup::with_name("ALA");
    let mut prev = AtomGroup::with_name("GLY");

    // Case A: Missing N in current residue
    curr.set_atom(
        "CA",
        Atom::new_with_pos("C", Position::new(1.0, 0.0, 0.0)).unwrap(),
    );
    prev.set_atom(
        "C",
        Atom::new_with_pos("C", Position::new(-1.0, 0.0, 0.0)).unwrap(),
    );
    assert!(add_backbone_hydrogens_to_residue_in_place(&mut curr, Some(&prev)).is_err());

    // Case B: Missing C in preceding residue
    curr.set_atom(
        "N",
        Atom::new_with_pos("N", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    let empty_prev = AtomGroup::with_name("GLY");
    assert!(add_backbone_hydrogens_to_residue_in_place(&mut curr, Some(&empty_prev)).is_err());

    // Case C: Coincident/degenerate coordinates (C_prev == N)
    prev.set_atom(
        "C",
        Atom::new_with_pos("C", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    assert!(add_backbone_hydrogens_to_residue_in_place(&mut curr, Some(&prev)).is_err());
}

// 7. Immutable API and build_nterm_hydrogens direct function test
#[test]
fn test_add_backbone_hydrogens_immutable_and_build_nterm() {
    let n = Position::new(0.0, 0.0, 0.0);
    let ca = Position::new(1.4, 0.0, 0.0);

    let mut ala = AtomGroup::with_name("ALA");
    ala.set_atom("N", Atom::new_with_pos("N", n).unwrap());
    ala.set_atom("CA", Atom::new_with_pos("C", ca).unwrap());

    // Direct build_nterm_hydrogens: 3 hydrogens for standard residue
    let h_vec = build_nterm_hydrogens(&ala).expect("build_nterm_hydrogens should succeed");
    assert_eq!(h_vec.len(), 3);
    for (_name, h) in &h_vec {
        let bond_len = distance(&n, &h.xyz);
        assert!((bond_len - 1.0).abs() < 1e-4);
    }

    // Direct build_nterm_hydrogens: 2 hydrogens for secondary amine (proline)
    let mut pro = AtomGroup::with_name("PRO");
    pro.set_atom("N", Atom::new_with_pos("N", n).unwrap());
    pro.set_atom("CA", Atom::new_with_pos("C", ca).unwrap());
    let h_vec_pro = build_nterm_hydrogens(&pro).expect("build_nterm_hydrogens should succeed");
    assert_eq!(h_vec_pro.len(), 2);

    // Immutable add_backbone_hydrogens_to_residue
    let hydrated = add_backbone_hydrogens_to_residue(&ala, None)
        .expect("immutable hydrogenation should succeed");
    assert!(!ala.has_atom("H1")); // original unmodified
    assert!(hydrated.has_atom("H1"));
    assert!(hydrated.has_atom("H2"));
    assert!(hydrated.has_atom("H3"));
}
