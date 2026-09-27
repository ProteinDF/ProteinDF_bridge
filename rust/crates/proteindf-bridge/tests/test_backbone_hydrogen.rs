// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::path::PathBuf;

use proteindf_bridge::atom::Atom;
use proteindf_bridge::atom_group::AtomGroup;
use proteindf_bridge::backbone_hydrogen::{
    add_backbone_hydrogens_to_residue, add_backbone_hydrogens_to_residue_in_place,
    build_backbone_amide_hydrogen, build_nterm_hydrogens, NH3_TETRAHEDRAL_HALF_ANGLE,
    STANDARD_AMIDE_NH_BOND_LENGTH,
};
use proteindf_bridge::format::Pdb;
use proteindf_bridge::modeling::Modeling;
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

    let mut checked_count = 0;

    // Test Chain A residues 2..21 (insulin A chain)
    // Predecessor of A/2 is A/1
    for resid in 2..=21 {
        let prev_key = (resid - 1).to_string();
        let curr_key = resid.to_string();

        let prev_res = get_residue(&protein, "A", &prev_key);
        let curr_res = get_residue(&protein, "A", &curr_key);

        // PRO residues have no backbone amide H
        if curr_res.name == "PRO" {
            continue;
        }

        let c_prev = &prev_res.get_atom("C").expect("prev C exists").xyz;
        let n_curr = &curr_res.get_atom("N").expect("curr N exists").xyz;
        let ca_curr = &curr_res.get_atom("CA").expect("curr CA exists").xyz;

        let calc_h =
            build_backbone_amide_hydrogen(c_prev, n_curr, ca_curr).expect("calc_h should succeed");

        // If experimental H is present in 1hls, check distance and angle consistency
        if let Some(exp_h) = curr_res.get_atom("H") {
            checked_count += 1;
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

    // Ensure that experimental comparisons were actually executed and not skipped
    assert!(
        checked_count >= 15,
        "Expected at least 15 residues with experimental H compared, but got {checked_count}"
    );
}

// 3. Proline exclusion test:
// In an internal chain, adding backbone hydrogen to a PRO residue must do nothing (return 0 added),
// and any spurious existing H atom (even with PDB serial_name key format like "1234_H") must be removed.
#[test]
fn test_backbone_hydrogen_proline_skipped() {
    let pdb_path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&pdb_path, None).expect("failed to load 1hls.pdb");
    let protein = pdb.get_atomgroup(None, None).expect("failed to get ag");

    // In 1hls, Chain B residue 28 is PRO
    let prev_res = get_residue(&protein, "B", "27"); // THR
    let pro_res = get_residue(&protein, "B", "28"); // PRO

    let mut pro_copy = pro_res.clone();
    // Add a spurious H with PDB-style "{serial}_{name}" storage key
    let mut spurious_h = Atom::new_with_pos("H", Position::new(0.0, 0.0, 0.0)).unwrap();
    spurious_h.name = "H".to_string();
    pro_copy.set_atom("1234_H", spurious_h);
    assert!(pro_copy.has_atom("H"));

    let report = add_backbone_hydrogens_to_residue_in_place(&mut pro_copy, Some(prev_res))
        .expect("PRO backbone hydrogenation should succeed");

    assert_eq!(
        report.added_hydrogens, 0,
        "Proline must NOT receive a backbone amide hydrogen"
    );
    assert!(
        !pro_copy.has_atom("H"),
        "Spurious H with key '1234_H' must be purged on internal PRO"
    );
}

// 4. N-terminal hydrogenation test:
// When prev_residue is None, the residue is treated as an N-terminus and receives
// ammonium hydrogens (H1, H2, H3).
// Verifies that:
// - 3 hydrogens are added with length ~1.0 A
// - Each H-N-CA angle matches tetrahedral angle arccos(-1/3) ~ 109.47 deg.
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

    // Check bond lengths N-H and tetrahedral angle H-N-CA (ideal: arccos(-1/3) ~ 109.47 deg)
    let n = gly_stripped.get_atom("N").unwrap();
    let ca = gly_stripped.get_atom("CA").unwrap();
    let ideal_tet_angle = (-1.0_f64 / 3.0).acos(); // ~1.9106 rad

    for h_name in ["H1", "H2", "H3"] {
        let h = gly_stripped.get_atom(h_name).unwrap();
        let bond_len = distance(&n.xyz, &h.xyz);
        assert!(
            (bond_len - 1.0).abs() < 1e-4,
            "N-{h_name} bond length was {bond_len}, expected 1.0 A"
        );

        let ang = angle_rad(&h.xyz, &n.xyz, &ca.xyz);
        assert!(
            (ang - ideal_tet_angle).abs() < 0.02,
            "H-N-CA angle for {h_name} was {:.2} deg ({:.4} rad), expected {:.2} deg ({:.4} rad)",
            ang.to_degrees(),
            ang,
            ideal_tet_angle.to_degrees(),
            ideal_tet_angle
        );
    }
}

// 5. N-terminal PRO test:
// When an N-terminal residue is PRO, it receives 2 hydrogens (H1, H2) corresponding to secondary amine,
// and any spurious existing H3 (even with key "{serial}_{name}") is removed by name.
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

    // Spurious existing H3 with PDB-style "{serial}_{name}" storage key (e.g. "5678_H3")
    let mut spurious_h3 = Atom::new_with_pos("H", Position::new(0.0, 1.0, 0.0)).unwrap();
    spurious_h3.name = "H3".to_string();
    pro.set_atom("5678_H3", spurious_h3);
    assert!(pro.has_atom("H3"));

    let report = add_backbone_hydrogens_to_residue_in_place(&mut pro, None)
        .expect("N-terminal PRO hydrogenation should succeed");

    assert_eq!(
        report.added_hydrogens, 2,
        "N-terminal PRO should receive 2 hydrogens"
    );
    assert!(pro.has_atom("H1"));
    assert!(pro.has_atom("H2"));
    assert!(
        !pro.has_atom("H3"),
        "Spurious H3 keyed as '5678_H3' must be purged on N-terminal PRO"
    );

    // Atomicity regression test: If build_nterm_hydrogens fails (e.g. missing CA),
    // spurious H3 must NOT be removed from the caller's AtomGroup.
    let mut pro_err = AtomGroup::with_name("PRO");
    pro_err.set_atom(
        "N",
        Atom::new_with_pos("N", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    let mut spurious_h3_err = Atom::new_with_pos("H", Position::new(0.0, 1.0, 0.0)).unwrap();
    spurious_h3_err.name = "H3".to_string();
    pro_err.set_atom("5678_H3", spurious_h3_err);
    assert!(pro_err.has_atom("H3"));

    // Missing CA causes build_nterm_hydrogens to fail
    let err_result = add_backbone_hydrogens_to_residue_in_place(&mut pro_err, None);
    assert!(err_result.is_err(), "Must fail when CA is missing");
    assert!(
        pro_err.has_atom("H3"),
        "Spurious H3 must NOT be removed if hydrogenation fails (atomicity violation)"
    );
}

// 6. Degenerate, missing coordinate, and distance error handling:
// Verifies proper BridgeError propagation when required atoms are missing, coincident, or disconnected.
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

    // Case C: Coincident coordinates (C_prev == N, caught by distance check < 0.8 A)
    prev.set_atom(
        "C",
        Atom::new_with_pos("C", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    assert!(add_backbone_hydrogens_to_residue_in_place(&mut curr, Some(&prev)).is_err());

    // Case D: Collinear coordinates with valid peptide bond distance (1.33 A):
    // C_prev at (-1.33, 0, 0), N at (0, 0, 0), CA at (1.40, 0, 0) -> strictly collinear.
    // Distance C_prev-N is 1.33 A (passes distance guard 0.8-2.5 A), but vectors are opposite,
    // so vec_cn + vec_can = (0, 0, 0) and build_backbone_amide_hydrogen returns None.
    let mut collinear_curr = AtomGroup::with_name("ALA");
    collinear_curr.set_atom(
        "N",
        Atom::new_with_pos("N", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    collinear_curr.set_atom(
        "CA",
        Atom::new_with_pos("C", Position::new(1.40, 0.0, 0.0)).unwrap(),
    );
    let mut collinear_prev = AtomGroup::with_name("GLY");
    collinear_prev.set_atom(
        "C",
        Atom::new_with_pos("C", Position::new(-1.33, 0.0, 0.0)).unwrap(),
    );
    let collinear_err =
        add_backbone_hydrogens_to_residue_in_place(&mut collinear_curr, Some(&collinear_prev));
    assert!(
        collinear_err.is_err(),
        "Must reject collinear C_prev, N, CA geometry"
    );
    let err_msg = collinear_err.unwrap_err().to_string();
    assert!(
        err_msg.contains("degenerate backbone geometry"),
        "Error message should mention degenerate geometry, got: {err_msg}"
    );

    // Case E: Disconnected / out-of-bounds peptide bond distance (> 2.5 A)
    let mut far_prev = AtomGroup::with_name("GLY");
    far_prev.set_atom(
        "C",
        Atom::new_with_pos("C", Position::new(-5.0, 0.0, 0.0)).unwrap(),
    );
    assert!(add_backbone_hydrogens_to_residue_in_place(&mut curr, Some(&far_prev)).is_err());
}

// 7. Immutable API and build_nterm_hydrogens direct function test
#[test]
fn test_add_backbone_hydrogens_immutable_and_build_nterm() {
    let n = Position::new(0.0, 0.0, 0.0);
    let ca = Position::new(1.4, 0.0, 0.0);

    let mut ala = AtomGroup::with_name("ALA");
    ala.set_atom("N", Atom::new_with_pos("N", n).unwrap());
    ala.set_atom("CA", Atom::new_with_pos("C", ca).unwrap());

    // Verify NH3_TETRAHEDRAL_HALF_ANGLE constant value (arccos(1/3) ~ 70.53 deg)
    assert!(
        (NH3_TETRAHEDRAL_HALF_ANGLE - (1.0_f64 / 3.0).acos()).abs() < 1e-12,
        "NH3_TETRAHEDRAL_HALF_ANGLE must equal arccos(1/3)"
    );

    // Direct build_nterm_hydrogens: 3 hydrogens for standard residue
    let h_vec = build_nterm_hydrogens(&ala).expect("build_nterm_hydrogens should succeed");
    assert_eq!(h_vec.len(), 3);
    let ideal_tet_angle = (-1.0_f64 / 3.0).acos(); // ~109.47 deg
    for (_name, h) in &h_vec {
        let bond_len = distance(&n, &h.xyz);
        assert!((bond_len - 1.0).abs() < 1e-4);
        let ang = angle_rad(&h.xyz, &n, &ca);
        assert!(
            (ang - ideal_tet_angle).abs() < 0.02,
            "Direct build_nterm_hydrogens H-N-CA angle was {:.2} deg, expected 109.47 deg",
            ang.to_degrees()
        );
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

// 8. Regression test for arbitary_rotate_matrix (Rodrigues rotation formula, (1, 2) entry ny*nz)
#[test]
fn test_arbitary_rotate_matrix_rodrigues_formula() {
    let modeling = Modeling::new().expect("Modeling should initialize");

    // Choose two general non-axial 3D vectors:
    let v_src = Position::new(1.0, 2.0, 3.0);
    let v_dst = Position::new(3.0, -1.0, 2.0);

    let rot = modeling
        .arbitary_rotate_matrix(v_src, v_dst)
        .expect("arbitary_rotate_matrix should succeed for non-parallel vectors");

    // 1. Matrix orthogonality: R * R^T = I
    for i in 0..3 {
        for j in 0..3 {
            let mut dot = 0.0;
            for k in 0..3 {
                dot += rot.get(i, k).unwrap() * rot.get(j, k).unwrap();
            }
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!(
                (dot - expected).abs() < 1e-6,
                "R * R^T ({i}, {j}) was {dot}, expected {expected}"
            );
        }
    }

    // 2. Determinant should be +1 (proper rotation)
    let g = |r: usize, c: usize| rot.get(r, c).unwrap();
    let det = g(0, 0) * (g(1, 1) * g(2, 2) - g(1, 2) * g(2, 1))
        - g(0, 1) * (g(1, 0) * g(2, 2) - g(1, 2) * g(2, 0))
        + g(0, 2) * (g(1, 0) * g(2, 1) - g(1, 1) * g(2, 0));
    assert!(
        (det - 1.0).abs() < 1e-6,
        "Determinant of rotation matrix was {det}, expected 1.0"
    );

    // 3. arbitary_rotate_matrix(a, b) generates the rotation matrix aligning b to a
    // Rotating normalized v_dst must yield normalized v_src
    let mut rotated = v_dst;
    rotated.norm().unwrap();
    rotated.rotate(&rot).unwrap();

    let mut expected_src = v_src;
    expected_src.norm().unwrap();

    let diff = distance(&rotated, &expected_src);
    assert!(
        diff < 1e-6,
        "Rotated vector deviated from expected direction: diff = {diff}"
    );
}
