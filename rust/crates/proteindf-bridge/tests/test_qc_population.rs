// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

//! Tests for orbital-resolved Mulliken atomic population calculations.
//!
//! Verifies:
//! 1. Synthetic 2-atom, 2-AO system with analytically verified expected values.
//! 2. Error handling for dimension mismatches (overlap matrix, AO count, atom index).
//! 3. Molecule with dummy atom without basis functions (contributes 0.0).

use proteindf_bridge::matrix::Matrix;
use proteindf_bridge::qc_result::{
    calc_orbital_mulliken_contributions, AtomicOrbital, MoCoefficients, OrbitalInfo, RunType,
};
use proteindf_bridge::SymmetricMatrix;

/// Tests orbital-resolved Mulliken population calculation on a synthetic 2-atom, 2-AO system.
///
/// # Analytic Derivation of Expected Values:
/// Let atom 0 have AO 0, and atom 1 have AO 1.
/// Overlap matrix S:
///   S = [[1.0, 0.2],
///        [0.2, 1.0]]
///
/// MO coefficients C (columns are MO 0 and MO 1):
///   C_0 = [c_00, c_10]^T = [0.8, 0.5]^T
///   C_1 = [c_01, c_11]^T = [-0.6, 0.7]^T
///
/// Matrix product SC = S * C:
///   (SC)_0 = [1.0*0.8 + 0.2*0.5, 0.2*0.8 + 1.0*0.5]^T = [0.90, 0.66]^T
///   (SC)_1 = [1.0*(-0.6) + 0.2*0.7, 0.2*(-0.6) + 1.0*0.7]^T = [-0.46, 0.58]^T
///
/// MO normalization (C_i^T S C_i):
///   MO 0: 0.8 * 0.90 + 0.5 * 0.66 = 0.72 + 0.33 = 1.05
///   MO 1: (-0.6) * (-0.46) + 0.7 * 0.58 = 0.276 + 0.406 = 0.682
///
/// Atomic contributions q_A^{(i)} = sum_{mu in A} C_{mu, i} * (SC)_{mu, i}:
/// For MO 0:
///   Atom 0 (AO 0): C_{0,0} * (SC)_{0,0} = 0.8 * 0.90 = 0.72
///   Atom 1 (AO 1): C_{1,0} * (SC)_{1,0} = 0.5 * 0.66 = 0.33
///   Sum over atoms = 0.72 + 0.33 = 1.05
///
/// For MO 1:
///   Atom 0 (AO 0): C_{0,1} * (SC)_{0,1} = (-0.6) * (-0.46) = 0.276
///   Atom 1 (AO 1): C_{1,1} * (SC)_{1,1} = 0.7 * 0.58 = 0.406
///   Sum over atoms = 0.276 + 0.406 = 0.682
///
/// If C is normalized so that C_i^T S C_i = 1.0:
/// Normalized MO 0:
///   scale_0 = 1.0 / sqrt(1.05) ≈ 0.9759000729485332
///   C_0_norm = C_0 * scale_0
///   q_0^{(0)} = 0.72 / 1.05 ≈ 0.6857142857142857
///   q_1^{(0)} = 0.33 / 1.05 ≈ 0.3142857142857143
///   q_0^{(0)} + q_1^{(0)} = 1.0
///
/// Normalized MO 1:
///   scale_1 = 1.0 / sqrt(0.682) ≈ 1.210887361718021
///   C_1_norm = C_1 * scale_1
///   q_0^{(1)} = 0.276 / 0.682 ≈ 0.40469208211143697
///   q_1^{(1)} = 0.406 / 0.682 ≈ 0.595307917888563
///   q_0^{(1)} + q_1^{(1)} = 1.0
#[test]
fn test_calc_orbital_mulliken_contributions_synthetic() {
    let mut s = SymmetricMatrix::new(2);
    s.set(0, 0, 1.0);
    s.set(0, 1, 0.2);
    s.set(1, 1, 1.0);

    let scale_0 = 1.0 / 1.05_f64.sqrt();
    let scale_1 = 1.0 / 0.682_f64.sqrt();

    let mut c_mat = Matrix::new(2, 2);
    c_mat.set(0, 0, 0.8 * scale_0);
    c_mat.set(1, 0, 0.5 * scale_0);
    c_mat.set(0, 1, -0.6 * scale_1);
    c_mat.set(1, 1, 0.7 * scale_1);

    let mo_coeffs = MoCoefficients::new(RunType::Rks, 1, c_mat);

    let orb_info = OrbitalInfo {
        num_of_aos: 2,
        orbitals: vec![
            AtomicOrbital {
                ao_index: 0,
                atom_index: 0,
                atom_symbol: "H".to_string(),
                basisset_name: "test".to_string(),
                cgto_index: 0,
                shell_type: "s".to_string(),
                basis_type: "s".to_string(),
            },
            AtomicOrbital {
                ao_index: 1,
                atom_index: 1,
                atom_symbol: "H".to_string(),
                basisset_name: "test".to_string(),
                cgto_index: 1,
                shell_type: "s".to_string(),
                basis_type: "s".to_string(),
            },
        ],
    };

    let contribs = calc_orbital_mulliken_contributions(&mo_coeffs, &s, &orb_info, 2)
        .expect("calculation failed");

    assert_eq!(contribs.len(), 2);

    // MO 0
    assert_eq!(contribs[0].orbital_index, 0);
    assert_eq!(contribs[0].atom_contributions.len(), 2);
    assert_eq!(contribs[0].atom_contributions[0].0, 0);
    assert_eq!(contribs[0].atom_contributions[1].0, 1);

    let q0_mo0 = contribs[0].atom_contributions[0].1;
    let q1_mo0 = contribs[0].atom_contributions[1].1;
    let expected_q0_mo0 = 0.72 / 1.05;
    let expected_q1_mo0 = 0.33 / 1.05;

    assert!(
        (q0_mo0 - expected_q0_mo0).abs() < 1e-12,
        "MO 0 Atom 0 contribution diff: {q0_mo0} vs {expected_q0_mo0}"
    );
    assert!(
        (q1_mo0 - expected_q1_mo0).abs() < 1e-12,
        "MO 0 Atom 1 contribution diff: {q1_mo0} vs {expected_q1_mo0}"
    );
    assert!(
        ((q0_mo0 + q1_mo0) - 1.0).abs() < 1e-12,
        "MO 0 sum across atoms must be 1.0"
    );

    // MO 1
    assert_eq!(contribs[1].orbital_index, 1);
    assert_eq!(contribs[1].atom_contributions.len(), 2);
    let q0_mo1 = contribs[1].atom_contributions[0].1;
    let q1_mo1 = contribs[1].atom_contributions[1].1;
    let expected_q0_mo1 = 0.276 / 0.682;
    let expected_q1_mo1 = 0.406 / 0.682;

    assert!(
        (q0_mo1 - expected_q0_mo1).abs() < 1e-12,
        "MO 1 Atom 0 contribution diff: {q0_mo1} vs {expected_q0_mo1}"
    );
    assert!(
        (q1_mo1 - expected_q1_mo1).abs() < 1e-12,
        "MO 1 Atom 1 contribution diff: {q1_mo1} vs {expected_q1_mo1}"
    );
    assert!(
        ((q0_mo1 + q1_mo1) - 1.0).abs() < 1e-12,
        "MO 1 sum across atoms must be 1.0"
    );
}

/// Tests that a dummy atom (e.g. 'X' with atom index 1) without basis functions
/// correctly receives 0.0 contribution in all molecular orbitals.
#[test]
fn test_calc_orbital_mulliken_dummy_atom() {
    let mut s = SymmetricMatrix::new(2);
    s.set(0, 0, 1.0);
    s.set(0, 1, 0.0);
    s.set(1, 1, 1.0);

    let mut c_mat = Matrix::new(2, 2);
    c_mat.set(0, 0, 1.0);
    c_mat.set(1, 0, 0.0);
    c_mat.set(0, 1, 0.0);
    c_mat.set(1, 1, 1.0);

    let mo_coeffs = MoCoefficients::new(RunType::Rks, 1, c_mat);

    // 3 atoms: Atom 0 (AO 0), Atom 1 (dummy atom, no AO), Atom 2 (AO 1)
    let orb_info = OrbitalInfo {
        num_of_aos: 2,
        orbitals: vec![
            AtomicOrbital {
                ao_index: 0,
                atom_index: 0,
                atom_symbol: "C".to_string(),
                basisset_name: "test".to_string(),
                cgto_index: 0,
                shell_type: "s".to_string(),
                basis_type: "s".to_string(),
            },
            AtomicOrbital {
                ao_index: 1,
                atom_index: 2,
                atom_symbol: "O".to_string(),
                basisset_name: "test".to_string(),
                cgto_index: 1,
                shell_type: "s".to_string(),
                basis_type: "s".to_string(),
            },
        ],
    };

    let contribs = calc_orbital_mulliken_contributions(&mo_coeffs, &s, &orb_info, 3)
        .expect("calculation failed");

    assert_eq!(contribs.len(), 2);
    for contrib in &contribs {
        assert_eq!(contrib.atom_contributions.len(), 3);
        // Atom 1 (dummy atom) contribution must be exactly 0.0
        assert_eq!(contrib.atom_contributions[1], (1, 0.0));
        let sum: f64 = contrib.atom_contributions.iter().map(|(_, val)| *val).sum();
        assert!((sum - 1.0).abs() < 1e-12);
    }
}

/// Tests error handling for dimension mismatches:
/// - S dimension does not match C rows
/// - OrbitalInfo AO count does not match C rows
/// - AO atom_index >= num_atoms
#[test]
fn test_calc_orbital_mulliken_dimension_mismatches() {
    let s_2x2 = SymmetricMatrix::new(2);
    let s_3x3 = SymmetricMatrix::new(3);

    let c_mat = Matrix::new(2, 2);
    let mo_coeffs = MoCoefficients::new(RunType::Rks, 1, c_mat);

    let orb_info_2 = OrbitalInfo {
        num_of_aos: 2,
        orbitals: vec![
            AtomicOrbital {
                ao_index: 0,
                atom_index: 0,
                atom_symbol: "H".to_string(),
                basisset_name: "test".to_string(),
                cgto_index: 0,
                shell_type: "s".to_string(),
                basis_type: "s".to_string(),
            },
            AtomicOrbital {
                ao_index: 1,
                atom_index: 1,
                atom_symbol: "H".to_string(),
                basisset_name: "test".to_string(),
                cgto_index: 1,
                shell_type: "s".to_string(),
                basis_type: "s".to_string(),
            },
        ],
    };

    let orb_info_3 = OrbitalInfo {
        num_of_aos: 3,
        orbitals: vec![],
    };

    // 1. S dimension mismatch
    let err1 = calc_orbital_mulliken_contributions(&mo_coeffs, &s_3x3, &orb_info_2, 2);
    assert!(err1.is_err(), "should fail on S dimension mismatch");

    // 2. OrbitalInfo AO count mismatch
    let err2 = calc_orbital_mulliken_contributions(&mo_coeffs, &s_2x2, &orb_info_3, 2);
    assert!(
        err2.is_err(),
        "should fail on OrbitalInfo AO count mismatch"
    );

    // 3. atom_index >= num_atoms
    let err3 = calc_orbital_mulliken_contributions(&mo_coeffs, &s_2x2, &orb_info_2, 1);
    assert!(err3.is_err(), "should fail when atom_index >= num_atoms");
}
