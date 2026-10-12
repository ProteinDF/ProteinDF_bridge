// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

//! Tests for ProteinDF HDF5 reader (`qc_result::PdfParamHdf5`).
//!
//! Expected values were determined by inspecting the actual test fixtures
//! (`rust/crates/proteindf-bridge/tests/data/*_reference.h5`) using `h5py` and
//! `ProteinDF_pytools` (`PdfParam_H5`, `OrbInfo`, `basisset.py`) in an isolated
//! temporary Python virtual environment outside the repository (`/tmp/.../vh5`).

#![cfg(feature = "hdf5")]

use std::path::PathBuf;

use proteindf_bridge::qc_result::{PdfParamHdf5, RunType};
use proteindf_bridge::BOHR_TO_ANGSTROM;

fn fixture_path(filename: &str) -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("tests");
    path.push("data");
    path.push(filename);
    path
}

/// Tests reading `H2O_HF_reference.h5`:
/// - RKS method, 3 atoms, 24 AOs, 24 MOs, 49 iterations, converged = true.
/// - Molecule: O, H, H with coordinates in Bohr.
/// - Basis set: O-DZVP2.H and O-DZVP2.O, 24 expanded AOs according to ProteinDF rules.
/// - Molecular orbital coefficients, orbital energies, occupations, Mulliken charges.
/// - Verifies optional fields not present in `reference.h5` (`s`, `P`, `F`, `h`, `h2`) are None / empty.
#[test]
fn test_h2o_hf_reference() {
    let path = fixture_path("H2O_HF_reference.h5");
    let result = PdfParamHdf5::load(&path).expect("failed to load H2O_HF_reference.h5");

    // 1. Parameters (verified against root attributes in h5py)
    assert_eq!(result.params.method, "rks");
    assert_eq!(result.params.xc_functional, "HF");
    assert_eq!(result.params.num_of_atoms, 3);
    assert_eq!(result.params.num_of_aos, 24);
    assert_eq!(result.params.num_of_mos, 24);
    assert_eq!(result.params.iterations, 49);
    assert_eq!(result.params.max_iterations, 100);
    assert!(result.params.scf_converged);
    assert_eq!(
        result.params.step_control,
        "create integral guess scf force"
    );
    assert_eq!(result.params.guess, "core");
    assert_eq!(result.params.cut_value, "1.0e-10");
    assert_eq!(result.params.cdam_tau, "1.0E-10");
    assert_eq!(result.params.cd_epsilon, "1.0E-6");
    assert_eq!(result.params.scf_acceleration, "damping");
    assert_eq!(result.params.damping_factor, Some(0.67));

    // 2. Molecule (verified against molecule/xyzc in h5py: O at origin, H at x=1.8897..., H at y=1.8897...)
    assert_eq!(result.molecule.len(), 3);
    let o = &result.molecule.atoms[0];
    assert_eq!(o.symbol, "O");
    assert_eq!(o.charge, 8.0);
    assert_eq!(o.pos_bohr.x, 0.0);
    assert_eq!(o.pos_bohr.y, 0.0);
    assert_eq!(o.pos_bohr.z, 0.0);

    let h1 = &result.molecule.atoms[1];
    assert_eq!(h1.symbol, "H");
    assert_eq!(h1.charge, 1.0);
    assert!((h1.pos_bohr.x - 1.889_726_125_457_828_1).abs() < 1e-12);
    assert_eq!(h1.pos_bohr.y, 0.0);
    assert_eq!(h1.pos_bohr.z, 0.0);

    let h2 = &result.molecule.atoms[2];
    assert_eq!(h2.symbol, "H");
    assert_eq!(h2.charge, 1.0);
    assert_eq!(h2.pos_bohr.x, 0.0);
    assert!((h2.pos_bohr.y - 1.889_726_125_457_828_1).abs() < 1e-12);
    assert_eq!(h2.pos_bohr.z, 0.0);

    // AtomGroup conversion: Bohr to Angstrom conversion check (1.8897261254578281 Bohr * BOHR_TO_ANGSTROM ≈ 1.0 Å)
    let ag = result
        .molecule
        .get_atomgroup(Some("water"))
        .expect("failed to convert to AtomGroup");
    assert_eq!(ag.name, "water");
    assert_eq!(ag.get_number_of_atoms(), 3);
    let atom_h1 = ag.get_atom("1").expect("atom 1 not found");
    assert!((atom_h1.xyz.x - 1.889_726_125_457_828_1 * BOHR_TO_ANGSTROM).abs() < 1e-10);
    assert!((atom_h1.xyz.x - 1.0).abs() < 1e-5);

    // 3. Basis set and OrbitalInfo (verified against ProteinDF OrbInfo rules: 24 AOs)
    assert_eq!(result.basisset.definitions.len(), 2);
    assert_eq!(result.orbital_info.num_of_aos, 24);
    assert_eq!(result.orbital_info.len(), 24);

    // Atom 0 (O): CGTOs: 3 s-shells (s), 2 p-shells (px,py,pz), 1 d-shell (dxy,dyz,dzx,dxx-yy,dzz) -> 3 + 6 + 5 = 14 AOs
    assert_eq!(result.orbital_info.atom_id(0), Some(0));
    assert_eq!(result.orbital_info.orbitals[0].shell_type, "s");
    assert_eq!(result.orbital_info.orbitals[0].basis_type, "s");
    assert_eq!(result.orbital_info.orbitals[3].shell_type, "p");
    assert_eq!(result.orbital_info.orbitals[3].basis_type, "px");
    assert_eq!(result.orbital_info.orbitals[4].basis_type, "py");
    assert_eq!(result.orbital_info.orbitals[5].basis_type, "pz");
    assert_eq!(result.orbital_info.orbitals[9].shell_type, "d");
    assert_eq!(result.orbital_info.orbitals[9].basis_type, "dxy");
    assert_eq!(result.orbital_info.orbitals[13].basis_type, "dzz");

    // Atom 1 (H): CGTOs: 2 s-shells (s), 1 p-shell (px,py,pz) -> 2 + 3 = 5 AOs (AOs 14..18)
    assert_eq!(result.orbital_info.atom_id(14), Some(1));
    assert_eq!(result.orbital_info.orbitals[14].atom_symbol, "H");
    assert_eq!(result.orbital_info.orbitals[14].shell_type, "s");
    assert_eq!(result.orbital_info.orbitals[16].shell_type, "p");
    assert_eq!(result.orbital_info.orbitals[16].basis_type, "px");

    // Atom 2 (H): CGTOs: 2 s-shells, 1 p-shell -> 5 AOs (AOs 19..23)
    assert_eq!(result.orbital_info.atom_id(19), Some(2));
    assert_eq!(result.orbital_info.atom_id(23), Some(2));

    // 4. Total energies across iterations (/TEs: length 49)
    assert_eq!(result.total_energies.len(), 49);
    assert!((result.total_energies[0] - (-60.633_567_801_727_23)).abs() < 1e-10);
    assert!((result.total_energies[48] - (-76.032_302_389_239_29)).abs() < 1e-10);

    // 5. MO coefficients (/C/rks_49: shape=(24, 24))
    let c = result
        .get_mo_coefficients(&RunType::Rks, Some(49))
        .expect("missing MO coefficients for iteration 49");
    assert_eq!(c.n_ao(), 24);
    assert_eq!(c.n_mo(), 24);
    assert_eq!(c.matrix.rows(), 24);
    assert_eq!(c.matrix.cols(), 24);
    assert_eq!(c.iteration, 49);
    // C[0, 0] = 1.00053064, C[0, 1] = -0.228648715, C[1, 0] = -0.00684479965 (from h5py)
    assert!((c.get(0, 0).unwrap() - 1.000_530_64).abs() < 1e-6);
    assert!((c.get(0, 1).unwrap() - (-0.228_648_715)).abs() < 1e-6);
    assert!((c.get(1, 0).unwrap() - (-0.006_844_799_65)).abs() < 1e-6);

    let mo0 = c.mo_vector(0).expect("failed to extract MO 0 vector");
    assert_eq!(mo0.len(), 24);
    assert!((mo0[0] - 1.000_530_64).abs() < 1e-6);

    // 6. Orbital energies (/energy_level: 49 iterations present, final iteration 49)
    assert_eq!(result.orbital_energies.len(), 49);
    let oe49 = result
        .get_orbital_energy(&RunType::Rks, None)
        .expect("missing final orbital energy");
    assert_eq!(oe49.iteration, 49);
    assert_eq!(oe49.len(), 24);
    // Energy of orbital 0 (O 1s): -20.56028171 Hartree, HOMO (orb 4): -0.51133860 Hartree
    assert!((oe49.energies[0] - (-20.560_281_71)).abs() < 1e-6);
    assert!((oe49.energies[4] - (-0.511_338_60)).abs() < 1e-6);

    // Check specific earlier iteration (e.g. iteration 1)
    let oe1 = result
        .get_orbital_energy(&RunType::Rks, Some(1))
        .expect("missing iteration 1 orbital energy");
    assert_eq!(oe1.iteration, 1);
    assert!((oe1.energies[0] - (-17.120_732_02)).abs() < 1e-6);

    // 7. Occupation (/occ/rks: 10 electrons in closed shell -> 5 doubly occupied orbitals)
    let occ = result
        .get_orbital_occupation(&RunType::Rks)
        .expect("missing occupation for rks");
    assert_eq!(occ.occupations.len(), 24);
    for i in 0..5 {
        assert_eq!(occ.occupations[i], 2.0);
    }
    for i in 5..24 {
        assert_eq!(occ.occupations[i], 0.0);
    }

    // 8. Mulliken population (/pop/mulliken/rks_49: O=-0.59339845, H1=+0.29669921, H2=+0.29669924)
    let pop = result
        .get_mulliken_population(&RunType::Rks, None)
        .expect("missing mulliken population");
    assert_eq!(pop.iteration, 49);
    assert_eq!(pop.charges.len(), 3);
    assert!((pop.charges[0] - (-0.593_398_45)).abs() < 1e-6);
    assert!((pop.charges[1] - 0.296_699_21).abs() < 1e-6);
    assert!((pop.charges[2] - 0.296_699_24).abs() < 1e-6);

    // 9. Verify optional matrices that are not in reference.h5 are None / empty
    assert!(result.s_matrix.is_none());
    assert!(result.density_matrices.is_empty());
    assert!(result.fock_matrices.is_empty());
    assert!(result.h_matrix.is_none());
    assert!(result.h2_matrix.is_none());
}

/// Tests reading `Ne_HF_reference.h5`:
/// - Single-atom system (1 atom, 14 AOs, 14 MOs, 30 iterations).
/// - Net Mulliken charge on single neutral atom is 0.0.
#[test]
fn test_ne_hf_reference() {
    let path = fixture_path("Ne_HF_reference.h5");
    let result = PdfParamHdf5::load(&path).expect("failed to load Ne_HF_reference.h5");

    assert_eq!(result.params.method, "rks");
    assert_eq!(result.params.num_of_atoms, 1);
    assert_eq!(result.params.num_of_aos, 14);
    assert_eq!(result.params.num_of_mos, 14);
    assert_eq!(result.params.iterations, 30);
    assert!(result.params.scf_converged);

    // Molecule
    assert_eq!(result.molecule.len(), 1);
    assert_eq!(result.molecule.atoms[0].symbol, "Ne");
    assert_eq!(result.molecule.atoms[0].charge, 10.0);

    // Orbital info: 14 AOs for single Ne
    assert_eq!(result.orbital_info.len(), 14);
    assert_eq!(result.orbital_info.atom_id(0), Some(0));
    assert_eq!(result.orbital_info.atom_id(13), Some(0));

    // MO coefficients: 14x14 matrix
    let c = result
        .get_mo_coefficients(&RunType::Rks, None)
        .expect("missing MO coefficients");
    assert_eq!(c.n_ao(), 14);
    assert_eq!(c.n_mo(), 14);
    assert_eq!(c.matrix.rows(), 14);
    assert_eq!(c.matrix.cols(), 14);
    assert_eq!(c.iteration, 30);
    assert!((c.get(0, 0).unwrap() - 0.998_877_889).abs() < 1e-6);

    // Orbital energy at iteration 30: 1s energy = -32.74420438 Hartree
    let oe = result
        .get_orbital_energy(&RunType::Rks, None)
        .expect("missing orbital energy");
    assert_eq!(oe.iteration, 30);
    assert!((oe.energies[0] - (-32.744_204_38)).abs() < 1e-6);

    // Mulliken charge: [0.0]
    let pop = result
        .get_mulliken_population(&RunType::Rks, None)
        .expect("missing mulliken population");
    assert_eq!(pop.charges.len(), 1);
    assert_eq!(pop.charges[0], 0.0);

    // Non-existent datasets
    assert!(result.s_matrix.is_none());
}

/// Tests reading `N2_UHF_reference.h5`:
/// - UKS method with missing groups (`C`, `energy_level`, `occ`, `pop` are absent in this minimal reference file).
/// - Confirms defensive handling: loader successfully returns without error, missing fields are safely empty/None.
#[test]
fn test_n2_uhf_reference() {
    let path = fixture_path("N2_UHF_reference.h5");
    let result = PdfParamHdf5::load(&path).expect("failed to load N2_UHF_reference.h5");

    assert_eq!(result.params.method, "uks");
    assert_eq!(result.params.xc_functional, "HF");
    assert_eq!(result.params.num_of_atoms, 2);
    assert_eq!(result.params.num_of_aos, 28);
    assert_eq!(result.params.num_of_mos, 28);
    assert_eq!(result.params.iterations, 36);
    assert!(result.params.scf_converged);

    // 2 Nitrogen atoms: N at 0.0, N at 1.8897261254578281 Bohr
    assert_eq!(result.molecule.len(), 2);
    assert_eq!(result.molecule.atoms[0].symbol, "N");
    assert_eq!(result.molecule.atoms[0].charge, 7.0);
    assert_eq!(result.molecule.atoms[1].symbol, "N");
    assert_eq!(result.molecule.atoms[1].charge, 7.0);

    // Basis set: 28 AOs (14 AOs per N atom)
    assert_eq!(result.orbital_info.len(), 28);
    assert_eq!(result.orbital_info.atom_id(0), Some(0));
    assert_eq!(result.orbital_info.atom_id(13), Some(0));
    assert_eq!(result.orbital_info.atom_id(14), Some(1));
    assert_eq!(result.orbital_info.atom_id(27), Some(1));

    // TEs: 36 iterations
    assert_eq!(result.total_energies.len(), 36);
    assert!((result.total_energies[0] - (-91.423_143_074_964_59)).abs() < 1e-10);
    assert!((result.total_energies[35] - (-108.925_183_195_859_99)).abs() < 1e-10);

    // Absent groups in this file should be safely empty / None
    assert!(result.mo_coefficients.is_empty());
    assert!(result
        .get_mo_coefficients(&RunType::UksAlpha, None)
        .is_none());
    assert!(result
        .get_mo_coefficients(&RunType::UksBeta, None)
        .is_none());
    assert!(result.orbital_energies.is_empty());
    assert!(result.orbital_occupations.is_empty());
    assert!(result.mulliken_charges.is_empty());
    assert!(result.s_matrix.is_none());
    assert!(result.density_matrices.is_empty());
    assert!(result.fock_matrices.is_empty());
}
