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

use proteindf_bridge::qc_result::{calc_orbital_mulliken_contributions, PdfParamHdf5, RunType};
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

// -----------------------------------------------------------------------------
// Synthetic HDF5 fixtures for UKS, ROKS, and optional matrices (PR#60 Review Request 2)
//
// The datasets and attributes follow ProteinDF_pytools (`pdfparam_hdf5.py`):
// - Matrix I/O `_set_matrix` (599-611行目): attributes `row`, `col`, `type` ("SY" or "GE")
// - Vector I/O `_set_vector` (628-638行目): attribute `dim`
// - Path generators (643-686行目):
//   - `_get_h5path_occ` (643-645行目): "occ/{run_type}"
//   - `_get_h5path_s` (647-649行目): "s"
//   - `_get_h5path_h` (651-653行目): "h"
//   - `_get_h5path_h2` (655-657行目): "h2"
//   - `_get_h5path_c_matrix` (659-662行目): "C/{run_type}_{iteration}"
//   - `_get_h5path_density_matrix` (664-667行目): "P/{run_type}_{iteration}"
//   - `_get_h5path_f_matrix` (669-672行目): "F/{run_type}_{iteration}"
//   - `_get_h5path_TEs` (674-676行目): "TEs"
//   - `_get_h5path_energy_level` (678-681行目): "energy_level/{run_type}_{iteration}"
//   - `_get_h5path_pop_mulliken_atom` (683-686行目): "pop/mulliken/{run_type}_{iteration}"
// -----------------------------------------------------------------------------

use hdf5_metno::types::FixedAscii;
use std::path::Path;

fn str_ascii<const N: usize>(s: &str) -> FixedAscii<N> {
    FixedAscii::<N>::from_ascii(s).expect("valid ascii")
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SyStorageMode {
    FullSymmetric,
    AsymmetricUpperNonzero,
    LowerTriangleOnly,
}

fn write_synthetic_qc_h5(
    file_path: &Path,
    method: &str,
    sy_mode: SyStorageMode,
) -> Result<(), Box<dyn std::error::Error>> {
    let file = hdf5_metno::File::create(file_path)?;

    // 1. Root attributes
    file.new_attr::<FixedAscii<32>>()
        .create("method")?
        .write_scalar(&str_ascii::<32>(method))?;
    file.new_attr::<FixedAscii<32>>()
        .create("xc_functional")?
        .write_scalar(&str_ascii::<32>("B3LYP"))?;
    file.new_attr::<i64>()
        .create("num_of_atoms")?
        .write_scalar(&2i64)?;
    file.new_attr::<i64>()
        .create("num_of_AOs")?
        .write_scalar(&2i64)?;
    file.new_attr::<i64>()
        .create("num_of_MOs")?
        .write_scalar(&2i64)?;
    file.new_attr::<i64>()
        .create("iterations")?
        .write_scalar(&2i64)?;
    file.new_attr::<bool>()
        .create("scf_converged")?
        .write_scalar(&true)?;

    // 2. Molecule (/molecule/xyzc: 2 atoms, shape (2, 6))
    let mol_data = vec![
        str_ascii::<32>("H"),
        str_ascii::<32>(""),
        str_ascii::<32>("0.000000"),
        str_ascii::<32>("0.000000"),
        str_ascii::<32>("0.000000"),
        str_ascii::<32>("1.000000"),
        str_ascii::<32>("H"),
        str_ascii::<32>(""),
        str_ascii::<32>("0.000000"),
        str_ascii::<32>("0.000000"),
        str_ascii::<32>("1.400000"),
        str_ascii::<32>("1.000000"),
    ];
    let mol_arr = ndarray::Array2::from_shape_vec((2, 6), mol_data)?;
    file.new_dataset_builder()
        .with_data(&mol_arr)
        .create("molecule/xyzc")?;

    // 3. Basis set (/basisset/{name, cgto, pgto})
    let name_data = vec![str_ascii::<32>("H-basis"), str_ascii::<32>("H")];
    let name_arr = ndarray::Array2::from_shape_vec((1, 2), name_data)?;
    file.new_dataset_builder()
        .with_data(&name_arr)
        .create("basisset/name")?;

    let cgto_data = vec![
        str_ascii::<32>("H-basis"),
        str_ascii::<32>("0"),
        str_ascii::<32>("s"),
        str_ascii::<32>("1.0"),
    ];
    let cgto_arr = ndarray::Array2::from_shape_vec((1, 4), cgto_data)?;
    file.new_dataset_builder()
        .with_data(&cgto_arr)
        .create("basisset/cgto")?;

    let pgto_data = vec![
        str_ascii::<32>("H-basis"),
        str_ascii::<32>("0"),
        str_ascii::<32>("0"),
        str_ascii::<32>("1.0"),
        str_ascii::<32>("1.0"),
    ];
    let pgto_arr = ndarray::Array2::from_shape_vec((1, 5), pgto_data)?;
    file.new_dataset_builder()
        .with_data(&pgto_arr)
        .create("basisset/pgto")?;

    // 4. Total energies (/TEs)
    let tes_arr = ndarray::Array1::from(vec![-1.10, -1.15]);
    let ds_tes = file
        .new_dataset_builder()
        .with_data(&tes_arr)
        .create("TEs")?;
    ds_tes
        .new_attr::<i64>()
        .create("dim")?
        .write_scalar(&2i64)?;

    // Helper closure to write a 2D matrix with pdfparam_hdf5 attributes (_set_matrix 599-611行目)
    let write_mat =
        |path: &str, mat_type: &str, data: Vec<f64>| -> Result<(), Box<dyn std::error::Error>> {
            let arr = ndarray::Array2::from_shape_vec((2, 2), data)?;
            let ds = file.new_dataset_builder().with_data(&arr).create(path)?;
            ds.new_attr::<i64>().create("row")?.write_scalar(&2i64)?;
            ds.new_attr::<i64>().create("col")?.write_scalar(&2i64)?;
            ds.new_attr::<FixedAscii<32>>()
                .create("type")?
                .write_scalar(&str_ascii::<32>(mat_type))?;
            Ok(())
        };

    // Helper closure to write a 1D vector with pdfparam_hdf5 attributes (_set_vector 628-638行目)
    let write_vec = |path: &str, data: Vec<f64>| -> Result<(), Box<dyn std::error::Error>> {
        let len = data.len() as i64;
        let arr = ndarray::Array1::from(data);
        let ds = file.new_dataset_builder().with_data(&arr).create(path)?;
        ds.new_attr::<i64>().create("dim")?.write_scalar(&len)?;
        Ok(())
    };

    // 5. Symmetric matrices (s, h, h2)
    let (s_data, h_data, h2_data) = match sy_mode {
        SyStorageMode::AsymmetricUpperNonzero => (
            // Asymmetric matrix: element (0, 1) = 0.5 != element (1, 0) = 0.2
            vec![1.0, 0.5, 0.2, 1.0],
            vec![-1.2, -0.3, -0.3, -1.2],
            vec![0.5, 0.1, 0.1, 0.5],
        ),
        SyStorageMode::LowerTriangleOnly => (
            // Lower triangle only: upper triangle element (0, 1) = 0.0, lower (1, 0) = 0.35
            vec![1.0, 0.0, 0.35, 1.0],
            vec![-1.2, 0.0, -0.35, -1.2],
            vec![0.5, 0.0, 0.25, 0.5],
        ),
        SyStorageMode::FullSymmetric => (
            // Fully populated symmetric matrix
            vec![1.0, 0.2, 0.2, 1.0],
            vec![-1.2, -0.3, -0.3, -1.2],
            vec![0.5, 0.1, 0.1, 0.5],
        ),
    };
    write_mat("s", "SY", s_data)?;
    write_mat("h", "SY", h_data)?;
    write_mat("h2", "SY", h2_data)?;

    // 6. Spin-specific orbital datasets
    let run_types = match method {
        "uks" => vec!["uks_alpha", "uks_beta"],
        "roks" => vec!["roks_close", "roks_open"],
        _ => vec!["rks"],
    };

    for (idx, rt) in run_types.iter().enumerate() {
        // Different values for each spin component (alpha vs beta, close vs open)
        let delta = idx as f64 * 0.1;

        // C matrix: path "C/{run_type}_2" (_get_h5path_c_matrix 659-662行目)
        let c_data = vec![1.0 - delta, 0.0 + delta, 0.0 + delta, 1.0 - delta];
        write_mat(&format!("C/{}_{}", rt, 2), "GE", c_data)?;

        // Energy levels: path "energy_level/{run_type}_2" (_get_h5path_energy_level 678-681行目)
        let oe_data = vec![-0.50 + delta, 0.20 + delta];
        write_vec(&format!("energy_level/{}_{}", rt, 2), oe_data)?;

        // Occupations: path "occ/{run_type}" (_get_h5path_occ 643-645行目)
        let occ_data = if idx == 0 {
            vec![1.0, 0.0]
        } else {
            vec![0.0, 0.0]
        };
        write_vec(&format!("occ/{}", rt), occ_data)?;

        // Mulliken charges: path "pop/mulliken/{run_type}_2" (_get_h5path_pop_mulliken_atom 683-686行目)
        let pop_data = vec![0.50 - delta, 0.50 + delta];
        write_vec(&format!("pop/mulliken/{}_{}", rt, 2), pop_data)?;

        // P and F matrices across iterations 1 and 2 (_get_h5path_density_matrix 664-667行目, _get_h5path_f_matrix 669-672行目)
        for itr in 1..=2 {
            let itr_offset = (itr - 1) as f64 * 0.05;
            let p_data = match sy_mode {
                SyStorageMode::LowerTriangleOnly => {
                    vec![
                        1.0 - delta + itr_offset,
                        0.0,
                        0.15 + itr_offset,
                        delta + itr_offset,
                    ]
                }
                _ => vec![1.0 - delta + itr_offset, 0.0, 0.0, delta + itr_offset],
            };
            write_mat(&format!("P/{}_{}", rt, itr), "SY", p_data)?;

            let f_data = match sy_mode {
                SyStorageMode::LowerTriangleOnly => {
                    vec![
                        -0.50 + delta + itr_offset,
                        0.0,
                        -0.25 + itr_offset,
                        0.20 + delta + itr_offset,
                    ]
                }
                _ => vec![
                    -0.50 + delta + itr_offset,
                    0.05 + itr_offset,
                    0.05 + itr_offset,
                    0.20 + delta + itr_offset,
                ],
            };
            write_mat(&format!("F/{}_{}", rt, itr), "SY", f_data)?;
        }
    }

    Ok(())
}

/// Tests UKS orbital data and optional matrices (`s`, `P`, `F`, `h`, `h2`) from a synthetic HDF5 file.
///
/// Verifies:
/// - UKS alpha and beta spin components have different coefficients, energies, occupations, charges.
/// - Top-level symmetric matrices `s`, `h`, `h2` are successfully read as `SymmetricMatrix`.
/// - Group-level symmetric matrices `P` and `F` are successfully read.
/// - Format follows ProteinDF `pdfparam_hdf5.py` (_get_h5path_* 643-686行目, _set_matrix 599-611行目).
#[test]
fn test_synthetic_uks_hdf5() {
    let tmp_dir = std::env::temp_dir().join(format!("test_qc_hdf5_uks_{}", std::process::id()));
    std::fs::create_dir_all(&tmp_dir).expect("failed to create temp dir");
    let h5_path = tmp_dir.join("synthetic_uks.h5");

    write_synthetic_qc_h5(&h5_path, "uks", SyStorageMode::FullSymmetric)
        .expect("failed to write synthetic UKS HDF5");

    let result = PdfParamHdf5::load(&h5_path).expect("failed to load synthetic UKS HDF5");

    // 1. Root parameters
    assert_eq!(result.params.method, "uks");
    assert_eq!(result.params.xc_functional, "B3LYP");
    assert_eq!(result.params.num_of_atoms, 2);
    assert_eq!(result.params.num_of_aos, 2);
    assert_eq!(result.params.num_of_mos, 2);
    assert_eq!(result.params.iterations, 2);
    assert!(result.params.scf_converged);

    // 2. Symmetric matrices (s, h, h2)
    let s = result
        .s_matrix
        .as_ref()
        .expect("s_matrix should be present");
    assert_eq!(s.dim(), 2);
    assert!((s.get(0, 0).unwrap() - 1.0).abs() < 1e-10);
    assert!((s.get(0, 1).unwrap() - 0.2).abs() < 1e-10);
    assert!((s.get(1, 0).unwrap() - 0.2).abs() < 1e-10);

    let h = result
        .h_matrix
        .as_ref()
        .expect("h_matrix should be present");
    assert_eq!(h.dim(), 2);
    assert!((h.get(0, 0).unwrap() - (-1.2)).abs() < 1e-10);
    assert!((h.get(0, 1).unwrap() - (-0.3)).abs() < 1e-10);

    let h2 = result
        .h2_matrix
        .as_ref()
        .expect("h2_matrix should be present");
    assert_eq!(h2.dim(), 2);
    assert!((h2.get(0, 0).unwrap() - 0.5).abs() < 1e-10);
    assert!((h2.get(0, 1).unwrap() - 0.1).abs() < 1e-10);

    // 3. Group-level matrices (P, F): verify retrieval methods, spin separation, and iteration tracking
    assert_eq!(result.density_matrices.len(), 4);
    assert_eq!(result.fock_matrices.len(), 4);

    // 3a. Density matrix: retrieval by run_type and iteration
    let p_alpha_1 = result
        .get_density_matrix(&RunType::UksAlpha, Some(1))
        .expect("missing UKS alpha P matrix at iteration 1");
    let p_alpha_2 = result
        .get_density_matrix(&RunType::UksAlpha, Some(2))
        .expect("missing UKS alpha P matrix at iteration 2");
    let p_alpha_latest = result
        .get_density_matrix(&RunType::UksAlpha, None)
        .expect("missing latest UKS alpha P matrix");
    assert_eq!(p_alpha_latest, p_alpha_2);

    let p_beta_1 = result
        .get_density_matrix(&RunType::UksBeta, Some(1))
        .expect("missing UKS beta P matrix at iteration 1");
    let p_beta_2 = result
        .get_density_matrix(&RunType::UksBeta, Some(2))
        .expect("missing UKS beta P matrix at iteration 2");
    let p_beta_latest = result
        .get_density_matrix(&RunType::UksBeta, None)
        .expect("missing latest UKS beta P matrix");
    assert_eq!(p_beta_latest, p_beta_2);

    // Alpha vs Beta must have distinct values (no confusion between spins)
    assert_ne!(p_alpha_2.get(0, 0).unwrap(), p_beta_2.get(0, 0).unwrap());
    assert!((p_alpha_2.get(0, 0).unwrap() - 1.05).abs() < 1e-10);
    assert!((p_beta_2.get(0, 0).unwrap() - 0.95).abs() < 1e-10);

    // Iteration 1 vs Iteration 2 must have distinct values (no confusion between iterations)
    assert_ne!(p_alpha_1.get(0, 0).unwrap(), p_alpha_2.get(0, 0).unwrap());
    assert_ne!(p_beta_1.get(0, 0).unwrap(), p_beta_2.get(0, 0).unwrap());
    assert!((p_alpha_1.get(0, 0).unwrap() - 1.00).abs() < 1e-10);
    assert!((p_beta_1.get(0, 0).unwrap() - 0.90).abs() < 1e-10);

    // Non-existent searches
    assert!(result.get_density_matrix(&RunType::Rks, None).is_none());
    assert!(result
        .get_density_matrix(&RunType::UksAlpha, Some(99))
        .is_none());

    // 3b. Fock matrix: retrieval by run_type and iteration
    let f_alpha_1 = result
        .get_fock_matrix(&RunType::UksAlpha, Some(1))
        .expect("missing UKS alpha F matrix at iteration 1");
    let f_alpha_2 = result
        .get_fock_matrix(&RunType::UksAlpha, Some(2))
        .expect("missing UKS alpha F matrix at iteration 2");
    let f_alpha_latest = result
        .get_fock_matrix(&RunType::UksAlpha, None)
        .expect("missing latest UKS alpha F matrix");
    assert_eq!(f_alpha_latest, f_alpha_2);

    let f_beta_1 = result
        .get_fock_matrix(&RunType::UksBeta, Some(1))
        .expect("missing UKS beta F matrix at iteration 1");
    let f_beta_2 = result
        .get_fock_matrix(&RunType::UksBeta, Some(2))
        .expect("missing UKS beta F matrix at iteration 2");
    let f_beta_latest = result
        .get_fock_matrix(&RunType::UksBeta, None)
        .expect("missing latest UKS beta F matrix");
    assert_eq!(f_beta_latest, f_beta_2);

    // Alpha vs Beta must have distinct values
    assert_ne!(f_alpha_2.get(0, 0).unwrap(), f_beta_2.get(0, 0).unwrap());
    // Iteration 1 vs Iteration 2 must have distinct values
    assert_ne!(f_alpha_1.get(0, 0).unwrap(), f_alpha_2.get(0, 0).unwrap());
    assert_ne!(f_beta_1.get(0, 0).unwrap(), f_beta_2.get(0, 0).unwrap());

    // Non-existent searches
    assert!(result.get_fock_matrix(&RunType::Rks, None).is_none());
    assert!(result
        .get_fock_matrix(&RunType::UksAlpha, Some(99))
        .is_none());

    // 4. UKS MO coefficients: alpha vs beta must have distinct values
    let c_alpha = result
        .get_mo_coefficients(&RunType::UksAlpha, Some(2))
        .expect("missing UKS alpha MO coefficients");
    let c_beta = result
        .get_mo_coefficients(&RunType::UksBeta, Some(2))
        .expect("missing UKS beta MO coefficients");
    assert_eq!(c_alpha.n_ao(), 2);
    assert_eq!(c_alpha.n_mo(), 2);
    assert_eq!(c_beta.n_ao(), 2);
    assert_eq!(c_beta.n_mo(), 2);

    // Alpha C[0, 0] = 1.0, Beta C[0, 0] = 0.9 (distinct!)
    assert!((c_alpha.get(0, 0).unwrap() - 1.0).abs() < 1e-10);
    assert!((c_beta.get(0, 0).unwrap() - 0.9).abs() < 1e-10);
    assert_ne!(c_alpha.get(0, 0).unwrap(), c_beta.get(0, 0).unwrap());

    // 5. UKS Orbital energies: alpha vs beta must have distinct values
    let oe_alpha = result
        .get_orbital_energy(&RunType::UksAlpha, Some(2))
        .expect("missing UKS alpha orbital energy");
    let oe_beta = result
        .get_orbital_energy(&RunType::UksBeta, Some(2))
        .expect("missing UKS beta orbital energy");
    assert!((oe_alpha.energies[0] - (-0.50)).abs() < 1e-10);
    assert!((oe_beta.energies[0] - (-0.40)).abs() < 1e-10);
    assert_ne!(oe_alpha.energies[0], oe_beta.energies[0]);

    // 6. UKS Orbital occupations: alpha vs beta must have distinct values
    let occ_alpha = result
        .get_orbital_occupation(&RunType::UksAlpha)
        .expect("missing UKS alpha occupation");
    let occ_beta = result
        .get_orbital_occupation(&RunType::UksBeta)
        .expect("missing UKS beta occupation");
    assert_eq!(occ_alpha.occupations[0], 1.0);
    assert_eq!(occ_beta.occupations[0], 0.0);
    assert_ne!(occ_alpha.occupations[0], occ_beta.occupations[0]);

    // 7. UKS Mulliken populations: alpha vs beta must have distinct values
    let pop_alpha = result
        .get_mulliken_population(&RunType::UksAlpha, Some(2))
        .expect("missing UKS alpha Mulliken population");
    let pop_beta = result
        .get_mulliken_population(&RunType::UksBeta, Some(2))
        .expect("missing UKS beta Mulliken population");
    assert!((pop_alpha.charges[0] - 0.50).abs() < 1e-10);
    assert!((pop_beta.charges[0] - 0.40).abs() < 1e-10);
    assert_ne!(pop_alpha.charges[0], pop_beta.charges[0]);

    // Cleanup
    let _ = std::fs::remove_dir_all(&tmp_dir);
}

/// Tests ROKS orbital data reading from a synthetic HDF5 file.
///
/// Verifies:
/// - ROKS closed-shell (`roks_close`) and open-shell (`roks_open`) run types are supported.
/// - MO coefficients, energies, occupations, and Mulliken charges are correctly read.
/// - Format follows ProteinDF `pdfparam_hdf5.py` (_get_h5path_* 643-686行目).
#[test]
fn test_synthetic_roks_hdf5() {
    let tmp_dir = std::env::temp_dir().join(format!("test_qc_hdf5_roks_{}", std::process::id()));
    std::fs::create_dir_all(&tmp_dir).expect("failed to create temp dir");
    let h5_path = tmp_dir.join("synthetic_roks.h5");

    write_synthetic_qc_h5(&h5_path, "roks", SyStorageMode::FullSymmetric)
        .expect("failed to write synthetic ROKS HDF5");

    let result = PdfParamHdf5::load(&h5_path).expect("failed to load synthetic ROKS HDF5");

    assert_eq!(result.params.method, "roks");

    // 1. ROKS MO coefficients (roks_close and roks_open)
    let c_close = result
        .get_mo_coefficients(&RunType::RoksClose, Some(2))
        .expect("missing ROKS close MO coefficients");
    let c_open = result
        .get_mo_coefficients(&RunType::RoksOpen, Some(2))
        .expect("missing ROKS open MO coefficients");
    assert_eq!(c_close.n_ao(), 2);
    assert_eq!(c_open.n_ao(), 2);
    assert!((c_close.get(0, 0).unwrap() - 1.0).abs() < 1e-10);
    assert!((c_open.get(0, 0).unwrap() - 0.9).abs() < 1e-10);

    // 2. ROKS orbital energies
    let oe_close = result
        .get_orbital_energy(&RunType::RoksClose, Some(2))
        .expect("missing ROKS close orbital energy");
    let oe_open = result
        .get_orbital_energy(&RunType::RoksOpen, Some(2))
        .expect("missing ROKS open orbital energy");
    assert!((oe_close.energies[0] - (-0.50)).abs() < 1e-10);
    assert!((oe_open.energies[0] - (-0.40)).abs() < 1e-10);

    // 3. ROKS occupations
    let occ_close = result
        .get_orbital_occupation(&RunType::RoksClose)
        .expect("missing ROKS close occupation");
    let occ_open = result
        .get_orbital_occupation(&RunType::RoksOpen)
        .expect("missing ROKS open occupation");
    assert_eq!(occ_close.occupations[0], 1.0);
    assert_eq!(occ_open.occupations[0], 0.0);

    // 4. ROKS Mulliken population
    let pop_close = result
        .get_mulliken_population(&RunType::RoksClose, Some(2))
        .expect("missing ROKS close Mulliken population");
    let pop_open = result
        .get_mulliken_population(&RunType::RoksOpen, Some(2))
        .expect("missing ROKS open Mulliken population");
    assert!((pop_close.charges[0] - 0.50).abs() < 1e-10);
    assert!((pop_open.charges[0] - 0.40).abs() < 1e-10);

    // Cleanup
    let _ = std::fs::remove_dir_all(&tmp_dir);
}

/// Tests that reading a non-symmetric matrix in HDF5 (`s` dataset) triggers a `BridgeError`.
///
/// Verifies the symmetry validation decision from PR#60 / PR#61:
/// When the upper triangle contains non-zero elements, each element must match its
/// lower-triangle counterpart within [`MATRIX_SYMMETRY_TOLERANCE`] ($10^{-8}$).
/// If an asymmetric non-zero element is found, loading fails with a [`BridgeError`].
#[test]
fn test_synthetic_asymmetric_s_matrix_error() {
    let tmp_dir = std::env::temp_dir().join(format!("test_qc_hdf5_asymm_{}", std::process::id()));
    std::fs::create_dir_all(&tmp_dir).expect("failed to create temp dir");
    let h5_path = tmp_dir.join("synthetic_asymm.h5");

    // Write with AsymmetricUpperNonzero (element (0, 1) = 0.5, element (1, 0) = 0.2, diff = 0.3)
    write_synthetic_qc_h5(&h5_path, "rks", SyStorageMode::AsymmetricUpperNonzero)
        .expect("failed to write synthetic HDF5 with asymmetric S");

    let err = PdfParamHdf5::load(&h5_path).expect_err("loading asymmetric matrix must fail");
    let err_msg = err.to_string();
    assert!(
        err_msg.contains("not symmetric"),
        "error message must mention symmetry failure: {err_msg}"
    );
    assert!(
        err_msg.contains("matrix 's'"),
        "error message must specify matrix 's': {err_msg}"
    );

    // Cleanup
    let _ = std::fs::remove_dir_all(&tmp_dir);
}

/// Tests reading symmetric matrices (`s`, `P`, `F`, `h`, `h2`) stored with lower-triangle
/// elements only (upper triangle all 0.0), matching `pdf-archive-h5.py` output.
///
/// Verifies:
/// - Upper triangle containing all 0.0 is accepted without error.
/// - Lower triangle values are mirrored / folded so that symmetric access
///   (`matrix.get(0, 1) == matrix.get(1, 0)`) returns the correct lower-triangle value.
#[test]
fn test_synthetic_lower_triangle_sy_matrix() {
    let tmp_dir =
        std::env::temp_dir().join(format!("test_qc_hdf5_lower_tri_{}", std::process::id()));
    std::fs::create_dir_all(&tmp_dir).expect("failed to create temp dir");
    let h5_path = tmp_dir.join("synthetic_lower_tri.h5");

    // Write with LowerTriangleOnly:
    // s: (0, 0)=1.0, (1, 0)=0.35, (0, 1)=0.0, (1, 1)=1.0
    // h: (0, 0)=-1.2, (1, 0)=-0.35, (0, 1)=0.0, (1, 1)=-1.2
    // h2: (0, 0)=0.5, (1, 0)=0.25, (0, 1)=0.0, (1, 1)=0.5
    // P: (0, 0)=1.0, (1, 0)=0.15, (0, 1)=0.0, (1, 1)=0.0
    // F: (0, 0)=-0.50, (1, 0)=-0.25, (0, 1)=0.0, (1, 1)=0.20
    write_synthetic_qc_h5(&h5_path, "rks", SyStorageMode::LowerTriangleOnly)
        .expect("failed to write synthetic lower-triangle HDF5");

    let result =
        PdfParamHdf5::load(&h5_path).expect("failed to load synthetic lower-triangle HDF5");

    // 1. s matrix: verify folded value (0, 1) equals (1, 0) == 0.35
    let s = result.s_matrix.as_ref().expect("missing s matrix");
    assert_eq!(s.dim(), 2);
    assert!((s.get(0, 0).unwrap() - 1.0).abs() < 1e-10);
    assert!((s.get(1, 1).unwrap() - 1.0).abs() < 1e-10);
    assert!((s.get(1, 0).unwrap() - 0.35).abs() < 1e-10);
    assert!((s.get(0, 1).unwrap() - 0.35).abs() < 1e-10);

    // 2. h matrix: verify folded value
    let h = result.h_matrix.as_ref().expect("missing h matrix");
    assert!((h.get(1, 0).unwrap() - (-0.35)).abs() < 1e-10);
    assert!((h.get(0, 1).unwrap() - (-0.35)).abs() < 1e-10);

    // 3. h2 matrix: verify folded value
    let h2 = result.h2_matrix.as_ref().expect("missing h2 matrix");
    assert!((h2.get(1, 0).unwrap() - 0.25).abs() < 1e-10);
    assert!((h2.get(0, 1).unwrap() - 0.25).abs() < 1e-10);

    // 4. P matrix: verify folded value
    assert_eq!(result.density_matrices.len(), 2);
    let p = result
        .get_density_matrix(&RunType::Rks, Some(1))
        .expect("missing RKS density matrix at iteration 1");
    assert!((p.get(1, 0).unwrap() - 0.15).abs() < 1e-10);
    assert!((p.get(0, 1).unwrap() - 0.15).abs() < 1e-10);

    // 5. F matrix: verify folded value
    assert_eq!(result.fock_matrices.len(), 2);
    let f = result
        .get_fock_matrix(&RunType::Rks, Some(1))
        .expect("missing RKS Fock matrix at iteration 1");
    assert!((f.get(1, 0).unwrap() - (-0.25)).abs() < 1e-10);
    assert!((f.get(0, 1).unwrap() - (-0.25)).abs() < 1e-10);

    // Cleanup
    let _ = std::fs::remove_dir_all(&tmp_dir);
}

/// Tests that atoms without basis set definitions (dummy atoms like 'X')
/// have 0 atomic orbitals and do not disrupt the atom index mapping of subsequent atoms.
///
/// Verifies:
/// - Synthetic molecule: Atom 0 = C (has basis set), Atom 1 = X (dummy, no basis set), Atom 2 = H (has basis set).
/// - Atom 1 ('X') contributes 0 AOs.
/// - Subsequent atom 2 ('H') correctly maps its AO to `atom_index = 2`, not shifted to 1.
#[test]
fn test_synthetic_dummy_atom_orbital_mapping() {
    let tmp_dir =
        std::env::temp_dir().join(format!("test_qc_hdf5_dummy_atom_{}", std::process::id()));
    std::fs::create_dir_all(&tmp_dir).expect("failed to create temp dir");
    let h5_path = tmp_dir.join("synthetic_dummy_atom.h5");

    let file = hdf5_metno::File::create(&h5_path).expect("failed to create HDF5 file");

    // 1. Root attributes
    let write_str_attr = |name: &str, val: &str| {
        file.new_attr::<FixedAscii<32>>()
            .create(name)
            .unwrap()
            .write_scalar(&str_ascii::<32>(val))
            .unwrap();
    };
    let write_int_attr = |name: &str, val: i64| {
        file.new_attr::<i64>()
            .create(name)
            .unwrap()
            .write_scalar(&val)
            .unwrap();
    };
    write_str_attr("method", "rks");
    write_str_attr("xc_functional", "B3LYP");
    write_int_attr("num_of_atoms", 3);
    write_int_attr("num_of_AOs", 2);
    write_int_attr("num_of_MOs", 2);
    write_int_attr("iterations", 1);
    file.new_attr::<bool>()
        .create("scf_converged")
        .unwrap()
        .write_scalar(&true)
        .unwrap();

    // 2. Molecule with 3 atoms: C (index 0), X (index 1, dummy), H (index 2)
    let xyzc_data = vec![
        str_ascii::<32>("C"),
        str_ascii::<32>(""),
        str_ascii::<32>("0.0"),
        str_ascii::<32>("0.0"),
        str_ascii::<32>("0.0"),
        str_ascii::<32>("6.0"),
        str_ascii::<32>("X"),
        str_ascii::<32>(""),
        str_ascii::<32>("1.0"),
        str_ascii::<32>("0.0"),
        str_ascii::<32>("0.0"),
        str_ascii::<32>("1.0"),
        str_ascii::<32>("H"),
        str_ascii::<32>(""),
        str_ascii::<32>("2.0"),
        str_ascii::<32>("0.0"),
        str_ascii::<32>("0.0"),
        str_ascii::<32>("1.0"),
    ];
    let xyzc_arr = ndarray::Array2::from_shape_vec((3, 6), xyzc_data).unwrap();
    file.new_dataset_builder()
        .with_data(&xyzc_arr)
        .create("molecule/xyzc")
        .unwrap();

    // 3. Basis set definitions: only C and H have basis sets. 'X' is omitted.
    let name_data = vec![
        str_ascii::<32>("C-basis"),
        str_ascii::<32>("C"),
        str_ascii::<32>("H-basis"),
        str_ascii::<32>("H"),
    ];
    let name_arr = ndarray::Array2::from_shape_vec((2, 2), name_data).unwrap();
    file.new_dataset_builder()
        .with_data(&name_arr)
        .create("basisset/name")
        .unwrap();

    let cgto_data = vec![
        str_ascii::<32>("C-basis"),
        str_ascii::<32>("0"),
        str_ascii::<32>("s"),
        str_ascii::<32>("1.0"),
        str_ascii::<32>("H-basis"),
        str_ascii::<32>("0"),
        str_ascii::<32>("s"),
        str_ascii::<32>("1.0"),
    ];
    let cgto_arr = ndarray::Array2::from_shape_vec((2, 4), cgto_data).unwrap();
    file.new_dataset_builder()
        .with_data(&cgto_arr)
        .create("basisset/cgto")
        .unwrap();

    let pgto_data = vec![
        str_ascii::<32>("C-basis"),
        str_ascii::<32>("0"),
        str_ascii::<32>("0"),
        str_ascii::<32>("1.0"),
        str_ascii::<32>("1.0"),
        str_ascii::<32>("H-basis"),
        str_ascii::<32>("0"),
        str_ascii::<32>("0"),
        str_ascii::<32>("1.0"),
        str_ascii::<32>("1.0"),
    ];
    let pgto_arr = ndarray::Array2::from_shape_vec((2, 5), pgto_data).unwrap();
    file.new_dataset_builder()
        .with_data(&pgto_arr)
        .create("basisset/pgto")
        .unwrap();

    drop(file);

    let result = PdfParamHdf5::load(&h5_path).expect("failed to load synthetic dummy atom HDF5");

    assert_eq!(result.molecule.atoms.len(), 3);
    assert_eq!(result.molecule.atoms[0].symbol, "C");
    assert_eq!(result.molecule.atoms[1].symbol, "X");
    assert_eq!(result.molecule.atoms[2].symbol, "H");

    let orb_info = &result.orbital_info;
    assert_eq!(orb_info.num_of_aos, 2);
    assert_eq!(orb_info.orbitals.len(), 2);

    // AO 0 corresponds to Atom 0 (C)
    assert_eq!(orb_info.orbitals[0].ao_index, 0);
    assert_eq!(orb_info.orbitals[0].atom_index, 0);
    assert_eq!(orb_info.orbitals[0].atom_symbol, "C");

    // AO 1 corresponds to Atom 2 (H). It must NOT be shifted to atom_index 1!
    assert_eq!(orb_info.orbitals[1].ao_index, 1);
    assert_eq!(orb_info.orbitals[1].atom_index, 2);
    assert_eq!(orb_info.orbitals[1].atom_symbol, "H");

    let _ = std::fs::remove_dir_all(&tmp_dir);
}

/// Helper function to verify mathematical and physical consistency of real QC calculation data:
/// 1. Structural dimensions (num_of_atoms, num_of_aos, num_of_mos, iterations)
/// 2. Final SCF energy matches reference within 1e-6 Hartree
/// 3. MO coefficients satisfy generalized orthogonality condition: $C^T S C = I$ within 1e-8
/// 4. Density matrix matches orbital-constructed density $P = \sum_i n_i C_i C_i^T$ within 1e-8
/// 5. $\text{tr}(PS) = N_{\text{electrons}}$ within 1e-8
#[allow(clippy::needless_range_loop, clippy::too_many_arguments)]
fn verify_real_qc_result(
    result: &proteindf_bridge::qc_result::QcResult,
    expected_atoms: usize,
    expected_aos: usize,
    expected_mos: usize,
    expected_iterations: usize,
    expected_final_energy: f64,
    expected_electrons: f64,
    expected_mulliken_charges: &[f64],
) {
    // 1. Structural dimensions
    assert_eq!(result.params.num_of_atoms, expected_atoms);
    assert_eq!(result.params.num_of_aos, expected_aos);
    assert_eq!(result.params.num_of_mos, expected_mos);
    assert_eq!(result.params.iterations, expected_iterations);
    assert!(result.params.scf_converged);

    // 2. Presence of required matrices
    let s = result.s_matrix.as_ref().expect("s matrix must be present");
    assert_eq!(s.dim(), expected_aos);
    assert!(
        !result.density_matrices.is_empty(),
        "density_matrices must not be empty"
    );
    assert!(
        !result.fock_matrices.is_empty(),
        "fock_matrices must not be empty"
    );
    assert!(result.h_matrix.is_some(), "h matrix must be present");

    // 3. Final total energy
    let final_energy = *result
        .total_energies
        .last()
        .expect("total_energies must not be empty");
    let energy_diff = (final_energy - expected_final_energy).abs();
    assert!(
        energy_diff < 1e-6,
        "final energy diff {energy_diff} exceeds 1e-6 (got {final_energy}, expected {expected_final_energy})"
    );

    // 4. C^T S C == I (orthonormality of MOs)
    let c = result
        .get_mo_coefficients(&RunType::Rks, Some(expected_iterations))
        .expect("missing final MO coefficients");
    assert_eq!(c.n_ao(), expected_aos);
    assert_eq!(c.n_mo(), expected_mos);

    // Compute SC = S * C (shape: n_ao x n_mo)
    let mut sc = vec![vec![0.0; expected_mos]; expected_aos];
    for p in 0..expected_aos {
        for i in 0..expected_mos {
            let mut sum = 0.0;
            for q in 0..expected_aos {
                sum += s.get(p, q).unwrap() * c.get(q, i).unwrap();
            }
            sc[p][i] = sum;
        }
    }

    // Compute C^T * (SC) (shape: n_mo x n_mo) and check |(C^T S C)_ij - delta_ij| < 1e-8
    for i in 0..expected_mos {
        for j in 0..expected_mos {
            let mut ct_sc_ij = 0.0;
            for p in 0..expected_aos {
                ct_sc_ij += c.get(p, i).unwrap() * sc[p][j];
            }
            let expected_val = if i == j { 1.0 } else { 0.0 };
            let diff = (ct_sc_ij - expected_val).abs();
            assert!(
                diff < 1e-8,
                "C^T S C element ({i}, {j}) diff {diff} exceeds 1e-8 (got {ct_sc_ij}, expected {expected_val})"
            );
        }
    }

    // 5. P_file == P_calc from C and occ: P_pq = sum_i occ_i * C_pi * C_qi
    let occ = result
        .get_orbital_occupation(&RunType::Rks)
        .expect("missing orbital occupations");
    assert_eq!(occ.occupations.len(), expected_mos);

    let p_entry = result
        .get_density_matrix(&RunType::Rks, Some(expected_iterations))
        .expect("missing final density matrix");
    let p_file = &p_entry.matrix;
    assert_eq!(p_file.dim(), expected_aos);

    for p in 0..expected_aos {
        for q in 0..=p {
            let mut p_calc_pq = 0.0;
            for i in 0..expected_mos {
                p_calc_pq += occ.occupations[i] * c.get(p, i).unwrap() * c.get(q, i).unwrap();
            }
            let p_file_pq = p_file.get(p, q).unwrap();
            let diff = (p_file_pq - p_calc_pq).abs();
            assert!(
                diff < 1e-8,
                "P matrix element ({p}, {q}) diff {diff} exceeds 1e-8 (file: {p_file_pq}, calc: {p_calc_pq})"
            );
        }
    }

    // 6. tr(PS) == N_electrons: tr(PS) = sum_p sum_q P_pq * S_qp
    let mut tr_ps = 0.0;
    for p in 0..expected_aos {
        for q in 0..expected_aos {
            tr_ps += p_file.get(p, q).unwrap() * s.get(q, p).unwrap();
        }
    }
    let elec_diff = (tr_ps - expected_electrons).abs();
    assert!(
        elec_diff < 1e-8,
        "tr(PS) diff {elec_diff} exceeds 1e-8 (got {tr_ps}, expected {expected_electrons})"
    );

    // 7. Orbital-resolved Mulliken atomic population verification:
    // (1) For every MO, sum of atomic contributions across all atoms equals 1.0 within 1e-10.
    // (2) For every atom A, net charge Z_A - sum_i (occ_i * q_A^{(i)}) matches ProteinDF_test
    //     pdfresults.std.h5 dataset pop/mulliken/rks_{final_iteration} within 1e-6.
    let orb_contribs =
        calc_orbital_mulliken_contributions(c, s, &result.orbital_info, expected_atoms)
            .expect("failed to calculate orbital Mulliken contributions");

    assert_eq!(orb_contribs.len(), expected_mos);

    // (1) Check MO normalization: sum across atoms == 1.0 within 1e-10 for every MO
    for (i, contrib) in orb_contribs.iter().enumerate() {
        assert_eq!(contrib.orbital_index, i);
        assert_eq!(contrib.atom_contributions.len(), expected_atoms);
        let sum_atoms: f64 = contrib.atom_contributions.iter().map(|(_, val)| *val).sum();
        let norm_diff = (sum_atoms - 1.0).abs();
        assert!(
            norm_diff < 1e-10,
            "MO {i} atom contributions sum diff {norm_diff} exceeds 1e-10 (got {sum_atoms}, expected 1.0)"
        );
    }

    // (2) Check net Mulliken charge consistency against expected values:
    assert_eq!(
        expected_mulliken_charges.len(),
        expected_atoms,
        "expected_mulliken_charges length must match expected_atoms"
    );
    assert_eq!(
        result.molecule.atoms.len(),
        expected_atoms,
        "molecule atoms count must match expected_atoms"
    );

    for (atom_idx, expected_charge) in expected_mulliken_charges.iter().enumerate() {
        let z_atom = result.molecule.atoms[atom_idx].charge;
        let mut gross_pop = 0.0;
        for i in 0..expected_mos {
            let occ_i = occ.occupations[i];
            let q_a_i = orb_contribs[i].atom_contributions[atom_idx].1;
            gross_pop += occ_i * q_a_i;
        }
        let net_charge = z_atom - gross_pop;
        let charge_diff = (net_charge - expected_charge).abs();
        assert!(
            charge_diff < 1e-6,
            "Atom {atom_idx} net Mulliken charge diff {charge_diff} exceeds 1e-6 (got {net_charge}, expected {expected_charge})"
        );
    }
}

/// Tests real calculation result `CH3COOX_HF.h5` regenerated by ProteinDF.
///
/// Marked `#[ignore]`. Runs when `--run-ignored` is specified.
/// Reads from directory specified by environment variable `PROTEINDF_BRIDGE_QC_DATA`.
/// Fails with a panic if the environment variable is not set or the file is missing.
///
/// # Determination of Expected Values and Verification Process
/// - System: CH3COOX_HF (acetic acid with dummy charge X).
/// - Dimensions:
///   - 8 atoms (including dummy atom X), 42 AOs, 42 MOs, 88 iterations.
///   - Verified by inspecting root attributes `/num_of_atoms`, `/num_of_AOs`, `/num_of_MOs`,
///     and `/iterations` in both `ProteinDF_test/serial_dev/CH3COOX_HF.mmap/pdfresults.std.h5`
///     and the regenerated file using h5py.
/// - Final total energy:
///   - ProteinDF_test baseline (`pdfresults.std.h5` dataset `/TEs`[-1]): `-227.32196883081747` Hartree.
///   - Regenerated file via ProteinDF (`bf581a8` build, dataset `/TEs`[-1]): `-227.32196883082707` Hartree.
///   - The expected value in this test is set to the ProteinDF_test reference value (`-227.32196883081747`),
///     verifying that the regenerated calculation matches it within `1e-6` Hartree (actual difference: 9.6e-12).
/// - Overlap matrix $S$ and MO orthonormality $C^T S C = I$:
///   - Read dataset `/s` (shape 42x42, type SY, lower-triangle storage) and unfolded symmetric entries.
///   - Read dataset `/C/rks_88` (shape 42x42, type GE).
///   - Evaluated matrix product $C^T S C$ via numpy; confirmed max absolute deviation from identity matrix $< 1e-8$.
/// - Density matrix $P$ consistency and total electrons:
///   - Read dataset `/occ/rks` (occupation vector, 42 elements; 16 orbitals doubly occupied, 26 unoccupied).
///   - Constructed density matrix $P_{\text{calc}} = \sum_i \mathrm{occ}_i C_i C_i^T$.
///   - Verified that $P_{\text{calc}}$ matches dataset `/P/rks_88` (lower-triangle storage) within $1e-8$.
///   - Computed $\mathrm{tr}(PS) = \sum_{pq} P_{pq} S_{qp} = 32.0$ electrons (neutral acetic acid $CH_3COOH$
///     has 32 electrons; dummy atom X has no basis functions or electrons), verified within $1e-8$.
/// - Orbital-resolved Mulliken atomic population and gross population charges:
///   - Calculated via `calc_orbital_mulliken_contributions` using $C$, $S$, and OrbitalInfo.
///   - Verified (1) sum of atomic contributions across atoms equals 1.0 within 1e-10 for all 42 MOs.
///   - Verified (2) atomic net Mulliken charge $Z_A - \sum_i \mathrm{occ}_i q_A^{(i)}$ matches baseline
///     `ProteinDF_test/serial_dev/CH3COOX_HF.mmap/pdfresults.std.h5` dataset `/pop/mulliken/rks_88`
///     within 1e-6 (actual max difference: < 3e-13).
///   - Baseline charges read from `/pop/mulliken/rks_88` via h5py:
///     `[-0.44462725, 0.21815275, 0.17203644, 0.17315344, 0.57957603, -0.80335457, -0.89493684, 1.0]`
///     (Full precision values: `[-0.4446272473755002, 0.21815275205681606, 0.17203643804703905,
///       0.1731534405514823, 0.5795760274037338, -0.8033545702707343, -0.894936840412818, 1.0]`).
#[test]
#[ignore]
fn test_real_qc_data_ch3coox_hf() {
    let dir_str = std::env::var("PROTEINDF_BRIDGE_QC_DATA").expect(
        "PROTEINDF_BRIDGE_QC_DATA environment variable not set. \
         Must point to directory containing regenerated CH3COOX_HF.h5 and GLY.h5",
    );
    let dir = std::path::PathBuf::from(dir_str);
    assert!(
        dir.is_dir(),
        "PROTEINDF_BRIDGE_QC_DATA is not a directory: {}",
        dir.display()
    );
    let file_path = dir.join("CH3COOX_HF.h5");
    assert!(
        file_path.exists(),
        "CH3COOX_HF.h5 not found in: {}",
        file_path.display()
    );

    let result = PdfParamHdf5::load(&file_path).expect("failed to load CH3COOX_HF.h5");
    let expected_mulliken_charges = [
        -0.444_627_247_375_500_2,
        0.218_152_752_056_816_06,
        0.172_036_438_047_039_05,
        0.173_153_440_551_482_3,
        0.579_576_027_403_733_8,
        -0.803_354_570_270_734_3,
        -0.894_936_840_412_818,
        1.0,
    ];
    verify_real_qc_result(
        &result,
        8,                       // atoms (including dummy atom X)
        42,                      // AOs
        42,                      // MOs
        88,                      // iterations
        -227.321_968_830_817_47, // final energy from ProteinDF_test pdfresults.std.h5
        32.0,                    // electrons
        &expected_mulliken_charges,
    );
}

/// Tests real calculation result `GLY.h5` regenerated by ProteinDF.
///
/// Marked `#[ignore]`. Runs when `--run-ignored` is specified.
/// Reads from directory specified by environment variable `PROTEINDF_BRIDGE_QC_DATA`.
/// Fails with a panic if the environment variable is not set or the file is missing.
///
/// # Determination of Expected Values and Verification Process
/// - System: GLY (glycine molecule).
/// - Dimensions:
///   - 10 atoms, 95 AOs, 95 MOs, 49 iterations.
///   - Verified by inspecting root attributes `/num_of_atoms`, `/num_of_AOs`, `/num_of_MOs`,
///     and `/iterations` in both `ProteinDF_test/parallel/GLY.omp/pdfresults.std.h5`
///     and the regenerated file using h5py.
/// - Final total energy:
///   - ProteinDF_test baseline (`pdfresults.std.h5` dataset `/TEs`[-1]): `-282.86416358371673` Hartree.
///   - Regenerated file via ProteinDF (`bf581a8` build, dataset `/TEs`[-1]): `-282.86416358371594` Hartree.
///   - The expected value in this test is set to the ProteinDF_test reference value (`-282.86416358371673`),
///     verifying that the regenerated calculation matches it within `1e-6` Hartree (actual difference: 7.9e-13).
/// - Overlap matrix $S$ and MO orthonormality $C^T S C = I$:
///   - Read dataset `/s` (shape 95x95, type SY, lower-triangle storage) and unfolded symmetric entries.
///   - Read dataset `/C/rks_49` (shape 95x95, type GE).
///   - Evaluated matrix product $C^T S C$ via numpy; confirmed max absolute deviation from identity matrix $< 1e-8$.
/// - Density matrix $P$ consistency and total electrons:
///   - Read dataset `/occ/rks` (occupation vector, 95 elements; 20 orbitals doubly occupied, 75 unoccupied).
///   - Constructed density matrix $P_{\text{calc}} = \sum_i \mathrm{occ}_i C_i C_i^T$.
///   - Verified that $P_{\text{calc}}$ matches dataset `/P/rks_49` (lower-triangle storage) within $1e-8$.
/// - Orbital-resolved Mulliken atomic population and gross population charges:
///   - Calculated via `calc_orbital_mulliken_contributions` using $C$, $S$, and OrbitalInfo.
///   - Verified (1) sum of atomic contributions across atoms equals 1.0 within 1e-10 for all 95 MOs.
///   - Verified (2) atomic net Mulliken charge $Z_A - \sum_i \mathrm{occ}_i q_A^{(i)}$ matches baseline
///     `ProteinDF_test/parallel/GLY.omp/pdfresults.std.h5` dataset `/pop/mulliken/rks_49`
///     within 1e-6 (actual max difference: < 2e-13).
///   - Baseline charges read from `/pop/mulliken/rks_49` via h5py:
///     `[-0.45308302, 0.23025299, 0.19545874, -0.07892775, 0.11297073, 0.07262471, 0.47490416, -0.46885192, -0.3762826, 0.29093396]`
///     (Full precision values: `[-0.45308301703768894, 0.23025299491437712, 0.19545873920798773,
///       -0.07892775457186652, 0.11297072892275017, 0.07262470962875667, 0.4749041564437064,
///       -0.4688519230365511, -0.3762825986824083, 0.2909339642109071]`).
#[test]
#[ignore]
fn test_real_qc_data_gly() {
    let dir_str = std::env::var("PROTEINDF_BRIDGE_QC_DATA").expect(
        "PROTEINDF_BRIDGE_QC_DATA environment variable not set. \
         Must point to directory containing regenerated CH3COOX_HF.h5 and GLY.h5",
    );
    let dir = std::path::PathBuf::from(dir_str);
    assert!(
        dir.is_dir(),
        "PROTEINDF_BRIDGE_QC_DATA is not a directory: {}",
        dir.display()
    );
    let file_path = dir.join("GLY.h5");
    assert!(
        file_path.exists(),
        "GLY.h5 not found in: {}",
        file_path.display()
    );

    let result = PdfParamHdf5::load(&file_path).expect("failed to load GLY.h5");
    let expected_mulliken_charges = [
        -0.453_083_017_037_688_94,
        0.230_252_994_914_377_12,
        0.195_458_739_207_987_73,
        -0.078_927_754_571_866_52,
        0.112_970_728_922_750_17,
        0.072_624_709_628_756_67,
        0.474_904_156_443_706_4,
        -0.468_851_923_036_551_1,
        -0.376_282_598_682_408_3,
        0.290_933_964_210_907_1,
    ];
    verify_real_qc_result(
        &result,
        10,                      // atoms
        95,                      // AOs
        95,                      // MOs
        49,                      // iterations
        -282.864_163_583_716_73, // final energy from ProteinDF_test pdfresults.std.h5
        40.0,                    // electrons
        &expected_mulliken_charges,
    );
}
