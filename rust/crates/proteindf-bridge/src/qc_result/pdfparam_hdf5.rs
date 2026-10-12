// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

//! ProteinDF HDF5 (`reference.h5` / `pdfresults.h5`) loader.

use std::path::Path;

use hdf5_metno::types::{FixedAscii, VarLenAscii, VarLenUnicode};
use hdf5_metno::File;


use crate::error::{BridgeError, Result};
use crate::position::Position;
use crate::qc_result::types::*;

/// HDF5 reader for ProteinDF calculation results.
pub struct PdfParamHdf5;

impl PdfParamHdf5 {
    /// Loads a [`QcResult`] from the given HDF5 file path.
    pub fn load(path: impl AsRef<Path>) -> Result<QcResult> {
        let p = path.as_ref();
        let file = File::open(p).map_err(|e| {
            BridgeError::input_error(
                p.display().to_string(),
                format!("failed to open HDF5 file: {}", e),
            )
        })?;

        Self::read_from_file(&file, p.display().to_string())
    }

    /// Internal function reading from an open HDF5 [`File`].
    fn read_from_file(file: &File, source_name: String) -> Result<QcResult> {
        let params = Self::read_parameters(file, &source_name)?;
        let molecule = Self::read_molecule(file, &source_name)?;
        let basisset = Self::read_basisset(file, &source_name)?;
        let orbital_info = OrbitalInfo::build(&molecule, &basisset)?;

        let total_energies = Self::read_total_energies(file).unwrap_or_default();
        let orbital_energies = Self::read_orbital_energies(file)?;
        let orbital_occupations = Self::read_orbital_occupations(file)?;
        let mo_coefficients = Self::read_mo_coefficients(file)?;
        let s_matrix = Self::read_symmetric_matrix(file, "s")?;
        let density_matrices = Self::read_group_matrices(file, "P")?;
        let fock_matrices = Self::read_group_matrices(file, "F")?;
        let h_matrix = Self::read_symmetric_matrix(file, "h")?;
        let h2_matrix = Self::read_symmetric_matrix(file, "h2")?;
        let mulliken_charges = Self::read_mulliken(file)?;

        Ok(QcResult {
            params,
            molecule,
            basisset,
            orbital_info,
            total_energies,
            orbital_energies,
            orbital_occupations,
            mo_coefficients,
            s_matrix,
            density_matrices,
            fock_matrices,
            h_matrix,
            h2_matrix,
            mulliken_charges,
        })
    }

    /// Reads root-level calculation parameters.
    fn read_parameters(file: &File, source: &str) -> Result<QcParameters> {
        let method = Self::read_str_attr(file, "method")?.unwrap_or_default();
        let xc_functional = Self::read_str_attr(file, "xc_functional")?.unwrap_or_default();
        let num_of_atoms = Self::read_int_attr(file, "num_of_atoms")?
            .ok_or_else(|| BridgeError::input_error(source, "missing root attribute 'num_of_atoms'"))?;
        let num_of_aos = Self::read_int_attr(file, "num_of_AOs")?
            .ok_or_else(|| BridgeError::input_error(source, "missing root attribute 'num_of_AOs'"))?;
        let num_of_mos = Self::read_int_attr(file, "num_of_MOs")?
            .ok_or_else(|| BridgeError::input_error(source, "missing root attribute 'num_of_MOs'"))?;
        let iterations = Self::read_int_attr(file, "iterations")?.unwrap_or(0);
        let max_iterations = Self::read_int_attr(file, "max_iterations")?.unwrap_or(0);
        let scf_converged = Self::read_bool_attr(file, "scf_converged")?.unwrap_or(false);
        let step_control = Self::read_str_attr(file, "step_control")?.unwrap_or_default();
        let guess = Self::read_str_attr(file, "guess")?.unwrap_or_default();
        let cut_value = Self::read_str_attr(file, "cut_value")?.unwrap_or_default();
        let cdam_tau = Self::read_str_attr(file, "CDAM_tau")?.unwrap_or_default();
        let cd_epsilon = Self::read_str_attr(file, "CD_epsilon")?.unwrap_or_default();
        let scf_acceleration = Self::read_str_attr(file, "scf_acceleration")?.unwrap_or_default();
        let damping_factor = Self::read_f64_attr(file, "scf_acceleration_damping_damping_factor")?;
        let comment = Self::read_str_attr(file, "comment")?.unwrap_or_default();

        Ok(QcParameters {
            method,
            xc_functional,
            num_of_atoms,
            num_of_aos,
            num_of_mos,
            iterations,
            max_iterations,
            scf_converged,
            step_control,
            guess,
            cut_value,
            cdam_tau,
            cd_epsilon,
            scf_acceleration,
            damping_factor,
            comment,
        })
    }

    /// Reads molecular structure from `/molecule/xyzc`.
    fn read_molecule(file: &File, source: &str) -> Result<QcMolecule> {
        let ds = file.dataset("molecule/xyzc").map_err(|e| {
            BridgeError::input_error(source, format!("dataset 'molecule/xyzc' not found: {}", e))
        })?;

        let shape = ds.shape();
        if shape.len() != 2 || shape[1] < 6 {
            return Err(BridgeError::input_error(
                "molecule/xyzc",
                format!("expected 2D array with >=6 columns, found {:?}", shape),
            ));
        }

        let n_atoms = shape[0];
        let raw: ndarray::Array2<FixedAscii<32>> = ds.read_2d().map_err(|e| {
            BridgeError::input_error("molecule/xyzc", format!("failed to read FixedAscii: {}", e))
        })?;

        let mut atoms = Vec::with_capacity(n_atoms);
        for i in 0..n_atoms {
            let symbol = raw[[i, 0]].as_str().trim().to_string();
            let label = raw[[i, 1]].as_str().trim().to_string();
            let x_str = raw[[i, 2]].as_str().trim();
            let y_str = raw[[i, 3]].as_str().trim();
            let z_str = raw[[i, 4]].as_str().trim();
            let charge_str = raw[[i, 5]].as_str().trim();

            let x = x_str.parse::<f64>().map_err(|_| {
                BridgeError::input_error(x_str, format!("invalid x coordinate at atom {}", i))
            })?;
            let y = y_str.parse::<f64>().map_err(|_| {
                BridgeError::input_error(y_str, format!("invalid y coordinate at atom {}", i))
            })?;
            let z = z_str.parse::<f64>().map_err(|_| {
                BridgeError::input_error(z_str, format!("invalid z coordinate at atom {}", i))
            })?;
            let charge = charge_str.parse::<f64>().map_err(|_| {
                BridgeError::input_error(charge_str, format!("invalid charge/atomic number at atom {}", i))
            })?;

            atoms.push(QcAtom {
                symbol,
                label,
                pos_bohr: Position::new(x, y, z),
                charge,
            });
        }

        Ok(QcMolecule { atoms })
    }

    /// Reads basis set definition from `/basisset/{name, cgto, pgto}`.
    fn read_basisset(file: &File, source: &str) -> Result<QcBasisSet> {
        // 1. Read basisset/name
        let ds_name = file.dataset("basisset/name").map_err(|e| {
            BridgeError::input_error(source, format!("dataset 'basisset/name' not found: {}", e))
        })?;
        let name_shape = ds_name.shape();
        let n_defs = name_shape[0];

        // Read name pairs: try VarLenAscii first, fall back to FixedAscii<32>
        let mut def_names = Vec::with_capacity(n_defs);
        if let Ok(raw_name) = ds_name.read_2d::<VarLenAscii>() {
            for i in 0..n_defs {
                let name = raw_name[[i, 0]].as_str().trim().to_string();
                let atom_label = raw_name[[i, 1]].as_str().trim().to_string();
                def_names.push((name, atom_label));
            }
        } else if let Ok(raw_name) = ds_name.read_2d::<FixedAscii<32>>() {
            for i in 0..n_defs {
                let name = raw_name[[i, 0]].as_str().trim().to_string();
                let atom_label = raw_name[[i, 1]].as_str().trim().to_string();
                def_names.push((name, atom_label));
            }
        } else {
            return Err(BridgeError::input_error(
                "basisset/name",
                "unsupported string encoding for 'basisset/name'",
            ));
        }

        // 2. Read basisset/cgto
        let ds_cgto = file.dataset("basisset/cgto").map_err(|e| {
            BridgeError::input_error(source, format!("dataset 'basisset/cgto' not found: {}", e))
        })?;
        let cgto_shape = ds_cgto.shape();
        let n_cgtos = cgto_shape[0];
        let raw_cgto: ndarray::Array2<FixedAscii<32>> = ds_cgto.read_2d().map_err(|e| {
            BridgeError::input_error("basisset/cgto", format!("failed to read cgto: {}", e))
        })?;

        // 3. Read basisset/pgto
        let ds_pgto = file.dataset("basisset/pgto").map_err(|e| {
            BridgeError::input_error(source, format!("dataset 'basisset/pgto' not found: {}", e))
        })?;
        let pgto_shape = ds_pgto.shape();
        let n_pgtos = pgto_shape[0];
        let raw_pgto: ndarray::Array2<FixedAscii<32>> = ds_pgto.read_2d().map_err(|e| {
            BridgeError::input_error("basisset/pgto", format!("failed to read pgto: {}", e))
        })?;

        // Parse PGTO entries: group by (name, cgto_id)
        let mut pgtos_map: std::collections::HashMap<(String, usize), Vec<QcPgto>> =
            std::collections::HashMap::new();

        for i in 0..n_pgtos {
            let name = raw_pgto[[i, 0]].as_str().trim().to_string();
            let cgto_id_str = raw_pgto[[i, 1]].as_str().trim();
            let pgto_id_str = raw_pgto[[i, 2]].as_str().trim();
            let coef_str = raw_pgto[[i, 3]].as_str().trim();
            let exp_str = raw_pgto[[i, 4]].as_str().trim();

            let cgto_id = cgto_id_str.parse::<usize>().map_err(|_| {
                BridgeError::input_error(cgto_id_str, format!("invalid cgto_id at pgto {}", i))
            })?;
            let pgto_id = pgto_id_str.parse::<usize>().map_err(|_| {
                BridgeError::input_error(pgto_id_str, format!("invalid pgto_id at pgto {}", i))
            })?;
            let coef = coef_str.parse::<f64>().map_err(|_| {
                BridgeError::input_error(coef_str, format!("invalid coef at pgto {}", i))
            })?;
            let exp = exp_str.parse::<f64>().map_err(|_| {
                BridgeError::input_error(exp_str, format!("invalid exp at pgto {}", i))
            })?;

            pgtos_map
                .entry((name, cgto_id))
                .or_default()
                .push(QcPgto {
                    pgto_id,
                    coef,
                    exp,
                });
        }

        // Parse CGTO entries: group by name
        let mut cgtos_map: std::collections::HashMap<String, Vec<QcCgto>> =
            std::collections::HashMap::new();

        for i in 0..n_cgtos {
            let name = raw_cgto[[i, 0]].as_str().trim().to_string();
            let cgto_id_str = raw_cgto[[i, 1]].as_str().trim();
            let shell_type = raw_cgto[[i, 2]].as_str().trim().to_lowercase();
            let scale_str = raw_cgto[[i, 3]].as_str().trim();

            let cgto_id = cgto_id_str.parse::<usize>().map_err(|_| {
                BridgeError::input_error(cgto_id_str, format!("invalid cgto_id at cgto {}", i))
            })?;
            let scale_factor = scale_str.parse::<f64>().map_err(|_| {
                BridgeError::input_error(scale_str, format!("invalid scale_factor at cgto {}", i))
            })?;

            let mut pgtos = pgtos_map.remove(&(name.clone(), cgto_id)).unwrap_or_default();
            pgtos.sort_by_key(|p| p.pgto_id);

            cgtos_map.entry(name).or_default().push(QcCgto {
                cgto_id,
                shell_type,
                scale_factor,
                pgtos,
            });
        }

        // Assemble definitions
        let mut definitions = Vec::with_capacity(n_defs);
        for (name, atom_label) in def_names {
            let mut cgtos = cgtos_map.remove(&name).unwrap_or_default();
            cgtos.sort_by_key(|c| c.cgto_id);
            definitions.push(QcBasisSetDefinition {
                name,
                atom_label,
                cgtos,
            });
        }

        Ok(QcBasisSet { definitions })
    }

    /// Reads total energy history from `/TEs`.
    fn read_total_energies(file: &File) -> Option<Vec<f64>> {
        let ds = file.dataset("TEs").ok()?;
        ds.read_1d::<f64>().ok().map(|a| a.to_vec())
    }

    /// Reads orbital energies from `/energy_level`.
    fn read_orbital_energies(file: &File) -> Result<Vec<OrbitalEnergy>> {
        let group = match file.group("energy_level") {
            Ok(g) => g,
            Err(_) => return Ok(Vec::new()),
        };

        let mut list = Vec::new();
        for member_name in group.member_names().unwrap_or_default() {
            if let Some((run_type, iteration)) = Self::parse_member_name(&member_name) {
                if let Ok(ds) = group.dataset(&member_name) {
                    if let Ok(arr) = ds.read_1d::<f64>() {
                        list.push(OrbitalEnergy {
                            run_type,
                            iteration,
                            energies: arr.to_vec(),
                        });
                    }
                }
            }
        }

        list.sort_by(|a, b| {
            a.run_type
                .as_str()
                .cmp(b.run_type.as_str())
                .then(a.iteration.cmp(&b.iteration))
        });
        Ok(list)
    }

    /// Reads orbital occupations from `/occ`.
    fn read_orbital_occupations(file: &File) -> Result<Vec<OrbitalOccupation>> {
        let group = match file.group("occ") {
            Ok(g) => g,
            Err(_) => return Ok(Vec::new()),
        };

        let mut list = Vec::new();
        for member_name in group.member_names().unwrap_or_default() {
            let run_type = RunType::parse(&member_name);
            if let Ok(ds) = group.dataset(&member_name) {
                if let Ok(arr) = ds.read_1d::<f64>() {
                    list.push(OrbitalOccupation {
                        run_type,
                        occupations: arr.to_vec(),
                    });
                }
            }
        }

        list.sort_by(|a, b| a.run_type.as_str().cmp(b.run_type.as_str()));
        Ok(list)
    }

    /// Reads MO coefficients from `/C`.
    fn read_mo_coefficients(file: &File) -> Result<Vec<MoCoefficients>> {
        let group = match file.group("C") {
            Ok(g) => g,
            Err(_) => return Ok(Vec::new()),
        };

        let mut list = Vec::new();
        for member_name in group.member_names().unwrap_or_default() {
            if let Some((run_type, iteration)) = Self::parse_member_name(&member_name) {
                if let Ok(ds) = group.dataset(&member_name) {
                    let shape = ds.shape();
                    if shape.len() == 2 {
                        let n_ao = shape[0];
                        let n_mo = shape[1];
                        if let Ok(arr) = ds.read_2d::<f64>() {
                            let (data, _) = arr.into_raw_vec_and_offset();
                            list.push(MoCoefficients::new(run_type, iteration, n_ao, n_mo, data)?);
                        }

                    }
                }
            }
        }

        list.sort_by(|a, b| {
            a.run_type
                .as_str()
                .cmp(b.run_type.as_str())
                .then(a.iteration.cmp(&b.iteration))
        });
        Ok(list)
    }

    /// Reads Mulliken charges from `/pop/mulliken`.
    fn read_mulliken(file: &File) -> Result<Vec<MullikenPopulation>> {
        let group = match file.group("pop/mulliken") {
            Ok(g) => g,
            Err(_) => return Ok(Vec::new()),
        };

        let mut list = Vec::new();
        for member_name in group.member_names().unwrap_or_default() {
            if let Some((run_type, iteration)) = Self::parse_member_name(&member_name) {
                if let Ok(ds) = group.dataset(&member_name) {
                    if let Ok(arr) = ds.read_1d::<f64>() {
                        list.push(MullikenPopulation {
                            run_type,
                            iteration,
                            charges: arr.to_vec(),
                        });
                    }
                }
            }
        }

        list.sort_by(|a, b| {
            a.run_type
                .as_str()
                .cmp(b.run_type.as_str())
                .then(a.iteration.cmp(&b.iteration))
        });
        Ok(list)
    }

    /// Reads a 2D symmetric matrix dataset at a top-level path (e.g. `s`, `h`, `h2`).
    fn read_symmetric_matrix(file: &File, path: &str) -> Result<Option<SymmetricMatrix>> {
        let ds = match file.dataset(path) {
            Ok(d) => d,
            Err(_) => return Ok(None),
        };

        let shape = ds.shape();
        if shape.len() != 2 || shape[0] != shape[1] {
            return Err(BridgeError::input_error(
                path,
                format!("expected square matrix, found shape {:?}", shape),
            ));
        }

        let dim = shape[0];
        let arr = ds.read_2d::<f64>().map_err(|e| {
            BridgeError::input_error(path, format!("failed to read matrix: {}", e))
        })?;

        let (data, _) = arr.into_raw_vec_and_offset();
        Ok(Some(SymmetricMatrix::from_full(dim, data)?))
    }

    /// Reads symmetric matrices under a group (e.g. `P`, `F`).
    fn read_group_matrices(file: &File, group_name: &str) -> Result<Vec<SymmetricMatrix>> {
        let group = match file.group(group_name) {
            Ok(g) => g,
            Err(_) => return Ok(Vec::new()),
        };

        let mut list = Vec::new();
        for member_name in group.member_names().unwrap_or_default() {
            if let Ok(ds) = group.dataset(&member_name) {
                let shape = ds.shape();
                if shape.len() == 2 && shape[0] == shape[1] {
                    let dim = shape[0];
                    if let Ok(arr) = ds.read_2d::<f64>() {
                        let (data, _) = arr.into_raw_vec_and_offset();
                        list.push(SymmetricMatrix::from_full(dim, data)?);
                    }
                }
            }
        }


        Ok(list)
    }

    /// Parses dataset name into `(RunType, iteration)`:
    /// - `"rks_49"` -> `(RunType::Rks, 49)`
    /// - `"uks_alpha_39"` -> `(RunType::UksAlpha, 39)`
    /// - `"uks_beta_39"` -> `(RunType::UksBeta, 39)`
    fn parse_member_name(name: &str) -> Option<(RunType, usize)> {
        let last_underscore = name.rfind('_')?;
        let run_str = &name[..last_underscore];
        let itr_str = &name[last_underscore + 1..];
        let iteration = itr_str.parse::<usize>().ok()?;
        let run_type = RunType::parse(run_str);
        Some((run_type, iteration))
    }

    // Helper functions for reading root attributes

    fn read_str_attr(file: &File, name: &str) -> Result<Option<String>> {
        let attr = match file.attr(name) {
            Ok(a) => a,
            Err(_) => return Ok(None),
        };

        if let Ok(u) = attr.read_scalar::<VarLenUnicode>() {
            return Ok(Some(u.as_str().to_string()));
        }
        if let Ok(a) = attr.read_scalar::<VarLenAscii>() {
            return Ok(Some(a.as_str().to_string()));
        }
        if let Ok(f) = attr.read_scalar::<FixedAscii<32>>() {
            return Ok(Some(f.as_str().trim().to_string()));
        }

        Ok(None)
    }

    fn read_int_attr(file: &File, name: &str) -> Result<Option<usize>> {
        let attr = match file.attr(name) {
            Ok(a) => a,
            Err(_) => return Ok(None),
        };

        if let Ok(v) = attr.read_scalar::<i64>() {
            return Ok(Some(v as usize));
        }
        if let Ok(v) = attr.read_scalar::<u64>() {
            return Ok(Some(v as usize));
        }
        if let Ok(v) = attr.read_scalar::<i32>() {
            return Ok(Some(v as usize));
        }
        if let Ok(v) = attr.read_scalar::<u32>() {
            return Ok(Some(v as usize));
        }

        Ok(None)
    }

    fn read_f64_attr(file: &File, name: &str) -> Result<Option<f64>> {
        let attr = match file.attr(name) {
            Ok(a) => a,
            Err(_) => return Ok(None),
        };

        if let Ok(v) = attr.read_scalar::<f64>() {
            return Ok(Some(v));
        }
        if let Ok(v) = attr.read_scalar::<f32>() {
            return Ok(Some(v as f64));
        }

        Ok(None)
    }

    fn read_bool_attr(file: &File, name: &str) -> Result<Option<bool>> {
        let attr = match file.attr(name) {
            Ok(a) => a,
            Err(_) => return Ok(None),
        };

        if let Ok(v) = attr.read_scalar::<bool>() {
            return Ok(Some(v));
        }

        Ok(None)
    }
}
