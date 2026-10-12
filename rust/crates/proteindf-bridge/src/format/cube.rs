// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

//! Gaussian CUBE format parser.
//!
//! Specification reference:
//! - Gaussian cubegen utility documentation (<https://gaussian.com/cubegen/>)
//! - Paul Bourke, Gaussian Cube File Format (<https://paulbourke.net/dataformats/cube/>)

use std::fmt;
use std::fs;
use std::path::Path;

use crate::atom::Atom;
use crate::atom_group::AtomGroup;
use crate::error::{BridgeError, Result};
use crate::periodic_table::PeriodicTable;
use crate::position::Position;

/// Conversion factor from Bohr (atomic units of length) to Angstroms.
///
/// Reference: CODATA 2018 / 2022 recommended value for Bohr radius
/// $a_0 = 0.529177210903 \times 10^{-10}$ m = 0.529177210903 Å.
pub const BOHR_TO_ANGSTROM: f64 = 0.529_177_210_903;

/// Length unit of coordinates and grid parameters in a CUBE file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CubeUnit {
    /// Atomic units (bohr). Default when grid point counts are positive.
    Bohr,
    /// Angstroms. Indicated when grid point counts are negative.
    Angstrom,
}

impl fmt::Display for CubeUnit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CubeUnit::Bohr => write!(f, "bohr"),
            CubeUnit::Angstrom => write!(f, "angstrom"),
        }
    }
}

/// A single atom record in a Gaussian CUBE file.
#[derive(Debug, Clone, PartialEq)]
pub struct CubeAtom {
    /// Atomic number (Z).
    pub atomic_number: usize,
    /// Nuclear charge.
    pub charge: f64,
    /// X coordinate (in file unit).
    pub x: f64,
    /// Y coordinate (in file unit).
    pub y: f64,
    /// Z coordinate (in file unit).
    pub z: f64,
}

/// Gaussian CUBE file data container and parser.
#[derive(Debug, Clone, PartialEq)]
pub struct Cube {
    comment1: String,
    comment2: String,
    origin: [f64; 3],
    unit: CubeUnit,
    n_grid: [usize; 3],
    spacing: [[f64; 3]; 3],
    atoms: Vec<CubeAtom>,
    orbitals: Option<Vec<i32>>,
    data: Vec<Vec<f64>>,
}

impl Default for Cube {
    fn default() -> Self {
        Self::new()
    }
}

impl Cube {
    /// Creates a new empty `Cube` instance.
    pub fn new() -> Self {
        Self {
            comment1: String::new(),
            comment2: String::new(),
            origin: [0.0, 0.0, 0.0],
            unit: CubeUnit::Bohr,
            n_grid: [0, 0, 0],
            spacing: [[0.0; 3]; 3],
            atoms: Vec::new(),
            orbitals: None,
            data: Vec::new(),
        }
    }

    /// Loads a Gaussian CUBE file from the given path.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        let mut cube = Self::new();
        cube.load(path)?;
        Ok(cube)
    }

    /// Parses a Gaussian CUBE format string.
    pub fn parse_str(content: &str) -> Result<Self> {
        let mut cube = Self::new();
        cube.parse(content)?;
        Ok(cube)
    }

    /// Loads Gaussian CUBE data from a file path.
    pub fn load(&mut self, path: impl AsRef<Path>) -> Result<()> {
        let p = path.as_ref();
        let content = fs::read_to_string(p).map_err(|e| {
            BridgeError::input_error(
                p.display().to_string(),
                format!("failed to read file: {}", e),
            )
        })?;
        self.parse(&content)
    }

    /// First comment line.
    pub fn comment1(&self) -> &str {
        &self.comment1
    }

    /// Sets first comment line.
    pub fn set_comment1(&mut self, val: impl Into<String>) {
        self.comment1 = val.into();
    }

    /// Second comment line.
    pub fn comment2(&self) -> &str {
        &self.comment2
    }

    /// Sets second comment line.
    pub fn set_comment2(&mut self, val: impl Into<String>) {
        self.comment2 = val.into();
    }

    /// Volumetric data origin [x, y, z] (in file unit).
    pub fn origin(&self) -> [f64; 3] {
        self.origin
    }

    /// Length unit of coordinates and grid vectors.
    pub fn unit(&self) -> CubeUnit {
        self.unit
    }

    /// Number of grid points along each axis [nx, ny, nz].
    pub fn n_grid(&self) -> [usize; 3] {
        self.n_grid
    }

    /// Step vectors along each axis [[dx_x, dx_y, dx_z], [dy_x, dy_y, dy_z], [dz_x, dz_y, dz_z]].
    pub fn spacing(&self) -> [[f64; 3]; 3] {
        self.spacing
    }

    /// List of atoms in the file.
    pub fn atoms(&self) -> &[CubeAtom] {
        &self.atoms
    }

    /// Number of atoms.
    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }

    /// Returns true if this CUBE file contains molecular orbital data (`NAtoms < 0`).
    pub fn is_molecular_orbital(&self) -> bool {
        self.orbitals.is_some()
    }

    /// Orbital numbers if molecular orbital CUBE, or None.
    pub fn orbitals(&self) -> Option<&[i32]> {
        self.orbitals.as_deref()
    }

    /// Number of datasets (number of orbitals if MO CUBE, or 1 for standard scalar field).
    pub fn num_datasets(&self) -> usize {
        self.data.len()
    }

    /// Computes the coordinate of a grid point `(ix, iy, iz)` in file units.
    ///
    /// Formula: `origin + ix * d1 + iy * d2 + iz * d3`.
    pub fn grid_point_position(&self, ix: usize, iy: usize, iz: usize) -> Result<[f64; 3]> {
        if ix >= self.n_grid[0] || iy >= self.n_grid[1] || iz >= self.n_grid[2] {
            return Err(BridgeError::value_error(
                format!("index ({}, {}, {})", ix, iy, iz),
                format!(
                    "grid index out of bounds for grid size ({}, {}, {})",
                    self.n_grid[0], self.n_grid[1], self.n_grid[2]
                ),
            ));
        }

        let pos = [
            self.origin[0]
                + ix as f64 * self.spacing[0][0]
                + iy as f64 * self.spacing[1][0]
                + iz as f64 * self.spacing[2][0],
            self.origin[1]
                + ix as f64 * self.spacing[0][1]
                + iy as f64 * self.spacing[1][1]
                + iz as f64 * self.spacing[2][1],
            self.origin[2]
                + ix as f64 * self.spacing[0][2]
                + iy as f64 * self.spacing[1][2]
                + iz as f64 * self.spacing[2][2],
        ];
        Ok(pos)
    }

    /// Retrieves the scalar value at grid point `(ix, iy, iz)` for dataset 0.
    pub fn value(&self, ix: usize, iy: usize, iz: usize) -> Result<f64> {
        self.orbital_value(0, ix, iy, iz)
    }

    /// Retrieves the scalar value at grid point `(ix, iy, iz)` for the specified dataset / orbital index.
    pub fn orbital_value(
        &self,
        dataset_idx: usize,
        ix: usize,
        iy: usize,
        iz: usize,
    ) -> Result<f64> {
        if dataset_idx >= self.data.len() {
            return Err(BridgeError::value_error(
                format!("dataset_idx {}", dataset_idx),
                format!(
                    "dataset index out of bounds, total datasets: {}",
                    self.data.len()
                ),
            ));
        }
        if ix >= self.n_grid[0] || iy >= self.n_grid[1] || iz >= self.n_grid[2] {
            return Err(BridgeError::value_error(
                format!("index ({}, {}, {})", ix, iy, iz),
                format!(
                    "grid index out of bounds for grid size ({}, {}, {})",
                    self.n_grid[0], self.n_grid[1], self.n_grid[2]
                ),
            ));
        }

        let flat_idx = (ix * self.n_grid[1] + iy) * self.n_grid[2] + iz;
        Ok(self.data[dataset_idx][flat_idx])
    }

    /// Returns the flat slice of values for the given dataset index (length `nx * ny * nz`).
    pub fn data(&self, dataset_idx: usize) -> Result<&[f64]> {
        self.data
            .get(dataset_idx)
            .map(|v| v.as_slice())
            .ok_or_else(|| {
                BridgeError::value_error(
                    format!("dataset_idx {}", dataset_idx),
                    format!(
                        "dataset index out of bounds, total datasets: {}",
                        self.data.len()
                    ),
                )
            })
    }

    /// Returns all datasets.
    pub fn all_data(&self) -> &[Vec<f64>] {
        &self.data
    }

    /// Builds and returns an [`AtomGroup`] from the atoms in this CUBE file.
    ///
    /// Atomic coordinates are converted to Angstroms:
    /// - If the file unit is Bohr, coordinates are multiplied by [`BOHR_TO_ANGSTROM`].
    /// - If the file unit is Angstrom, coordinates are preserved as is.
    pub fn get_atomgroup(&self) -> Result<AtomGroup> {
        let name = if !self.comment1.trim().is_empty() {
            self.comment1.trim()
        } else if !self.comment2.trim().is_empty() {
            self.comment2.trim()
        } else {
            "cube_structure"
        };

        let mut ag = AtomGroup::with_name(name);

        for (i, atom_data) in self.atoms.iter().enumerate() {
            let symbol = PeriodicTable::get_symbol(atom_data.atomic_number)?;
            let pos_angstrom = match self.unit {
                CubeUnit::Bohr => Position::new(
                    atom_data.x * BOHR_TO_ANGSTROM,
                    atom_data.y * BOHR_TO_ANGSTROM,
                    atom_data.z * BOHR_TO_ANGSTROM,
                ),
                CubeUnit::Angstrom => Position::new(atom_data.x, atom_data.y, atom_data.z),
            };

            let mut atom = Atom::new_with_pos(symbol, pos_angstrom)?;
            atom.charge = atom_data.charge;
            atom.name = format!("{}{}", symbol, i + 1);
            ag.set_atom(&i.to_string(), atom);
        }

        Ok(ag)
    }

    /// Alias for [`get_atomgroup`].
    pub fn get_atom_group(&self) -> Result<AtomGroup> {
        self.get_atomgroup()
    }

    /// Parses Gaussian CUBE format string into this instance.
    pub fn parse(&mut self, content: &str) -> Result<()> {
        let mut lines = content.lines().enumerate();

        // Line 1: Comment 1
        let (_, l1) = lines.next().ok_or_else(|| {
            BridgeError::input_error("line 1", "unexpected EOF at line 1 (comment 1)")
        })?;
        let comment1 = l1.trim().to_string();

        // Line 2: Comment 2
        let (_, l2) = lines.next().ok_or_else(|| {
            BridgeError::input_error("line 2", "unexpected EOF at line 2 (comment 2)")
        })?;
        let comment2 = l2.trim().to_string();

        // Line 3: Natoms, Origin (x, y, z), [NVal]
        let (line_idx_3, l3) = lines.next().ok_or_else(|| {
            BridgeError::input_error("line 3", "unexpected EOF at line 3 (atoms/origin)")
        })?;
        let line_num_3 = line_idx_3 + 1;
        let tokens_3: Vec<&str> = l3.split_whitespace().collect();
        if tokens_3.len() < 4 {
            return Err(BridgeError::input_error(
                format!("line {}: '{}'", line_num_3, l3),
                "insufficient fields for atoms count and origin (expected at least 4)",
            ));
        }

        let natoms_signed = tokens_3[0].parse::<i64>().map_err(|_| {
            BridgeError::input_error(
                format!("line {}: '{}'", line_num_3, tokens_3[0]),
                "invalid integer for atom count",
            )
        })?;

        let ox = tokens_3[1].parse::<f64>().map_err(|_| {
            BridgeError::input_error(
                format!("line {}: '{}'", line_num_3, tokens_3[1]),
                "invalid float for origin X",
            )
        })?;
        let oy = tokens_3[2].parse::<f64>().map_err(|_| {
            BridgeError::input_error(
                format!("line {}: '{}'", line_num_3, tokens_3[2]),
                "invalid float for origin Y",
            )
        })?;
        let oz = tokens_3[3].parse::<f64>().map_err(|_| {
            BridgeError::input_error(
                format!("line {}: '{}'", line_num_3, tokens_3[3]),
                "invalid float for origin Z",
            )
        })?;

        let is_mo = natoms_signed < 0;
        let num_atoms = natoms_signed.unsigned_abs() as usize;

        // Lines 4-6: 3 grid axis specifications
        let mut n_grid = [0usize; 3];
        let mut spacing = [[0.0f64; 3]; 3];
        let mut unit = CubeUnit::Bohr;

        for axis in 0..3 {
            let (line_idx, l) = lines.next().ok_or_else(|| {
                BridgeError::input_error(
                    format!("line {}", axis + 4),
                    format!("unexpected EOF at grid axis line {}", axis + 4),
                )
            })?;
            let line_num = line_idx + 1;
            let tokens: Vec<&str> = l.split_whitespace().collect();
            if tokens.len() < 4 {
                return Err(BridgeError::input_error(
                    format!("line {}: '{}'", line_num, l),
                    "insufficient fields for grid axis (expected at least 4: N, dx, dy, dz)",
                ));
            }

            let n_axis_signed = tokens[0].parse::<i64>().map_err(|_| {
                BridgeError::input_error(
                    format!("line {}: '{}'", line_num, tokens[0]),
                    "invalid integer for grid point count",
                )
            })?;

            if n_axis_signed == 0 {
                return Err(BridgeError::input_error(
                    format!("line {}: '{}'", line_num, tokens[0]),
                    "grid point count must not be 0",
                ));
            }

            if n_axis_signed < 0 {
                unit = CubeUnit::Angstrom;
            }
            n_grid[axis] = n_axis_signed.unsigned_abs() as usize;

            for coord in 0..3 {
                spacing[axis][coord] = tokens[coord + 1].parse::<f64>().map_err(|_| {
                    BridgeError::input_error(
                        format!("line {}: '{}'", line_num, tokens[coord + 1]),
                        format!("invalid float for grid step coordinate {}", coord + 1),
                    )
                })?;
            }
        }

        // Lines 7 .. 7 + num_atoms: Atom lines
        let mut atoms = Vec::with_capacity(num_atoms);
        for atom_idx in 0..num_atoms {
            let (line_idx, l) = lines.next().ok_or_else(|| {
                BridgeError::input_error(
                    format!("line {}", atom_idx + 7),
                    format!(
                        "unexpected EOF at atom {} (line {})",
                        atom_idx + 1,
                        atom_idx + 7
                    ),
                )
            })?;
            let line_num = line_idx + 1;
            let tokens: Vec<&str> = l.split_whitespace().collect();
            if tokens.len() < 5 {
                return Err(BridgeError::input_error(
                    format!("line {}: '{}'", line_num, l),
                    "insufficient fields for atom record (expected at least 5: Z, charge, x, y, z)",
                ));
            }

            let atomic_number = tokens[0].parse::<usize>().map_err(|_| {
                BridgeError::input_error(
                    format!("line {}: '{}'", line_num, tokens[0]),
                    "invalid atomic number in atom record",
                )
            })?;

            let charge = tokens[1].parse::<f64>().map_err(|_| {
                BridgeError::input_error(
                    format!("line {}: '{}'", line_num, tokens[1]),
                    "invalid charge in atom record",
                )
            })?;
            let x = tokens[2].parse::<f64>().map_err(|_| {
                BridgeError::input_error(
                    format!("line {}: '{}'", line_num, tokens[2]),
                    "invalid x coordinate in atom record",
                )
            })?;
            let y = tokens[3].parse::<f64>().map_err(|_| {
                BridgeError::input_error(
                    format!("line {}: '{}'", line_num, tokens[3]),
                    "invalid y coordinate in atom record",
                )
            })?;
            let z = tokens[4].parse::<f64>().map_err(|_| {
                BridgeError::input_error(
                    format!("line {}: '{}'", line_num, tokens[4]),
                    "invalid z coordinate in atom record",
                )
            })?;

            atoms.push(CubeAtom {
                atomic_number,
                charge,
                x,
                y,
                z,
            });
        }

        // If MO CUBE, parse orbital line(s)
        let (orbitals, num_datasets) = if is_mo {
            let mut mo_tokens: Vec<(usize, String)> = Vec::new();
            // We need 1 number for NMO, then NMO orbital numbers
            let mut expected_nmo: Option<usize> = None;

            while expected_nmo.is_none() || mo_tokens.len() < 1 + expected_nmo.unwrap() {
                let (line_idx, l) = lines.next().ok_or_else(|| {
                    BridgeError::input_error(
                        "MO header",
                        "unexpected EOF while reading molecular orbital numbers",
                    )
                })?;
                let line_num = line_idx + 1;
                for tok in l.split_whitespace() {
                    mo_tokens.push((line_num, tok.to_string()));
                    if expected_nmo.is_none() {
                        let nmo = tok.parse::<usize>().map_err(|_| {
                            BridgeError::input_error(
                                format!("line {}: '{}'", line_num, tok),
                                "invalid integer for number of molecular orbitals",
                            )
                        })?;
                        if nmo == 0 {
                            return Err(BridgeError::input_error(
                                format!("line {}: '{}'", line_num, tok),
                                "number of molecular orbitals must be greater than 0",
                            ));
                        }
                        expected_nmo = Some(nmo);
                    }
                }
            }

            let nmo = expected_nmo.unwrap();
            let mut orb_list = Vec::with_capacity(nmo);
            for (line_num, tok) in mo_tokens.into_iter().skip(1).take(nmo) {
                let orb_id = tok.parse::<i32>().map_err(|_| {
                    BridgeError::input_error(
                        format!("line {}: '{}'", line_num, tok),
                        "invalid integer for orbital number",
                    )
                })?;
                orb_list.push(orb_id);
            }

            (Some(orb_list), nmo)
        } else {
            (None, 1)
        };

        // Parse volumetric grid data
        let total_grid_points = n_grid[0] * n_grid[1] * n_grid[2];
        let total_values_expected = total_grid_points * num_datasets;

        let mut datasets: Vec<Vec<f64>> = (0..num_datasets)
            .map(|_| Vec::with_capacity(total_grid_points))
            .collect();

        let mut read_values = 0usize;
        let mut current_point_idx = 0usize;
        let mut current_dataset_idx = 0usize;

        for (line_idx, l) in lines {
            let line_num = line_idx + 1;
            let tokens = l.split_whitespace();
            for tok in tokens {
                if read_values >= total_values_expected {
                    return Err(BridgeError::input_error(
                        format!("line {}: '{}'", line_num, tok),
                        format!(
                            "too many volumetric values in CUBE file (expected {}, found more)",
                            total_values_expected
                        ),
                    ));
                }

                let val = tok.parse::<f64>().map_err(|_| {
                    BridgeError::input_error(
                        format!("line {}: '{}'", line_num, tok),
                        "invalid float value in volumetric data",
                    )
                })?;

                datasets[current_dataset_idx].push(val);
                read_values += 1;

                current_dataset_idx += 1;
                if current_dataset_idx == num_datasets {
                    current_dataset_idx = 0;
                    current_point_idx += 1;
                }
            }
        }

        let _ = current_point_idx;

        if read_values < total_values_expected {
            return Err(BridgeError::input_error(
                "volumetric data",
                format!(
                    "insufficient volumetric values in CUBE file (expected {}, found {})",
                    total_values_expected, read_values
                ),
            ));
        }

        self.comment1 = comment1;
        self.comment2 = comment2;
        self.origin = [ox, oy, oz];
        self.unit = unit;
        self.n_grid = n_grid;
        self.spacing = spacing;
        self.atoms = atoms;
        self.orbitals = orbitals;
        self.data = datasets;

        Ok(())
    }
}
