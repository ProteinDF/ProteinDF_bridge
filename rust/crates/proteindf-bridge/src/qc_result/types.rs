// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

//! Data types representing ProteinDF quantum chemistry calculation results.

use std::fmt;

use crate::atom::Atom;
use crate::atom_group::AtomGroup;
use crate::error::{BridgeError, Result};
use crate::format::cube::BOHR_TO_ANGSTROM;
use crate::matrix::{Matrix, SymmetricMatrix};
use crate::position::Position;

/// Spin calculation method.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QcMethod {
    /// Restricted Kohn-Sham / Hartree-Fock (closed shell).
    Rks,
    /// Unrestricted Kohn-Sham / Hartree-Fock (open shell, alpha/beta spin).
    Uks,
    /// Restricted open-shell Kohn-Sham.
    Roks,
}

impl QcMethod {
    /// Parses a method string (e.g. `"rks"`, `"uks"`, `"roks"`).
    pub fn parse(s: &str) -> Option<Self> {
        match s.to_ascii_lowercase().as_str() {
            "rks" => Some(QcMethod::Rks),
            "uks" => Some(QcMethod::Uks),
            "roks" => Some(QcMethod::Roks),
            _ => None,
        }
    }

    /// Returns the canonical method name.
    pub fn as_str(&self) -> &'static str {
        match self {
            QcMethod::Rks => "rks",
            QcMethod::Uks => "uks",
            QcMethod::Roks => "roks",
        }
    }
}

impl fmt::Display for QcMethod {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

/// Run type specifying spin component or calculation type.
///
/// In ProteinDF (`pdfparam_object.py` 61-68行目):
/// - RKS: `"rks"`
/// - UKS: `"uks_alpha"`, `"uks_beta"`
/// - ROKS: `"roks_close"`, `"roks_open"`
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum RunType {
    /// RKS single component.
    Rks,
    /// UKS Alpha spin component (`"uks_alpha"`).
    UksAlpha,
    /// UKS Beta spin component (`"uks_beta"`).
    UksBeta,
    /// ROKS closed shell component (`"roks_close"`).
    RoksClose,
    /// ROKS open shell component (`"roks_open"`).
    RoksOpen,
    /// Custom / arbitrary run type identifier.
    Other(String),
}

impl RunType {
    /// Parses a run type string.
    pub fn parse(s: &str) -> Self {
        match s.to_ascii_lowercase().as_str() {
            "rks" => RunType::Rks,
            "uks_alpha" => RunType::UksAlpha,
            "uks_beta" => RunType::UksBeta,
            "roks_close" => RunType::RoksClose,
            "roks_open" => RunType::RoksOpen,
            other => RunType::Other(other.to_string()),
        }
    }

    /// Returns the string representation.
    pub fn as_str(&self) -> &str {
        match self {
            RunType::Rks => "rks",
            RunType::UksAlpha => "uks_alpha",
            RunType::UksBeta => "uks_beta",
            RunType::RoksClose => "roks_close",
            RunType::RoksOpen => "roks_open",
            RunType::Other(s) => s.as_str(),
        }
    }
}

impl fmt::Display for RunType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

/// Root-level quantum chemistry calculation parameters and metadata.
#[derive(Debug, Clone, PartialEq)]
pub struct QcParameters {
    /// Calculation method (`"rks"`, `"uks"`, etc.).
    pub method: String,
    /// XC functional or method name (e.g. `"HF"`, `"B3LYP"`).
    pub xc_functional: String,
    /// Number of atoms ($N_{\text{atom}}$).
    pub num_of_atoms: usize,
    /// Number of atomic basis functions / AOs ($N_{\text{AO}}$).
    pub num_of_aos: usize,
    /// Number of molecular orbitals ($N_{\text{MO}}$).
    pub num_of_mos: usize,
    /// Executed SCF iterations.
    pub iterations: usize,
    /// Maximum SCF iterations setting.
    pub max_iterations: usize,
    /// Whether the SCF calculation converged.
    pub scf_converged: bool,
    /// Step control string (e.g. `"create integral guess scf force"`).
    pub step_control: String,
    /// Initial guess type (e.g. `"core"`).
    pub guess: String,
    /// Cutoff threshold (string representation, e.g. `"1.0e-10"`).
    pub cut_value: String,
    /// Cholesky decomposition tau parameter (e.g. `"1.0E-10"`).
    pub cdam_tau: String,
    /// Cholesky decomposition epsilon parameter (e.g. `"1.0E-6"`).
    pub cd_epsilon: String,
    /// SCF acceleration method (e.g. `"damping"`).
    pub scf_acceleration: String,
    /// Damping factor if damping acceleration is used.
    pub damping_factor: Option<f64>,
    /// User comment.
    pub comment: String,
}

/// Single atom in a quantum chemistry calculation result.
#[derive(Debug, Clone, PartialEq)]
pub struct QcAtom {
    /// Element symbol (e.g. `"O"`, `"H"`, `"N"`).
    pub symbol: String,
    /// Atom label (e.g. `""`, or residue/atom modifier).
    pub label: String,
    /// Coordinates in atomic units / Bohr ($x$, $y$, $z$).
    pub pos_bohr: Position,
    /// Nuclear charge / atomic number $Z$.
    pub charge: f64,
}

impl QcAtom {
    /// Formatted atom label (`symbol` or `symbol@label`).
    ///
    /// Follows ProteinDF convention (`orbinfo.py` 99-103行目).
    pub fn atom_label(&self) -> String {
        if self.label.is_empty() {
            self.symbol.clone()
        } else {
            format!("{}@{}", self.symbol, self.label)
        }
    }

    /// Coordinates converted to Angstroms ($x, y, z$).
    pub fn pos_angstrom(&self) -> Position {
        Position::new(
            self.pos_bohr.x * BOHR_TO_ANGSTROM,
            self.pos_bohr.y * BOHR_TO_ANGSTROM,
            self.pos_bohr.z * BOHR_TO_ANGSTROM,
        )
    }
}

/// Molecular geometry in quantum chemistry calculation result.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct QcMolecule {
    /// List of atoms.
    pub atoms: Vec<QcAtom>,
}

impl QcMolecule {
    /// Number of atoms.
    pub fn len(&self) -> usize {
        self.atoms.len()
    }

    /// Whether there are no atoms.
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    /// Slice of atoms.
    pub fn atoms(&self) -> &[QcAtom] {
        &self.atoms
    }

    /// Converts molecular geometry to an [`AtomGroup`].
    ///
    /// Atomic coordinates are converted from Bohr to Angstroms using
    /// [`BOHR_TO_ANGSTROM`].
    pub fn get_atomgroup(&self, name: Option<&str>) -> Result<AtomGroup> {
        let group_name = name.unwrap_or("qc_molecule");
        let mut ag = AtomGroup::with_name(group_name);

        for (i, atom_data) in self.atoms.iter().enumerate() {
            let symbol = &atom_data.symbol;
            let pos_angstrom = atom_data.pos_angstrom();
            let mut atom = Atom::new_with_pos(symbol, pos_angstrom)?;
            atom.charge = atom_data.charge;
            if !atom_data.label.is_empty() {
                atom.name = atom_data.label.clone();
            } else {
                atom.name = format!("{}{}", symbol, i + 1);
            }
            ag.set_atom(&i.to_string(), atom);
        }

        Ok(ag)
    }
}

/// Primitive Gaussian-Type Orbital (PGTO).
#[derive(Debug, Clone, PartialEq)]
pub struct QcPgto {
    /// PGTO index within CGTO.
    pub pgto_id: usize,
    /// Contraction coefficient.
    pub coef: f64,
    /// Gaussian exponent $\alpha$.
    pub exp: f64,
}

/// Contracted Gaussian-Type Orbital (CGTO).
#[derive(Debug, Clone, PartialEq)]
pub struct QcCgto {
    /// CGTO index within basis set definition.
    pub cgto_id: usize,
    /// Shell type string (e.g. `"s"`, `"p"`, `"d"`, `"f"`).
    pub shell_type: String,
    /// Scale factor.
    pub scale_factor: f64,
    /// Primitive GTOs belonging to this CGTO.
    pub pgtos: Vec<QcPgto>,
}

/// Basis set definition for an atom type / label.
#[derive(Debug, Clone, PartialEq)]
pub struct QcBasisSetDefinition {
    /// Basis set name (e.g. `"O-DZVP2.H"`).
    pub name: String,
    /// Corresponding atom label (e.g. `"H"`).
    pub atom_label: String,
    /// Contracted GTOs in this basis set.
    pub cgtos: Vec<QcCgto>,
}

/// Basis set collection for the system.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct QcBasisSet {
    /// Mapping of atom_label to basis set definition.
    pub definitions: Vec<QcBasisSetDefinition>,
}

impl QcBasisSet {
    /// Finds basis set definition for an atom label.
    pub fn find_by_atom_label(&self, atom_label: &str) -> Option<&QcBasisSetDefinition> {
        self.definitions.iter().find(|d| d.atom_label == atom_label)
    }

    /// Finds basis set definition by basis set name.
    pub fn find_by_name(&self, name: &str) -> Option<&QcBasisSetDefinition> {
        self.definitions.iter().find(|d| d.name == name)
    }
}

/// Information about a single atomic orbital (basis function).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomicOrbital {
    /// Global AO index $0 \le p < N_{\text{AO}}$.
    pub ao_index: usize,
    /// Index of the atom this AO belongs to ($0 \le A < N_{\text{atom}}$).
    pub atom_index: usize,
    /// Element symbol of the parent atom.
    pub atom_symbol: String,
    /// Basis set name.
    pub basisset_name: String,
    /// CGTO index within the basis set.
    pub cgto_index: usize,
    /// Shell type (`"s"`, `"p"`, `"d"`, `"f"`).
    pub shell_type: String,
    /// Basis type / component name:
    /// - `s`: `"s"`
    /// - `p`: `"px"`, `"py"`, `"pz"`
    /// - `d`: `"dxy"`, `"dyz"`, `"dzx"`, `"dxx-yy"`, `"dzz"`
    /// - `f`: `"z3"`, `"xz2"`, `"yz2"`, `"3x2y-y3"`, `"x3-3xy2"`, `"xyz"`, `"x2z-y2z"`
    pub basis_type: String,
}

/// Orbital metadata and AO-to-atom mapping table.
///
/// Follows ProteinDF `orbinfo.py` and `basisset.py`:
/// - Iterates through each atom $A$ in order ($0 .. N_{\text{atom}}-1$).
/// - Finds the basis set for the atom's label (`symbol` or `symbol@label`).
/// - Iterates through each CGTO in the basis set.
/// - Expands the CGTO into $2l + 1$ basis functions according to `shell_type`:
///   - s: 1 (`"s"`)
///   - p: 3 (`"px"`, `"py"`, `"pz"`)
///   - d: 5 (`"dxy"`, `"dyz"`, `"dzx"`, `"dxx-yy"`, `"dzz"`)
///   - f: 7 (`"z3"`, `"xz2"`, `"yz2"`, `"3x2y-y3"`, `"x3-3xy2"`, `"xyz"`, `"x2z-y2z"`)
#[derive(Debug, Clone, PartialEq, Default)]
pub struct OrbitalInfo {
    /// Total number of atomic orbitals.
    pub num_of_aos: usize,
    /// List of atomic orbital definitions.
    pub orbitals: Vec<AtomicOrbital>,
}

impl OrbitalInfo {
    /// Builds `OrbitalInfo` from `QcMolecule` and `QcBasisSet`.
    pub fn build(molecule: &QcMolecule, basisset: &QcBasisSet) -> Result<Self> {
        let mut orbitals = Vec::new();
        let mut ao_index = 0;

        for (atom_idx, atom) in molecule.atoms.iter().enumerate() {
            let atom_label = atom.atom_label();
            // In ProteinDF (pdfparam_object.py lines 661-670 `_get_basisset_common`), atoms without
            // a basis set definition (e.g. dummy atoms like 'X') return an empty BasisSet (0 AOs).
            let Some(def) = basisset.find_by_atom_label(&atom_label) else {
                continue;
            };

            for cgto in &def.cgtos {
                let basis_components =
                    Self::shell_components(&cgto.shell_type).ok_or_else(|| {
                        BridgeError::input_error(
                            &cgto.shell_type,
                            format!(
                                "unsupported shell type '{}' in basis set '{}'",
                                cgto.shell_type, def.name
                            ),
                        )
                    })?;

                for b_type in basis_components {
                    orbitals.push(AtomicOrbital {
                        ao_index,
                        atom_index: atom_idx,
                        atom_symbol: atom.symbol.clone(),
                        basisset_name: def.name.clone(),
                        cgto_index: cgto.cgto_id,
                        shell_type: cgto.shell_type.to_lowercase(),
                        basis_type: b_type.to_string(),
                    });
                    ao_index += 1;
                }
            }
        }

        Ok(OrbitalInfo {
            num_of_aos: orbitals.len(),
            orbitals,
        })
    }

    /// Returns basis function component names for a shell type.
    ///
    /// Corresponds to ProteinDF `basisset.py` `_orb_types` and `get_basis_type`.
    pub fn shell_components(shell_type: &str) -> Option<&'static [&'static str]> {
        match shell_type.to_ascii_lowercase().as_str() {
            "s" => Some(&["s"]),
            "p" => Some(&["px", "py", "pz"]),
            "d" => Some(&["dxy", "dyz", "dzx", "dxx-yy", "dzz"]),
            "f" => Some(&["z3", "xz2", "yz2", "3x2y-y3", "x3-3xy2", "xyz", "x2z-y2z"]),
            _ => None,
        }
    }

    /// Number of orbitals.
    pub fn len(&self) -> usize {
        self.orbitals.len()
    }

    /// Whether empty.
    pub fn is_empty(&self) -> bool {
        self.orbitals.is_empty()
    }

    /// Returns atom index for an AO index.
    pub fn atom_id(&self, ao_idx: usize) -> Option<usize> {
        self.orbitals.get(ao_idx).map(|ao| ao.atom_index)
    }

    /// Returns AO definition for an AO index.
    pub fn get(&self, ao_idx: usize) -> Option<&AtomicOrbital> {
        self.orbitals.get(ao_idx)
    }
}

/// Orbital energies for a spin component or calculation.
#[derive(Debug, Clone, PartialEq)]
pub struct OrbitalEnergy {
    /// Run type (e.g. `RunType::Rks`, `RunType::UksAlpha`).
    pub run_type: RunType,
    /// Iteration index.
    pub iteration: usize,
    /// Energy levels in Hartree (a.u.), length $N_{\text{MO}}$.
    pub energies: Vec<f64>,
}

impl OrbitalEnergy {
    /// Number of molecular orbitals.
    pub fn len(&self) -> usize {
        self.energies.len()
    }

    /// Whether empty.
    pub fn is_empty(&self) -> bool {
        self.energies.is_empty()
    }
}

/// Orbital occupations.
#[derive(Debug, Clone, PartialEq)]
pub struct OrbitalOccupation {
    /// Run type (e.g. `RunType::Rks`, `RunType::UksAlpha`).
    pub run_type: RunType,
    /// Occupation numbers, length $N_{\text{MO}}$.
    pub occupations: Vec<f64>,
}

/// Molecular orbital coefficient matrix $C$ ($N_{\text{AO}} \times N_{\text{MO}}$).
///
/// In ProteinDF HDF5 (`docs/qc-result-format.md` §1.3 (3)):
/// - Shape: `(N_AO, N_MO)`
/// - Row $p$: AO index ($0 \le p < N_{\text{AO}}$)
/// - Column $i$: MO index ($0 \le i < N_{\text{MO}}$)
/// - Value: $C_{pi}$ (expansion coefficient of AO $p$ in MO $i$)
///
/// Represented internally using [`Matrix`].
#[derive(Debug, Clone, PartialEq)]
pub struct MoCoefficients {
    /// Run type.
    pub run_type: RunType,
    /// Iteration index.
    pub iteration: usize,
    /// MO coefficient matrix ($N_{\text{AO}} \times N_{\text{MO}}$).
    pub matrix: Matrix,
}

impl MoCoefficients {
    /// Creates a new `MoCoefficients` container.
    pub fn new(run_type: RunType, iteration: usize, matrix: Matrix) -> Self {
        Self {
            run_type,
            iteration,
            matrix,
        }
    }

    /// Number of rows ($N_{\text{AO}}$).
    pub fn n_ao(&self) -> usize {
        self.matrix.rows()
    }

    /// Number of columns ($N_{\text{MO}}$).
    pub fn n_mo(&self) -> usize {
        self.matrix.cols()
    }

    /// Coefficient $C_{pi}$ for AO $p$ and MO $i$.
    pub fn get(&self, ao_idx: usize, mo_idx: usize) -> Option<f64> {
        self.matrix.get(ao_idx, mo_idx).ok()
    }

    /// Column vector $\mathbf{C}_i$ for MO $i$ (length $N_{\text{AO}}$).
    pub fn mo_vector(&self, mo_idx: usize) -> Result<Vec<f64>> {
        if mo_idx >= self.matrix.cols() {
            return Err(BridgeError::value_error(
                format!("mo_idx {}", mo_idx),
                format!("MO index out of range (n_mo: {})", self.matrix.cols()),
            ));
        }
        let mut vec = Vec::with_capacity(self.matrix.rows());
        for p in 0..self.matrix.rows() {
            vec.push(self.matrix.get(p, mo_idx)?);
        }
        Ok(vec)
    }
}

/// Tolerance for verifying matrix symmetry ($|A_{ij} - A_{ji}| \le 10^{-8}$)
/// when reading full $N \times N$ symmetric matrix datasets from HDF5.
pub const MATRIX_SYMMETRY_TOLERANCE: f64 = 1.0e-8;

/// Mulliken population analysis result (net atomic charges).
#[derive(Debug, Clone, PartialEq)]
pub struct MullikenPopulation {
    /// Run type.
    pub run_type: RunType,
    /// Iteration index.
    pub iteration: usize,
    /// Net Mulliken charges for each atom, length $N_{\text{atom}}$.
    pub charges: Vec<f64>,
}

/// Comprehensive quantum chemistry calculation results container.
#[derive(Debug, Clone, PartialEq)]
pub struct QcResult {
    /// Calculation parameters and metadata.
    pub params: QcParameters,
    /// Molecular structure.
    pub molecule: QcMolecule,
    /// Basis set definitions.
    pub basisset: QcBasisSet,
    /// Atomic orbital info and mapping.
    pub orbital_info: OrbitalInfo,
    /// Total SCF energies across iterations (`/TEs`), length `iterations`.
    pub total_energies: Vec<f64>,
    /// Orbital energies by run type (`/energy_level`).
    pub orbital_energies: Vec<OrbitalEnergy>,
    /// Orbital occupations by run type (`/occ`).
    pub orbital_occupations: Vec<OrbitalOccupation>,
    /// MO coefficient matrices by run type (`/C`).
    pub mo_coefficients: Vec<MoCoefficients>,
    /// Overlap matrix $S$ (`/s`), present in `pdfresults.h5`.
    ///
    /// Accepts lower-triangle-only storage (upper triangle all 0.0) or full symmetric data
    /// verified within [`MATRIX_SYMMETRY_TOLERANCE`] ($10^{-8}$).
    pub s_matrix: Option<SymmetricMatrix>,
    /// Density matrices by run type (`/P`).
    ///
    /// Accepts lower-triangle-only storage or full symmetric data.
    pub density_matrices: Vec<SymmetricMatrix>,
    /// Fock matrices by run type (`/F`).
    ///
    /// Accepts lower-triangle-only storage or full symmetric data.
    pub fock_matrices: Vec<SymmetricMatrix>,
    /// One-electron Hamiltonian $h$ (`/h`).
    ///
    /// Accepts lower-triangle-only storage or full symmetric data.
    pub h_matrix: Option<SymmetricMatrix>,
    /// Two-electron Hamiltonian component $h_2$ (`/h2`).
    ///
    /// Accepts lower-triangle-only storage or full symmetric data.
    pub h2_matrix: Option<SymmetricMatrix>,
    /// Mulliken atomic charges (`/pop/mulliken`).
    pub mulliken_charges: Vec<MullikenPopulation>,
}

impl QcResult {
    /// Returns the MO coefficients for a run type at a given iteration.
    ///
    /// If iteration is `None`, the latest iteration is returned.
    pub fn get_mo_coefficients(
        &self,
        run_type: &RunType,
        iteration: Option<usize>,
    ) -> Option<&MoCoefficients> {
        let matching: Vec<&MoCoefficients> = self
            .mo_coefficients
            .iter()
            .filter(|c| &c.run_type == run_type)
            .collect();
        if let Some(target_itr) = iteration {
            matching.into_iter().find(|c| c.iteration == target_itr)
        } else {
            matching.into_iter().max_by_key(|c| c.iteration)
        }
    }

    /// Returns the orbital energies for a run type at a given iteration.
    ///
    /// If iteration is `None`, the latest iteration is returned.
    pub fn get_orbital_energy(
        &self,
        run_type: &RunType,
        iteration: Option<usize>,
    ) -> Option<&OrbitalEnergy> {
        let matching: Vec<&OrbitalEnergy> = self
            .orbital_energies
            .iter()
            .filter(|e| &e.run_type == run_type)
            .collect();
        if let Some(target_itr) = iteration {
            matching.into_iter().find(|e| e.iteration == target_itr)
        } else {
            matching.into_iter().max_by_key(|e| e.iteration)
        }
    }

    /// Returns the orbital occupation for a run type.
    pub fn get_orbital_occupation(&self, run_type: &RunType) -> Option<&OrbitalOccupation> {
        self.orbital_occupations
            .iter()
            .find(|occ| &occ.run_type == run_type)
    }

    /// Returns the Mulliken atomic charges for a run type.
    pub fn get_mulliken_population(
        &self,
        run_type: &RunType,
        iteration: Option<usize>,
    ) -> Option<&MullikenPopulation> {
        let matching: Vec<&MullikenPopulation> = self
            .mulliken_charges
            .iter()
            .filter(|p| &p.run_type == run_type)
            .collect();
        if let Some(target_itr) = iteration {
            matching.into_iter().find(|p| p.iteration == target_itr)
        } else {
            matching.into_iter().max_by_key(|p| p.iteration)
        }
    }
}
