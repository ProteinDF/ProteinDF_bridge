// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

//! Quantum chemistry calculation result I/O and data structures.

pub mod population;
pub mod types;

#[cfg(feature = "hdf5")]
pub mod pdfparam_hdf5;

pub use crate::matrix::{Matrix, SymmetricMatrix};
#[cfg(feature = "hdf5")]
pub use pdfparam_hdf5::PdfParamHdf5;
pub use population::{calc_orbital_mulliken_contributions, OrbitalContribution};
pub use types::{
    AtomicOrbital, DensityMatrix, FockMatrix, MoCoefficients, MullikenPopulation, OrbitalEnergy,
    OrbitalInfo, OrbitalOccupation, QcAtom, QcBasisSet, QcBasisSetDefinition, QcCgto, QcMethod,
    QcMolecule, QcParameters, QcPgto, QcResult, RunType, MATRIX_SYMMETRY_TOLERANCE,
};
