// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom_group::PyAtomGroup;
use crate::error::to_py_err;
use crate::matrix::PySymmetricMatrix;
use proteindf_bridge::bond::Bond as CoreBond;
use pyo3::prelude::*;

/// Heuristic bond detector based on covalent radii.
///
/// Bonds are established between pairs of atoms (p, q) satisfying:
///
///     distance(p, q) <= covalent_radius(p) + covalent_radius(q) + 0.45 Å
///
/// where 0.45 Å is the covalent bond tolerance.
///
/// Key behaviors:
///     - No bonds across models: When the structure contains multiple models
///       (satisfying Format.is_models), bonds are not formed between atoms
///       belonging to different models.
///     - Preservation of existing bonds: Pre-existing bonds and their bond
///       orders in the AtomGroup are preserved and not duplicated.
///     - Dense matrix threshold: If the number of atoms exceeds 2000,
///       the pairwise matrices (distmat and bondmat) are None to prevent
///       excessive memory usage, while bond connectivity is still recorded
///       in the AtomGroup in O(N) time.
#[pyclass(name = "Bond", module = "proteindf_bridge.rs")]
#[derive(Default)]
pub struct PyBond {
    pub(crate) inner: CoreBond,
}

#[pymethods]
impl PyBond {
    /// Initializes a new Bond detector instance.
    #[new]
    pub fn new() -> Self {
        Self {
            inner: CoreBond::new(),
        }
    }

    /// Detects and sets up covalent bonds for the given molecular structure.
    ///
    /// Computes interatomic distances and adds detected covalent bonds to `mol`.
    /// Existing bonds and their bond orders are preserved without duplicate
    /// registration. If `mol` represents multiple models, no bonds are formed
    /// across different models.
    ///
    /// Args:
    ///     mol (AtomGroup): Molecular structure to analyze and update.
    ///
    /// Returns:
    ///     None
    ///
    /// Raises:
    ///     TypeError: If `mol` is not an AtomGroup instance.
    ///     BrInputError: If calculating covalent radii fails (e.g., unknown element).
    pub fn setup(&mut self, py: Python<'_>, mol: &mut PyAtomGroup) -> PyResult<()> {
        let mut core = mol.to_core(py);
        self.inner.setup_heuristic(&mut core).map_err(to_py_err)?;
        mol.update_from_core(py, core)?;
        Ok(())
    }

    /// Pairwise atomic distance matrix, or None if atom count exceeds 2000.
    ///
    /// Returns a SymmetricMatrix when the number of atoms is 2000 or fewer,
    /// or None if there are more than 2000 atoms to protect against O(N^2)
    /// memory consumption.
    #[getter]
    pub fn distmat(&self) -> Option<PySymmetricMatrix> {
        self.inner
            .distmat
            .as_ref()
            .map(|m| PySymmetricMatrix::from_core(m.clone()))
    }

    /// Pairwise bond connectivity matrix, or None if atom count exceeds 2000.
    ///
    /// Returns a SymmetricMatrix (1.0 for bonded, 0.0 otherwise) when the number
    /// of atoms is 2000 or fewer, or None if there are more than 2000 atoms to
    /// protect against O(N^2) memory consumption.
    #[getter]
    pub fn bondmat(&self) -> Option<PySymmetricMatrix> {
        self.inner
            .bondmat
            .as_ref()
            .map(|m| PySymmetricMatrix::from_core(m.clone()))
    }
}
