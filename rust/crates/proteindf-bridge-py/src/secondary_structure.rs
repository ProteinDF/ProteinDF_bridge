// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use pyo3::prelude::*;

use proteindf_bridge::secondary_structure::{
    apply_secondary_structure as core_apply_secondary_structure,
    calc_secondary_structure as core_calc_secondary_structure,
    SecondaryStructure as CoreSecondaryStructure,
};

use crate::atom_group::PyAtomGroup;

/// Secondary structure assignment for a residue.
#[pyclass(
    name = "SecondaryStructure",
    module = "proteindf_bridge_rs",
    skip_from_py_object
)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PySecondaryStructure {
    inner: CoreSecondaryStructure,
}

impl From<CoreSecondaryStructure> for PySecondaryStructure {
    fn from(inner: CoreSecondaryStructure) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PySecondaryStructure {
    /// Residue key in the chain (e.g. "1").
    #[getter]
    pub fn residue_key(&self) -> String {
        self.inner.residue_key.clone()
    }

    /// Residue name (e.g. "ALA").
    #[getter]
    pub fn residue_name(&self) -> String {
        self.inner.residue_name.clone()
    }

    /// 3-state secondary structure code: "H" (Helix), "E" (Strand), or "-" (Loop).
    #[getter]
    pub fn code(&self) -> String {
        self.inner.code.to_string()
    }

    fn __repr__(&self) -> String {
        format!(
            "SecondaryStructure(residue_key='{}', residue_name='{}', code='{}')",
            self.inner.residue_key,
            self.inner.residue_name,
            self.inner.code.as_char()
        )
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}

/// Calculates 3-state secondary structure assignments for all residues in a chain.
#[pyfunction]
pub fn calc_secondary_structure(chain: &PyAtomGroup) -> Vec<PySecondaryStructure> {
    core_calc_secondary_structure(&chain.inner)
        .into_iter()
        .map(PySecondaryStructure::from)
        .collect()
}

/// Applies 3-state secondary structure assignments ('H', 'E', '-') directly to the residue groups of the given chain in-place.
///
/// Note: This modifies the AtomGroup in-place. If called on a copy or a sub-tree copy (such as one returned by indexing or get_group()), the original root tree will remain unchanged.
#[pyfunction]
pub fn apply_secondary_structure(mut chain: PyRefMut<'_, PyAtomGroup>) {
    core_apply_secondary_structure(&mut chain.inner);
}
