// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use proteindf_bridge::atom_group::SchemaViolation as CoreSchemaViolation;
use pyo3::prelude::*;

/// A violation of the standard protein schema (`/model_N/chain_id/res_key/atom_key`).
#[pyclass(
    name = "SchemaViolation",
    module = "proteindf_bridge_rs",
    skip_from_py_object
)]
#[derive(Clone)]
pub struct PySchemaViolation {
    pub(crate) inner: CoreSchemaViolation,
}

impl From<CoreSchemaViolation> for PySchemaViolation {
    fn from(inner: CoreSchemaViolation) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PySchemaViolation {
    /// Returns the violation kind as a string:
    /// - `"DirectAtomsAtNonResidueLevel"`
    /// - `"SubgroupsInResidue"`
    /// - `"ExcessiveDepth"`
    #[getter]
    pub fn violation_type(&self) -> String {
        match &self.inner {
            CoreSchemaViolation::DirectAtomsAtNonResidueLevel { .. } => {
                "DirectAtomsAtNonResidueLevel".to_string()
            }
            CoreSchemaViolation::SubgroupsInResidue { .. } => "SubgroupsInResidue".to_string(),
            CoreSchemaViolation::ExcessiveDepth { .. } => "ExcessiveDepth".to_string(),
        }
    }

    /// Returns the hierarchy path where the violation was detected.
    #[getter]
    pub fn path(&self) -> String {
        match &self.inner {
            CoreSchemaViolation::DirectAtomsAtNonResidueLevel { path, .. } => path.clone(),
            CoreSchemaViolation::SubgroupsInResidue { path, .. } => path.clone(),
            CoreSchemaViolation::ExcessiveDepth { path, .. } => path.clone(),
        }
    }

    /// Returns the path depth if applicable.
    #[getter]
    pub fn depth(&self) -> Option<usize> {
        match &self.inner {
            CoreSchemaViolation::DirectAtomsAtNonResidueLevel { depth, .. } => Some(*depth),
            CoreSchemaViolation::SubgroupsInResidue { .. } => None,
            CoreSchemaViolation::ExcessiveDepth { depth, .. } => Some(*depth),
        }
    }

    /// Returns the atom keys if direct atoms were found at non-residue level.
    #[getter]
    pub fn atom_keys(&self) -> Option<Vec<String>> {
        match &self.inner {
            CoreSchemaViolation::DirectAtomsAtNonResidueLevel { atom_keys, .. } => {
                Some(atom_keys.clone())
            }
            _ => None,
        }
    }

    /// Returns the subgroup keys if subgroups were found inside a residue group.
    #[getter]
    pub fn group_keys(&self) -> Option<Vec<String>> {
        match &self.inner {
            CoreSchemaViolation::SubgroupsInResidue { group_keys, .. } => Some(group_keys.clone()),
            _ => None,
        }
    }

    /// Returns a human-readable explanation of the schema violation.
    #[getter]
    pub fn description(&self) -> String {
        format!("{}", self.inner)
    }

    pub fn __str__(&self) -> String {
        format!("{}", self.inner)
    }

    pub fn __repr__(&self) -> String {
        format!(
            "SchemaViolation(type='{}', path='{}', description='{}')",
            self.violation_type(),
            self.path(),
            self.inner
        )
    }
}
