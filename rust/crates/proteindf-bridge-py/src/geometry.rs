// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use proteindf_bridge::geometry::{
    GeometryReport as CoreGeometryReport, GeometryViolation as CoreGeometryViolation,
    GeometryViolationKind as CoreGeometryViolationKind,
};
use pyo3::prelude::*;
use pyo3::types::PyList;

/// A single geometry violation (bond length, bond angle or planarity).
#[pyclass(
    name = "GeometryViolation",
    module = "proteindf_bridge.rs",
    skip_from_py_object
)]
#[derive(Clone, Debug)]
pub struct PyGeometryViolation {
    pub(crate) inner: CoreGeometryViolation,
}

#[pymethods]
impl PyGeometryViolation {
    /// Violation kind: `"bond_length"`, `"bond_angle"` or `"planarity"`.
    #[getter]
    pub fn kind(&self) -> &'static str {
        match self.inner.kind {
            CoreGeometryViolationKind::BondLength => "bond_length",
            CoreGeometryViolationKind::BondAngle => "bond_angle",
            CoreGeometryViolationKind::Planarity => "planarity",
        }
    }

    /// Hierarchical path of the residue group.
    #[getter]
    pub fn path(&self) -> String {
        self.inner.path.clone()
    }

    /// Component identifier of the residue (e.g. `"ALA"`).
    #[getter]
    pub fn comp_id(&self) -> String {
        self.inner.comp_id.clone()
    }

    /// Atom names involved in the violation.
    #[getter]
    pub fn atoms(&self) -> Vec<String> {
        self.inner.atoms.clone()
    }

    /// Ideal value (Angstrom for lengths, degrees for angles, 0.0 for planarity).
    #[getter]
    pub fn ideal(&self) -> f64 {
        self.inner.ideal
    }

    /// Observed value (same unit as `ideal`).
    #[getter]
    pub fn observed(&self) -> f64 {
        self.inner.observed
    }

    /// Absolute deviation `|observed - ideal|`.
    #[getter]
    pub fn deviation(&self) -> f64 {
        self.inner.deviation
    }

    pub fn __repr__(&self) -> String {
        format!(
            "GeometryViolation(kind={:?}, path={:?}, comp_id={:?}, atoms={:?}, ideal={}, observed={}, deviation={})",
            self.kind(),
            self.inner.path,
            self.inner.comp_id,
            self.inner.atoms,
            self.inner.ideal,
            self.inner.observed,
            self.inner.deviation
        )
    }
}

/// List of geometry violations found by `AtomGroup.check_geometry`.
#[pyclass(
    name = "GeometryReport",
    module = "proteindf_bridge.rs",
    skip_from_py_object
)]
#[derive(Debug)]
pub struct PyGeometryReport {
    pub(crate) violations: Vec<Py<PyGeometryViolation>>,
}

impl PyGeometryReport {
    pub fn new(py: Python<'_>, report: CoreGeometryReport) -> PyResult<Self> {
        let violations = report
            .violations
            .into_iter()
            .map(|v| Py::new(py, PyGeometryViolation { inner: v }))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(Self { violations })
    }
}

#[pymethods]
impl PyGeometryReport {
    /// All violations in traversal order (bond lengths, bond angles, planarity).
    #[getter]
    pub fn violations<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        PyList::new(py, &self.violations)
    }

    pub fn __len__(&self) -> usize {
        self.violations.len()
    }

    pub fn __bool__(&self) -> bool {
        !self.violations.is_empty()
    }

    pub fn __repr__(&self) -> String {
        format!("GeometryReport(violations={})", self.violations.len())
    }
}
