// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use proteindf_bridge::hydrogenation::HydrogenationReport as CoreHydrogenationReport;
use proteindf_bridge::orchestrator::OverallHydrogenationReport as CoreOverallHydrogenationReport;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

/// Per-residue or per-component report of hydrogen addition and removal.
#[pyclass(
    name = "HydrogenationReport",
    module = "proteindf_bridge_rs",
    skip_from_py_object
)]
#[derive(Clone, Debug)]
pub struct PyHydrogenationReport {
    pub(crate) inner: CoreHydrogenationReport,
}

impl PyHydrogenationReport {
    pub fn from_core(report: CoreHydrogenationReport) -> Self {
        Self { inner: report }
    }
}

#[pymethods]
impl PyHydrogenationReport {
    /// Number of hydrogen atoms added to this component.
    #[getter]
    pub fn added_hydrogens(&self) -> usize {
        self.inner.added_hydrogens
    }

    /// Names of the hydrogen atoms added to this component.
    #[getter]
    pub fn added_atom_names(&self) -> Vec<String> {
        self.inner.added_atom_names.clone()
    }

    /// Number of hydrogen (or spurious) atoms removed from this component.
    #[getter]
    pub fn removed_hydrogens(&self) -> usize {
        self.inner.removed_hydrogens
    }

    /// Names of the hydrogen (or spurious) atoms removed from this component.
    #[getter]
    pub fn removed_atom_names(&self) -> Vec<String> {
        self.inner.removed_atom_names.clone()
    }

    pub fn __repr__(&self) -> String {
        format!(
            "HydrogenationReport(added={}, removed={})",
            self.inner.added_hydrogens, self.inner.removed_hydrogens
        )
    }
}

/// Overall summary report of the hydrogenation process across an AtomGroup structure.
#[pyclass(
    name = "OverallHydrogenationReport",
    module = "proteindf_bridge_rs",
    skip_from_py_object
)]
pub struct PyOverallHydrogenationReport {
    total_added_hydrogens: usize,
    total_removed_hydrogens: usize,
    hydrogenated_residues: usize,
    skipped_residues: Py<PyList>,
    step_errors: Py<PyList>,
    residue_reports: Py<PyDict>,
}

impl PyOverallHydrogenationReport {
    pub fn new(py: Python<'_>, report: CoreOverallHydrogenationReport) -> PyResult<Self> {
        let skipped_list = PyList::new(py, report.skipped_residues)?;
        let errors_list = PyList::new(py, report.step_errors)?;

        let reports_dict = PyDict::new(py);
        for (path, res_rep) in report.residue_reports {
            let py_rep = Py::new(py, PyHydrogenationReport::from_core(res_rep))?;
            reports_dict.set_item(path, py_rep)?;
        }

        Ok(Self {
            total_added_hydrogens: report.total_added_hydrogens,
            total_removed_hydrogens: report.total_removed_hydrogens,
            hydrogenated_residues: report.hydrogenated_residues,
            skipped_residues: skipped_list.unbind(),
            step_errors: errors_list.unbind(),
            residue_reports: reports_dict.unbind(),
        })
    }
}

#[pymethods]
impl PyOverallHydrogenationReport {
    /// Total number of hydrogens added across all residues/components.
    #[getter]
    pub fn total_added_hydrogens(&self) -> usize {
        self.total_added_hydrogens
    }

    /// Total number of hydrogens (or spurious atoms) removed across all residues/components.
    #[getter]
    pub fn total_removed_hydrogens(&self) -> usize {
        self.total_removed_hydrogens
    }

    /// Number of residues/components that were modified (at least one hydrogen added or removed).
    #[getter]
    pub fn hydrogenated_residues(&self) -> usize {
        self.hydrogenated_residues
    }

    /// Residues/components where NO changes were made (completely skipped due to missing CCD
    /// template, insufficient heavy atoms such as HOH water, or geometric failure),
    /// returned as a list of `(path, reason)` tuples.
    #[getter]
    pub fn skipped_residues(&self, py: Python<'_>) -> Py<PyList> {
        self.skipped_residues.clone_ref(py)
    }

    /// Processing failures, geometric errors, or missing sidechain templates on partially modified residues,
    /// returned as a list of `(path, error_message)` tuples.
    #[getter]
    pub fn step_errors(&self, py: Python<'_>) -> Py<PyList> {
        self.step_errors.clone_ref(py)
    }

    /// Detailed per-residue hydrogenation reports, keyed by residue path.
    #[getter]
    pub fn residue_reports(&self, py: Python<'_>) -> Py<PyDict> {
        self.residue_reports.clone_ref(py)
    }

    pub fn __repr__(&self, py: Python<'_>) -> String {
        let skipped_len = self.skipped_residues.bind(py).len();
        let errors_len = self.step_errors.bind(py).len();
        format!(
            "OverallHydrogenationReport(added={}, removed={}, residues={}, skipped={}, errors={})",
            self.total_added_hydrogens,
            self.total_removed_hydrogens,
            self.hydrogenated_residues,
            skipped_len,
            errors_len
        )
    }
}
