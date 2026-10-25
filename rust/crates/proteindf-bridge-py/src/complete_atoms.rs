// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use proteindf_bridge::complete_atoms::{
    CompletionReport as CoreCompletionReport, CompletionStatus as CoreCompletionStatus,
    ResidueCompletion as CoreResidueCompletion,
};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

/// Per-residue result of missing heavy-atom completion.
#[pyclass(
    name = "ResidueCompletion",
    module = "proteindf_bridge.rs",
    skip_from_py_object
)]
#[derive(Clone, Debug)]
pub struct PyResidueCompletion {
    pub(crate) inner: CoreResidueCompletion,
}

impl PyResidueCompletion {
    pub fn from_core(inner: CoreResidueCompletion) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PyResidueCompletion {
    /// Hierarchical path of the residue group (e.g. `"/model_1/A/1/"`).
    #[getter]
    pub fn path(&self) -> &str {
        &self.inner.path
    }

    /// Chemical component identifier (residue group's `name`, e.g. `"ALA"`).
    #[getter]
    pub fn comp_id(&self) -> &str {
        &self.inner.comp_id
    }

    /// Completion outcome status (`"completed"`, `"unchanged"`, `"no_template"`,
    /// `"insufficient_anchors"`, or `"poor_fit"`).
    #[getter]
    pub fn status(&self) -> &'static str {
        match self.inner.status {
            CoreCompletionStatus::Completed => "completed",
            CoreCompletionStatus::Unchanged => "unchanged",
            CoreCompletionStatus::NoTemplate => "no_template",
            CoreCompletionStatus::InsufficientAnchors => "insufficient_anchors",
            CoreCompletionStatus::PoorFit => "poor_fit",
        }
    }

    /// Names of heavy atoms added to the residue, in template order.
    #[getter]
    pub fn added(&self) -> Vec<String> {
        self.inner.added.clone()
    }

    /// Number of anchor heavy atoms used for superposition.
    #[getter]
    pub fn anchors(&self) -> usize {
        self.inner.anchors
    }

    /// Anchor superposition RMSD in Angstroms, or `None` if superposition was not computed.
    #[getter]
    pub fn rmsd(&self) -> Option<f64> {
        self.inner.rmsd
    }

    /// Serializes this residue completion result to a JSON-compatible Python dictionary.
    pub fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("path", &self.inner.path)?;
        dict.set_item("comp_id", &self.inner.comp_id)?;
        dict.set_item("status", self.status())?;
        dict.set_item("added", &self.inner.added)?;
        dict.set_item("anchors", self.inner.anchors)?;
        dict.set_item("rmsd", self.inner.rmsd)?;
        Ok(dict)
    }

    pub fn __repr__(&self) -> String {
        format!(
            "ResidueCompletion(path={:?}, comp_id={:?}, status={:?}, added={:?}, anchors={}, rmsd={:?})",
            self.inner.path,
            self.inner.comp_id,
            self.status(),
            self.inner.added,
            self.inner.anchors,
            self.inner.rmsd,
        )
    }
}

/// Structure-wide report of missing heavy-atom completion across residues.
#[pyclass(
    name = "CompletionReport",
    module = "proteindf_bridge.rs",
    skip_from_py_object
)]
#[derive(Debug)]
pub struct PyCompletionReport {
    pub(crate) residues: Vec<Py<PyResidueCompletion>>,
    pub(crate) completed: Vec<Py<PyResidueCompletion>>,
    pub(crate) skipped: Vec<Py<PyResidueCompletion>>,
    pub(crate) total_added: usize,
}

impl PyCompletionReport {
    pub fn new(py: Python<'_>, report: CoreCompletionReport) -> PyResult<Self> {
        let total_added = report.total_added();
        let mut residues = Vec::with_capacity(report.residues.len());
        let mut completed = Vec::new();
        let mut skipped = Vec::new();

        for res in report.residues {
            let status = res.status;
            let py_res = Py::new(py, PyResidueCompletion::from_core(res))?;
            match status {
                CoreCompletionStatus::Completed => {
                    completed.push(py_res.clone_ref(py));
                }
                CoreCompletionStatus::NoTemplate
                | CoreCompletionStatus::InsufficientAnchors
                | CoreCompletionStatus::PoorFit => {
                    skipped.push(py_res.clone_ref(py));
                }
                CoreCompletionStatus::Unchanged => {}
            }
            residues.push(py_res);
        }

        Ok(Self {
            residues,
            completed,
            skipped,
            total_added,
        })
    }
}

#[pymethods]
impl PyCompletionReport {
    /// Completion results for all visited residues in traversal order.
    #[getter]
    pub fn residues<'py>(&self, py: Python<'py>) -> Bound<'py, PyList> {
        PyList::new(py, &self.residues).expect("build residues list")
    }

    /// All residues with status `"completed"`.
    #[getter]
    pub fn completed<'py>(&self, py: Python<'py>) -> Bound<'py, PyList> {
        PyList::new(py, &self.completed).expect("build completed list")
    }

    /// All skipped residues (`"no_template"`, `"insufficient_anchors"`, or `"poor_fit"`).
    #[getter]
    pub fn skipped<'py>(&self, py: Python<'py>) -> Bound<'py, PyList> {
        PyList::new(py, &self.skipped).expect("build skipped list")
    }

    /// Total number of heavy atoms added across all residues.
    #[getter]
    pub fn total_added(&self) -> usize {
        self.total_added
    }

    /// Serializes this report to a JSON-compatible Python dictionary.
    pub fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("total_added", self.total_added)?;

        let res_list = PyList::empty(py);
        for r in &self.residues {
            res_list.append(r.borrow(py).to_dict(py)?)?;
        }
        dict.set_item("residues", res_list)?;

        let comp_list = PyList::empty(py);
        for r in &self.completed {
            comp_list.append(r.borrow(py).to_dict(py)?)?;
        }
        dict.set_item("completed", comp_list)?;

        let skip_list = PyList::empty(py);
        for r in &self.skipped {
            skip_list.append(r.borrow(py).to_dict(py)?)?;
        }
        dict.set_item("skipped", skip_list)?;

        Ok(dict)
    }

    pub fn __str__(&self, py: Python<'_>) -> PyResult<String> {
        let mut lines = Vec::new();
        lines.push(format!(
            "CompletionReport: {} heavy atom(s) added across {} completed residue(s), {} skipped residue(s)",
            self.total_added,
            self.completed.len(),
            self.skipped.len()
        ));

        if !self.completed.is_empty() {
            lines.push("Completed residues:".to_string());
            for r_py in &self.completed {
                let r = r_py.borrow(py);
                lines.push(format!(
                    "  {} ({}): added {:?}",
                    r.inner.path, r.inner.comp_id, r.inner.added
                ));
            }
        }

        if !self.skipped.is_empty() {
            lines.push("Skipped residues:".to_string());
            for r_py in &self.skipped {
                let r = r_py.borrow(py);
                lines.push(format!(
                    "  {} ({}): status {}",
                    r.inner.path,
                    r.inner.comp_id,
                    r.status()
                ));
            }
        }

        Ok(lines.join("\n"))
    }

    pub fn __repr__(&self) -> String {
        format!(
            "CompletionReport(total_added={}, completed={}, skipped={}, residues={})",
            self.total_added,
            self.completed.len(),
            self.skipped.len(),
            self.residues.len()
        )
    }
}
