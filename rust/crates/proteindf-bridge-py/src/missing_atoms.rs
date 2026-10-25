// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::collections::BTreeMap;

use proteindf_bridge::missing_atoms::{
    MissingAtomsReport as CoreMissingAtomsReport, ResidueAtomCheck as CoreResidueAtomCheck,
    ResidueStatus as CoreResidueStatus,
};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

/// Per-residue result of missing-atom inspection.
#[pyclass(
    name = "ResidueAtomCheck",
    module = "proteindf_bridge.rs",
    skip_from_py_object
)]
#[derive(Clone, Debug)]
pub struct PyResidueAtomCheck {
    pub(crate) inner: CoreResidueAtomCheck,
}

impl PyResidueAtomCheck {
    pub fn from_core(inner: CoreResidueAtomCheck) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PyResidueAtomCheck {
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

    /// Inspection status (`"complete"`, `"incomplete"`, or `"no_template"`).
    #[getter]
    pub fn status(&self) -> &'static str {
        match self.inner.status {
            CoreResidueStatus::Complete => "complete",
            CoreResidueStatus::Incomplete => "incomplete",
            CoreResidueStatus::NoTemplate => "no_template",
        }
    }

    /// Number of expected atoms (after filtering by element and removing optional atoms).
    #[getter]
    pub fn expected(&self) -> usize {
        self.inner.expected
    }

    /// Number of expected atoms that are present in the residue.
    #[getter]
    pub fn present(&self) -> usize {
        self.inner.present
    }

    /// Missing expected atom names, in template order.
    #[getter]
    pub fn missing(&self) -> Vec<String> {
        self.inner.missing.clone()
    }

    /// Atom names present in the residue that are neither expected nor optional,
    /// in tree order with duplicates deduplicated.
    #[getter]
    pub fn unexpected(&self) -> Vec<String> {
        self.inner.unexpected.clone()
    }

    /// Serializes this residue check result to a JSON-compatible Python dictionary.
    pub fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("path", &self.inner.path)?;
        dict.set_item("comp_id", &self.inner.comp_id)?;
        dict.set_item("status", self.status())?;
        dict.set_item("expected", self.inner.expected)?;
        dict.set_item("present", self.inner.present)?;
        dict.set_item("missing", &self.inner.missing)?;
        dict.set_item("unexpected", &self.inner.unexpected)?;
        Ok(dict)
    }

    pub fn __repr__(&self) -> String {
        format!(
            "ResidueAtomCheck(path={:?}, comp_id={:?}, status={:?}, missing={:?}, unexpected={:?})",
            self.inner.path,
            self.inner.comp_id,
            self.status(),
            self.inner.missing,
            self.inner.unexpected
        )
    }
}

/// Structure-wide report of missing-atom inspection across residues.
#[pyclass(
    name = "MissingAtomsReport",
    module = "proteindf_bridge.rs",
    skip_from_py_object
)]
#[derive(Debug)]
pub struct PyMissingAtomsReport {
    pub(crate) residues: Vec<Py<PyResidueAtomCheck>>,
    pub(crate) incomplete: Vec<Py<PyResidueAtomCheck>>,
    pub(crate) no_template: Vec<Py<PyResidueAtomCheck>>,
    pub(crate) total_missing: usize,
    pub(crate) is_complete: bool,
    pub(crate) no_template_components: BTreeMap<String, usize>,
}

impl PyMissingAtomsReport {
    pub fn new(py: Python<'_>, report: CoreMissingAtomsReport) -> PyResult<Self> {
        let total_missing = report.total_missing();
        let is_complete = report.is_complete();
        let no_template_components = report.no_template_components();

        let mut residues = Vec::with_capacity(report.residues.len());
        let mut incomplete = Vec::new();
        let mut no_template = Vec::new();

        for res in report.residues {
            let status = res.status;
            let py_res = Py::new(py, PyResidueAtomCheck::from_core(res))?;
            if status == CoreResidueStatus::Incomplete {
                incomplete.push(py_res.clone_ref(py));
            } else if status == CoreResidueStatus::NoTemplate {
                no_template.push(py_res.clone_ref(py));
            }
            residues.push(py_res);
        }

        Ok(Self {
            residues,
            incomplete,
            no_template,
            total_missing,
            is_complete,
            no_template_components,
        })
    }
}

#[pymethods]
impl PyMissingAtomsReport {
    /// Inspection results for all visited residues in traversal order.
    #[getter]
    pub fn residues<'py>(&self, py: Python<'py>) -> Bound<'py, PyList> {
        PyList::new(py, &self.residues).expect("build residues list")
    }

    /// All residues with status `"incomplete"`.
    #[getter]
    pub fn incomplete<'py>(&self, py: Python<'py>) -> Bound<'py, PyList> {
        PyList::new(py, &self.incomplete).expect("build incomplete list")
    }

    /// All residues with status `"no_template"`.
    #[getter]
    pub fn no_template<'py>(&self, py: Python<'py>) -> Bound<'py, PyList> {
        PyList::new(py, &self.no_template).expect("build no_template list")
    }

    /// Total number of missing atoms across all residues.
    #[getter]
    pub fn total_missing(&self) -> usize {
        self.total_missing
    }

    /// True if there are no incomplete residues.
    #[getter]
    pub fn is_complete(&self) -> bool {
        self.is_complete
    }

    /// Map from component ID to residue count for all residues without templates.
    #[getter]
    pub fn no_template_components<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        for (k, v) in &self.no_template_components {
            dict.set_item(k, v)?;
        }
        Ok(dict)
    }

    /// Serializes this report to a JSON-compatible Python dictionary.
    pub fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("total_missing", self.total_missing)?;
        dict.set_item("is_complete", self.is_complete)?;

        let comp_dict = PyDict::new(py);
        for (k, v) in &self.no_template_components {
            comp_dict.set_item(k, v)?;
        }
        dict.set_item("no_template_components", comp_dict)?;

        let res_list = PyList::empty(py);
        for r in &self.residues {
            let r_dict = r.borrow(py).to_dict(py)?;
            res_list.append(r_dict)?;
        }
        dict.set_item("residues", res_list)?;

        let inc_list = PyList::empty(py);
        for r in &self.incomplete {
            let r_dict = r.borrow(py).to_dict(py)?;
            inc_list.append(r_dict)?;
        }
        dict.set_item("incomplete", inc_list)?;

        let notemp_list = PyList::empty(py);
        for r in &self.no_template {
            let r_dict = r.borrow(py).to_dict(py)?;
            notemp_list.append(r_dict)?;
        }
        dict.set_item("no_template", notemp_list)?;

        Ok(dict)
    }

    pub fn __str__(&self, py: Python<'_>) -> PyResult<String> {
        let mut lines = Vec::new();
        if self.is_complete && self.no_template_components.is_empty() {
            lines.push("MissingAtomsReport: complete (0 missing atoms)".to_string());
        } else {
            lines.push(format!(
                "MissingAtomsReport: {} missing atom(s) across {} incomplete residue(s), {} missing template component(s)",
                self.total_missing,
                self.incomplete.len(),
                self.no_template_components.len()
            ));
        }

        if !self.incomplete.is_empty() {
            lines.push("Incomplete residues:".to_string());
            for r_py in &self.incomplete {
                let r = r_py.borrow(py);
                lines.push(format!(
                    "  {} ({}): missing {}",
                    r.inner.path,
                    r.inner.comp_id,
                    r.inner.missing.join(", ")
                ));
            }
        }

        if !self.no_template_components.is_empty() {
            lines.push("Components without template:".to_string());
            for (comp, count) in &self.no_template_components {
                lines.push(format!("  {}: {} residue(s)", comp, count));
            }
        }

        Ok(lines.join("\n"))
    }

    pub fn __repr__(&self) -> String {
        format!(
            "MissingAtomsReport(complete={}, total_missing={}, incomplete={}, no_template={})",
            if self.is_complete { "True" } else { "False" },
            self.total_missing,
            self.incomplete.len(),
            self.no_template.len()
        )
    }
}
