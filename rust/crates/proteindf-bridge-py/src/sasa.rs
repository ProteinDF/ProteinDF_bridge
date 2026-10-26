// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use proteindf_bridge::sasa::{
    AtomSasa as CoreAtomSasa, ResidueSasa as CoreResidueSasa, SasaReport as CoreSasaReport,
};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

/// Solvent-accessible surface area for a single heavy atom.
#[pyclass(name = "AtomSasa", module = "proteindf_bridge.rs", skip_from_py_object)]
#[derive(Clone, Debug)]
pub struct PyAtomSasa {
    pub(crate) inner: CoreAtomSasa,
}

impl PyAtomSasa {
    pub fn from_core(inner: CoreAtomSasa) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PyAtomSasa {
    /// Hierarchical path of the atom (e.g. `"/model_1/A/1/CA"`).
    #[getter]
    pub fn path(&self) -> &str {
        &self.inner.path
    }

    /// Atom name (e.g. `"CA"`).
    #[getter]
    pub fn name(&self) -> &str {
        &self.inner.name
    }

    /// Solvent-accessible surface area in Å².
    #[getter]
    pub fn area(&self) -> f64 {
        self.inner.area
    }

    /// Serializes this atom SASA result to a JSON-compatible Python dictionary.
    pub fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("path", &self.inner.path)?;
        dict.set_item("name", &self.inner.name)?;
        dict.set_item("area", self.inner.area)?;
        Ok(dict)
    }

    pub fn __repr__(&self) -> String {
        format!(
            "AtomSasa(path={:?}, name={:?}, area={})",
            self.inner.path, self.inner.name, self.inner.area
        )
    }

    pub fn __str__(&self) -> String {
        self.__repr__()
    }
}

/// Solvent-accessible surface area aggregated for a single residue/component.
#[pyclass(
    name = "ResidueSasa",
    module = "proteindf_bridge.rs",
    skip_from_py_object
)]
#[derive(Clone, Debug)]
pub struct PyResidueSasa {
    pub(crate) inner: CoreResidueSasa,
}

impl PyResidueSasa {
    pub fn from_core(inner: CoreResidueSasa) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PyResidueSasa {
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

    /// Sum of heavy-atom solvent-accessible surface areas in this residue in Å².
    #[getter]
    pub fn area(&self) -> f64 {
        self.inner.area
    }

    /// Serializes this residue SASA result to a JSON-compatible Python dictionary.
    pub fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("path", &self.inner.path)?;
        dict.set_item("comp_id", &self.inner.comp_id)?;
        dict.set_item("area", self.inner.area)?;
        Ok(dict)
    }

    pub fn __repr__(&self) -> String {
        format!(
            "ResidueSasa(path={:?}, comp_id={:?}, area={})",
            self.inner.path, self.inner.comp_id, self.inner.area
        )
    }

    pub fn __str__(&self) -> String {
        self.__repr__()
    }
}

/// Structure-wide report of per-atom and per-residue solvent-accessible surface areas.
#[pyclass(
    name = "SasaReport",
    module = "proteindf_bridge.rs",
    skip_from_py_object
)]
#[derive(Debug)]
pub struct PySasaReport {
    pub(crate) atoms: Vec<Py<PyAtomSasa>>,
    pub(crate) residues: Vec<Py<PyResidueSasa>>,
    pub(crate) total_area: f64,
}

impl PySasaReport {
    pub fn new(py: Python<'_>, report: CoreSasaReport) -> PyResult<Self> {
        let total_area: f64 = report.atoms.iter().map(|a| a.area).sum();

        let mut atoms = Vec::with_capacity(report.atoms.len());
        for atom in report.atoms {
            atoms.push(Py::new(py, PyAtomSasa::from_core(atom))?);
        }

        let mut residues = Vec::with_capacity(report.residues.len());
        for res in report.residues {
            residues.push(Py::new(py, PyResidueSasa::from_core(res))?);
        }

        Ok(Self {
            atoms,
            residues,
            total_area,
        })
    }
}

#[pymethods]
impl PySasaReport {
    /// Per-atom SASA results for all heavy atoms in traversal order.
    #[getter]
    pub fn atoms<'py>(&self, py: Python<'py>) -> Bound<'py, PyList> {
        PyList::new(py, &self.atoms).expect("build atoms list")
    }

    /// Per-residue SASA results for all leaf residue groups in traversal order.
    #[getter]
    pub fn residues<'py>(&self, py: Python<'py>) -> Bound<'py, PyList> {
        PyList::new(py, &self.residues).expect("build residues list")
    }

    /// Sum of heavy-atom solvent-accessible surface areas in Å².
    #[getter]
    pub fn total_area(&self) -> f64 {
        self.total_area
    }

    /// Serializes this SASA report to a JSON-compatible Python dictionary.
    pub fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("total_area", self.total_area)?;

        let atom_list = PyList::empty(py);
        for a in &self.atoms {
            let a_dict = a.borrow(py).to_dict(py)?;
            atom_list.append(a_dict)?;
        }
        dict.set_item("atoms", atom_list)?;

        let res_list = PyList::empty(py);
        for r in &self.residues {
            let r_dict = r.borrow(py).to_dict(py)?;
            res_list.append(r_dict)?;
        }
        dict.set_item("residues", res_list)?;

        Ok(dict)
    }

    pub fn __str__(&self) -> String {
        format!(
            "SasaReport(residues={}, atoms={}, total_area={:.2})",
            self.residues.len(),
            self.atoms.len(),
            self.total_area
        )
    }

    pub fn __repr__(&self) -> String {
        format!(
            "SasaReport(residues={}, atoms={}, total_area={})",
            self.residues.len(),
            self.atoms.len(),
            self.total_area
        )
    }
}
