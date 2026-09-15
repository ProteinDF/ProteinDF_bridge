// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom_group::PyAtomGroup;
use crate::error::to_py_err;
use proteindf_bridge::format::pdb::Pdb as CorePdb;
use pyo3::prelude::*;

#[pyclass(name = "Pdb", module = "proteindf_bridge_rs")]
#[derive(Clone)]
pub struct PyPdb {
    pub(crate) inner: CorePdb,
}

#[pymethods]
impl PyPdb {
    #[new]
    #[pyo3(signature = (file_path=None, mode=None))]
    pub fn new(file_path: Option<&str>, mode: Option<&str>) -> PyResult<Self> {
        let mut inner = CorePdb::new(mode);
        if let Some(path) = file_path {
            inner.load(path).map_err(to_py_err)?;
        }
        Ok(Self { inner })
    }

    pub fn load(&mut self, file_path: &str) -> PyResult<()> {
        self.inner.load(file_path).map_err(to_py_err)
    }

    pub fn save(&self, file_path: &str) -> PyResult<()> {
        self.inner.save(file_path).map_err(to_py_err)
    }

    pub fn get_text(&self) -> String {
        self.inner.get_text()
    }

    pub fn __str__(&self) -> String {
        self.inner.get_text()
    }

    pub fn renumber(&mut self) {
        self.inner.renumber();
    }

    #[getter]
    pub fn mode(&self) -> Option<String> {
        self.inner.mode().map(|s| s.to_string())
    }

    #[setter]
    pub fn set_mode(&mut self, mode: Option<String>) {
        self.inner.set_mode(mode.as_deref());
    }

    #[pyo3(signature = (select_model=None, select_altloc=Some("A")))]
    pub fn get_atomgroup(
        &self,
        select_model: Option<usize>,
        select_altloc: Option<&str>,
    ) -> PyResult<PyAtomGroup> {
        let ag = self
            .inner
            .get_atomgroup(select_model, select_altloc)
            .map_err(to_py_err)?;
        Ok(PyAtomGroup::from_core(ag))
    }

    #[pyo3(signature = (atomgroup, is_charge2tempfactor=false))]
    pub fn set_by_atomgroup(
        &mut self,
        atomgroup: &PyAtomGroup,
        is_charge2tempfactor: bool,
    ) -> PyResult<()> {
        self.inner
            .set_by_atomgroup(&atomgroup.inner, is_charge2tempfactor)
            .map_err(to_py_err)
    }

    pub fn get_modpdb_atomgroup(&self, ag_protein: &PyAtomGroup) -> PyAtomGroup {
        let ag = self.inner.get_modpdb_atomgroup(&ag_protein.inner);
        PyAtomGroup::from_core(ag)
    }
}
