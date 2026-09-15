// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom_group::PyAtomGroup;
use crate::error::to_py_err;
use proteindf_bridge::error::BridgeError;
use proteindf_bridge::format::gro::SimpleGro as CoreSimpleGro;
use pyo3::prelude::*;
use std::fs;

#[pyclass(name = "SimpleGro", module = "proteindf_bridge_rs")]
#[derive(Clone, Default)]
pub struct PySimpleGro {
    pub(crate) inner: CoreSimpleGro,
}

#[pymethods]
impl PySimpleGro {
    #[new]
    #[pyo3(signature = (file_path=None))]
    pub fn new(file_path: Option<&str>) -> PyResult<Self> {
        let mut inner = CoreSimpleGro::new();
        if let Some(path) = file_path {
            inner.load(path).map_err(to_py_err)?;
        }
        Ok(Self { inner })
    }

    #[getter]
    pub fn title(&self) -> String {
        self.inner.title.clone()
    }

    #[setter]
    pub fn set_title(&mut self, val: String) {
        self.inner.title = val;
    }

    #[getter]
    pub fn num_of_atoms(&self) -> usize {
        self.inner.num_of_atoms
    }

    #[getter]
    pub fn box_vectors(&self) -> Vec<f64> {
        self.inner.box_vectors.clone()
    }

    pub fn load(&mut self, file_path: &str) -> PyResult<()> {
        self.inner.load(file_path).map_err(to_py_err)
    }

    pub fn save(&self, file_path: &str) -> PyResult<()> {
        fs::write(file_path, self.inner.get_text()).map_err(|e| {
            to_py_err(BridgeError::general(format!(
                "failed to write GRO to {}: {}",
                file_path, e
            )))
        })
    }

    pub fn get_atomgroup(&self) -> PyResult<PyAtomGroup> {
        let ag = self.inner.get_atomgroup().map_err(to_py_err)?;
        Ok(PyAtomGroup::from_core(ag))
    }

    pub fn set_by_atomgroup(&mut self, atomgroup: &PyAtomGroup) {
        self.inner.set_by_atomgroup(&atomgroup.inner);
    }

    pub fn get_text(&self) -> String {
        self.inner.get_text()
    }

    pub fn __str__(&self) -> String {
        self.inner.get_text()
    }
}
