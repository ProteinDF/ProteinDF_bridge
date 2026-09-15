// Copyright (C) 2014 The ProteinDF development team.
// see also AUTHORS and README if provided.
//
// This file is a part of the ProteinDF software package.

use crate::atom_group::PyAtomGroup;
use crate::error::to_py_err;
use proteindf_bridge::format::mol2::SimpleMol2 as CoreSimpleMol2;
use pyo3::prelude::*;

#[pyclass(name = "SimpleMol2", module = "proteindf_bridge_rs")]
#[derive(Clone, Default)]
pub struct PySimpleMol2 {
    pub(crate) inner: CoreSimpleMol2,
}

#[pymethods]
impl PySimpleMol2 {
    #[new]
    #[pyo3(signature = (atomgroup=None))]
    pub fn new(atomgroup: Option<&PyAtomGroup>) -> PyResult<Self> {
        let mut mol2 = CoreSimpleMol2::new();
        if let Some(ag) = atomgroup {
            mol2.set_by_atomgroup(&ag.inner);
        }
        Ok(Self { inner: mol2 })
    }

    pub fn set_by_atomgroup(&mut self, atomgroup: &PyAtomGroup) {
        self.inner.set_by_atomgroup(&atomgroup.inner);
    }

    pub fn save(&self, file_path: &str) -> PyResult<()> {
        self.inner.save(file_path).map_err(to_py_err)
    }

    pub fn get_text(&self) -> PyResult<String> {
        self.inner.get_text().map_err(to_py_err)
    }

    pub fn __str__(&self) -> PyResult<String> {
        self.inner.get_text().map_err(to_py_err)
    }
}
