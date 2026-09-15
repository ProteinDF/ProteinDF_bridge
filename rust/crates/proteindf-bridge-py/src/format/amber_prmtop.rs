// Copyright (C) 2014 The ProteinDF development team.
// see also AUTHORS and README if provided.
//
// This file is a part of the ProteinDF software package.

use crate::atom_group::PyAtomGroup;
use crate::error::to_py_err;
use crate::position::PyPosition;
use proteindf_bridge::format::amber_prmtop::AmberPrmtop as CoreAmberPrmtop;
use pyo3::prelude::*;

#[pyclass(name = "AmberPrmtop", module = "proteindf_bridge_rs")]
#[derive(Clone, Default)]
pub struct PyAmberPrmtop {
    pub(crate) inner: CoreAmberPrmtop,
}

#[pymethods]
impl PyAmberPrmtop {
    #[new]
    #[pyo3(signature = (prmtop_path=None, inpcrd_path=None))]
    pub fn new(prmtop_path: Option<&str>, inpcrd_path: Option<&str>) -> PyResult<Self> {
        let mut inner = CoreAmberPrmtop::new();
        if let (Some(top), Some(crd)) = (prmtop_path, inpcrd_path) {
            inner.load(top, crd).map_err(to_py_err)?;
        }
        Ok(Self { inner })
    }

    pub fn load(&mut self, prmtop_path: &str, inpcrd_path: &str) -> PyResult<()> {
        self.inner.load(prmtop_path, inpcrd_path).map_err(to_py_err)
    }

    #[getter]
    pub fn atom_names(&self) -> Vec<String> {
        self.inner.atom_names().to_vec()
    }

    #[getter]
    pub fn charges(&self) -> Vec<f64> {
        self.inner.charges().to_vec()
    }

    #[getter]
    pub fn atomic_numbers(&self) -> Vec<usize> {
        self.inner.atomic_numbers().to_vec()
    }

    #[getter]
    pub fn xyz(&self) -> Vec<PyPosition> {
        self.inner
            .xyz()
            .iter()
            .map(|p| PyPosition::from_core(*p))
            .collect()
    }

    pub fn get_atomgroup(&self) -> PyResult<PyAtomGroup> {
        let ag = self.inner.get_atomgroup().map_err(to_py_err)?;
        Ok(PyAtomGroup::from_core(ag))
    }
}
