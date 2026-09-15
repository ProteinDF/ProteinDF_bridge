// Copyright (C) 2014 The ProteinDF development team.
// see also AUTHORS and README if provided.
//
// This file is a part of the ProteinDF software package.

use crate::atom_group::PyAtomGroup;
use crate::error::to_py_err;
use proteindf_bridge::error::BridgeError;
use proteindf_bridge::format::mmcif::SimpleMmcif as CoreSimpleMmcif;
use pyo3::prelude::*;

#[pyclass(name = "SimpleMmcif", module = "proteindf_bridge_rs")]
#[derive(Clone, Default)]
pub struct PySimpleMmcif {
    pub(crate) inner: CoreSimpleMmcif,
}

#[pymethods]
impl PySimpleMmcif {
    #[new]
    #[pyo3(signature = (file_path=None))]
    pub fn new(file_path: Option<&str>) -> PyResult<Self> {
        let mut inner = CoreSimpleMmcif::new();
        if let Some(path) = file_path {
            inner.load(path).map_err(to_py_err)?;
        }
        Ok(Self { inner })
    }

    pub fn load(&mut self, file_path: &str) -> PyResult<()> {
        self.inner.load(file_path).map_err(to_py_err)
    }

    pub fn get_molecule_names(&self) -> Vec<String> {
        self.inner.get_molecule_names()
    }

    #[pyo3(signature = (name=None))]
    pub fn get_atomgroup(&self, name: Option<&str>) -> PyResult<PyAtomGroup> {
        let ag = match name {
            Some(n) => self.inner.get_atomgroup(n).map_err(to_py_err)?,
            None => {
                let names = self.inner.get_molecule_names();
                let first = names.first().ok_or_else(|| {
                    to_py_err(BridgeError::input_error(
                        "mmCIF",
                        "No data blocks found in file",
                    ))
                })?;
                self.inner.get_atomgroup(first).map_err(to_py_err)?
            }
        };
        Ok(PyAtomGroup::from_core(ag))
    }

    #[pyo3(signature = (select_model=None, select_altloc=Some("A")))]
    pub fn get_structure_atomgroup(
        &self,
        select_model: Option<usize>,
        select_altloc: Option<&str>,
    ) -> PyResult<PyAtomGroup> {
        let ag = self
            .inner
            .get_structure_atomgroup(select_model, select_altloc)
            .map_err(to_py_err)?;
        Ok(PyAtomGroup::from_core(ag))
    }

    #[pyo3(signature = (block_name, select_model=None, select_altloc=Some("A")))]
    pub fn get_structure_atomgroup_for_block(
        &self,
        block_name: &str,
        select_model: Option<usize>,
        select_altloc: Option<&str>,
    ) -> PyResult<PyAtomGroup> {
        let ag = self
            .inner
            .get_structure_atomgroup_for_block(block_name, select_model, select_altloc)
            .map_err(to_py_err)?;
        Ok(PyAtomGroup::from_core(ag))
    }
}
