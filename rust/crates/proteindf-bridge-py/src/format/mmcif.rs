// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom_group::PyAtomGroup;
use crate::error::to_py_err;
use proteindf_bridge::error::BridgeError;
use proteindf_bridge::format::mmcif::SimpleMmcif as CoreSimpleMmcif;
use proteindf_bridge::format::mmcif_writer::{validate_for_mmcif_write, MmcifWriteOptions};
use pyo3::prelude::*;

#[pyclass(
    name = "SimpleMmcif",
    module = "proteindf_bridge_rs",
    skip_from_py_object
)]
#[derive(Clone)]
pub struct PySimpleMmcif {
    pub(crate) inner: CoreSimpleMmcif,
    atomgroup: Option<PyAtomGroup>,
    data_block_name: String,
    charge_to_b_factor: bool,
}

impl Default for PySimpleMmcif {
    fn default() -> Self {
        Self {
            inner: CoreSimpleMmcif::default(),
            atomgroup: None,
            data_block_name: "structure".to_string(),
            charge_to_b_factor: false,
        }
    }
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
        Ok(Self {
            inner,
            atomgroup: None,
            data_block_name: "structure".to_string(),
            charge_to_b_factor: false,
        })
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

    #[pyo3(signature = (
        atomgroup,
        data_block_name = "structure",
        charge_to_b_factor = false,
        is_charge2tempfactor = false
    ))]
    pub fn set_by_atomgroup(
        &mut self,
        atomgroup: &PyAtomGroup,
        data_block_name: &str,
        charge_to_b_factor: bool,
        is_charge2tempfactor: bool,
    ) -> PyResult<()> {
        let opts = MmcifWriteOptions {
            data_block_name: data_block_name.to_string(),
            charge_to_b_factor: charge_to_b_factor || is_charge2tempfactor,
        };
        validate_for_mmcif_write(&atomgroup.inner, &opts).map_err(to_py_err)?;
        self.atomgroup = Some(atomgroup.clone());
        self.data_block_name = data_block_name.to_string();
        self.charge_to_b_factor = charge_to_b_factor || is_charge2tempfactor;
        Ok(())
    }

    #[pyo3(signature = (
        file_path,
        atomgroup = None,
        data_block_name = None,
        charge_to_b_factor = None,
        is_charge2tempfactor = None
    ))]
    pub fn save(
        &self,
        file_path: &str,
        atomgroup: Option<&PyAtomGroup>,
        data_block_name: Option<&str>,
        charge_to_b_factor: Option<bool>,
        is_charge2tempfactor: Option<bool>,
    ) -> PyResult<()> {
        let ag = match (atomgroup, &self.atomgroup) {
            (Some(ag), _) => &ag.inner,
            (None, Some(ag)) => &ag.inner,
            (None, None) => {
                return Err(to_py_err(BridgeError::input_error(
                    "mmCIF",
                    "No AtomGroup specified for save(); call set_by_atomgroup() first or pass atomgroup",
                )));
            }
        };
        let c2b = charge_to_b_factor
            .or(is_charge2tempfactor)
            .unwrap_or(self.charge_to_b_factor);
        let block_name = data_block_name
            .map(|s| s.to_string())
            .unwrap_or_else(|| self.data_block_name.clone());
        let opts = MmcifWriteOptions {
            data_block_name: block_name,
            charge_to_b_factor: c2b,
        };
        CoreSimpleMmcif::save_structure(ag, file_path, &opts).map_err(to_py_err)
    }

    #[staticmethod]
    #[pyo3(signature = (
        atomgroup,
        file_path,
        data_block_name = "structure",
        charge_to_b_factor = false,
        is_charge2tempfactor = false
    ))]
    pub fn save_structure(
        atomgroup: &PyAtomGroup,
        file_path: &str,
        data_block_name: &str,
        charge_to_b_factor: bool,
        is_charge2tempfactor: bool,
    ) -> PyResult<()> {
        let opts = MmcifWriteOptions {
            data_block_name: data_block_name.to_string(),
            charge_to_b_factor: charge_to_b_factor || is_charge2tempfactor,
        };
        CoreSimpleMmcif::save_structure(&atomgroup.inner, file_path, &opts).map_err(to_py_err)
    }

    #[staticmethod]
    #[pyo3(signature = (
        atomgroup,
        data_block_name = "structure",
        charge_to_b_factor = false,
        is_charge2tempfactor = false
    ))]
    pub fn write_structure(
        atomgroup: &PyAtomGroup,
        data_block_name: &str,
        charge_to_b_factor: bool,
        is_charge2tempfactor: bool,
    ) -> PyResult<String> {
        let opts = MmcifWriteOptions {
            data_block_name: data_block_name.to_string(),
            charge_to_b_factor: charge_to_b_factor || is_charge2tempfactor,
        };
        let mut buf = Vec::new();
        CoreSimpleMmcif::write_structure(&atomgroup.inner, &mut buf, &opts).map_err(to_py_err)?;
        String::from_utf8(buf)
            .map_err(|e| to_py_err(BridgeError::input_error("write", e.to_string())))
    }

    pub fn get_text(&self) -> PyResult<String> {
        if let Some(ref ag) = self.atomgroup {
            let opts = MmcifWriteOptions {
                data_block_name: self.data_block_name.clone(),
                charge_to_b_factor: self.charge_to_b_factor,
            };
            let mut buf = Vec::new();
            CoreSimpleMmcif::write_structure(&ag.inner, &mut buf, &opts).map_err(to_py_err)?;
            String::from_utf8(buf)
                .map_err(|e| to_py_err(BridgeError::input_error("write", e.to_string())))
        } else {
            Err(to_py_err(BridgeError::input_error(
                "mmCIF",
                "No AtomGroup specified; call set_by_atomgroup() first",
            )))
        }
    }

    pub fn __str__(&self) -> PyResult<String> {
        self.get_text()
    }
}
