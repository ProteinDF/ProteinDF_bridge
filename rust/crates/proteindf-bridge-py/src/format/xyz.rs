// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom_group::PyAtomGroup;
use crate::error::to_py_err;
use proteindf_bridge::format::xyz::Xyz as CoreXyz;
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

#[pyclass(name = "Xyz", module = "proteindf_bridge_rs")]
#[derive(Clone, Default)]
pub struct PyXyz {
    pub(crate) inner: CoreXyz,
}

#[pymethods]
impl PyXyz {
    #[new]
    #[pyo3(signature = (*args))]
    pub fn new(args: &Bound<'_, pyo3::types::PyTuple>) -> PyResult<Self> {
        let mut xyz = CoreXyz::new();
        if args.len() == 1 {
            let first = args.get_item(0)?;
            if let Ok(path) = first.extract::<String>() {
                xyz.load(&path).map_err(to_py_err)?;
            } else if let Ok(ag) = first.extract::<PyRef<PyAtomGroup>>() {
                xyz.set_by_atomgroup(&ag.inner);
            } else {
                return Err(PyTypeError::new_err(
                    "Xyz.__init__: expected filepath (str) or AtomGroup",
                ));
            }
        } else if args.len() > 1 {
            return Err(PyTypeError::new_err(
                "Xyz.__init__: illegal number of arguments",
            ));
        }
        Ok(Self { inner: xyz })
    }

    #[getter]
    pub fn comment(&self) -> String {
        self.inner.comment().to_string()
    }

    #[setter]
    pub fn set_comment(&mut self, val: String) {
        self.inner.set_comment(val);
    }

    pub fn load(&mut self, file_path: &str) -> PyResult<()> {
        self.inner.load(file_path).map_err(to_py_err)
    }

    pub fn save(&self, file_path: &str) -> PyResult<()> {
        self.inner.save(file_path).map_err(to_py_err)
    }

    pub fn get_atomgroup(&self) -> PyResult<PyAtomGroup> {
        let ag = self.inner.get_atom_group().map_err(to_py_err)?;
        Ok(PyAtomGroup::from_core(ag))
    }

    pub fn get_atom_group(&self) -> PyResult<PyAtomGroup> {
        self.get_atomgroup()
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
