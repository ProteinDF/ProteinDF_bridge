// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom_group::PyAtomGroup;
use proteindf_bridge::ssbond::SSBond as CoreSSBond;
use pyo3::prelude::*;

#[pyclass(name = "SSBond", module = "proteindf_bridge.rs", skip_from_py_object)]
#[derive(Clone)]
pub struct PySSBond {
    pub(crate) inner: CoreSSBond,
}

#[pymethods]
impl PySSBond {
    #[new]
    pub fn new(py: Python<'_>, model: &PyAtomGroup) -> Self {
        Self {
            inner: CoreSSBond::new(&model.to_core(py)),
        }
    }

    pub fn get_bonds(&mut self) -> Vec<(String, String)> {
        self.inner.get_bonds().to_vec()
    }

    #[staticmethod]
    pub fn find_bonds(py: Python<'_>, model: &PyAtomGroup) -> Vec<(String, String)> {
        CoreSSBond::find_bonds(&model.to_core(py))
    }
}
