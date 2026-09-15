// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom_group::PyAtomGroup;
use proteindf_bridge::ssbond::SSBond as CoreSSBond;
use pyo3::prelude::*;

#[pyclass(name = "SSBond", module = "proteindf_bridge_rs")]
#[derive(Clone)]
pub struct PySSBond {
    pub(crate) inner: CoreSSBond,
}

#[pymethods]
impl PySSBond {
    #[new]
    pub fn new(model: &PyAtomGroup) -> Self {
        Self {
            inner: CoreSSBond::new(&model.inner),
        }
    }

    pub fn get_bonds(&mut self) -> Vec<(String, String)> {
        self.inner.get_bonds().to_vec()
    }

    #[staticmethod]
    pub fn find_bonds(model: &PyAtomGroup) -> Vec<(String, String)> {
        CoreSSBond::find_bonds(&model.inner)
    }
}
