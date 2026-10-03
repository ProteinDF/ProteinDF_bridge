// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom_group::PyAtomGroup;
use crate::error::to_py_err;
use crate::matrix::PySymmetricMatrix;
use proteindf_bridge::bond::Bond as CoreBond;
use pyo3::prelude::*;

#[pyclass(name = "Bond", module = "proteindf_bridge_rs")]
#[derive(Default)]
pub struct PyBond {
    pub(crate) inner: CoreBond,
}

#[pymethods]
impl PyBond {
    #[new]
    pub fn new() -> Self {
        Self {
            inner: CoreBond::new(),
        }
    }

    pub fn setup(&mut self, mol: &mut PyAtomGroup) -> PyResult<()> {
        self.inner
            .setup_heuristic(&mut mol.inner)
            .map_err(to_py_err)
    }

    #[getter]
    pub fn distmat(&self) -> Option<PySymmetricMatrix> {
        self.inner
            .distmat
            .as_ref()
            .map(|m| PySymmetricMatrix::from_core(m.clone()))
    }

    #[getter]
    pub fn bondmat(&self) -> Option<PySymmetricMatrix> {
        self.inner
            .bondmat
            .as_ref()
            .map(|m| PySymmetricMatrix::from_core(m.clone()))
    }
}
