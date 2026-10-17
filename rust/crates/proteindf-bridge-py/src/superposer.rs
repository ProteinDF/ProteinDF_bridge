// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom_group::PyAtomGroup;
use crate::error::to_py_err;
use crate::matrix::PyMatrix;
use crate::position::PyPosition;
use proteindf_bridge::superposer::Superposer as CoreSuperposer;
use pyo3::prelude::*;

#[pyclass(
    name = "Superposer",
    module = "proteindf_bridge.rs",
    skip_from_py_object
)]
#[derive(Clone)]
pub struct PySuperposer {
    pub(crate) inner: CoreSuperposer,
}

#[pymethods]
impl PySuperposer {
    #[new]
    pub fn new(
        py: Python<'_>,
        atom_group1: &PyAtomGroup,
        atom_group2: &PyAtomGroup,
    ) -> PyResult<Self> {
        let inner = CoreSuperposer::new(&atom_group1.to_core(py), &atom_group2.to_core(py))
            .map_err(to_py_err)?;
        Ok(Self { inner })
    }

    #[getter]
    pub fn num_of_positions(&self) -> usize {
        self.inner.num_of_positions()
    }

    #[getter]
    pub fn positions1(&self) -> Vec<PyPosition> {
        self.inner
            .positions1()
            .iter()
            .map(|p| PyPosition::from_core(*p))
            .collect()
    }

    #[getter]
    pub fn positions2(&self) -> Vec<PyPosition> {
        self.inner
            .positions2()
            .iter()
            .map(|p| PyPosition::from_core(*p))
            .collect()
    }

    #[getter]
    pub fn center1(&self) -> PyPosition {
        PyPosition::from_core(self.inner.center1())
    }

    #[getter]
    pub fn center2(&self) -> PyPosition {
        PyPosition::from_core(self.inner.center2())
    }

    #[getter]
    pub fn shift_positions1(&self) -> Vec<PyPosition> {
        self.inner
            .shift_positions1()
            .iter()
            .map(|p| PyPosition::from_core(*p))
            .collect()
    }

    #[getter]
    pub fn shift_positions2(&self) -> Vec<PyPosition> {
        self.inner
            .shift_positions2()
            .iter()
            .map(|p| PyPosition::from_core(*p))
            .collect()
    }

    #[getter]
    pub fn rotation_mat(&self) -> PyMatrix {
        PyMatrix::from_core(self.inner.rotation_mat().clone())
    }

    #[getter]
    pub fn update_positions1(&self) -> Vec<PyPosition> {
        self.inner
            .update_positions1()
            .iter()
            .map(|p| PyPosition::from_core(*p))
            .collect()
    }

    #[getter]
    pub fn update_positions2(&self) -> Vec<PyPosition> {
        self.inner
            .update_positions2()
            .iter()
            .map(|p| PyPosition::from_core(*p))
            .collect()
    }

    #[getter]
    pub fn rmsd(&self) -> f64 {
        self.inner.rmsd()
    }

    pub fn superimpose(&self, py: Python<'_>, atomgroup: &PyAtomGroup) -> PyResult<PyAtomGroup> {
        let res = self
            .inner
            .superimpose(&atomgroup.to_core(py))
            .map_err(to_py_err)?;
        PyAtomGroup::from_core(py, res)
    }
}
