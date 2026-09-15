// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom_group::PyAtomGroup;
use crate::error::to_py_err;
use crate::matrix::{PyMatrix, PySymmetricMatrix};
use crate::position::PyPosition;
use crate::vector::PyVector;
use proteindf_bridge::superposer_quaternion::SuperposerQuaternion as CoreSuperposerQuaternion;
use pyo3::prelude::*;

#[pyclass(name = "Superposer_quaternion", module = "proteindf_bridge_rs")]
#[derive(Clone)]
pub struct PySuperposerQuaternion {
    pub(crate) inner: CoreSuperposerQuaternion,
}

#[pymethods]
impl PySuperposerQuaternion {
    #[new]
    pub fn new(atomgroup1: &PyAtomGroup, atomgroup2: &PyAtomGroup) -> PyResult<Self> {
        let inner = CoreSuperposerQuaternion::new(&atomgroup1.inner, &atomgroup2.inner)
            .map_err(to_py_err)?;
        Ok(Self { inner })
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
    #[allow(non_snake_case)]
    pub fn r_A(&self) -> Vec<PyPosition> {
        self.inner
            .r_a()
            .iter()
            .map(|p| PyPosition::from_core(*p))
            .collect()
    }

    #[getter]
    pub fn r_a(&self) -> Vec<PyPosition> {
        self.r_A()
    }

    #[getter]
    #[allow(non_snake_case)]
    pub fn r_B(&self) -> Vec<PyPosition> {
        self.inner
            .r_b()
            .iter()
            .map(|p| PyPosition::from_core(*p))
            .collect()
    }

    #[getter]
    pub fn r_b(&self) -> Vec<PyPosition> {
        self.r_B()
    }

    #[getter]
    pub fn va(&self) -> Vec<PyPosition> {
        self.inner
            .va()
            .iter()
            .map(|p| PyPosition::from_core(*p))
            .collect()
    }

    #[getter]
    pub fn vb(&self) -> Vec<PyPosition> {
        self.inner
            .vb()
            .iter()
            .map(|p| PyPosition::from_core(*p))
            .collect()
    }

    #[getter]
    #[allow(non_snake_case)]
    pub fn matB(&self) -> PySymmetricMatrix {
        PySymmetricMatrix::from_core(self.inner.mat_b().clone())
    }

    #[getter]
    pub fn mat_b(&self) -> PySymmetricMatrix {
        self.matB()
    }

    #[getter]
    pub fn eigval(&self) -> PyVector {
        PyVector::from_core(self.inner.eigval().clone())
    }

    #[getter]
    pub fn eigvec(&self) -> PyMatrix {
        PyMatrix::from_core(self.inner.eigvec().clone())
    }

    #[getter]
    #[allow(non_snake_case)]
    pub fn matR(&self) -> PyMatrix {
        PyMatrix::from_core(self.inner.mat_r().clone())
    }

    #[getter]
    pub fn mat_r(&self) -> PyMatrix {
        self.matR()
    }

    #[getter]
    pub fn rotation_mat(&self) -> PyMatrix {
        self.matR()
    }

    #[getter]
    pub fn rmsd(&self) -> f64 {
        self.inner.rmsd()
    }

    pub fn calc(&self) -> f64 {
        self.inner.rmsd()
    }

    pub fn superimpose(&self, atomgroup: &PyAtomGroup) -> PyResult<PyAtomGroup> {
        let res = self
            .inner
            .superimpose(&atomgroup.inner)
            .map_err(to_py_err)?;
        Ok(PyAtomGroup::from_core(res))
    }
}
