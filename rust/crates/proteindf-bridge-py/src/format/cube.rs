// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom_group::PyAtomGroup;
use crate::error::to_py_err;
use proteindf_bridge::format::cube::{Cube as CoreCube, CubeUnit};
use pyo3::prelude::*;

#[pyclass(name = "Cube", module = "proteindf_bridge.rs", skip_from_py_object)]
#[derive(Clone, Default)]
pub struct PyCube {
    pub(crate) inner: CoreCube,
}

type Spacing3D = ((f64, f64, f64), (f64, f64, f64), (f64, f64, f64));

#[pymethods]
impl PyCube {
    #[new]
    #[pyo3(signature = (file_path=None))]
    pub fn new(file_path: Option<&str>) -> PyResult<Self> {
        let mut inner = CoreCube::new();
        if let Some(path) = file_path {
            inner.load(path).map_err(to_py_err)?;
        }
        Ok(Self { inner })
    }

    #[getter]
    pub fn comment1(&self) -> String {
        self.inner.comment1().to_string()
    }

    #[setter]
    pub fn set_comment1(&mut self, val: String) {
        self.inner.set_comment1(val);
    }

    #[getter]
    pub fn comment2(&self) -> String {
        self.inner.comment2().to_string()
    }

    #[setter]
    pub fn set_comment2(&mut self, val: String) {
        self.inner.set_comment2(val);
    }

    #[getter]
    pub fn origin(&self) -> (f64, f64, f64) {
        let o = self.inner.origin();
        (o[0], o[1], o[2])
    }

    #[getter]
    pub fn unit(&self) -> String {
        match self.inner.unit() {
            CubeUnit::Bohr => "bohr".to_string(),
            CubeUnit::Angstrom => "angstrom".to_string(),
        }
    }

    #[getter]
    pub fn n_grid(&self) -> (usize, usize, usize) {
        let n = self.inner.n_grid();
        (n[0], n[1], n[2])
    }

    #[getter]
    pub fn spacing(&self) -> Spacing3D {
        let s = self.inner.spacing();
        (
            (s[0][0], s[0][1], s[0][2]),
            (s[1][0], s[1][1], s[1][2]),
            (s[2][0], s[2][1], s[2][2]),
        )
    }

    #[getter]
    pub fn num_of_atoms(&self) -> usize {
        self.inner.num_atoms()
    }

    #[getter]
    pub fn is_molecular_orbital(&self) -> bool {
        self.inner.is_molecular_orbital()
    }

    #[getter]
    pub fn orbitals(&self) -> Option<Vec<i32>> {
        self.inner.orbitals().map(|s| s.to_vec())
    }

    #[getter]
    pub fn num_datasets(&self) -> usize {
        self.inner.num_datasets()
    }

    pub fn load(&mut self, file_path: &str) -> PyResult<()> {
        self.inner.load(file_path).map_err(to_py_err)
    }

    pub fn parse_str(&mut self, content: &str) -> PyResult<()> {
        self.inner.parse(content).map_err(to_py_err)
    }

    pub fn grid_point_position(
        &self,
        ix: usize,
        iy: usize,
        iz: usize,
    ) -> PyResult<(f64, f64, f64)> {
        let p = self
            .inner
            .grid_point_position(ix, iy, iz)
            .map_err(to_py_err)?;
        Ok((p[0], p[1], p[2]))
    }

    #[pyo3(signature = (ix, iy, iz, dataset_idx=0))]
    pub fn value(&self, ix: usize, iy: usize, iz: usize, dataset_idx: usize) -> PyResult<f64> {
        self.inner
            .orbital_value(dataset_idx, ix, iy, iz)
            .map_err(to_py_err)
    }

    #[pyo3(signature = (dataset_idx=0))]
    pub fn data(&self, dataset_idx: usize) -> PyResult<Vec<f64>> {
        let slice = self.inner.data(dataset_idx).map_err(to_py_err)?;
        Ok(slice.to_vec())
    }

    pub fn get_atomgroup(&self) -> PyResult<PyAtomGroup> {
        let ag = self.inner.get_atomgroup().map_err(to_py_err)?;
        Ok(PyAtomGroup::from_core(ag))
    }
}
