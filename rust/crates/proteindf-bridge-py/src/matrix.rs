// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::error::to_py_err;
use crate::vector::PyVector;
use proteindf_bridge::matrix::{Matrix as CoreMatrix, SymmetricMatrix as CoreSymmetricMatrix};
use pyo3::exceptions::{PyIndexError, PyTypeError};
use pyo3::prelude::*;

#[pyclass(name = "Matrix", module = "proteindf_bridge_rs")]
#[derive(Clone)]
pub struct PyMatrix {
    pub(crate) inner: CoreMatrix,
}

impl PyMatrix {
    pub fn from_core(mat: CoreMatrix) -> Self {
        Self { inner: mat }
    }
}

#[pymethods]
impl PyMatrix {
    #[new]
    #[pyo3(signature = (*args))]
    pub fn new(args: &Bound<'_, pyo3::types::PyTuple>) -> PyResult<Self> {
        let len = args.len();
        if len == 0 {
            Ok(Self {
                inner: CoreMatrix::default_1x1(),
            })
        } else if len == 1 {
            let arg = args.get_item(0)?;
            if let Ok(m) = arg.extract::<PyRef<PyMatrix>>() {
                Ok(Self {
                    inner: m.inner.clone(),
                })
            } else if let Ok(sm) = arg.extract::<PyRef<PySymmetricMatrix>>() {
                Ok(Self {
                    inner: sm.inner.get_general_matrix(),
                })
            } else if let Ok(list2d) = arg.extract::<Vec<Vec<f64>>>() {
                if list2d.is_empty() || list2d[0].is_empty() {
                    return Err(PyTypeError::new_err("Matrix cannot be empty"));
                }
                let rows = list2d.len();
                let cols = list2d[0].len();
                let mut data = Vec::with_capacity(rows * cols);
                for row in list2d {
                    if row.len() != cols {
                        return Err(PyTypeError::new_err(
                            "All rows must have the same number of columns",
                        ));
                    }
                    data.extend(row);
                }
                Ok(Self {
                    inner: CoreMatrix::from_vec(rows, cols, data),
                })
            } else {
                Err(PyTypeError::new_err(format!(
                    "Unsupported argument for Matrix: {}",
                    arg.get_type()
                )))
            }
        } else if len == 2 {
            let rows = args.get_item(0)?.extract::<usize>()?;
            let cols = args.get_item(1)?.extract::<usize>()?;
            Ok(Self {
                inner: CoreMatrix::new(rows, cols),
            })
        } else {
            Err(PyTypeError::new_err("Matrix takes at most 2 arguments"))
        }
    }

    #[getter]
    pub fn rows(&self) -> usize {
        self.inner.rows()
    }

    #[getter]
    pub fn cols(&self) -> usize {
        self.inner.cols()
    }

    #[getter]
    pub fn max(&self) -> f64 {
        self.inner.max()
    }

    #[getter]
    pub fn min(&self) -> f64 {
        self.inner.min()
    }

    #[getter]
    pub fn data<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let np = py.import("numpy")?;
        let list = self.to_list();
        np.call_method1("array", (list,))
    }

    pub fn get(&self, row: usize, col: usize) -> PyResult<f64> {
        self.inner.get(row, col).map_err(to_py_err)
    }

    pub fn set(&mut self, row: usize, col: usize, value: f64) -> PyResult<()> {
        if row >= self.inner.rows() || col >= self.inner.cols() {
            return Err(PyIndexError::new_err("Matrix index out of bounds"));
        }
        self.inner.set(row, col, value);
        Ok(())
    }

    pub fn add(&mut self, row: usize, col: usize, value: f64) -> PyResult<()> {
        if row >= self.inner.rows() || col >= self.inner.cols() {
            return Err(PyIndexError::new_err("Matrix index out of bounds"));
        }
        self.inner.add(row, col, value);
        Ok(())
    }

    pub fn transpose(slf: Py<Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let mut borrow = slf.borrow_mut(py);
        borrow.inner = borrow.inner.transpose();
        drop(borrow);
        Ok(slf)
    }

    pub fn select(
        &self,
        start_row: usize,
        start_col: usize,
        end_row: usize,
        end_col: usize,
    ) -> PyResult<Self> {
        self.inner
            .select(start_row, start_col, end_row, end_col)
            .map(|m| Self { inner: m })
            .map_err(to_py_err)
    }

    pub fn get_row_vector(&self, row: usize) -> PyResult<PyVector> {
        self.inner
            .get_row_vector(row)
            .map(PyVector::from_core)
            .map_err(to_py_err)
    }

    pub fn get_col_vector(&self, col: usize) -> PyResult<PyVector> {
        self.inner
            .get_col_vector(col)
            .map(PyVector::from_core)
            .map_err(to_py_err)
    }

    pub fn resize(&mut self, new_rows: usize, new_cols: usize) {
        self.inner.resize(new_rows, new_cols);
    }

    pub fn inverse(&self) -> PyResult<Self> {
        self.inner
            .inverse()
            .map(|m| Self { inner: m })
            .map_err(to_py_err)
    }

    pub fn to_list(&self) -> Vec<Vec<f64>> {
        let mut res = Vec::with_capacity(self.inner.rows());
        for r in 0..self.inner.rows() {
            let mut row = Vec::with_capacity(self.inner.cols());
            for c in 0..self.inner.cols() {
                row.push(self.inner.get(r, c).unwrap_or(0.0));
            }
            res.push(row);
        }
        res
    }

    pub fn __add__(&self, other: &PyMatrix) -> PyResult<Self> {
        if self.inner.rows() != other.inner.rows() || self.inner.cols() != other.inner.cols() {
            return Err(PyTypeError::new_err(
                "Matrix dimensions must match for addition",
            ));
        }
        Ok(Self {
            inner: &self.inner + &other.inner,
        })
    }

    pub fn __sub__(&self, other: &PyMatrix) -> PyResult<Self> {
        if self.inner.rows() != other.inner.rows() || self.inner.cols() != other.inner.cols() {
            return Err(PyTypeError::new_err(
                "Matrix dimensions must match for subtraction",
            ));
        }
        Ok(Self {
            inner: &self.inner - &other.inner,
        })
    }

    pub fn __mul__<'py>(
        &self,
        py: Python<'py>,
        other: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyAny>> {
        if let Ok(scalar) = other.extract::<f64>() {
            let res = Self {
                inner: &self.inner * scalar,
            };
            Ok(res.into_pyobject(py)?.into_any())
        } else if let Ok(other_mat) = other.extract::<PyRef<PyMatrix>>() {
            if self.inner.cols() != other_mat.inner.rows() {
                return Err(PyTypeError::new_err(format!(
                    "Inner dimensions must match: {} != {}",
                    self.inner.cols(),
                    other_mat.inner.rows()
                )));
            }
            let res = Self {
                inner: &self.inner * &other_mat.inner,
            };
            Ok(res.into_pyobject(py)?.into_any())
        } else if let Ok(other_vec) = other.extract::<PyRef<PyVector>>() {
            if self.inner.cols() != other_vec.inner.len() {
                return Err(PyTypeError::new_err(format!(
                    "Matrix cols must match Vector length: {} != {}",
                    self.inner.cols(),
                    other_vec.inner.len()
                )));
            }
            let res = PyVector::from_core(&self.inner * &other_vec.inner);
            Ok(res.into_pyobject(py)?.into_any())
        } else {
            Err(PyTypeError::new_err(
                "Unsupported type for Matrix multiplication",
            ))
        }
    }

    pub fn __rmul__(&self, scalar: f64) -> Self {
        Self {
            inner: &self.inner * scalar,
        }
    }

    pub fn __repr__(&self) -> String {
        format!("Matrix({}x{})", self.inner.rows(), self.inner.cols())
    }

    pub fn __str__(&self) -> String {
        format!("{}", self.inner)
    }

    pub fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        if let Ok(other_mat) = other.extract::<PyRef<PyMatrix>>() {
            if self.inner.rows() != other_mat.inner.rows()
                || self.inner.cols() != other_mat.inner.cols()
            {
                return false;
            }
            for r in 0..self.inner.rows() {
                for c in 0..self.inner.cols() {
                    let a = self.inner.get(r, c).unwrap_or(0.0);
                    let b = other_mat.inner.get(r, c).unwrap_or(0.0);
                    if (a - b).abs() > 1.0e-5 {
                        return false;
                    }
                }
            }
            true
        } else {
            false
        }
    }
}

#[pyclass(name = "SymmetricMatrix", module = "proteindf_bridge_rs")]
#[derive(Clone)]
pub struct PySymmetricMatrix {
    pub(crate) inner: CoreSymmetricMatrix,
}

impl PySymmetricMatrix {
    pub fn from_core(mat: CoreSymmetricMatrix) -> Self {
        Self { inner: mat }
    }
}

#[pymethods]
impl PySymmetricMatrix {
    #[new]
    #[pyo3(signature = (*args))]
    pub fn new(args: &Bound<'_, pyo3::types::PyTuple>) -> PyResult<Self> {
        let len = args.len();
        if len == 0 {
            Ok(Self {
                inner: CoreSymmetricMatrix::new(1),
            })
        } else if len == 1 {
            let arg = args.get_item(0)?;
            if let Ok(sm) = arg.extract::<PyRef<PySymmetricMatrix>>() {
                Ok(Self {
                    inner: sm.inner.clone(),
                })
            } else if let Ok(dim) = arg.extract::<usize>() {
                Ok(Self {
                    inner: CoreSymmetricMatrix::new(dim),
                })
            } else if let Ok(list2d) = arg.extract::<Vec<Vec<f64>>>() {
                let dim = list2d.len();
                let mut sm = CoreSymmetricMatrix::new(dim);
                for (r, row) in list2d.into_iter().enumerate() {
                    if row.len() != dim {
                        return Err(PyTypeError::new_err("SymmetricMatrix must be square"));
                    }
                    for (c, val) in row.into_iter().enumerate() {
                        sm.set(r, c, val);
                    }
                }
                Ok(Self { inner: sm })
            } else {
                Err(PyTypeError::new_err(format!(
                    "Unsupported argument for SymmetricMatrix: {}",
                    arg.get_type()
                )))
            }
        } else {
            Err(PyTypeError::new_err(
                "SymmetricMatrix takes at most 1 argument",
            ))
        }
    }

    #[getter]
    pub fn dim(&self) -> usize {
        self.inner.dim()
    }

    #[getter]
    pub fn rows(&self) -> usize {
        self.inner.dim()
    }

    #[getter]
    pub fn cols(&self) -> usize {
        self.inner.dim()
    }

    pub fn get(&self, row: usize, col: usize) -> PyResult<f64> {
        self.inner.get(row, col).map_err(to_py_err)
    }

    pub fn set(&mut self, row: usize, col: usize, value: f64) -> PyResult<()> {
        if row >= self.inner.dim() || col >= self.inner.dim() {
            return Err(PyIndexError::new_err("Index out of bounds"));
        }
        self.inner.set(row, col, value);
        Ok(())
    }

    pub fn add(&mut self, row: usize, col: usize, value: f64) -> PyResult<()> {
        if row >= self.inner.dim() || col >= self.inner.dim() {
            return Err(PyIndexError::new_err("Index out of bounds"));
        }
        self.inner.add(row, col, value);
        Ok(())
    }

    pub fn get_general_matrix(&self) -> PyMatrix {
        PyMatrix::from_core(self.inner.get_general_matrix())
    }

    pub fn resize(&mut self, new_dim: usize) {
        self.inner.resize(new_dim);
    }

    pub fn get_raw_data(&self) -> Vec<f64> {
        self.inner.get_raw_data()
    }

    pub fn eig(&self) -> PyResult<(PyVector, PyMatrix)> {
        let (w, v) = self.inner.eig().map_err(to_py_err)?;
        Ok((PyVector::from_core(w), PyMatrix::from_core(v)))
    }

    pub fn __repr__(&self) -> String {
        format!("SymmetricMatrix({}x{})", self.inner.dim(), self.inner.dim())
    }

    pub fn __str__(&self) -> String {
        format!("{}", self.inner.get_general_matrix())
    }

    pub fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        if let Ok(other_sm) = other.extract::<PyRef<PySymmetricMatrix>>() {
            if self.inner.dim() != other_sm.inner.dim() {
                return false;
            }
            for r in 0..self.inner.dim() {
                for c in 0..=r {
                    let a = self.inner.get(r, c).unwrap_or(0.0);
                    let b = other_sm.inner.get(r, c).unwrap_or(0.0);
                    if (a - b).abs() > 1.0e-5 {
                        return false;
                    }
                }
            }
            true
        } else {
            false
        }
    }
}
