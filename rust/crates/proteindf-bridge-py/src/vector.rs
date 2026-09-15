// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use proteindf_bridge::vector::Vector as CoreVector;
use pyo3::exceptions::{PyIndexError, PyTypeError};
use pyo3::prelude::*;

#[pyclass(name = "Vector", module = "proteindf_bridge_rs")]
#[derive(Clone)]
pub struct PyVector {
    pub(crate) inner: CoreVector,
}

impl PyVector {
    pub fn from_core(vec: CoreVector) -> Self {
        Self { inner: vec }
    }
}

#[pymethods]
impl PyVector {
    #[new]
    #[pyo3(signature = (obj=None))]
    pub fn new(obj: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        match obj {
            None => Ok(Self {
                inner: CoreVector::new(0),
            }),
            Some(any) => {
                if let Ok(vec) = any.extract::<PyRef<PyVector>>() {
                    Ok(Self {
                        inner: vec.inner.clone(),
                    })
                } else if let Ok(size) = any.extract::<usize>() {
                    Ok(Self {
                        inner: CoreVector::new(size),
                    })
                } else if let Ok(list) = any.extract::<Vec<f64>>() {
                    Ok(Self {
                        inner: CoreVector::from_vec(list),
                    })
                } else {
                    Err(PyTypeError::new_err(format!(
                        "Unsupported type for Vector: {}",
                        any.get_type()
                    )))
                }
            }
        }
    }

    #[getter]
    pub fn max(&self) -> PyResult<f64> {
        if self.inner.is_empty() {
            return Err(PyTypeError::new_err("cannot compute max of empty Vector"));
        }
        Ok(self.inner.max())
    }

    #[getter]
    pub fn min(&self) -> PyResult<f64> {
        if self.inner.is_empty() {
            return Err(PyTypeError::new_err("cannot compute min of empty Vector"));
        }
        Ok(self.inner.min())
    }

    pub fn abs(&self) -> Self {
        Self {
            inner: self.inner.abs(),
        }
    }

    #[getter]
    pub fn data<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let np = py.import("numpy")?;
        np.call_method1("array", (self.inner.to_vec(),))
    }

    pub fn size(&self) -> usize {
        self.inner.len()
    }

    pub fn resize(&mut self, new_size: usize) {
        self.inner.resize(new_size);
    }

    pub fn get(&self, index: usize) -> PyResult<f64> {
        self.inner
            .get(index)
            .map_err(|_| PyIndexError::new_err("index out of range"))
    }

    pub fn set(&mut self, index: usize, value: f64) -> PyResult<()> {
        if index >= self.inner.len() {
            return Err(PyIndexError::new_err("index out of range"));
        }
        self.inner.set(index, value);
        Ok(())
    }

    pub fn to_list(&self) -> Vec<f64> {
        self.inner.to_vec()
    }

    pub fn get_ndarray<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.data(py)
    }

    pub fn argsort(&self) -> Self {
        let indices = self.inner.argsort();
        let f64_indices: Vec<f64> = indices.into_iter().map(|x| x as f64).collect();
        Self {
            inner: CoreVector::from_vec(f64_indices),
        }
    }

    pub fn flip(&self) -> Self {
        Self {
            inner: self.inner.flip(),
        }
    }

    pub fn dot(&self, other: &PyVector) -> PyResult<f64> {
        if self.inner.len() != other.inner.len() {
            return Err(PyTypeError::new_err(format!(
                "Vector lengths must match for dot product: {} != {}",
                self.inner.len(),
                other.inner.len()
            )));
        }
        Ok(self.inner.dot(&other.inner))
    }

    pub fn __len__(&self) -> usize {
        self.inner.len()
    }

    pub fn __getitem__(&self, index: isize) -> PyResult<f64> {
        let len = self.inner.len() as isize;
        let idx = if index < 0 { len + index } else { index };
        if idx < 0 || idx >= len {
            return Err(PyIndexError::new_err("index out of range"));
        }
        Ok(self.inner[idx as usize])
    }

    pub fn __setitem__(&mut self, index: isize, value: f64) -> PyResult<()> {
        let len = self.inner.len() as isize;
        let idx = if index < 0 { len + index } else { index };
        if idx < 0 || idx >= len {
            return Err(PyIndexError::new_err("index out of range"));
        }
        self.inner[idx as usize] = value;
        Ok(())
    }

    pub fn __add__(&self, other: &PyVector) -> PyResult<Self> {
        if self.inner.len() != other.inner.len() {
            return Err(PyTypeError::new_err(
                "Vector lengths must match for addition",
            ));
        }
        Ok(Self {
            inner: &self.inner + &other.inner,
        })
    }

    pub fn __sub__(&self, other: &PyVector) -> PyResult<Self> {
        if self.inner.len() != other.inner.len() {
            return Err(PyTypeError::new_err(
                "Vector lengths must match for subtraction",
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
        if let Ok(val) = other.extract::<f64>() {
            let res = Self {
                inner: &self.inner * val,
            };
            Ok(res.into_pyobject(py)?.into_any())
        } else if let Ok(other_vec) = other.extract::<PyRef<PyVector>>() {
            let dot_prod = self.dot(&other_vec)?;
            Ok(dot_prod.into_pyobject(py)?.into_any())
        } else {
            Err(PyTypeError::new_err("Unsupported type for multiplication"))
        }
    }

    pub fn __rmul__(&self, scalar: f64) -> Self {
        Self {
            inner: &self.inner * scalar,
        }
    }

    pub fn __neg__(&self) -> Self {
        Self {
            inner: -&self.inner,
        }
    }

    pub fn __repr__(&self) -> String {
        format!("Vector({:?})", self.inner.to_vec())
    }

    pub fn __str__(&self) -> String {
        let mut output = String::new();
        let len = self.inner.len();
        for order in (0..len).step_by(10) {
            output.push('\n');
            for j in order..std::cmp::min(order + 10, len) {
                output.push_str(&format!("   {:5} th", j + 1));
            }
            output.push('\n');
            for _ in order..std::cmp::min(order + 10, len) {
                output.push_str("-----------");
            }
            output.push_str("----\n\n");
            for j in order..std::cmp::min(order + 10, len) {
                output.push_str(&format!(" {:10.6}", self.inner[j]));
            }
            output.push('\n');
        }
        output
    }

    pub fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        if let Ok(other_vec) = other.extract::<PyRef<PyVector>>() {
            self.inner == other_vec.inner
        } else {
            false
        }
    }
}
