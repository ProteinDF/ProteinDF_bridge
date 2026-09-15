// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::error::to_py_err;
use crate::matrix::PyMatrix;
use proteindf_bridge::position::Position as CorePosition;
use pyo3::exceptions::{PyIndexError, PyTypeError};
use pyo3::prelude::*;
use std::str::FromStr;

#[pyclass(name = "Position", module = "proteindf_bridge_rs")]
#[derive(Clone)]
pub struct PyPosition {
    pub(crate) inner: CorePosition,
}

impl PyPosition {
    pub fn from_core(pos: CorePosition) -> Self {
        Self { inner: pos }
    }
}

fn extract_position(arg: &Bound<'_, PyAny>) -> PyResult<CorePosition> {
    if let Ok(pos) = arg.extract::<PyRef<PyPosition>>() {
        Ok(pos.inner)
    } else if let Ok(s) = arg.extract::<String>() {
        CorePosition::from_str(&s).map_err(to_py_err)
    } else if let Ok(list) = arg.extract::<Vec<f64>>() {
        CorePosition::from_slice(&list).map_err(to_py_err)
    } else {
        Err(PyTypeError::new_err(format!(
            "Unsupported type for Position: {}",
            arg.get_type()
        )))
    }
}

#[pymethods]
impl PyPosition {
    #[new]
    #[pyo3(signature = (*args))]
    pub fn new(args: &Bound<'_, pyo3::types::PyTuple>) -> PyResult<Self> {
        let len = args.len();
        if len == 0 {
            Ok(Self {
                inner: CorePosition::default(),
            })
        } else if len == 1 {
            let arg = args.get_item(0)?;
            let core_pos = extract_position(&arg)?;
            Ok(Self { inner: core_pos })
        } else if len == 3 {
            let x = args.get_item(0)?.extract::<f64>()?;
            let y = args.get_item(1)?.extract::<f64>()?;
            let z = args.get_item(2)?.extract::<f64>()?;
            Ok(Self {
                inner: CorePosition::new(x, y, z),
            })
        } else {
            Err(PyTypeError::new_err("Position takes 0, 1, or 3 arguments"))
        }
    }

    #[getter]
    pub fn x(&self) -> f64 {
        self.inner.x
    }

    #[setter]
    pub fn set_x(&mut self, val: f64) {
        self.inner.x = val;
    }

    #[getter]
    pub fn y(&self) -> f64 {
        self.inner.y
    }

    #[setter]
    pub fn set_y(&mut self, val: f64) {
        self.inner.y = val;
    }

    #[getter]
    pub fn z(&self) -> f64 {
        self.inner.z
    }

    #[setter]
    pub fn set_z(&mut self, val: f64) {
        self.inner.z = val;
    }

    #[getter]
    pub fn xyz(&self) -> [f64; 3] {
        self.inner.xyz()
    }

    #[setter]
    pub fn set_xyz(&mut self, arg: &Bound<'_, PyAny>) -> PyResult<()> {
        let new_pos = extract_position(arg)?;
        self.inner.move_to(new_pos);
        Ok(())
    }

    #[getter]
    pub fn position(&self) -> [f64; 3] {
        self.xyz()
    }

    #[setter]
    pub fn set_position(&mut self, arg: &Bound<'_, PyAny>) -> PyResult<()> {
        self.set_xyz(arg)
    }

    #[getter]
    pub fn epsilon(&self) -> f64 {
        self.inner.epsilon
    }

    #[setter]
    pub fn set_epsilon(&mut self, val: f64) {
        self.inner.epsilon = val;
    }

    pub fn move_to(&mut self, arg: &Bound<'_, PyAny>) -> PyResult<()> {
        let new_pos = extract_position(arg)?;
        self.inner.move_to(new_pos);
        Ok(())
    }

    #[pyo3(signature = (other=None))]
    pub fn square_distance_from(&self, other: Option<&Bound<'_, PyAny>>) -> PyResult<f64> {
        let target = match other {
            Some(o) => extract_position(o)?,
            None => CorePosition::default(),
        };
        Ok(self.inner.square_distance_from(&target))
    }

    #[pyo3(signature = (other=None))]
    pub fn distance_from(&self, other: Option<&Bound<'_, PyAny>>) -> PyResult<f64> {
        let target = match other {
            Some(o) => extract_position(o)?,
            None => CorePosition::default(),
        };
        Ok(self.inner.distance_from(&target))
    }

    pub fn norm(slf: Py<Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let mut borrow = slf.borrow_mut(py);
        borrow.inner.norm().map_err(to_py_err)?;
        drop(borrow);
        Ok(slf)
    }

    pub fn rotate(&mut self, mat: &PyMatrix) -> PyResult<()> {
        self.inner.rotate(&mat.inner).map_err(to_py_err)
    }

    pub fn dot(&self, other: &Bound<'_, PyAny>) -> PyResult<f64> {
        let target = extract_position(other)?;
        Ok(self.inner.dot(&target))
    }

    pub fn cross(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
        let target = extract_position(other)?;
        Ok(Self {
            inner: self.inner.cross(&target),
        })
    }

    pub fn get_raw_data(&self) -> [f64; 3] {
        self.inner.get_raw_data()
    }

    pub fn __getitem__(&self, index: isize) -> PyResult<f64> {
        match index {
            0 | -3 => Ok(self.inner.x),
            1 | -2 => Ok(self.inner.y),
            2 | -1 => Ok(self.inner.z),
            _ => Err(PyIndexError::new_err("index out of range")),
        }
    }

    pub fn __setitem__(&mut self, index: isize, value: f64) -> PyResult<()> {
        match index {
            0 | -3 => {
                self.inner.x = value;
                Ok(())
            }
            1 | -2 => {
                self.inner.y = value;
                Ok(())
            }
            2 | -1 => {
                self.inner.z = value;
                Ok(())
            }
            _ => Err(PyIndexError::new_err("index out of range")),
        }
    }

    pub fn __abs__(&self) -> f64 {
        self.inner.length()
    }

    pub fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
        let target = extract_position(other)?;
        Ok(Self {
            inner: self.inner + target,
        })
    }

    pub fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
        let target = extract_position(other)?;
        Ok(Self {
            inner: self.inner - target,
        })
    }

    pub fn __mul__<'py>(
        &self,
        py: Python<'py>,
        other: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyAny>> {
        if let Ok(val) = other.extract::<f64>() {
            let res = Self {
                inner: self.inner * val,
            };
            Ok(res.into_pyobject(py)?.into_any())
        } else if let Ok(target) = extract_position(other) {
            let dot_prod = self.inner.dot(&target);
            Ok(dot_prod.into_pyobject(py)?.into_any())
        } else {
            Err(PyTypeError::new_err(
                "Unsupported type for Position multiplication",
            ))
        }
    }

    pub fn __rmul__(&self, scalar: f64) -> Self {
        Self {
            inner: self.inner * scalar,
        }
    }

    pub fn __truediv__(&self, scalar: f64) -> PyResult<Self> {
        if scalar.abs() < 1.0e-15 {
            return Err(pyo3::exceptions::PyZeroDivisionError::new_err(
                "division by zero",
            ));
        }
        Ok(Self {
            inner: self.inner / scalar,
        })
    }

    pub fn __neg__(&self) -> Self {
        Self { inner: -self.inner }
    }

    pub fn __repr__(&self) -> String {
        format!(
            "Position([{:.6}, {:.6}, {:.6}])",
            self.inner.x, self.inner.y, self.inner.z
        )
    }

    pub fn __str__(&self) -> String {
        format!("{}", self.inner)
    }

    pub fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        if let Ok(pos) = extract_position(other) {
            self.inner.distance_from(&pos) < self.inner.epsilon
        } else {
            false
        }
    }
}
