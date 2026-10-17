// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::error::to_py_err;
use proteindf_bridge::position::Position as CorePosition;
use pyo3::exceptions::{PyIndexError, PyTypeError};
use pyo3::prelude::*;
use std::str::FromStr;

/// Represents a 3D Cartesian position vector (x, y, z).
///
/// Supports basic vector arithmetic (addition, subtraction, scalar multiplication,
/// dot product, cross product, and normalization).
///
/// Args:
///     *args: 0, 1, or 3 positional arguments.
///         - If 0 arguments: Initializes to origin (0.0, 0.0, 0.0).
///         - If 1 argument: Another Position, a sequence of 3 numbers, or a coordinate string.
///         - If 3 arguments: Three numbers (x, y, z).
///
/// Raises:
///     TypeError: If the arguments are invalid in count or type.
///     BrInputError: If parsing from a string or slice fails.
#[pyclass(name = "Position", module = "proteindf_bridge.rs", from_py_object)]
#[derive(Clone)]
pub struct PyPosition {
    pub(crate) inner: CorePosition,
}

impl PyPosition {
    pub fn from_core(pos: CorePosition) -> Self {
        Self { inner: pos }
    }
}

pub(crate) fn extract_position(arg: &Bound<'_, PyAny>) -> PyResult<CorePosition> {
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

static NUMBERS_REAL: pyo3::sync::PyOnceLock<Py<PyAny>> = pyo3::sync::PyOnceLock::new();

fn extract_scalar(py: Python<'_>, obj: &Bound<'_, PyAny>) -> Option<f64> {
    if let Ok(val) = obj.extract::<f64>() {
        return Some(val);
    }
    let real_cls = NUMBERS_REAL
        .get_or_try_init(py, || -> PyResult<Py<PyAny>> {
            let m = py.import("numbers")?;
            let cls = m.getattr("Real")?;
            Ok(cls.unbind())
        })
        .ok()?;

    if obj.is_instance(real_cls.bind(py)).unwrap_or(false) {
        if let Ok(f) = obj.call_method0("__float__") {
            if let Ok(val) = f.extract::<f64>() {
                return Some(val);
            }
        }
    }
    None
}

#[pymethods]
impl PyPosition {
    /// Initializes a new Position instance.
    ///
    /// Args:
    ///     *args: 0, 1, or 3 positional arguments.
    ///         - If 0 arguments: Initializes to (0.0, 0.0, 0.0).
    ///         - If 1 argument: Another Position, a sequence of 3 numbers, or a coordinate string.
    ///         - If 3 arguments: Three numbers (x, y, z).
    ///
    /// Raises:
    ///     TypeError: If the arguments are invalid in count or type.
    ///     BrInputError: If parsing from a string or slice fails.
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

    /// The x-coordinate value.
    #[getter]
    pub fn x(&self) -> f64 {
        self.inner.x
    }

    #[setter]
    pub fn set_x(&mut self, val: f64) {
        self.inner.x = val;
    }

    /// The y-coordinate value.
    #[getter]
    pub fn y(&self) -> f64 {
        self.inner.y
    }

    #[setter]
    pub fn set_y(&mut self, val: f64) {
        self.inner.y = val;
    }

    /// The z-coordinate value.
    #[getter]
    pub fn z(&self) -> f64 {
        self.inner.z
    }

    #[setter]
    pub fn set_z(&mut self, val: f64) {
        self.inner.z = val;
    }

    /// The coordinates as a list of three floats [x, y, z].
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

    /// Alias for xyz. The coordinates as a list of three floats [x, y, z].
    #[getter]
    pub fn position(&self) -> [f64; 3] {
        self.xyz()
    }

    #[setter]
    pub fn set_position(&mut self, arg: &Bound<'_, PyAny>) -> PyResult<()> {
        self.set_xyz(arg)
    }

    /// Distance tolerance threshold used for equality comparisons.
    #[getter]
    pub fn epsilon(&self) -> f64 {
        self.inner.epsilon
    }

    #[setter]
    pub fn set_epsilon(&mut self, val: f64) {
        self.inner.epsilon = val;
    }

    /// Sets the position to new coordinates in-place.
    ///
    /// Args:
    ///     arg: Target coordinates as a Position, a sequence of 3 floats, or a coordinate string.
    ///
    /// Raises:
    ///     TypeError: If arg is not a valid coordinate specification.
    ///     BrInputError: If parsing from a string or slice fails.
    pub fn move_to(&mut self, arg: &Bound<'_, PyAny>) -> PyResult<()> {
        let new_pos = extract_position(arg)?;
        self.inner.move_to(new_pos);
        Ok(())
    }

    /// Calculates the squared Euclidean distance to another position.
    ///
    /// Args:
    ///     other: The target position as a Position, a sequence of 3 floats, or a coordinate string.
    ///         Defaults to None (the origin).
    ///
    /// Returns:
    ///     The squared Euclidean distance as a float.
    ///
    /// Raises:
    ///     TypeError: If other is not a valid coordinate specification.
    ///     BrInputError: If parsing other fails.
    #[pyo3(signature = (other=None))]
    pub fn square_distance_from(&self, other: Option<&Bound<'_, PyAny>>) -> PyResult<f64> {
        let target = match other {
            Some(o) => extract_position(o)?,
            None => CorePosition::default(),
        };
        Ok(self.inner.square_distance_from(&target))
    }

    /// Calculates the Euclidean distance to another position.
    ///
    /// Args:
    ///     other: The target position as a Position, a sequence of 3 floats, or a coordinate string.
    ///         Defaults to None (the origin).
    ///
    /// Returns:
    ///     The Euclidean distance as a float.
    ///
    /// Raises:
    ///     TypeError: If other is not a valid coordinate specification.
    ///     BrInputError: If parsing other fails.
    #[pyo3(signature = (other=None))]
    pub fn distance_from(&self, other: Option<&Bound<'_, PyAny>>) -> PyResult<f64> {
        let target = match other {
            Some(o) => extract_position(o)?,
            None => CorePosition::default(),
        };
        Ok(self.inner.distance_from(&target))
    }

    /// Normalizes the position vector in-place to unit length.
    ///
    /// Returns:
    ///     The normalized Position object (self).
    ///
    /// Raises:
    ///     BrValueError: If the vector has zero length.
    pub fn norm(slf: Py<Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let mut borrow = slf.borrow_mut(py);
        borrow.inner.norm().map_err(to_py_err)?;
        drop(borrow);
        Ok(slf)
    }

    /// Rotates the position vector in-place by a 3x3 rotation matrix.
    ///
    /// Args:
    ///     mat: A 3x3 rotation matrix (Matrix or an object with rows, cols, and get(row, col)).
    ///
    /// Raises:
    ///     TypeError: If mat is not a Matrix or an object with rows, cols, and get(row, col).
    ///     ValueError: If matrix dimensions are not 3x3.
    pub fn rotate(&mut self, mat: &Bound<'_, PyAny>) -> PyResult<()> {
        let core_mat = crate::matrix::extract_rotation_matrix(mat)?;
        self.inner.rotate(&core_mat).map_err(to_py_err)
    }

    /// Computes the dot product with another position vector.
    ///
    /// Args:
    ///     other: The other position as a Position, sequence of 3 floats, or coordinate string.
    ///
    /// Returns:
    ///     The dot product as a float.
    ///
    /// Raises:
    ///     TypeError: If other cannot be parsed as a position.
    ///     BrInputError: If parsing other fails.
    pub fn dot(&self, other: &Bound<'_, PyAny>) -> PyResult<f64> {
        let target = extract_position(other)?;
        Ok(self.inner.dot(&target))
    }

    /// Computes the cross product with another position vector.
    ///
    /// Args:
    ///     other: The other position as a Position, sequence of 3 floats, or coordinate string.
    ///
    /// Returns:
    ///     A new Position instance representing the cross product vector.
    ///
    /// Raises:
    ///     TypeError: If other cannot be parsed as a position.
    ///     BrInputError: If parsing other fails.
    pub fn cross(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
        let target = extract_position(other)?;
        Ok(Self {
            inner: self.inner.cross(&target),
        })
    }

    /// Returns the coordinates as a list of three floats [x, y, z].
    ///
    /// Returns:
    ///     A list containing the x, y, and z coordinates.
    pub fn get_raw_data(&self) -> [f64; 3] {
        self.inner.get_raw_data()
    }

    /// Calculates the dihedral (torsion) angle in degrees formed by four points.
    ///
    /// The angle is measured around the axis formed by p2 and p3 for the sequence
    /// (self, p2, p3, p4).
    ///
    /// Args:
    ///     p2: The second position.
    ///     p3: The third position.
    ///     p4: The fourth position.
    ///
    /// Returns:
    ///     The dihedral angle in degrees (-180.0 to 180.0).
    ///
    /// Raises:
    ///     TypeError: If any point cannot be parsed as a position.
    ///     BrInputError: If parsing any point fails.
    pub fn dihedral(
        &self,
        p2: &Bound<'_, PyAny>,
        p3: &Bound<'_, PyAny>,
        p4: &Bound<'_, PyAny>,
    ) -> PyResult<f64> {
        let pos2 = extract_position(p2)?;
        let pos3 = extract_position(p3)?;
        let pos4 = extract_position(p4)?;
        Ok(self.inner.dihedral(&pos2, &pos3, &pos4))
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

    pub fn __len__(&self) -> usize {
        3
    }

    pub fn __iter__<'py>(slf: Bound<'py, Self>) -> PyResult<Bound<'py, pyo3::types::PyIterator>> {
        let borrow = slf.borrow();
        let py = slf.py();
        let tuple = (borrow.inner.x, borrow.inner.y, borrow.inner.z).into_pyobject(py)?;
        tuple.try_iter()
    }

    pub fn __mul__<'py>(
        &self,
        py: Python<'py>,
        other: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyAny>> {
        if let Some(val) = extract_scalar(py, other) {
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

    pub fn __rmul__<'py>(
        &self,
        py: Python<'py>,
        other: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyAny>> {
        if let Some(val) = extract_scalar(py, other) {
            let res = Self {
                inner: self.inner * val,
            };
            Ok(res.into_pyobject(py)?.into_any())
        } else {
            Ok(py.NotImplemented().into_bound(py))
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

    pub fn __copy__(&self) -> Self {
        self.clone()
    }

    pub fn __deepcopy__(&self, _memo: &Bound<'_, PyAny>) -> Self {
        self.clone()
    }

    pub fn __getstate__(&self) -> (f64, f64, f64) {
        (self.inner.x, self.inner.y, self.inner.z)
    }

    pub fn __setstate__(&mut self, state: &Bound<'_, PyAny>) -> PyResult<()> {
        if state.cast::<pyo3::types::PyList>().is_err()
            && state.cast::<pyo3::types::PyTuple>().is_err()
        {
            return Err(PyTypeError::new_err(
                "Position.__setstate__ requires a list or tuple",
            ));
        }
        let list = state.extract::<Vec<f64>>()?;
        self.inner = CorePosition::from_slice(&list).map_err(to_py_err)?;
        Ok(())
    }

    pub fn __reduce__<'py>(
        slf: Bound<'py, Self>,
    ) -> PyResult<(Bound<'py, PyAny>, (f64, f64, f64))> {
        let cls = slf.get_type().into_any();
        let borrow = slf.borrow();
        Ok((cls, (borrow.inner.x, borrow.inner.y, borrow.inner.z)))
    }

    pub fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        if let Ok(pos) = extract_position(other) {
            self.inner.distance_from(&pos) < self.inner.epsilon
        } else {
            false
        }
    }
}
