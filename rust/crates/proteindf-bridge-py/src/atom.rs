// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::error::to_py_err;
use crate::matrix::PyMatrix;
use crate::position::PyPosition;
use proteindf_bridge::atom::Atom as CoreAtom;
use proteindf_bridge::position::Position;
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::str::FromStr;

#[pyclass(name = "Atom", module = "proteindf_bridge_rs")]
#[derive(Clone)]
pub struct PyAtom {
    pub(crate) inner: CoreAtom,
}

impl PyAtom {
    pub fn from_core(atom: CoreAtom) -> Self {
        Self { inner: atom }
    }
}

fn extract_pos(arg: &Bound<'_, PyAny>) -> PyResult<Position> {
    if let Ok(pos) = arg.extract::<PyRef<PyPosition>>() {
        Ok(pos.inner)
    } else if let Ok(s) = arg.extract::<String>() {
        Position::from_str(&s).map_err(to_py_err)
    } else if let Ok(list) = arg.extract::<Vec<f64>>() {
        Position::from_slice(&list).map_err(to_py_err)
    } else {
        Err(PyTypeError::new_err(format!(
            "Unsupported type for Position: {}",
            arg.get_type()
        )))
    }
}

#[pymethods]
impl PyAtom {
    #[new]
    #[pyo3(signature = (*args, **kwargs))]
    pub fn new(
        args: &Bound<'_, pyo3::types::PyTuple>,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Self> {
        let mut atom = CoreAtom::new();

        if args.len() == 1 {
            let first = args.get_item(0)?;
            if let Ok(other) = first.extract::<PyRef<PyAtom>>() {
                atom = other.inner.clone();
            } else if let Ok(sym) = first.extract::<String>() {
                atom = CoreAtom::from_symbol(&sym).map_err(to_py_err)?;
            } else if let Ok(num) = first.extract::<usize>() {
                atom.set_atomic_number(num);
            } else if let Ok(dict) = first.downcast::<PyDict>() {
                if let Some(z) = dict.get_item("Z")? {
                    atom.set_atomic_number(z.extract::<usize>()?);
                }
                if let Some(name) = dict.get_item("name")? {
                    atom.name = name.extract::<String>()?;
                }
                if let Some(q) = dict.get_item("Q")? {
                    atom.charge = q.extract::<f64>()?;
                }
                if let Some(xyz) = dict.get_item("xyz")? {
                    atom.xyz = extract_pos(&xyz)?;
                }
                if let Some(force) = dict.get_item("force")? {
                    atom.force = extract_pos(&force)?;
                }
            } else {
                return Err(PyTypeError::new_err("Unsupported argument for Atom"));
            }
        } else if args.len() > 1 {
            return Err(PyTypeError::new_err(
                "Atom takes at most 1 positional argument",
            ));
        }

        if let Some(kw) = kwargs {
            if let Some(sym) = kw.get_item("symbol")? {
                let symbol_str = sym.extract::<String>()?;
                atom.set_symbol(&symbol_str).map_err(to_py_err)?;
            }
            if let Some(name) = kw.get_item("name")? {
                atom.name = name.extract::<String>()?;
            }
            if let Some(label) = kw.get_item("label")? {
                atom.label = label.extract::<String>()?;
            }
            if let Some(charge) = kw.get_item("charge")? {
                atom.charge = charge.extract::<f64>()?;
            }
            if let Some(path) = kw.get_item("path")? {
                atom.path = path.extract::<String>()?;
            }
            if let Some(pos) = kw.get_item("position")? {
                atom.xyz = extract_pos(&pos)?;
            }
            if let Some(xyz) = kw.get_item("xyz")? {
                atom.xyz = extract_pos(&xyz)?;
            }
            if let Some(force) = kw.get_item("force")? {
                atom.force = extract_pos(&force)?;
            }
        }

        Ok(Self { inner: atom })
    }

    #[getter]
    pub fn atomic_number(&self) -> usize {
        self.inner.atomic_number()
    }

    #[setter]
    pub fn set_atomic_number(&mut self, val: usize) {
        self.inner.set_atomic_number(val);
    }

    #[getter]
    pub fn symbol(&self) -> PyResult<String> {
        self.inner
            .symbol()
            .map(|s| s.to_string())
            .map_err(to_py_err)
    }

    #[setter]
    pub fn set_symbol(&mut self, sym: &str) -> PyResult<()> {
        self.inner.set_symbol(sym).map_err(to_py_err)
    }

    #[getter]
    pub fn name(&self) -> String {
        self.inner.name.clone()
    }

    #[setter]
    pub fn set_name(&mut self, val: String) {
        self.inner.name = val;
    }

    #[getter]
    pub fn label(&self) -> String {
        self.inner.label.clone()
    }

    #[setter]
    pub fn set_label(&mut self, val: String) {
        self.inner.label = val;
    }

    #[getter]
    pub fn charge(&self) -> f64 {
        self.inner.charge
    }

    #[setter]
    pub fn set_charge(&mut self, val: f64) {
        self.inner.charge = val;
    }

    #[getter]
    pub fn path(&self) -> String {
        self.inner.path.clone()
    }

    #[setter]
    pub fn set_path(&mut self, val: String) {
        self.inner.path = val;
    }

    #[getter]
    pub fn xyz(&self) -> PyPosition {
        PyPosition::from_core(self.inner.xyz)
    }

    #[setter]
    pub fn set_xyz(&mut self, val: &Bound<'_, PyAny>) -> PyResult<()> {
        self.inner.xyz = extract_pos(val)?;
        Ok(())
    }

    #[getter]
    pub fn position(&self) -> PyPosition {
        self.xyz()
    }

    #[setter]
    pub fn set_position(&mut self, val: &Bound<'_, PyAny>) -> PyResult<()> {
        self.set_xyz(val)
    }

    #[getter]
    pub fn force(&self) -> PyPosition {
        PyPosition::from_core(self.inner.force)
    }

    #[setter]
    pub fn set_force(&mut self, val: &Bound<'_, PyAny>) -> PyResult<()> {
        self.inner.force = extract_pos(val)?;
        Ok(())
    }

    #[getter]
    pub fn is_real(&self) -> bool {
        self.inner.is_real()
    }

    #[getter]
    pub fn vdw(&self) -> PyResult<f64> {
        self.inner.vdw().map_err(to_py_err)
    }

    pub fn weight(&self) -> PyResult<f64> {
        self.inner.weight().map_err(to_py_err)
    }

    pub fn move_to(slf: Py<Self>, pos: &Bound<'_, PyAny>, py: Python<'_>) -> PyResult<Py<Self>> {
        let p = extract_pos(pos)?;
        let mut borrow = slf.borrow_mut(py);
        borrow.inner.move_to(p);
        drop(borrow);
        Ok(slf)
    }

    pub fn shift_by(slf: Py<Self>, dir: &Bound<'_, PyAny>, py: Python<'_>) -> PyResult<Py<Self>> {
        let d = extract_pos(dir)?;
        let mut borrow = slf.borrow_mut(py);
        borrow.inner.shift_by(d);
        drop(borrow);
        Ok(slf)
    }

    pub fn rotate(&mut self, rotmat: &PyMatrix) -> PyResult<()> {
        self.inner.rotate(&rotmat.inner).map_err(to_py_err)
    }

    pub fn get_raw_data<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("Z", self.inner.atomic_number())?;
        dict.set_item("name", &self.inner.name)?;
        dict.set_item("Q", self.inner.charge)?;
        dict.set_item("xyz", self.inner.xyz.get_raw_data().to_vec())?;
        dict.set_item("force", self.inner.force.get_raw_data().to_vec())?;
        Ok(dict)
    }

    pub fn set_by_raw_data(&mut self, data: &Bound<'_, PyDict>) -> PyResult<()> {
        if let Some(z) = data.get_item("Z")? {
            self.inner.set_atomic_number(z.extract::<usize>()?);
        }
        if let Some(name) = data.get_item("name")? {
            self.inner.name = name.extract::<String>()?;
        }
        if let Some(q) = data.get_item("Q")? {
            self.inner.charge = q.extract::<f64>()?;
        }
        if let Some(xyz) = data.get_item("xyz")? {
            self.inner.xyz = extract_pos(&xyz)?;
        }
        if let Some(force) = data.get_item("force")? {
            self.inner.force = extract_pos(&force)?;
        }
        Ok(())
    }

    pub fn __imul__(&mut self, scalar: f64) {
        self.inner.xyz *= scalar;
    }

    pub fn __repr__(&self) -> String {
        let sym = self.inner.symbol().unwrap_or("?");
        format!(
            "Atom('{}', xyz=({:.3}, {:.3}, {:.3}))",
            sym, self.inner.xyz.x, self.inner.xyz.y, self.inner.xyz.z
        )
    }

    pub fn __str__(&self) -> String {
        format!("{}", self.inner)
    }

    pub fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        if let Ok(other_atom) = other.extract::<PyRef<PyAtom>>() {
            self.inner == other_atom.inner
        } else {
            false
        }
    }
}
