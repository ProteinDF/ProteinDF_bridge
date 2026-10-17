// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::error::to_py_err;
use crate::position::PyPosition;
use proteindf_bridge::atom::Atom as CoreAtom;
use proteindf_bridge::position::Position;
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::str::FromStr;

type AtomState = (
    usize,
    String,
    String,
    f64,
    String,
    (f64, f64, f64),
    (f64, f64, f64),
);

#[pyclass(name = "Atom", module = "proteindf_bridge.rs")]
pub struct PyAtom {
    pub(crate) atomic_number: usize,
    pub(crate) name: String,
    pub(crate) label: String,
    pub(crate) charge: f64,
    pub(crate) path: String,
    pub(crate) xyz: Py<PyPosition>,
    pub(crate) force: Py<PyPosition>,
}

impl PyAtom {
    pub fn from_core(py: Python<'_>, atom: CoreAtom) -> PyResult<Self> {
        let xyz = Py::new(py, PyPosition::from_core(atom.xyz))?;
        let force = Py::new(py, PyPosition::from_core(atom.force))?;
        Ok(Self {
            atomic_number: atom.atomic_number(),
            name: atom.name,
            label: atom.label,
            charge: atom.charge,
            path: atom.path,
            xyz,
            force,
        })
    }

    pub fn to_core(&self, py: Python<'_>) -> CoreAtom {
        let xyz_core = self.xyz.borrow(py).inner;
        let force_core = self.force.borrow(py).inner;
        let mut atom = CoreAtom::new();
        atom.set_atomic_number(self.atomic_number);
        atom.name = self.name.clone();
        atom.label = self.label.clone();
        atom.charge = self.charge;
        atom.path = self.path.clone();
        atom.xyz = xyz_core;
        atom.force = force_core;
        atom
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

fn to_charge(val: &Bound<'_, PyAny>) -> PyResult<f64> {
    let py = val.py();
    let float_type = py.get_type::<pyo3::types::PyFloat>();
    float_type.call1((val,))?.extract()
}

#[pymethods]
impl PyAtom {
    #[new]
    #[pyo3(signature = (*args, **kwargs))]
    pub fn new(
        py: Python<'_>,
        args: &Bound<'_, pyo3::types::PyTuple>,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Self> {
        let mut atom = CoreAtom::new();

        if args.len() == 1 {
            let first = args.get_item(0)?;
            if let Ok(other) = first.extract::<PyRef<PyAtom>>() {
                atom = other.to_core(py);
            } else if let Ok(sym) = first.extract::<String>() {
                atom = CoreAtom::from_symbol(&sym).map_err(to_py_err)?;
            } else if let Ok(num) = first.extract::<usize>() {
                atom.set_atomic_number(num);
            } else if let Ok(dict) = first.cast::<PyDict>() {
                if let Some(z) = dict.get_item("Z")? {
                    atom.set_atomic_number(z.extract::<usize>()?);
                }
                if let Some(name) = dict.get_item("name")? {
                    atom.name = name.extract::<String>()?;
                }
                if let Some(q) = dict.get_item("Q")? {
                    atom.charge = to_charge(&q)?;
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
                atom.charge = to_charge(&charge)?;
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

        Self::from_core(py, atom)
    }

    #[getter]
    pub fn atomic_number(&self) -> usize {
        self.atomic_number
    }

    #[setter]
    pub fn set_atomic_number(&mut self, val: usize) {
        self.atomic_number = val;
    }

    #[getter]
    pub fn symbol(&self) -> PyResult<String> {
        proteindf_bridge::periodic_table::PeriodicTable::get_symbol(self.atomic_number)
            .map(|s| s.to_string())
            .map_err(to_py_err)
    }

    #[setter]
    pub fn set_symbol(&mut self, sym: &str) -> PyResult<()> {
        let z = proteindf_bridge::periodic_table::PeriodicTable::get_atomic_number(sym)
            .map_err(to_py_err)?;
        self.atomic_number = z;
        Ok(())
    }

    #[getter]
    pub fn name(&self) -> String {
        self.name.clone()
    }

    #[setter]
    pub fn set_name(&mut self, val: String) {
        self.name = val;
    }

    #[getter]
    pub fn label(&self) -> String {
        self.label.clone()
    }

    #[setter]
    pub fn set_label(&mut self, val: String) {
        self.label = val;
    }

    #[getter]
    pub fn charge(&self) -> f64 {
        self.charge
    }

    #[setter]
    pub fn set_charge(&mut self, val: &Bound<'_, PyAny>) -> PyResult<()> {
        self.charge = to_charge(val)?;
        Ok(())
    }

    #[getter]
    pub fn path(&self) -> String {
        self.path.clone()
    }

    #[setter]
    pub fn set_path(&mut self, val: String) {
        self.path = val;
    }

    #[getter]
    pub fn xyz(&self, py: Python<'_>) -> Py<PyPosition> {
        self.xyz.clone_ref(py)
    }

    #[setter]
    pub fn set_xyz(&mut self, py: Python<'_>, val: &Bound<'_, PyAny>) -> PyResult<()> {
        let pos = extract_pos(val)?;
        self.xyz = Py::new(py, PyPosition::from_core(pos))?;
        Ok(())
    }

    #[getter]
    pub fn position(&self, py: Python<'_>) -> Py<PyPosition> {
        self.xyz(py)
    }

    #[setter]
    pub fn set_position(&mut self, py: Python<'_>, val: &Bound<'_, PyAny>) -> PyResult<()> {
        self.set_xyz(py, val)
    }

    #[getter]
    pub fn force(&self, py: Python<'_>) -> Py<PyPosition> {
        self.force.clone_ref(py)
    }

    #[setter]
    pub fn set_force(&mut self, py: Python<'_>, val: &Bound<'_, PyAny>) -> PyResult<()> {
        let pos = extract_pos(val)?;
        self.force = Py::new(py, PyPosition::from_core(pos))?;
        Ok(())
    }

    #[getter]
    pub fn is_real(&self) -> bool {
        self.atomic_number > 0
    }

    #[getter]
    pub fn vdw(&self) -> PyResult<f64> {
        proteindf_bridge::periodic_table::PeriodicTable::vdw(self.atomic_number).map_err(to_py_err)
    }

    pub fn weight(&self) -> PyResult<f64> {
        proteindf_bridge::periodic_table::PeriodicTable::atomic_weight(self.atomic_number)
            .map_err(to_py_err)
    }

    pub fn move_to(slf: Py<Self>, pos: &Bound<'_, PyAny>, py: Python<'_>) -> PyResult<Py<Self>> {
        let p = extract_pos(pos)?;
        let borrow = slf.borrow(py);
        borrow.xyz.borrow_mut(py).inner.move_to(p);
        drop(borrow);
        Ok(slf)
    }

    pub fn shift_by(slf: Py<Self>, dir: &Bound<'_, PyAny>, py: Python<'_>) -> PyResult<Py<Self>> {
        let d = extract_pos(dir)?;
        let borrow = slf.borrow(py);
        borrow.xyz.borrow_mut(py).inner += d;
        drop(borrow);
        Ok(slf)
    }

    pub fn rotate(&self, py: Python<'_>, rotmat: &Bound<'_, PyAny>) -> PyResult<()> {
        let core_mat = crate::matrix::extract_rotation_matrix(rotmat)?;
        self.xyz
            .borrow_mut(py)
            .inner
            .rotate(&core_mat)
            .map_err(to_py_err)
    }

    pub fn get_raw_data<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("Z", self.atomic_number)?;
        dict.set_item("name", &self.name)?;
        dict.set_item("Q", self.charge)?;
        let xyz = self.xyz.borrow(py);
        dict.set_item("xyz", xyz.inner.get_raw_data().to_vec())?;
        let force = self.force.borrow(py);
        dict.set_item("force", force.inner.get_raw_data().to_vec())?;
        Ok(dict)
    }

    pub fn set_by_raw_data(&mut self, py: Python<'_>, data: &Bound<'_, PyDict>) -> PyResult<()> {
        if let Some(z) = data.get_item("Z")? {
            self.atomic_number = z.extract::<usize>()?;
        }
        if let Some(name) = data.get_item("name")? {
            self.name = name.extract::<String>()?;
        }
        if let Some(q) = data.get_item("Q")? {
            self.charge = to_charge(&q)?;
        }
        if let Some(xyz) = data.get_item("xyz")? {
            let pos = extract_pos(&xyz)?;
            self.xyz = Py::new(py, PyPosition::from_core(pos))?;
        }
        if let Some(force) = data.get_item("force")? {
            let pos = extract_pos(&force)?;
            self.force = Py::new(py, PyPosition::from_core(pos))?;
        }
        Ok(())
    }

    pub fn __imul__(&mut self, py: Python<'_>, scalar: f64) {
        self.xyz.borrow_mut(py).inner *= scalar;
    }

    pub fn __repr__(&self, py: Python<'_>) -> String {
        let sym = proteindf_bridge::periodic_table::PeriodicTable::get_symbol(self.atomic_number)
            .unwrap_or("?");
        let xyz = self.xyz.borrow(py);
        format!(
            "Atom('{}', xyz=({:.3}, {:.3}, {:.3}))",
            sym, xyz.inner.x, xyz.inner.y, xyz.inner.z
        )
    }

    pub fn __str__(&self, py: Python<'_>) -> String {
        let core = self.to_core(py);
        format!("{}", core)
    }

    pub fn __copy__(&self, py: Python<'_>) -> PyResult<Self> {
        let core = self.to_core(py);
        Self::from_core(py, core)
    }

    pub fn __deepcopy__(&self, py: Python<'_>, _memo: &Bound<'_, PyAny>) -> PyResult<Self> {
        let core = self.to_core(py);
        Self::from_core(py, core)
    }

    pub fn __getstate__<'py>(&self, py: Python<'py>) -> AtomState {
        let xyz = self.xyz.borrow(py);
        let force = self.force.borrow(py);
        (
            self.atomic_number,
            self.name.clone(),
            self.label.clone(),
            self.charge,
            self.path.clone(),
            (xyz.inner.x, xyz.inner.y, xyz.inner.z),
            (force.inner.x, force.inner.y, force.inner.z),
        )
    }

    pub fn __setstate__(&mut self, py: Python<'_>, state: &Bound<'_, PyAny>) -> PyResult<()> {
        if let Ok(tup) = state.extract::<AtomState>() {
            self.atomic_number = tup.0;
            self.name = tup.1;
            self.label = tup.2;
            self.charge = tup.3;
            self.path = tup.4;
            self.xyz = Py::new(
                py,
                PyPosition::from_core(Position::new(tup.5 .0, tup.5 .1, tup.5 .2)),
            )?;
            self.force = Py::new(
                py,
                PyPosition::from_core(Position::new(tup.6 .0, tup.6 .1, tup.6 .2)),
            )?;
            Ok(())
        } else if let Ok(dict) = state.cast::<PyDict>() {
            self.set_by_raw_data(py, dict)
        } else {
            Err(PyTypeError::new_err("Invalid state for Atom"))
        }
    }

    pub fn __reduce__<'py>(
        slf: Bound<'py, Self>,
    ) -> PyResult<(
        Bound<'py, PyAny>,
        Bound<'py, pyo3::types::PyTuple>,
        AtomState,
    )> {
        let py = slf.py();
        let cls = slf.get_type().into_any();
        let empty_args = pyo3::types::PyTuple::empty(py);
        let state = slf.borrow().__getstate__(py);
        Ok((cls, empty_args, state))
    }

    pub fn __eq__(&self, other: &Bound<'_, PyAny>, py: Python<'_>) -> bool {
        if let Ok(other_atom) = other.extract::<PyRef<PyAtom>>() {
            let self_core = self.to_core(py);
            let other_core = other_atom.to_core(py);
            self_core == other_core
        } else {
            false
        }
    }
}
