// Copyright (C) 2014 The ProteinDF development team.
// see also AUTHORS and README if provided.
//
// This file is a part of the ProteinDF software package.

use crate::atom::PyAtom;
use crate::error::to_py_err;
use crate::matrix::PyMatrix;
use crate::position::PyPosition;
use proteindf_bridge::atom_group::AtomGroup as CoreAtomGroup;
use proteindf_bridge::position::Position;
use pyo3::exceptions::{PyKeyError, PyTypeError};
use pyo3::prelude::*;
use std::collections::HashMap;
use std::str::FromStr;

#[pyclass(name = "AtomGroup", module = "proteindf_bridge_rs")]
#[derive(Clone, Default)]
pub struct PyAtomGroup {
    pub(crate) inner: CoreAtomGroup,
}

impl PyAtomGroup {
    pub fn from_core(group: CoreAtomGroup) -> Self {
        Self { inner: group }
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
impl PyAtomGroup {
    #[new]
    #[pyo3(signature = (*args, **kwargs))]
    pub fn new(
        args: &Bound<'_, pyo3::types::PyTuple>,
        kwargs: Option<&Bound<'_, pyo3::types::PyDict>>,
    ) -> PyResult<Self> {
        let mut group = CoreAtomGroup::new();

        if args.len() == 1 {
            let first = args.get_item(0)?;
            if let Ok(other) = first.extract::<PyRef<PyAtomGroup>>() {
                group = other.inner.clone();
            } else if let Ok(name) = first.extract::<String>() {
                group.name = name;
            } else {
                return Err(PyTypeError::new_err("Unsupported argument for AtomGroup"));
            }
        } else if args.len() > 1 {
            return Err(PyTypeError::new_err(
                "AtomGroup takes at most 1 positional argument",
            ));
        }

        if let Some(kw) = kwargs {
            if let Some(name) = kw.get_item("name")? {
                group.name = name.extract::<String>()?;
            }
        }

        Ok(Self { inner: group })
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
    pub fn path(&self) -> String {
        self.inner.path().to_string()
    }

    #[setter]
    pub fn set_path(&mut self, val: String) {
        self.inner.set_path(val);
    }

    pub fn get_number_of_atoms(&self) -> usize {
        self.inner.get_number_of_atoms()
    }

    pub fn get_number_of_all_atoms(&self) -> usize {
        self.inner.get_number_of_all_atoms()
    }

    pub fn get_number_of_groups(&self) -> usize {
        self.inner.get_number_of_groups()
    }

    pub fn get_number_of_bonds(&self) -> usize {
        self.inner.get_number_of_bonds()
    }

    pub fn get_atom(&self, key: &str) -> Option<PyAtom> {
        self.inner
            .get_atom(key)
            .map(|a| PyAtom::from_core(a.clone()))
    }

    pub fn set_atom(&mut self, key: &str, atom: &PyAtom) {
        self.inner.set_atom(key, atom.inner.clone());
    }

    pub fn has_atom(&self, key: &str) -> bool {
        self.inner.has_atom(key)
    }

    pub fn del_atom(&mut self, key: &str) -> Option<PyAtom> {
        self.inner.remove_atom(key).map(PyAtom::from_core)
    }

    pub fn get_group(&self, key: &str) -> Option<PyAtomGroup> {
        self.inner
            .get_group(key)
            .map(|g| PyAtomGroup::from_core(g.clone()))
    }

    pub fn set_group(&mut self, key: &str, group: &PyAtomGroup) {
        self.inner.set_group(key, group.inner.clone());
    }

    pub fn has_group(&self, key: &str) -> bool {
        self.inner.has_group(key)
    }

    pub fn has_groupkey(&self, key: &str) -> bool {
        self.inner.has_groupkey(key)
    }

    pub fn del_group(&mut self, key: &str) -> Option<PyAtomGroup> {
        self.inner.remove_group(key).map(PyAtomGroup::from_core)
    }

    pub fn atoms(&self) -> Vec<(String, PyAtom)> {
        self.inner
            .atoms()
            .map(|(k, a)| (k.clone(), PyAtom::from_core(a.clone())))
            .collect()
    }

    pub fn groups(&self) -> Vec<(String, PyAtomGroup)> {
        self.inner
            .groups()
            .map(|(k, g)| (k.clone(), PyAtomGroup::from_core(g.clone())))
            .collect()
    }

    pub fn get_atom_list(&self) -> Vec<PyAtom> {
        self.inner
            .get_atom_list()
            .into_iter()
            .map(PyAtom::from_core)
            .collect()
    }

    pub fn get_path_list(&self) -> Vec<String> {
        self.inner.get_path_list()
    }

    pub fn get_atom_kinds_count(&self) -> HashMap<String, usize> {
        self.inner.get_atom_kinds_count()
    }

    pub fn get_formula(&self) -> String {
        self.inner.get_formula()
    }

    pub fn formula(&self) -> String {
        self.inner.get_formula()
    }

    #[getter]
    pub fn weight(&self) -> f64 {
        self.inner
            .get_atom_list()
            .iter()
            .map(|a| a.weight().unwrap_or(0.0))
            .sum()
    }

    pub fn center(&self) -> PyPosition {
        PyPosition::from_core(self.inner.center())
    }

    pub fn r#box(&self) -> (PyPosition, PyPosition) {
        let (bmin, bmax) = self.inner.r#box();
        (PyPosition::from_core(bmin), PyPosition::from_core(bmax))
    }

    pub fn get_box(&self) -> (PyPosition, PyPosition) {
        self.r#box()
    }

    pub fn merge(&mut self, other: &PyAtomGroup) {
        self.inner.merge(&other.inner);
    }

    #[pyo3(signature = (atom1, atom2, order=1))]
    pub fn add_bond(&mut self, atom1: &PyAtom, atom2: &PyAtom, order: usize) {
        self.inner.add_bond(&atom1.inner, &atom2.inner, order);
    }

    pub fn bonds(&self) -> Vec<(String, String, usize)> {
        self.inner
            .bonds()
            .iter()
            .map(|b| (b.atom1_path.clone(), b.atom2_path.clone(), b.order))
            .collect()
    }

    pub fn get_bond_list(&mut self) -> Vec<(String, String, usize)> {
        self.inner
            .get_bond_list()
            .into_iter()
            .map(|b| (b.atom1_path, b.atom2_path, b.order))
            .collect()
    }

    pub fn shift_by(&mut self, dir: &Bound<'_, PyAny>) -> PyResult<()> {
        let d = extract_pos(dir)?;
        self.inner.shift_by(d);
        Ok(())
    }

    pub fn rotate(&mut self, rotmat: &PyMatrix) -> PyResult<()> {
        self.inner.rotate(&rotmat.inner).map_err(to_py_err)
    }

    pub fn sum_of_atomic_number(&self) -> f64 {
        self.inner.sum_of_atomic_number()
    }

    pub fn __getitem__(&self, key: &str) -> PyResult<PyAtomGroup> {
        self.get_group(key)
            .ok_or_else(|| PyKeyError::new_err(format!("Group key not found: {}", key)))
    }

    pub fn __and__(&self, other: &PyAtomGroup) -> Self {
        Self {
            inner: &self.inner & &other.inner,
        }
    }

    pub fn __or__(&self, other: &PyAtomGroup) -> Self {
        Self {
            inner: &self.inner | &other.inner,
        }
    }

    pub fn __xor__(&self, other: &PyAtomGroup) -> Self {
        Self {
            inner: &self.inner ^ &other.inner,
        }
    }

    pub fn __len__(&self) -> usize {
        self.inner.get_number_of_atoms()
    }

    pub fn __repr__(&self) -> String {
        format!(
            "AtomGroup(name='{}', path='{}', atoms={}, groups={})",
            self.inner.name,
            self.inner.path(),
            self.inner.get_number_of_atoms(),
            self.inner.get_number_of_groups()
        )
    }

    pub fn __str__(&self) -> String {
        format!("{}", self.inner)
    }
}
