// Copyright (C) 2014 The ProteinDF development team.
// see also AUTHORS and README if provided.
//
// This file is a part of the ProteinDF software package.

use crate::error::to_py_err;
use proteindf_bridge::periodic_table::PeriodicTable as CorePeriodicTable;
use pyo3::prelude::*;

#[pyclass(name = "PeriodicTable", module = "proteindf_bridge_rs")]
pub struct PyPeriodicTable;

fn resolve_atomic_number(atom: &Bound<'_, PyAny>) -> PyResult<usize> {
    if let Ok(num) = atom.extract::<usize>() {
        Ok(num)
    } else if let Ok(sym) = atom.extract::<String>() {
        CorePeriodicTable::get_atomic_number(&sym).map_err(to_py_err)
    } else {
        Err(pyo3::exceptions::PyTypeError::new_err(
            "Expected int or str for atom",
        ))
    }
}

#[pymethods]
impl PyPeriodicTable {
    #[staticmethod]
    pub fn get_num_of_atoms() -> usize {
        CorePeriodicTable::get_num_of_atoms()
    }

    #[staticmethod]
    pub fn get_symbol(atomic_number: usize) -> PyResult<String> {
        CorePeriodicTable::get_symbol(atomic_number)
            .map(|s| s.to_string())
            .map_err(to_py_err)
    }

    #[staticmethod]
    pub fn get_atomic_number(symbol: &str) -> PyResult<usize> {
        CorePeriodicTable::get_atomic_number(symbol).map_err(to_py_err)
    }

    #[staticmethod]
    pub fn symbol(atomic_number: usize) -> PyResult<String> {
        Self::get_symbol(atomic_number)
    }

    #[staticmethod]
    pub fn atomic_number(symbol: &str) -> PyResult<usize> {
        Self::get_atomic_number(symbol)
    }

    #[staticmethod]
    pub fn vdw(atom: &Bound<'_, PyAny>) -> PyResult<f64> {
        let num = resolve_atomic_number(atom)?;
        CorePeriodicTable::vdw(num).map_err(to_py_err)
    }

    #[staticmethod]
    pub fn get_vdw(atom: &Bound<'_, PyAny>) -> PyResult<f64> {
        Self::vdw(atom)
    }

    #[staticmethod]
    pub fn atomic_weight(atom: &Bound<'_, PyAny>) -> PyResult<f64> {
        let num = resolve_atomic_number(atom)?;
        CorePeriodicTable::atomic_weight(num).map_err(to_py_err)
    }

    #[staticmethod]
    pub fn get_weight(atom: &Bound<'_, PyAny>) -> PyResult<f64> {
        Self::atomic_weight(atom)
    }
}
