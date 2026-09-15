// Copyright (C) 2014 The ProteinDF development team.
// see also AUTHORS and README if provided.
//
// This file is a part of the ProteinDF software package.

use crate::atom::PyAtom;
use crate::atom_group::PyAtomGroup;
use crate::error::to_py_err;
use crate::position::PyPosition;
use proteindf_bridge::atom_group::Selector;
use proteindf_bridge::position::Position;
use proteindf_bridge::selector::{
    SelectAtom as CoreSelectAtom, SelectAtomGroup as CoreSelectAtomGroup,
    SelectName as CoreSelectName, SelectPath as CoreSelectPath,
    SelectPathRegex as CoreSelectPathRegex, SelectPathSimple as CoreSelectPathSimple,
    SelectPathWildcard as CoreSelectPathWildcard, SelectRange as CoreSelectRange,
    SelectSymbol as CoreSelectSymbol,
};
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use std::str::FromStr;

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

// -----------------------------------------------------------------------------
// 1. Select_Symbol
// -----------------------------------------------------------------------------
#[pyclass(name = "Select_Symbol", module = "proteindf_bridge_rs")]
#[derive(Clone)]
pub struct PySelectSymbol {
    pub(crate) inner: CoreSelectSymbol,
}

#[pymethods]
impl PySelectSymbol {
    #[new]
    pub fn new(symbol: &str) -> Self {
        Self {
            inner: CoreSelectSymbol::new(symbol),
        }
    }

    pub fn is_match(&self, obj: &Bound<'_, PyAny>) -> bool {
        if let Ok(atom) = obj.extract::<PyRef<PyAtom>>() {
            self.inner.is_match_atom(&atom.inner)
        } else {
            false
        }
    }
}

// -----------------------------------------------------------------------------
// 2. Select_Name
// -----------------------------------------------------------------------------
#[pyclass(name = "Select_Name", module = "proteindf_bridge_rs")]
#[derive(Clone)]
pub struct PySelectName {
    pub(crate) inner: CoreSelectName,
}

#[pymethods]
impl PySelectName {
    #[new]
    pub fn new(query: &str) -> Self {
        Self {
            inner: CoreSelectName::new(query),
        }
    }

    pub fn is_match(&self, obj: &Bound<'_, PyAny>) -> bool {
        if let Ok(atom) = obj.extract::<PyRef<PyAtom>>() {
            self.inner.is_match_atom(&atom.inner)
        } else if let Ok(ag) = obj.extract::<PyRef<PyAtomGroup>>() {
            self.inner.is_match_group(&ag.inner)
        } else {
            false
        }
    }
}

// -----------------------------------------------------------------------------
// 3. Select_Path_simple
// -----------------------------------------------------------------------------
#[pyclass(name = "Select_Path_simple", module = "proteindf_bridge_rs")]
#[derive(Clone)]
pub struct PySelectPathSimple {
    pub(crate) inner: CoreSelectPathSimple,
}

#[pymethods]
impl PySelectPathSimple {
    #[new]
    pub fn new(query: &str) -> Self {
        Self {
            inner: CoreSelectPathSimple::new(query),
        }
    }

    pub fn is_match(&self, obj: &Bound<'_, PyAny>) -> bool {
        if let Ok(atom) = obj.extract::<PyRef<PyAtom>>() {
            self.inner.is_match_atom(&atom.inner)
        } else if let Ok(ag) = obj.extract::<PyRef<PyAtomGroup>>() {
            self.inner.is_match_group(&ag.inner)
        } else {
            false
        }
    }
}

// -----------------------------------------------------------------------------
// 4. Select_Path_wildcard
// -----------------------------------------------------------------------------
#[pyclass(name = "Select_Path_wildcard", module = "proteindf_bridge_rs")]
#[derive(Clone)]
pub struct PySelectPathWildcard {
    pub(crate) inner: CoreSelectPathWildcard,
}

#[pymethods]
impl PySelectPathWildcard {
    #[new]
    pub fn new(query: &str) -> PyResult<Self> {
        let inner = CoreSelectPathWildcard::new(query).map_err(to_py_err)?;
        Ok(Self { inner })
    }

    pub fn is_match(&self, obj: &Bound<'_, PyAny>) -> bool {
        if let Ok(atom) = obj.extract::<PyRef<PyAtom>>() {
            self.inner.is_match_atom(&atom.inner)
        } else if let Ok(ag) = obj.extract::<PyRef<PyAtomGroup>>() {
            self.inner.is_match_group(&ag.inner)
        } else {
            false
        }
    }
}

// -----------------------------------------------------------------------------
// 5. Select_PathRegex
// -----------------------------------------------------------------------------
#[pyclass(name = "Select_PathRegex", module = "proteindf_bridge_rs")]
#[derive(Clone)]
pub struct PySelectPathRegex {
    pub(crate) inner: CoreSelectPathRegex,
}

#[pymethods]
impl PySelectPathRegex {
    #[new]
    pub fn new(query: &str) -> PyResult<Self> {
        let inner = CoreSelectPathRegex::new(query).map_err(to_py_err)?;
        Ok(Self { inner })
    }

    pub fn is_match(&self, obj: &Bound<'_, PyAny>) -> bool {
        if let Ok(atom) = obj.extract::<PyRef<PyAtom>>() {
            self.inner.is_match_atom(&atom.inner)
        } else if let Ok(ag) = obj.extract::<PyRef<PyAtomGroup>>() {
            self.inner.is_match_group(&ag.inner)
        } else {
            false
        }
    }
}

// -----------------------------------------------------------------------------
// 6. Select_Path
// -----------------------------------------------------------------------------
#[pyclass(name = "Select_Path", module = "proteindf_bridge_rs")]
#[derive(Clone)]
pub struct PySelectPath {
    pub(crate) inner: CoreSelectPath,
}

#[pymethods]
impl PySelectPath {
    #[new]
    #[pyo3(signature = (query, use_wildcard=false))]
    pub fn new(query: &str, use_wildcard: bool) -> PyResult<Self> {
        let inner = CoreSelectPath::new(query, use_wildcard).map_err(to_py_err)?;
        Ok(Self { inner })
    }

    pub fn is_match(&self, obj: &Bound<'_, PyAny>) -> bool {
        if let Ok(atom) = obj.extract::<PyRef<PyAtom>>() {
            self.inner.is_match_atom(&atom.inner)
        } else if let Ok(ag) = obj.extract::<PyRef<PyAtomGroup>>() {
            self.inner.is_match_group(&ag.inner)
        } else {
            false
        }
    }
}

// -----------------------------------------------------------------------------
// 7. Select_Range
// -----------------------------------------------------------------------------
#[pyclass(name = "Select_Range", module = "proteindf_bridge_rs")]
#[derive(Clone)]
pub struct PySelectRange {
    pub(crate) inner: CoreSelectRange,
}

#[pymethods]
impl PySelectRange {
    #[new]
    pub fn new(pos: &Bound<'_, PyAny>, d: f64) -> PyResult<Self> {
        let p = extract_pos(pos)?;
        Ok(Self {
            inner: CoreSelectRange::new(p, d),
        })
    }

    pub fn is_match(&self, obj: &Bound<'_, PyAny>) -> bool {
        if let Ok(atom) = obj.extract::<PyRef<PyAtom>>() {
            self.inner.is_match_atom(&atom.inner)
        } else {
            false
        }
    }
}

// -----------------------------------------------------------------------------
// 8. Select_Atom
// -----------------------------------------------------------------------------
#[pyclass(name = "Select_Atom", module = "proteindf_bridge_rs")]
#[derive(Clone)]
pub struct PySelectAtom {
    pub(crate) inner: CoreSelectAtom,
}

#[pymethods]
impl PySelectAtom {
    #[new]
    pub fn new(atom: &PyAtom, distance: f64) -> Self {
        Self {
            inner: CoreSelectAtom::new(atom.inner.clone(), distance),
        }
    }

    pub fn is_match(&self, obj: &Bound<'_, PyAny>) -> bool {
        if let Ok(atom) = obj.extract::<PyRef<PyAtom>>() {
            self.inner.is_match_atom(&atom.inner)
        } else {
            false
        }
    }
}

// -----------------------------------------------------------------------------
// 9. Select_AtomGroup
// -----------------------------------------------------------------------------
#[pyclass(name = "Select_AtomGroup", module = "proteindf_bridge_rs")]
#[derive(Clone)]
pub struct PySelectAtomGroup {
    pub(crate) inner: CoreSelectAtomGroup,
}

#[pymethods]
impl PySelectAtomGroup {
    #[new]
    #[pyo3(signature = (ref_atomgroup, range=1.0e-5))]
    pub fn new(ref_atomgroup: &PyAtomGroup, range: f64) -> Self {
        Self {
            inner: CoreSelectAtomGroup::new(&ref_atomgroup.inner, range),
        }
    }

    pub fn is_match(&self, obj: &Bound<'_, PyAny>) -> bool {
        if let Ok(atom) = obj.extract::<PyRef<PyAtom>>() {
            self.inner.is_match_atom(&atom.inner)
        } else {
            false
        }
    }
}
