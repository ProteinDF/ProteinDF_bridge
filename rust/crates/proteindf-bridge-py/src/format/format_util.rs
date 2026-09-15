// Copyright (C) 2014 The ProteinDF development team.
// see also AUTHORS and README if provided.
//
// This file is a part of the ProteinDF software package.

use crate::atom_group::PyAtomGroup;
use proteindf_bridge::format::Format as CoreFormat;
use pyo3::prelude::*;

#[pyclass(name = "Format", module = "proteindf_bridge_rs")]
pub struct PyFormat;

#[pymethods]
impl PyFormat {
    #[staticmethod]
    pub fn is_residue(model: &Bound<'_, PyAny>) -> bool {
        if let Ok(ag) = model.extract::<PyRef<PyAtomGroup>>() {
            CoreFormat::is_residue(&ag.inner)
        } else {
            false
        }
    }

    #[staticmethod]
    pub fn is_chain(model: &Bound<'_, PyAny>) -> bool {
        if let Ok(ag) = model.extract::<PyRef<PyAtomGroup>>() {
            CoreFormat::is_chain(&ag.inner)
        } else {
            false
        }
    }

    #[staticmethod]
    pub fn is_protein(model: &Bound<'_, PyAny>) -> bool {
        if let Ok(ag) = model.extract::<PyRef<PyAtomGroup>>() {
            CoreFormat::is_protein(&ag.inner)
        } else {
            false
        }
    }

    #[staticmethod]
    pub fn is_models(model: &Bound<'_, PyAny>) -> bool {
        if let Ok(ag) = model.extract::<PyRef<PyAtomGroup>>() {
            CoreFormat::is_models(&ag.inner)
        } else {
            false
        }
    }
}
