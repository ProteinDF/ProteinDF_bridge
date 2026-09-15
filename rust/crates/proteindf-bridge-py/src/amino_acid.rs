// Copyright (C) 2014 The ProteinDF development team.
// see also AUTHORS and README if provided.
//
// This file is a part of the ProteinDF software package.

use crate::atom_group::PyAtomGroup;
use proteindf_bridge::amino_acid::AminoAcid as CoreAminoAcid;
use pyo3::prelude::*;

#[pyclass(name = "AminoAcid", module = "proteindf_bridge_rs")]
#[derive(Clone, Copy, Default)]
pub struct PyAminoAcid;

#[pymethods]
impl PyAminoAcid {
    #[new]
    pub fn new() -> Self {
        Self
    }

    #[staticmethod]
    pub fn is_aminoacid(atomgroup: &Bound<'_, PyAny>) -> bool {
        if let Ok(ag) = atomgroup.extract::<PyRef<PyAtomGroup>>() {
            CoreAminoAcid::is_aminoacid(&ag.inner)
        } else if let Ok(name_attr) = atomgroup.getattr("name") {
            if let Ok(name) = name_attr.extract::<String>() {
                CoreAminoAcid::is_aminoacid_name(&name)
            } else {
                false
            }
        } else if let Ok(name) = atomgroup.extract::<String>() {
            CoreAminoAcid::is_aminoacid_name(&name)
        } else {
            false
        }
    }

    #[staticmethod]
    pub fn is_aminoacid_name(name: &str) -> bool {
        CoreAminoAcid::is_aminoacid_name(name)
    }

    #[staticmethod]
    pub fn aa_list() -> Vec<String> {
        CoreAminoAcid::aa_list()
            .iter()
            .map(|s| s.to_string())
            .collect()
    }

    #[classattr]
    #[allow(non_snake_case)]
    fn _AA_list() -> Vec<String> {
        CoreAminoAcid::aa_list()
            .iter()
            .map(|s| s.to_string())
            .collect()
    }
}
