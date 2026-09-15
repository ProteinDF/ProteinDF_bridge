// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom_group::PyAtomGroup;
use proteindf_bridge::ion_pair::IonPair as CoreIonPair;
use pyo3::prelude::*;

#[pyclass(name = "IonPair", module = "proteindf_bridge_rs")]
#[derive(Clone)]
pub struct PyIonPair {
    pub(crate) inner: CoreIonPair,
}

#[pymethods]
impl PyIonPair {
    #[new]
    pub fn new(model: &PyAtomGroup) -> Self {
        Self {
            inner: CoreIonPair::new(&model.inner),
        }
    }

    pub fn get_ion_pairs(&self) -> Vec<(String, String, String, String)> {
        self.inner
            .get_ion_pairs()
            .into_iter()
            .map(|r| r.into_tuple())
            .collect()
    }

    #[staticmethod]
    pub fn find_ion_pairs(model: &PyAtomGroup) -> Vec<(String, String, String, String)> {
        CoreIonPair::find_ion_pairs(&model.inner)
            .into_iter()
            .map(|r| r.into_tuple())
            .collect()
    }
}
