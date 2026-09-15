// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom_group::PyAtomGroup;
use crate::position::extract_position;
use proteindf_bridge::position::dihedral_angle as core_dihedral_angle;
use proteindf_bridge::ramachandran::{
    calc_phi_psi as core_calc_phi_psi, RamachandranAngle as CoreRamachandranAngle,
};
use pyo3::prelude::*;

#[pyclass(name = "RamachandranAngle", module = "proteindf_bridge_rs")]
#[derive(Clone)]
pub struct PyRamachandranAngle {
    pub(crate) inner: CoreRamachandranAngle,
}

#[pymethods]
impl PyRamachandranAngle {
    #[getter]
    pub fn residue_key(&self) -> &str {
        &self.inner.residue_key
    }

    #[getter]
    pub fn residue_name(&self) -> &str {
        &self.inner.residue_name
    }

    #[getter]
    pub fn phi(&self) -> Option<f64> {
        self.inner.phi
    }

    #[getter]
    pub fn psi(&self) -> Option<f64> {
        self.inner.psi
    }

    pub fn __repr__(&self) -> String {
        format!(
            "RamachandranAngle(residue_key='{}', residue_name='{}', phi={:?}, psi={:?})",
            self.inner.residue_key, self.inner.residue_name, self.inner.phi, self.inner.psi
        )
    }
}

#[pyfunction]
#[pyo3(name = "calc_phi_psi")]
pub fn py_calc_phi_psi(chain: &PyAtomGroup) -> Vec<PyRamachandranAngle> {
    core_calc_phi_psi(&chain.inner)
        .into_iter()
        .map(|inner| PyRamachandranAngle { inner })
        .collect()
}

#[pyfunction]
#[pyo3(name = "dihedral_angle")]
pub fn py_dihedral_angle(
    p1: &Bound<'_, PyAny>,
    p2: &Bound<'_, PyAny>,
    p3: &Bound<'_, PyAny>,
    p4: &Bound<'_, PyAny>,
) -> PyResult<f64> {
    let pos1 = extract_position(p1)?;
    let pos2 = extract_position(p2)?;
    let pos3 = extract_position(p3)?;
    let pos4 = extract_position(p4)?;
    Ok(core_dihedral_angle(&pos1, &pos2, &pos3, &pos4))
}
