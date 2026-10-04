// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use pyo3::prelude::*;

use proteindf_bridge::ch_pi::{
    calc_ch_pi_interactions_with_thresholds as core_calc_ch_pi_interactions_with_thresholds,
    ChPiInteraction as CoreChPiInteraction, DEFAULT_MAX_ANGLE_DEG, DEFAULT_MAX_DISTANCE,
};

use crate::atom_group::PyAtomGroup;

/// Represents a detected CH-pi interaction between a carbon atom and an aromatic ring.
#[pyclass(
    name = "ChPiInteraction",
    module = "proteindf_bridge.rs",
    skip_from_py_object
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyChPiInteraction {
    inner: CoreChPiInteraction,
}

impl From<CoreChPiInteraction> for PyChPiInteraction {
    fn from(inner: CoreChPiInteraction) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PyChPiInteraction {
    /// Path of the residue containing the carbon atom (e.g. "/model_1/A/10/").
    #[getter]
    pub fn carbon_path(&self) -> String {
        self.inner.carbon_path.clone()
    }

    /// Name of the carbon atom (e.g. "CG1").
    #[getter]
    pub fn carbon_atom(&self) -> String {
        self.inner.carbon_atom.clone()
    }

    /// Path of the residue containing the aromatic ring (e.g. "/model_1/B/5/").
    #[getter]
    pub fn ring_path(&self) -> String {
        self.inner.ring_path.clone()
    }

    /// Residue name of the aromatic ring (e.g. "HIS", "PHE", "TYR", "TRP").
    #[getter]
    pub fn ring_residue(&self) -> String {
        self.inner.ring_residue.clone()
    }

    /// Distance between carbon atom and ring centroid in Angstroms.
    #[getter]
    pub fn distance(&self) -> f64 {
        self.inner.distance
    }

    /// Angle between centroid-carbon vector and ring unit normal in degrees.
    #[getter]
    pub fn angle(&self) -> f64 {
        self.inner.angle
    }

    fn __repr__(&self) -> String {
        format!(
            "ChPiInteraction(carbon='{}{}', ring='{}' ({}), distance={:.4}, angle={:.2})",
            self.inner.carbon_path,
            self.inner.carbon_atom,
            self.inner.ring_path,
            self.inner.ring_residue,
            self.inner.distance,
            self.inner.angle
        )
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}

/// Detects CH-pi interactions with optional custom distance and angle thresholds (default 4.5 A, 40.0 deg).
#[pyfunction]
#[pyo3(signature = (root, max_distance=None, max_angle_deg=None))]
pub fn calc_ch_pi_interactions(
    root: &PyAtomGroup,
    max_distance: Option<f64>,
    max_angle_deg: Option<f64>,
) -> Vec<PyChPiInteraction> {
    let d = max_distance.unwrap_or(DEFAULT_MAX_DISTANCE);
    let a = max_angle_deg.unwrap_or(DEFAULT_MAX_ANGLE_DEG);
    core_calc_ch_pi_interactions_with_thresholds(&root.inner, d, a)
        .into_iter()
        .map(PyChPiInteraction::from)
        .collect()
}

/// Detects CH-pi interactions with explicit distance and angle thresholds.
#[pyfunction]
pub fn calc_ch_pi_interactions_with_thresholds(
    root: &PyAtomGroup,
    max_distance: f64,
    max_angle_deg: f64,
) -> Vec<PyChPiInteraction> {
    core_calc_ch_pi_interactions_with_thresholds(&root.inner, max_distance, max_angle_deg)
        .into_iter()
        .map(PyChPiInteraction::from)
        .collect()
}
