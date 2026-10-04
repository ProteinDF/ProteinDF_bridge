// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use pyo3::prelude::*;

use proteindf_bridge::hydrogen_bond::{
    calc_backbone_hbonds as core_calc_backbone_hbonds,
    calc_sidechain_hbonds as core_calc_sidechain_hbonds,
    calc_sidechain_hbonds_with_options as core_calc_sidechain_hbonds_with_options,
    HydrogenBond as CoreHydrogenBond, SidechainHydrogenBond as CoreSidechainHydrogenBond,
};

use crate::atom_group::PyAtomGroup;

/// Represents a detected backbone hydrogen bond.
#[pyclass(
    name = "HydrogenBond",
    module = "proteindf_bridge.rs",
    skip_from_py_object
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyHydrogenBond {
    inner: CoreHydrogenBond,
}

impl From<CoreHydrogenBond> for PyHydrogenBond {
    fn from(inner: CoreHydrogenBond) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PyHydrogenBond {
    /// Donor residue key (e.g. "6").
    #[getter]
    pub fn donor_residue_key(&self) -> String {
        self.inner.donor_residue_key.clone()
    }

    /// Acceptor residue key (e.g. "2").
    #[getter]
    pub fn acceptor_residue_key(&self) -> String {
        self.inner.acceptor_residue_key.clone()
    }

    /// Interaction energy in kcal/mol (Kabsch-Sander model).
    #[getter]
    pub fn energy(&self) -> f64 {
        self.inner.energy
    }

    fn __repr__(&self) -> String {
        format!(
            "HydrogenBond(donor='{}', acceptor='{}', energy={:.4})",
            self.inner.donor_residue_key, self.inner.acceptor_residue_key, self.inner.energy
        )
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}

/// Represents a detected sidechain hydrogen bond.
#[pyclass(
    name = "SidechainHydrogenBond",
    module = "proteindf_bridge.rs",
    skip_from_py_object
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PySidechainHydrogenBond {
    inner: CoreSidechainHydrogenBond,
}

impl From<CoreSidechainHydrogenBond> for PySidechainHydrogenBond {
    fn from(inner: CoreSidechainHydrogenBond) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PySidechainHydrogenBond {
    /// Donor residue path (e.g. "/model_1/A/8/").
    #[getter]
    pub fn donor_path(&self) -> String {
        self.inner.donor_path.clone()
    }

    /// Donor atom name (e.g. "OG1").
    #[getter]
    pub fn donor_atom(&self) -> String {
        self.inner.donor_atom.clone()
    }

    /// Acceptor residue path (e.g. "/model_1/A/4/").
    #[getter]
    pub fn acceptor_path(&self) -> String {
        self.inner.acceptor_path.clone()
    }

    /// Acceptor atom name (e.g. "O").
    #[getter]
    pub fn acceptor_atom(&self) -> String {
        self.inner.acceptor_atom.clone()
    }

    /// Distance between donor and acceptor heavy atoms in Angstroms.
    #[getter]
    pub fn distance(&self) -> f64 {
        self.inner.distance
    }

    /// Hydrogen bond angle (D-H...A) in degrees, if explicit hydrogen was evaluated.
    #[getter]
    pub fn angle(&self) -> Option<f64> {
        self.inner.angle
    }

    fn __repr__(&self) -> String {
        match self.inner.angle {
            Some(ang) => format!(
                "SidechainHydrogenBond(donor='{}{}', acceptor='{}{}', distance={:.4}, angle={:.2})",
                self.inner.donor_path,
                self.inner.donor_atom,
                self.inner.acceptor_path,
                self.inner.acceptor_atom,
                self.inner.distance,
                ang
            ),
            None => format!(
                "SidechainHydrogenBond(donor='{}{}', acceptor='{}{}', distance={:.4}, angle=None)",
                self.inner.donor_path,
                self.inner.donor_atom,
                self.inner.acceptor_path,
                self.inner.acceptor_atom,
                self.inner.distance
            ),
        }
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}

/// Detects backbone hydrogen bonds within a protein chain using the Kabsch-Sander electrostatic model.
#[pyfunction]
pub fn calc_backbone_hbonds(chain: &PyAtomGroup) -> Vec<PyHydrogenBond> {
    core_calc_backbone_hbonds(&chain.inner)
        .into_iter()
        .map(PyHydrogenBond::from)
        .collect()
}

/// Detects sidechain hydrogen bonds in the structure using default heavy-atom-only mode (distance < 3.5 A).
#[pyfunction]
pub fn calc_sidechain_hbonds(root: &PyAtomGroup) -> Vec<PySidechainHydrogenBond> {
    core_calc_sidechain_hbonds(&root.inner)
        .into_iter()
        .map(PySidechainHydrogenBond::from)
        .collect()
}

/// Detects sidechain hydrogen bonds with optional explicit hydrogen angle evaluation.
#[pyfunction]
pub fn calc_sidechain_hbonds_with_options(
    root: &PyAtomGroup,
    explicit_h_mode: bool,
) -> Vec<PySidechainHydrogenBond> {
    core_calc_sidechain_hbonds_with_options(&root.inner, explicit_h_mode)
        .into_iter()
        .map(PySidechainHydrogenBond::from)
        .collect()
}
