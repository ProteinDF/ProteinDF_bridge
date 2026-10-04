// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyType};

use proteindf_bridge::interaction_set::{
    Interaction as CoreInteraction, InteractionKind, InteractionSet as CoreInteractionSet,
};

use crate::atom_group::PyAtomGroup;
use crate::error::{to_py_err, BrValueError};

fn kind_to_str(kind: InteractionKind) -> &'static str {
    match kind {
        InteractionKind::Disulfide => "disulfide",
        InteractionKind::SaltBridge => "salt_bridge",
        InteractionKind::HydrogenBond => "hydrogen_bond",
        InteractionKind::ChPi => "ch_pi",
    }
}

fn parse_kind(s: &str) -> PyResult<InteractionKind> {
    match s.trim().to_ascii_lowercase().as_str() {
        "disulfide" => Ok(InteractionKind::Disulfide),
        "salt_bridge" | "saltbridge" => Ok(InteractionKind::SaltBridge),
        "hydrogen_bond" | "hydrogenbond" => Ok(InteractionKind::HydrogenBond),
        "ch_pi" | "chpi" => Ok(InteractionKind::ChPi),
        other => Err(BrValueError::new_err(format!(
            "Unknown interaction kind '{other}', expected 'disulfide', 'salt_bridge', 'hydrogen_bond', or 'ch_pi'"
        ))),
    }
}

/// Represents a detected non-covalent or disulfide interaction.
#[pyclass(
    name = "Interaction",
    module = "proteindf_bridge.rs",
    skip_from_py_object
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyInteraction {
    inner: CoreInteraction,
}

impl From<CoreInteraction> for PyInteraction {
    fn from(inner: CoreInteraction) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PyInteraction {
    /// Interaction classification ("disulfide", "salt_bridge", "hydrogen_bond", or "ch_pi").
    #[getter]
    pub fn kind(&self) -> String {
        kind_to_str(self.inner.kind).to_string()
    }

    /// Participating atom or residue paths.
    ///
    /// Returns a new list on each access to protect against external mutation.
    #[getter]
    pub fn atoms(&self) -> Vec<String> {
        self.inner.atoms.clone()
    }

    /// Interaction distance in Angstroms, if measured.
    #[getter]
    pub fn distance(&self) -> Option<f64> {
        self.inner.distance
    }

    /// Interaction angle in degrees, if measured.
    #[getter]
    pub fn angle(&self) -> Option<f64> {
        self.inner.angle
    }

    /// Specific role or description (e.g. "backbone (E=-1.688 kcal/mol)", "sidechain", "disulfide", "CH-pi (HIS)").
    #[getter]
    pub fn donor_acceptor_role(&self) -> Option<String> {
        self.inner.donor_acceptor_role.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "Interaction(kind='{}', atoms={:?}, distance={:?}, angle={:?}, role={:?})",
            kind_to_str(self.inner.kind),
            self.inner.atoms,
            self.inner.distance,
            self.inner.angle,
            self.inner.donor_acceptor_role
        )
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}

/// Aggregation set of detected interactions in a molecular structure.
#[pyclass(
    name = "InteractionSet",
    module = "proteindf_bridge.rs",
    from_py_object
)]
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PyInteractionSet {
    pub(crate) inner: CoreInteractionSet,
}

impl From<CoreInteractionSet> for PyInteractionSet {
    fn from(inner: CoreInteractionSet) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PyInteractionSet {
    /// Creates a new empty `InteractionSet`.
    #[new]
    pub fn new() -> Self {
        Self {
            inner: CoreInteractionSet::new(),
        }
    }

    /// Detects all supported interactions in the given model:
    /// 1. Disulfide bonds
    /// 2. Salt bridges
    /// 3. Hydrogen bonds (backbone + sidechain)
    /// 4. CH-pi interactions
    ///
    /// Optional `ch_pi_max_distance` (default 4.5 A) and `ch_pi_max_angle_deg` (default 40.0 deg)
    /// can be specified to customize CH-pi thresholds.
    #[classmethod]
    #[pyo3(signature = (model, ch_pi_max_distance=None, ch_pi_max_angle_deg=None))]
    pub fn detect_all(
        _cls: &Bound<'_, PyType>,
        model: &PyAtomGroup,
        ch_pi_max_distance: Option<f64>,
        ch_pi_max_angle_deg: Option<f64>,
    ) -> Self {
        let set =
            CoreInteractionSet::detect_all(&model.inner, ch_pi_max_distance, ch_pi_max_angle_deg);
        Self { inner: set }
    }

    /// Total number of interactions in the set.
    pub fn __len__(&self) -> usize {
        self.inner.len()
    }

    /// Returns `True` if there are no interactions.
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// List of all detected interactions.
    ///
    /// Returns a new list on each access to protect against external mutation.
    #[getter]
    pub fn interactions(&self) -> Vec<PyInteraction> {
        self.inner
            .interactions
            .iter()
            .cloned()
            .map(PyInteraction::from)
            .collect()
    }

    /// Filters interactions matching a specific interaction kind.
    ///
    /// Returns a new list of matching interactions.
    pub fn filter_by_kind(&self, kind: &str) -> PyResult<Vec<PyInteraction>> {
        let k = parse_kind(kind)?;
        Ok(self
            .inner
            .filter_by_kind(k)
            .into_iter()
            .cloned()
            .map(PyInteraction::from)
            .collect())
    }

    /// Counts interactions matching a specific interaction kind.
    pub fn count_by_kind(&self, kind: &str) -> PyResult<usize> {
        let k = parse_kind(kind)?;
        Ok(self.inner.count_by_kind(k))
    }

    /// Serializes `InteractionSet` to MessagePack bytes.
    pub fn to_msgpack<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        let bytes = self.inner.to_msgpack().map_err(to_py_err)?;
        Ok(PyBytes::new(py, &bytes))
    }

    /// Deserializes `InteractionSet` from MessagePack bytes.
    #[classmethod]
    pub fn from_msgpack(_cls: &Bound<'_, PyType>, bytes: &[u8]) -> PyResult<Self> {
        let set = CoreInteractionSet::from_msgpack(bytes).map_err(to_py_err)?;
        Ok(Self { inner: set })
    }

    /// Saves `InteractionSet` as a MessagePack file.
    pub fn save_msgpack(&self, path: &str) -> PyResult<()> {
        self.inner.save_msgpack(path).map_err(to_py_err)
    }

    /// Loads `InteractionSet` from a MessagePack file.
    #[classmethod]
    pub fn load_msgpack(_cls: &Bound<'_, PyType>, path: &str) -> PyResult<Self> {
        let set = CoreInteractionSet::load_msgpack(path).map_err(to_py_err)?;
        Ok(Self { inner: set })
    }

    /// Serializes `InteractionSet` to a YAML string.
    pub fn to_yaml(&self) -> PyResult<String> {
        self.inner.to_yaml().map_err(to_py_err)
    }

    /// Deserializes `InteractionSet` from a YAML string.
    #[classmethod]
    pub fn from_yaml(_cls: &Bound<'_, PyType>, yaml_str: &str) -> PyResult<Self> {
        let set = CoreInteractionSet::from_yaml(yaml_str).map_err(to_py_err)?;
        Ok(Self { inner: set })
    }

    /// Saves `InteractionSet` as a YAML file.
    pub fn save_yaml(&self, path: &str) -> PyResult<()> {
        self.inner.save_yaml(path).map_err(to_py_err)
    }

    /// Loads `InteractionSet` from a YAML file.
    #[classmethod]
    pub fn load_yaml(_cls: &Bound<'_, PyType>, path: &str) -> PyResult<Self> {
        let set = CoreInteractionSet::load_yaml(path).map_err(to_py_err)?;
        Ok(Self { inner: set })
    }

    fn __repr__(&self) -> String {
        format!(
            "InteractionSet(total={}, disulfide={}, salt_bridge={}, hydrogen_bond={}, ch_pi={})",
            self.inner.len(),
            self.inner.count_by_kind(InteractionKind::Disulfide),
            self.inner.count_by_kind(InteractionKind::SaltBridge),
            self.inner.count_by_kind(InteractionKind::HydrogenBond),
            self.inner.count_by_kind(InteractionKind::ChPi),
        )
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}
