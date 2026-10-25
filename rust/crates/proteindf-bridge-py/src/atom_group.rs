// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom::PyAtom;
use crate::ccd_templates::PyCcdTemplateDb;
use crate::complete_atoms::PyCompletionReport;
use crate::error::to_py_err;
use crate::hydrogenation::PyOverallHydrogenationReport;
use crate::missing_atoms::PyMissingAtomsReport;
use crate::position::PyPosition;
use crate::schema::PySchemaViolation;
use indexmap::IndexMap;
use proteindf_bridge::atom_group::{AtomGroup as CoreAtomGroup, BondRecord, Selector};
use proteindf_bridge::ccd_templates::CcdTemplateDb as CoreCcdTemplateDb;
use proteindf_bridge::complete_atoms::CompleteOptions;
use proteindf_bridge::missing_atoms::MissingAtomsOptions;
use proteindf_bridge::position::Position;
use proteindf_bridge::secondary_structure::{
    apply_secondary_structure as core_apply_secondary_structure, SsCode,
};
use pyo3::exceptions::{PyDeprecationWarning, PyKeyError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PySet, PyString};
use std::collections::HashMap;
use std::str::FromStr;

/// A hierarchical tree representation of molecular structures.
///
/// An `AtomGroup` represents a node in a molecular hierarchy containing child groups
/// (`AtomGroup`) and atoms (`Atom`). Child groups typically represent structural units
/// such as models, chains, or residues, while atoms represent the individual atomic sites.
///
/// Objects retrieved from the tree (via indexing, :meth:`get_group`, or :meth:`get_atom`)
/// are direct references to the nodes within the tree; modifying them mutates the tree
/// in-place. Conversely, inserting objects via :meth:`set_group` or :meth:`set_atom`
/// deep-copies the provided objects into the tree.
///
/// The constructor accepts:
/// - A group name string (e.g. ``AtomGroup("model_1")`` or ``AtomGroup(name="model_1")``).
/// - Another `AtomGroup` to clone (e.g. ``AtomGroup(other_group)``).
/// - A dictionary representing raw tree data (e.g. ``AtomGroup(raw_dict)``).
/// - No arguments (creates an empty group with default settings).
#[pyclass(
    name = "AtomGroup",
    module = "proteindf_bridge.rs",
    weakref,
    from_py_object
)]
pub struct PyAtomGroup {
    pub(crate) name: String,
    pub(crate) atoms: IndexMap<String, Py<PyAtom>>,
    pub(crate) groups: IndexMap<String, Py<PyAtomGroup>>,
    pub(crate) bonds: Vec<BondRecord>,
    pub(crate) secondary_structure: Option<SsCode>,
    pub(crate) parent: Option<Py<PyAtomGroup>>,
    pub(crate) my_key: Option<String>,
    pub(crate) root_path: Option<String>,
    pub(crate) sort_atoms: Option<String>,
    pub(crate) sort_groups: Option<String>,
}

impl Default for PyAtomGroup {
    fn default() -> Self {
        Self {
            name: String::new(),
            atoms: IndexMap::new(),
            groups: IndexMap::new(),
            bonds: Vec::new(),
            secondary_structure: None,
            parent: None,
            my_key: None,
            root_path: Some("/".to_string()),
            sort_atoms: None,
            sort_groups: None,
        }
    }
}

impl Clone for PyAtomGroup {
    fn clone(&self) -> Self {
        Python::attach(|py| self.deep_clone(py).expect("failed to clone PyAtomGroup"))
    }
}

impl PyAtomGroup {
    pub(crate) fn get_parent(&self, py: Python<'_>) -> Option<Py<PyAtomGroup>> {
        self.parent.as_ref().map(|p| p.clone_ref(py))
    }

    pub(crate) fn compute_path(&self, py: Python<'_>) -> String {
        if let Some(parent_py) = self.get_parent(py) {
            let parent_path = parent_py.borrow(py).compute_path(py);
            if let Some(key) = &self.my_key {
                format!("{}{}/", parent_path, key)
            } else {
                format!("{}{}/", parent_path, self.name)
            }
        } else {
            let mut p = self.root_path.clone().unwrap_or_else(|| "/".to_string());
            if p.is_empty() || !p.ends_with('/') {
                p.push('/');
            }
            p
        }
    }

    pub(crate) fn compute_depth(&self, py: Python<'_>) -> usize {
        if let Some(parent_py) = self.get_parent(py) {
            parent_py.borrow(py).compute_depth(py) + 1
        } else {
            0
        }
    }

    pub(crate) fn for_each_child_group<F>(&self, py: Python<'_>, mut f: F)
    where
        F: FnMut(&Bound<'_, PyAtomGroup>),
    {
        if self.sort_groups.as_deref() == Some("nice") {
            let mut grp_keys: Vec<String> = self.groups.keys().cloned().collect();
            proteindf_bridge::brd::sort_nicely(&mut grp_keys);
            for k in grp_keys {
                if let Some(g) = self.groups.get(&k) {
                    f(g.bind(py));
                }
            }
        } else {
            for g in self.groups.values() {
                f(g.bind(py));
            }
        }
    }

    pub(crate) fn for_each_direct_atom<F>(&self, py: Python<'_>, mut f: F)
    where
        F: FnMut(&Bound<'_, PyAtom>),
    {
        if self.sort_atoms.as_deref() == Some("nice") {
            let mut atm_keys: Vec<String> = self.atoms.keys().cloned().collect();
            proteindf_bridge::brd::sort_nicely(&mut atm_keys);
            for k in atm_keys {
                if let Some(a) = self.atoms.get(&k) {
                    f(a.bind(py));
                }
            }
        } else {
            for a in self.atoms.values() {
                f(a.bind(py));
            }
        }
    }

    pub(crate) fn for_each_atom_recursive<F>(&self, py: Python<'_>, f: &mut F)
    where
        F: FnMut(&Bound<'_, PyAtom>),
    {
        self.for_each_child_group(py, |g| {
            g.borrow().for_each_atom_recursive(py, f);
        });
        self.for_each_direct_atom(py, |a| {
            f(a);
        });
    }

    pub(crate) fn get_number_of_all_atoms_py(&self, py: Python<'_>) -> usize {
        let mut count = self.atoms.len();
        self.for_each_child_group(py, |g| {
            count += g.borrow().get_number_of_all_atoms_py(py);
        });
        count
    }

    pub(crate) fn center_py(&self, py: Python<'_>) -> proteindf_bridge::position::Position {
        let mut sum = proteindf_bridge::position::Position::default();
        let mut total_atoms = 0;

        self.for_each_child_group(py, |g| {
            let g_ref = g.borrow();
            let n = g_ref.get_number_of_all_atoms_py(py);
            if n > 0 {
                sum += g_ref.center_py(py) * (n as f64);
                total_atoms += n;
            }
        });

        self.for_each_direct_atom(py, |a| {
            sum += a.borrow().xyz.borrow(py).inner;
            total_atoms += 1;
        });

        if total_atoms > 0 {
            sum /= total_atoms as f64;
        }
        sum
    }

    pub(crate) fn box_py(
        &self,
        py: Python<'_>,
    ) -> (
        proteindf_bridge::position::Position,
        proteindf_bridge::position::Position,
    ) {
        let mut box_min = self.center_py(py);
        let mut box_max = box_min;

        self.for_each_child_group(py, |g| {
            let (grp_min, grp_max) = g.borrow().box_py(py);
            box_min.x = box_min.x.min(grp_min.x);
            box_min.y = box_min.y.min(grp_min.y);
            box_min.z = box_min.z.min(grp_min.z);
            box_max.x = box_max.x.max(grp_max.x);
            box_max.y = box_max.y.max(grp_max.y);
            box_max.z = box_max.z.max(grp_max.z);
        });

        self.for_each_direct_atom(py, |a| {
            let xyz = a.borrow().xyz.borrow(py).inner;
            box_min.x = box_min.x.min(xyz.x);
            box_min.y = box_min.y.min(xyz.y);
            box_min.z = box_min.z.min(xyz.z);
            box_max.x = box_max.x.max(xyz.x);
            box_max.y = box_max.y.max(xyz.y);
            box_max.z = box_max.z.max(xyz.z);
        });

        (box_min, box_max)
    }

    pub(crate) fn to_core(&self, py: Python<'_>) -> CoreAtomGroup {
        let path = self.compute_path(py);
        let depth = self.compute_depth(py);
        self.to_core_with_path(py, path, depth)
    }

    pub(crate) fn rotate_core(
        &mut self,
        py: Python<'_>,
        rotmat: &proteindf_bridge::matrix::Matrix,
    ) -> PyResult<()> {
        for a in self.atoms.values() {
            a.borrow(py)
                .xyz
                .borrow_mut(py)
                .inner
                .rotate(rotmat)
                .map_err(to_py_err)?;
        }
        for g in self.groups.values() {
            g.borrow_mut(py).rotate_core(py, rotmat)?;
        }
        Ok(())
    }

    pub(crate) fn to_core_with_path(
        &self,
        py: Python<'_>,
        path: String,
        depth: usize,
    ) -> CoreAtomGroup {
        let mut core = CoreAtomGroup::with_name(&self.name);
        core.set_path_with_depth(path, depth);
        let current_path = core.path().to_string();
        for (k, g) in &self.groups {
            let child_path = format!("{}{}/", current_path, k);
            let child_core = g.borrow(py).to_core_with_path(py, child_path, depth + 1);
            core.insert_group_with_paths(k, child_core);
        }
        for (k, a) in &self.atoms {
            let mut atom_core = a.borrow(py).to_core(py);
            atom_core.path = format!("{}{}", current_path, k);
            core.insert_atom_with_path(k, atom_core);
        }
        core.set_bonds(self.bonds.clone());
        core.set_secondary_structure(self.secondary_structure);
        core
    }

    pub(crate) fn from_core(py: Python<'_>, ag: CoreAtomGroup) -> PyResult<Self> {
        let root_path = ag.path().to_string();
        let name = ag.name.clone();
        let bonds = ag.bonds().to_vec();
        let secondary_structure = ag.secondary_structure();

        let mut atoms = IndexMap::with_capacity(ag.get_number_of_atoms());
        for (k, a) in ag.atoms() {
            let py_atom = PyAtom::from_core(py, a.clone())?;
            atoms.insert(k.clone(), Py::new(py, py_atom)?);
        }

        let mut groups = IndexMap::with_capacity(ag.get_number_of_groups());
        for (k, child_core) in ag.groups() {
            let child_py = Self::from_core_py(py, child_core, None, Some(k.clone()), None)?;
            groups.insert(k.clone(), child_py);
        }

        Ok(Self {
            name,
            atoms,
            groups,
            bonds,
            secondary_structure,
            parent: None,
            my_key: None,
            root_path: Some(root_path),
            sort_atoms: None,
            sort_groups: None,
        })
    }

    pub(crate) fn from_core_py(
        py: Python<'_>,
        ag: &CoreAtomGroup,
        parent: Option<Py<PyAtomGroup>>,
        my_key: Option<String>,
        root_path: Option<String>,
    ) -> PyResult<Py<PyAtomGroup>> {
        let mut atoms = IndexMap::with_capacity(ag.get_number_of_atoms());
        for (k, a) in ag.atoms() {
            let py_atom = PyAtom::from_core(py, a.clone())?;
            atoms.insert(k.clone(), Py::new(py, py_atom)?);
        }

        let new_group = PyAtomGroup {
            name: ag.name.clone(),
            atoms,
            groups: IndexMap::new(),
            bonds: ag.bonds().to_vec(),
            secondary_structure: ag.secondary_structure(),
            parent,
            my_key,
            root_path,
            sort_atoms: None,
            sort_groups: None,
        };
        let new_group_py = Py::new(py, new_group)?;

        for (k, child_core) in ag.groups() {
            let child_py = Self::from_core_py(
                py,
                child_core,
                Some(new_group_py.clone_ref(py)),
                Some(k.clone()),
                None,
            )?;
            new_group_py
                .borrow_mut(py)
                .groups
                .insert(k.clone(), child_py);
        }

        Ok(new_group_py)
    }

    pub(crate) fn update_from_core(
        &mut self,
        py: Python<'_>,
        mut core: CoreAtomGroup,
    ) -> PyResult<()> {
        self.name = core.name.clone();
        self.bonds = core.bonds().to_vec();
        self.secondary_structure = core.secondary_structure();

        // 1. atoms
        let mut new_atoms = IndexMap::with_capacity(core.get_number_of_atoms());
        for (k, atom_core) in core.atoms_mut() {
            if let Some(existing_atom) = self.atoms.shift_remove(k) {
                let mut b = existing_atom.borrow_mut(py);
                b.atomic_number = atom_core.atomic_number();
                b.name = atom_core.name.clone();
                b.label = atom_core.label.clone();
                b.charge = atom_core.charge;
                b.path = atom_core.path.clone();
                b.xyz.borrow_mut(py).inner = atom_core.xyz;
                b.force.borrow_mut(py).inner = atom_core.force;
                drop(b);
                new_atoms.insert(k.clone(), existing_atom);
            } else {
                let atom = PyAtom::from_core(py, atom_core.clone())?;
                new_atoms.insert(k.clone(), Py::new(py, atom)?);
            }
        }
        self.atoms = new_atoms;

        // 2. groups
        let mut new_groups = IndexMap::with_capacity(core.get_number_of_groups());
        let mut to_recurse = Vec::new();
        for (k, child_core) in core.groups_mut() {
            if let Some(existing_child) = self.groups.shift_remove(k) {
                to_recurse.push((existing_child.clone_ref(py), std::mem::take(child_core)));
                new_groups.insert(k.clone(), existing_child);
            } else {
                let child = PyAtomGroup::from_core(py, std::mem::take(child_core))?;
                let child_py = Py::new(py, child)?;
                child_py.borrow_mut(py).my_key = Some(k.clone());
                new_groups.insert(k.clone(), child_py);
            }
        }
        let parent_path = core.path().to_string();
        for (k, removed_child) in self.groups.drain(..) {
            let mut b = removed_child.borrow_mut(py);
            let old_path = format!("{}{}/", parent_path, k);
            b.parent = None;
            b.my_key = None;
            b.root_path = Some(old_path);
        }
        self.groups = new_groups;

        for (child_py, child_core) in to_recurse {
            child_py.borrow_mut(py).update_from_core(py, child_core)?;
        }

        Ok(())
    }

    pub(crate) fn deep_clone(&self, py: Python<'_>) -> PyResult<Self> {
        let current_path = self.compute_path(py);
        let mut atoms = IndexMap::with_capacity(self.atoms.len());
        for (k, atom_py) in &self.atoms {
            let atom_copy = atom_py.borrow(py).__copy__(py)?;
            atoms.insert(k.clone(), Py::new(py, atom_copy)?);
        }

        let mut groups = IndexMap::with_capacity(self.groups.len());
        for (k, child_py) in &self.groups {
            let child_borrow = child_py.borrow(py);
            let child_clone = Self::deep_clone_py(py, &child_borrow, None, Some(k.clone()))?;
            groups.insert(k.clone(), child_clone);
        }

        Ok(Self {
            name: self.name.clone(),
            atoms,
            groups,
            bonds: self.bonds.clone(),
            secondary_structure: self.secondary_structure,
            parent: None,
            my_key: None,
            root_path: Some(current_path),
            sort_atoms: self.sort_atoms.clone(),
            sort_groups: self.sort_groups.clone(),
        })
    }

    pub(crate) fn deep_clone_py(
        py: Python<'_>,
        original: &PyAtomGroup,
        parent: Option<Py<PyAtomGroup>>,
        my_key: Option<String>,
    ) -> PyResult<Py<PyAtomGroup>> {
        let mut atoms = IndexMap::with_capacity(original.atoms.len());
        for (k, atom_py) in &original.atoms {
            let atom_copy = atom_py.borrow(py).__copy__(py)?;
            atoms.insert(k.clone(), Py::new(py, atom_copy)?);
        }

        let new_group = PyAtomGroup {
            name: original.name.clone(),
            atoms,
            groups: IndexMap::new(),
            bonds: original.bonds.clone(),
            secondary_structure: original.secondary_structure,
            parent,
            my_key,
            root_path: original.root_path.clone(),
            sort_atoms: original.sort_atoms.clone(),
            sort_groups: original.sort_groups.clone(),
        };
        let new_group_py = Py::new(py, new_group)?;

        for (k, child_py) in &original.groups {
            let child_borrow = child_py.borrow(py);
            let child_clone = Self::deep_clone_py(
                py,
                &child_borrow,
                Some(new_group_py.clone_ref(py)),
                Some(k.clone()),
            )?;
            new_group_py
                .borrow_mut(py)
                .groups
                .insert(k.clone(), child_clone);
        }

        Ok(new_group_py)
    }

    pub(crate) fn clone_for_select(&self, py: Python<'_>) -> PyResult<Self> {
        let current_path = self.compute_path(py);
        let mut atoms = IndexMap::with_capacity(self.atoms.len());
        for (k, atom_py) in &self.atoms {
            let atom_copy = atom_py.borrow(py).__copy__(py)?;
            atoms.insert(k.clone(), Py::new(py, atom_copy)?);
        }

        let mut groups = IndexMap::with_capacity(self.groups.len());
        for (k, child_py) in &self.groups {
            let child_borrow = child_py.borrow(py);
            let child_clone =
                Self::clone_subtree_for_select(py, &child_borrow, None, Some(k.clone()))?;
            groups.insert(k.clone(), child_clone);
        }

        Ok(Self {
            name: self.name.clone(),
            atoms,
            groups,
            bonds: self.bonds.clone(),
            secondary_structure: self.secondary_structure,
            parent: None,
            my_key: None,
            root_path: Some(current_path),
            sort_atoms: None,
            sort_groups: None,
        })
    }

    pub(crate) fn clone_subtree_for_select(
        py: Python<'_>,
        original: &PyAtomGroup,
        parent: Option<Py<PyAtomGroup>>,
        my_key: Option<String>,
    ) -> PyResult<Py<PyAtomGroup>> {
        let mut atoms = IndexMap::with_capacity(original.atoms.len());
        for (k, atom_py) in &original.atoms {
            let atom_copy = atom_py.borrow(py).__copy__(py)?;
            atoms.insert(k.clone(), Py::new(py, atom_copy)?);
        }

        let new_group = PyAtomGroup {
            name: original.name.clone(),
            atoms,
            groups: IndexMap::new(),
            bonds: original.bonds.clone(),
            secondary_structure: original.secondary_structure,
            parent,
            my_key,
            root_path: None,
            sort_atoms: None,
            sort_groups: None,
        };
        let new_group_py = Py::new(py, new_group)?;

        for (k, child_py) in &original.groups {
            let child_borrow = child_py.borrow(py);
            let child_clone = Self::clone_subtree_for_select(
                py,
                &child_borrow,
                Some(new_group_py.clone_ref(py)),
                Some(k.clone()),
            )?;
            new_group_py
                .borrow_mut(py)
                .groups
                .insert(k.clone(), child_clone);
        }

        Ok(new_group_py)
    }

    fn shift_by_position(&mut self, py: Python<'_>, pos: Position) {
        for a in self.atoms.values() {
            a.borrow(py).xyz.borrow_mut(py).inner += pos;
        }
        for g in self.groups.values() {
            g.borrow_mut(py).shift_by_position(py, pos);
        }
    }

    fn collect_atom_kinds(&self, py: Python<'_>, set: &Bound<'_, PySet>) -> PyResult<()> {
        for g in self.groups.values() {
            g.borrow(py).collect_atom_kinds(py, set)?;
        }
        for a in self.atoms.values() {
            let sym = a.borrow(py).symbol()?;
            set.add(sym)?;
        }
        Ok(())
    }

    fn collect_xyz(&self, py: Python<'_>, output: &mut String) -> PyResult<()> {
        let mut grp_keys: Vec<String> = self.groups.keys().cloned().collect();
        if self.sort_groups.as_deref() == Some("nice") {
            proteindf_bridge::brd::sort_nicely(&mut grp_keys);
        }
        for k in grp_keys {
            if let Some(g) = self.groups.get(&k) {
                g.borrow(py).collect_xyz(py, output)?;
            }
        }

        let mut atm_keys: Vec<String> = self.atoms.keys().cloned().collect();
        if self.sort_atoms.as_deref() == Some("nice") {
            proteindf_bridge::brd::sort_nicely(&mut atm_keys);
        }
        for k in atm_keys {
            if let Some(a) = self.atoms.get(&k) {
                let atom = a.borrow(py);
                let sym = atom.symbol()?;
                let pos = atom.xyz.borrow(py);
                output.push_str(&format!(
                    "{:>2} {:8.3} {:8.3} {:8.3}\n",
                    sym, pos.inner.x, pos.inner.y, pos.inner.z
                ));
            }
        }
        Ok(())
    }
}

fn warn_deprecation(py: Python<'_>, message: &str) -> PyResult<()> {
    let warnings = py.import("warnings")?;
    let category = py.get_type::<PyDeprecationWarning>();
    warnings.call_method1("warn", (message, category, 1))?;
    Ok(())
}

fn ensure_child_link(
    slf: &Bound<'_, PyAtomGroup>,
    key: &str,
    child: &Py<PyAtomGroup>,
) -> PyResult<()> {
    let py = slf.py();
    {
        let b = child.borrow(py);
        if b.parent.is_some() && b.my_key.as_deref() == Some(key) {
            return Ok(());
        }
    }
    let mut b = child.borrow_mut(py);
    if b.parent.is_none() {
        b.parent = Some(slf.clone().unbind());
    }
    b.my_key = Some(key.to_string());
    Ok(())
}

fn assign_charges_rec(
    slf: &Bound<'_, PyAtomGroup>,
    charges: &Bound<'_, PyAny>,
    charge_index: &mut usize,
) -> PyResult<()> {
    let py = slf.py();
    let (child_groups, my_atoms) = {
        let borrow = slf.borrow();
        let mut grp_keys: Vec<String> = borrow.groups.keys().cloned().collect();
        if borrow.sort_groups.as_deref() == Some("nice") {
            proteindf_bridge::brd::sort_nicely(&mut grp_keys);
        }
        let child_groups: Vec<(String, Py<PyAtomGroup>)> = grp_keys
            .into_iter()
            .filter_map(|k| borrow.groups.get(&k).map(|g| (k, g.clone_ref(py))))
            .collect();

        let mut atm_keys: Vec<String> = borrow.atoms.keys().cloned().collect();
        if borrow.sort_atoms.as_deref() == Some("nice") {
            proteindf_bridge::brd::sort_nicely(&mut atm_keys);
        }
        let my_atoms: Vec<Py<PyAtom>> = atm_keys
            .into_iter()
            .filter_map(|k| borrow.atoms.get(&k).map(|a| a.clone_ref(py)))
            .collect();
        (child_groups, my_atoms)
    };

    for (k, g) in &child_groups {
        ensure_child_link(slf, k, g)?;
        assign_charges_rec(g.bind(py), charges, charge_index)?;
    }

    for a in &my_atoms {
        let item = charges.get_item(*charge_index)?;
        let charge_val: f64 = item.extract()?;
        a.borrow_mut(py).charge = charge_val;
        *charge_index += 1;
    }

    Ok(())
}

fn extract_key(key: &Bound<'_, PyAny>) -> PyResult<String> {
    if let Ok(s) = key.extract::<String>() {
        Ok(s)
    } else if let Ok(i) = key.extract::<i64>() {
        Ok(i.to_string())
    } else {
        key.str()?.extract::<String>()
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
    fn __traverse__(&self, visit: pyo3::PyVisit<'_>) -> Result<(), pyo3::PyTraverseError> {
        if let Some(parent) = &self.parent {
            visit.call(parent)?;
        }
        for child in self.groups.values() {
            visit.call(child)?;
        }
        for atom in self.atoms.values() {
            visit.call(atom)?;
        }
        Ok(())
    }

    fn __clear__(&mut self) {
        self.parent = None;
        self.groups.clear();
        self.atoms.clear();
    }

    /// Creates a new `AtomGroup`.
    ///
    /// Args:
    ///     *args: Optional positional argument: a group name (str), another `AtomGroup`
    ///         to clone, or a raw data dictionary. If a dictionary is given, unrecognized
    ///         keys are ignored.
    ///     **kwargs: Optional keyword arguments: `name` (str) or `path` (str).
    ///
    /// Returns:
    ///     AtomGroup: The newly created atom group.
    ///
    /// Raises:
    ///     TypeError: If given more than one positional argument or an unsupported argument type.
    ///     BrValueError: If nested atom data within the dictionary is invalid.
    #[new]
    #[pyo3(signature = (*args, **kwargs))]
    pub fn new(
        py: Python<'_>,
        args: &Bound<'_, pyo3::types::PyTuple>,
        kwargs: Option<&Bound<'_, pyo3::types::PyDict>>,
    ) -> PyResult<Self> {
        let mut group = Self::default();

        if args.len() == 1 {
            let first = args.get_item(0)?;
            if let Ok(other) = first.extract::<PyRef<PyAtomGroup>>() {
                group = other.deep_clone(py)?;
            } else if let Ok(dict) = first.cast::<pyo3::types::PyDict>() {
                let val = crate::value::py_to_value(py, dict.as_any())?;
                let mut core = CoreAtomGroup::new();
                proteindf_bridge::brd::atomgroup_set_by_dict_data(&mut core, &val)
                    .map_err(to_py_err)?;
                group = PyAtomGroup::from_core(py, core)?;
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
            if let Some(path) = kw.get_item("path")? {
                let mut p = path.extract::<String>()?;
                if p.is_empty() || !p.ends_with('/') {
                    p.push('/');
                }
                group.root_path = Some(p);
            }
        }

        Ok(group)
    }

    /// Name of this atom group.
    ///
    /// Returns:
    ///     str: The group name.
    #[getter]
    pub fn name(&self) -> String {
        self.name.clone()
    }

    /// Sets the name of this atom group.
    ///
    /// Args:
    ///     val (str): The new group name.
    #[setter]
    pub fn set_name(&mut self, val: String) {
        self.name = val;
    }

    /// Parent `AtomGroup` in the hierarchy, or None if this is a root group.
    ///
    /// Returns:
    ///     Optional[AtomGroup]: The parent group reference, or None.
    #[getter]
    pub fn parent(&self, py: Python<'_>) -> Option<Py<PyAtomGroup>> {
        self.get_parent(py)
    }

    /// Hierarchical path of this group from the root (e.g. "/model_1/A/").
    ///
    /// Returns:
    ///     str: The hierarchical path.
    #[getter]
    pub fn path(&self) -> String {
        Python::attach(|py| self.compute_path(py))
    }

    /// Sets the root path of this group.
    ///
    /// Args:
    ///     val (str): The root path string.
    ///
    /// Raises:
    ///     ValueError: If this group is a child group with a parent.
    #[setter]
    pub fn set_path(&mut self, py: Python<'_>, val: String) -> PyResult<()> {
        if self.get_parent(py).is_some() {
            return Err(PyValueError::new_err(
                "Cannot set path on a child group with a parent",
            ));
        }
        let mut p = val;
        if p.is_empty() || !p.ends_with('/') {
            p.push('/');
        }
        self.root_path = Some(p);
        Ok(())
    }

    /// Returns the depth of this group in the hierarchy.
    ///
    /// Returns:
    ///     int: The path depth (0 for root, 1 for model, 2 for chain, 3 for residue).
    pub fn path_depth(&self) -> usize {
        Python::attach(|py| self.compute_depth(py))
    }

    /// Checks whether this group is at the model level (path depth 1).
    ///
    /// Returns:
    ///     bool: True if at model level, False otherwise.
    pub fn is_model_level(&self) -> bool {
        self.path_depth() == 1
    }

    /// Checks whether this group is at the chain level (path depth 2).
    ///
    /// Returns:
    ///     bool: True if at chain level, False otherwise.
    pub fn is_chain_level(&self) -> bool {
        self.path_depth() == 2
    }

    /// Checks whether this group is at the residue level (path depth 3).
    ///
    /// Returns:
    ///     bool: True if at residue level, False otherwise.
    pub fn is_residue_level(&self) -> bool {
        self.path_depth() == 3
    }

    /// Validates this group against the standard hierarchy schema.
    ///
    /// Returns:
    ///     List[SchemaViolation]: A list of detected schema violations.
    pub fn validate_schema(&self) -> Vec<PySchemaViolation> {
        Python::attach(|py| {
            self.to_core(py)
                .validate_schema()
                .into_iter()
                .map(PySchemaViolation::from)
                .collect()
        })
    }

    /// Secondary structure code of this residue group ('H', 'E', '-', or None).
    ///
    /// Returns:
    ///     Optional[str]: 'H' for helix, 'E' for strand, '-' for loop, or None.
    #[getter]
    pub fn secondary_structure(&self) -> Option<String> {
        self.secondary_structure.map(|ss| ss.to_string())
    }

    /// Sets the secondary structure code of this residue group.
    ///
    /// Args:
    ///     val (Optional[str]): 'H', 'E', '-', or None.
    ///
    /// Raises:
    ///     ValueError: If `val` is not one of 'H', 'E', '-', or None.
    #[setter]
    pub fn set_secondary_structure(&mut self, val: Option<&str>) -> PyResult<()> {
        match val {
            None => {
                self.secondary_structure = None;
                Ok(())
            }
            Some("H") => {
                self.secondary_structure = Some(SsCode::Helix);
                Ok(())
            }
            Some("E") => {
                self.secondary_structure = Some(SsCode::Strand);
                Ok(())
            }
            Some("-") => {
                self.secondary_structure = Some(SsCode::Loop);
                Ok(())
            }
            Some(other) => Err(PyValueError::new_err(format!(
                "Invalid secondary structure code '{other}', expected 'H', 'E', '-', or None"
            ))),
        }
    }

    /// Returns the number of direct atoms contained in this group.
    ///
    /// Returns:
    ///     int: Number of direct atoms.
    pub fn get_number_of_atoms(&self) -> usize {
        self.atoms.len()
    }

    /// Returns the total number of atoms in this group and all its descendant groups.
    ///
    /// Returns:
    ///     int: Total number of atoms across the entire sub-tree.
    pub fn get_number_of_all_atoms(&self) -> usize {
        Python::attach(|py| self.get_number_of_all_atoms_py(py))
    }

    /// Returns the number of direct child groups contained in this group.
    ///
    /// Returns:
    ///     int: Number of direct child groups.
    pub fn get_number_of_groups(&self) -> usize {
        self.groups.len()
    }

    /// Returns the number of bonds stored directly in this group.
    ///
    /// Returns:
    ///     int: Number of direct bonds.
    pub fn get_number_of_bonds(&self) -> usize {
        self.bonds.len()
    }

    /// Retrieves an atom matching the given key or name from this group.
    ///
    /// Args:
    ///     key (Union[str, int]): The atom key or name to search for.
    ///
    /// Returns:
    ///     Optional[Atom]: The matching atom reference, or None if not found.
    pub fn get_atom(slf: &Bound<'_, Self>, key: &Bound<'_, PyAny>) -> PyResult<Option<Py<PyAtom>>> {
        let k = extract_key(key)?;
        let py = slf.py();
        let borrow = slf.borrow();
        let atom_opt = if let Some(a) = borrow.atoms.get(&k) {
            Some((k.clone(), a.clone_ref(py)))
        } else {
            borrow
                .atoms
                .iter()
                .find(|(_, a)| a.borrow(py).name.trim() == k.trim())
                .map(|(key_name, a)| (key_name.clone(), a.clone_ref(py)))
        };
        if let Some((atom_key, atom_py)) = atom_opt {
            let group_path = borrow.compute_path(py);
            atom_py.borrow_mut(py).path = format!("{}{}", group_path, atom_key);
            Ok(Some(atom_py))
        } else {
            Ok(None)
        }
    }

    /// Inserts an atom into this group or at a specified hierarchical path.
    ///
    /// If `key` contains path separators (e.g. ``"subgroup/C1"``), intermediate
    /// child groups are created as needed and the atom is inserted at the leaf.
    /// The inserted atom is deep-copied into the tree.
    ///
    /// Args:
    ///     key (Union[str, int]): The atom key or relative/absolute path.
    ///     atom (Atom): The atom to insert.
    ///
    /// Raises:
    ///     TypeError: If `atom` is not an `Atom` instance.
    pub fn set_atom(
        slf: &Bound<'_, Self>,
        key: &Bound<'_, PyAny>,
        atom: &Bound<'_, PyAtom>,
    ) -> PyResult<()> {
        let py = slf.py();
        let key_str = extract_key(key)?;
        let mut cur_key = key_str.as_str();
        while cur_key.starts_with('/') {
            cur_key = &cur_key[1..];
        }

        if let Some((grp_key, rest)) = cur_key.split_once('/') {
            let has_grp = slf.borrow().groups.contains_key(grp_key);
            if !has_grp {
                let new_group = Py::new(py, PyAtomGroup::default())?;
                let grp_key_py = pyo3::types::PyString::new(py, grp_key);
                Self::set_group(slf, grp_key_py.as_any(), new_group.bind(py))?;
            }
            let child_group = {
                let borrow = slf.borrow();
                borrow.groups.get(grp_key).unwrap().clone_ref(py)
            };
            let rest_py = pyo3::types::PyString::new(py, rest);
            Self::set_atom(child_group.bind(py), rest_py.as_any(), atom)?;
        } else {
            let mut atom_copy = atom.borrow().__copy__(py)?;
            let group_path = slf.borrow().compute_path(py);
            atom_copy.path = format!("{}{}", group_path, cur_key);
            let atom_py = Py::new(py, atom_copy)?;
            slf.borrow_mut().atoms.insert(cur_key.to_string(), atom_py);
        }
        Ok(())
    }

    /// Checks whether an atom matching the given key or name exists in this group.
    ///
    /// Args:
    ///     key (Union[str, int]): The atom key or name to search for.
    ///
    /// Returns:
    ///     bool: True if a matching atom is found, False otherwise.
    pub fn has_atom(&self, py: Python<'_>, key: &Bound<'_, PyAny>) -> PyResult<bool> {
        let k = extract_key(key)?;
        Ok(self.atoms.contains_key(&k)
            || self
                .atoms
                .values()
                .any(|a| a.borrow(py).name.trim() == k.trim()))
    }

    /// Checks whether an atom with the given key exists directly in this group.
    ///
    /// Args:
    ///     key (Union[str, int]): The atom key.
    ///
    /// Returns:
    ///     bool: True if an atom with the key exists, False otherwise.
    pub fn has_atomkey(&self, key: &Bound<'_, PyAny>) -> PyResult<bool> {
        let k = extract_key(key)?;
        Ok(self.atoms.contains_key(&k))
    }

    /// Checks whether an atom with the given name exists directly in this group.
    ///
    /// Args:
    ///     name (Union[str, int]): The atom name.
    ///
    /// Returns:
    ///     bool: True if an atom with the name exists, False otherwise.
    pub fn has_atomname(&self, py: Python<'_>, name: &Bound<'_, PyAny>) -> PyResult<bool> {
        let name_str = extract_key(name)?;
        let target = name_str.trim();
        for a in self.atoms.values() {
            if a.borrow(py).name.trim() == target {
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// Removes an atom with the given key from this group.
    ///
    /// Args:
    ///     key (Union[str, int]): The atom key to remove.
    ///
    /// Returns:
    ///     Optional[Atom]: The removed atom, or None if not found.
    pub fn del_atom(
        &mut self,
        _py: Python<'_>,
        key: &Bound<'_, PyAny>,
    ) -> PyResult<Option<Py<PyAtom>>> {
        let k = extract_key(key)?;
        Ok(self.atoms.shift_remove(&k))
    }

    /// Removes an atom with the given key from this group.
    ///
    /// Args:
    ///     key (Union[str, int]): The atom key to remove.
    ///
    /// Returns:
    ///     Optional[Atom]: The removed atom, or None if not found.
    pub fn remove_atom(
        &mut self,
        py: Python<'_>,
        key: &Bound<'_, PyAny>,
    ) -> PyResult<Option<Py<PyAtom>>> {
        self.del_atom(py, key)
    }

    /// Removes an atom with the given key from this group.
    ///
    /// .. deprecated::
    ///     Use :meth:`remove_atom` instead.
    ///
    /// Args:
    ///     key (Union[str, int]): The atom key to remove.
    pub fn erase_atom(&mut self, py: Python<'_>, key: &Bound<'_, PyAny>) -> PyResult<()> {
        warn_deprecation(py, "erase_atom() is deprecated, use remove_atom() instead.")?;
        self.remove_atom(py, key)?;
        Ok(())
    }

    /// Retrieves a child group matching the given key or name from this group.
    ///
    /// Args:
    ///     key (Union[str, int]): The child group key or name.
    ///
    /// Returns:
    ///     Optional[AtomGroup]: The child group reference, or None if not found.
    pub fn get_group(
        slf: &Bound<'_, Self>,
        key: &Bound<'_, PyAny>,
    ) -> PyResult<Option<Py<PyAtomGroup>>> {
        let k = extract_key(key)?;
        let py = slf.py();
        let borrow = slf.borrow();
        let child_opt = if let Some(g) = borrow.groups.get(&k) {
            Some((k.clone(), g.clone_ref(py)))
        } else {
            borrow
                .groups
                .iter()
                .find(|(_, g)| g.borrow(py).name == k)
                .map(|(key_name, g)| (key_name.clone(), g.clone_ref(py)))
        };
        drop(borrow);

        if let Some((child_key, child)) = child_opt {
            ensure_child_link(slf, &child_key, &child)?;
            Ok(Some(child))
        } else {
            Ok(None)
        }
    }

    /// Inserts a child group into this group with the specified key.
    ///
    /// The inserted group is deep-copied into the tree.
    ///
    /// Args:
    ///     key (Union[str, int]): The child group key.
    ///     group (AtomGroup): The child group to insert.
    ///
    /// Raises:
    ///     TypeError: If `group` is not an `AtomGroup` instance.
    pub fn set_group(
        slf: &Bound<'_, Self>,
        key: &Bound<'_, PyAny>,
        group: &Bound<'_, PyAtomGroup>,
    ) -> PyResult<()> {
        let k = extract_key(key)?;
        let py = slf.py();
        let parent_py = slf.clone().unbind();

        let child_borrow = group.borrow();
        let copy_py = Self::deep_clone_py(py, &child_borrow, Some(parent_py), Some(k.clone()))?;
        drop(child_borrow);

        let parent_path = slf.borrow().compute_path(py);
        let old = slf.borrow_mut().groups.insert(k.clone(), copy_py);
        if let Some(old_group) = old {
            let mut b = old_group.borrow_mut(py);
            let old_path = format!("{}{}/", parent_path, k);
            b.parent = None;
            b.my_key = None;
            b.root_path = Some(old_path);
        }
        Ok(())
    }

    /// Checks whether a child group matching the given key or name exists in this group.
    ///
    /// Args:
    ///     key (Union[str, int]): The child group key or name.
    ///
    /// Returns:
    ///     bool: True if a matching child group is found, False otherwise.
    pub fn has_group(&self, py: Python<'_>, key: &Bound<'_, PyAny>) -> PyResult<bool> {
        let k = extract_key(key)?;
        Ok(self.groups.contains_key(&k) || self.groups.values().any(|g| g.borrow(py).name == k))
    }

    /// Checks whether a child group with the given key exists directly in this group.
    ///
    /// Args:
    ///     key (Union[str, int]): The child group key.
    ///
    /// Returns:
    ///     bool: True if a child group with the key exists, False otherwise.
    pub fn has_groupkey(&self, key: &Bound<'_, PyAny>) -> PyResult<bool> {
        let k = extract_key(key)?;
        Ok(self.groups.contains_key(&k))
    }

    /// Checks whether a child group with the given name exists directly in this group.
    ///
    /// Args:
    ///     name (Union[str, int]): The child group name.
    ///
    /// Returns:
    ///     bool: True if a child group with the name exists, False otherwise.
    pub fn has_groupname(&self, py: Python<'_>, name: &Bound<'_, PyAny>) -> PyResult<bool> {
        let target = extract_key(name)?;
        for g in self.groups.values() {
            if g.borrow(py).name == target {
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// Removes a child group with the given key from this group.
    ///
    /// Args:
    ///     key (Union[str, int]): The child group key to remove.
    ///
    /// Returns:
    ///     Optional[AtomGroup]: The removed child group, or None if not found.
    pub fn del_group(
        slf: &Bound<'_, Self>,
        key: &Bound<'_, PyAny>,
    ) -> PyResult<Option<Py<PyAtomGroup>>> {
        let k = extract_key(key)?;
        let py = slf.py();
        let parent_path = slf.borrow().compute_path(py);
        let old = slf.borrow_mut().groups.shift_remove(&k);
        if let Some(old_group) = old {
            let mut b = old_group.borrow_mut(py);
            let old_path = format!("{}{}/", parent_path, k);
            b.parent = None;
            b.my_key = None;
            b.root_path = Some(old_path);
            drop(b);
            Ok(Some(old_group))
        } else {
            Ok(None)
        }
    }

    /// Removes a child group with the given key from this group.
    ///
    /// Args:
    ///     key (Union[str, int]): The child group key to remove.
    ///
    /// Returns:
    ///     Optional[AtomGroup]: The removed child group, or None if not found.
    pub fn remove_group(
        slf: &Bound<'_, Self>,
        key: &Bound<'_, PyAny>,
    ) -> PyResult<Option<Py<PyAtomGroup>>> {
        Self::del_group(slf, key)
    }

    /// Removes a child group with the given key from this group.
    ///
    /// .. deprecated::
    ///     Use :meth:`remove_group` instead.
    ///
    /// Args:
    ///     key (Union[str, int]): The child group key to remove.
    pub fn erase_group(slf: &Bound<'_, Self>, key: &Bound<'_, PyAny>) -> PyResult<()> {
        let py = slf.py();
        warn_deprecation(
            py,
            "erase_group() is deprecated, use remove_group() instead.",
        )?;
        Self::remove_group(slf, key)?;
        Ok(())
    }

    /// Sorting mode for direct atoms.
    ///
    /// Setting to ``"nice"`` enables natural (alphanumeric) sorting of atom keys
    /// (e.g. "C1", "C2", "C10" instead of alphabetical "C1", "C10", "C2") when
    /// accessing atoms. Setting to None maintains insertion order.
    ///
    /// Returns:
    ///     Optional[str]: The sorting mode (``"nice"`` or None).
    #[getter]
    pub fn sort_atoms(&self) -> Option<String> {
        self.sort_atoms.clone()
    }

    /// Sets the sorting mode for direct atoms.
    ///
    /// Args:
    ///     val (Optional[str]): The sorting mode (``"nice"`` or None).
    #[setter]
    pub fn set_sort_atoms(&mut self, val: Option<String>) {
        self.sort_atoms = val;
    }

    /// Sorting mode for direct child groups.
    ///
    /// Setting to ``"nice"`` enables natural (alphanumeric) sorting of group keys
    /// (e.g. "res1", "res2", "res10" instead of alphabetical "res1", "res10", "res2")
    /// when accessing child groups. Setting to None maintains insertion order.
    ///
    /// Returns:
    ///     Optional[str]: The sorting mode (``"nice"`` or None).
    #[getter]
    pub fn sort_groups(&self) -> Option<String> {
        self.sort_groups.clone()
    }

    /// Sets the sorting mode for direct child groups.
    ///
    /// Args:
    ///     val (Optional[str]): The sorting mode (``"nice"`` or None).
    #[setter]
    pub fn set_sort_groups(&mut self, val: Option<String>) {
        self.sort_groups = val;
    }

    /// Returns a list of (key, atom) pairs for all atoms directly in this group.
    ///
    /// If :attr:`sort_atoms` is set to ``"nice"``, the items are sorted naturally.
    ///
    /// Returns:
    ///     List[Tuple[str, Atom]]: List of (key, atom) pairs.
    pub fn atoms(slf: &Bound<'_, Self>) -> PyResult<Vec<(String, Py<PyAtom>)>> {
        let py = slf.py();
        let borrow = slf.borrow();
        let group_path = borrow.compute_path(py);
        let mut keys: Vec<String> = borrow.atoms.keys().cloned().collect();
        if borrow.sort_atoms.as_deref() == Some("nice") {
            proteindf_bridge::brd::sort_nicely(&mut keys);
        }
        let mut result = Vec::with_capacity(keys.len());
        for k in keys {
            if let Some(a) = borrow.atoms.get(&k) {
                a.borrow_mut(py).path = format!("{}{}", group_path, k);
                result.push((k, a.clone_ref(py)));
            }
        }
        Ok(result)
    }

    /// Returns a list of (key, group) pairs for all child groups directly in this group.
    ///
    /// If :attr:`sort_groups` is set to ``"nice"``, the items are sorted naturally.
    ///
    /// Returns:
    ///     List[Tuple[str, AtomGroup]]: List of (key, group) pairs.
    pub fn groups(slf: &Bound<'_, Self>) -> PyResult<Vec<(String, Py<PyAtomGroup>)>> {
        let py = slf.py();
        let borrow = slf.borrow();
        let mut keys: Vec<String> = borrow.groups.keys().cloned().collect();
        if borrow.sort_groups.as_deref() == Some("nice") {
            proteindf_bridge::brd::sort_nicely(&mut keys);
        }
        let pairs: Vec<(String, Py<PyAtomGroup>)> = keys
            .into_iter()
            .filter_map(|k| borrow.groups.get(&k).map(|g| (k, g.clone_ref(py))))
            .collect();
        drop(borrow);
        for (k, g) in &pairs {
            ensure_child_link(slf, k, g)?;
        }
        Ok(pairs)
    }

    /// Returns a list of keys of all direct child groups.
    ///
    /// If :attr:`sort_groups` is set to ``"nice"``, the keys are sorted naturally.
    ///
    /// Returns:
    ///     List[str]: List of child group keys.
    pub fn get_group_list(&self) -> Vec<String> {
        let mut keys: Vec<String> = self.groups.keys().cloned().collect();
        if self.sort_groups.as_deref() == Some("nice") {
            proteindf_bridge::brd::sort_nicely(&mut keys);
        }
        keys
    }

    /// Returns a list of keys of all direct atoms.
    ///
    /// If :attr:`sort_atoms` is set to ``"nice"``, the keys are sorted naturally.
    ///
    /// Returns:
    ///     List[str]: List of atom keys.
    pub fn get_atom_keys(&self) -> Vec<String> {
        let mut keys: Vec<String> = self.atoms.keys().cloned().collect();
        if self.sort_atoms.as_deref() == Some("nice") {
            proteindf_bridge::brd::sort_nicely(&mut keys);
        }
        keys
    }

    /// Recursively collects and returns all atoms in this group and all descendant groups.
    ///
    /// Returns:
    ///     List[Atom]: Flat list of all atoms in the hierarchy.
    pub fn get_atom_list(slf: &Bound<'_, Self>) -> PyResult<Vec<Py<PyAtom>>> {
        let py = slf.py();
        let mut list = Vec::new();
        let (group_path, child_groups, my_atoms) = {
            let borrow = slf.borrow();
            let group_path = borrow.compute_path(py);

            let mut grp_keys: Vec<String> = borrow.groups.keys().cloned().collect();
            if borrow.sort_groups.as_deref() == Some("nice") {
                proteindf_bridge::brd::sort_nicely(&mut grp_keys);
            }
            let child_groups: Vec<(String, Py<PyAtomGroup>)> = grp_keys
                .into_iter()
                .filter_map(|k| borrow.groups.get(&k).map(|g| (k, g.clone_ref(py))))
                .collect();

            let mut atm_keys: Vec<String> = borrow.atoms.keys().cloned().collect();
            if borrow.sort_atoms.as_deref() == Some("nice") {
                proteindf_bridge::brd::sort_nicely(&mut atm_keys);
            }
            let my_atoms: Vec<(String, Py<PyAtom>)> = atm_keys
                .into_iter()
                .filter_map(|k| borrow.atoms.get(&k).map(|a| (k, a.clone_ref(py))))
                .collect();

            (group_path, child_groups, my_atoms)
        };

        for (k, g) in &child_groups {
            ensure_child_link(slf, k, g)?;
            let child_atoms = Self::get_atom_list(g.bind(py))?;
            list.extend(child_atoms);
        }
        for (k, a) in &my_atoms {
            a.borrow_mut(py).path = format!("{}{}", group_path, k);
            list.push(a.clone_ref(py));
        }
        Ok(list)
    }

    /// Recursively collects and returns the paths of all atoms in this group and descendant groups.
    ///
    /// Returns:
    ///     List[str]: List of full paths for all atoms in the hierarchy.
    pub fn get_path_list(slf: &Bound<'_, Self>) -> PyResult<Vec<String>> {
        let py = slf.py();
        let atoms = Self::get_atom_list(slf)?;
        Ok(atoms.iter().map(|a| a.borrow(py).path.clone()).collect())
    }

    /// Recursively searches and collects atoms whose key or name matches the target string.
    ///
    /// Args:
    ///     key_or_name (Union[str, int]): The atom key or name to match.
    ///
    /// Returns:
    ///     List[Atom]: List of matching atoms found across the hierarchy.
    pub fn pickup_atoms(
        slf: &Bound<'_, Self>,
        key_or_name: &Bound<'_, PyAny>,
    ) -> PyResult<Vec<Py<PyAtom>>> {
        let py = slf.py();
        let target = extract_key(key_or_name)?;
        let mut result = Vec::new();

        let (group_path, child_groups, my_atoms) = {
            let borrow = slf.borrow();
            let group_path = borrow.compute_path(py);

            let mut grp_keys: Vec<String> = borrow.groups.keys().cloned().collect();
            if borrow.sort_groups.as_deref() == Some("nice") {
                proteindf_bridge::brd::sort_nicely(&mut grp_keys);
            }
            let child_groups: Vec<(String, Py<PyAtomGroup>)> = grp_keys
                .into_iter()
                .filter_map(|k| borrow.groups.get(&k).map(|g| (k, g.clone_ref(py))))
                .collect();

            let mut atm_keys: Vec<String> = borrow.atoms.keys().cloned().collect();
            if borrow.sort_atoms.as_deref() == Some("nice") {
                proteindf_bridge::brd::sort_nicely(&mut atm_keys);
            }
            let my_atoms: Vec<(String, Py<PyAtom>)> = atm_keys
                .into_iter()
                .filter_map(|k| borrow.atoms.get(&k).map(|a| (k, a.clone_ref(py))))
                .collect();
            (group_path, child_groups, my_atoms)
        };

        for (k, g) in &child_groups {
            ensure_child_link(slf, k, g)?;
            let child_atoms = Self::pickup_atoms(g.bind(py), key_or_name)?;
            result.extend(child_atoms);
        }

        for (k, a) in &my_atoms {
            if k == &target || a.borrow(py).name == target {
                a.borrow_mut(py).path = format!("{}{}", group_path, k);
                result.push(a.clone_ref(py));
            }
        }

        Ok(result)
    }

    /// Sequentially assigns partial charges from a sequence to all atoms in the hierarchy.
    ///
    /// The number of values in `charges` must exactly match the total atom count.
    ///
    /// Args:
    ///     charges (Sequence[float]): An iterable of numerical charge values.
    ///
    /// Raises:
    ///     TypeError: If `charges` is not a valid sequence.
    ///     ValueError: If the length of `charges` does not match the total atom count.
    pub fn assign_charges(slf: &Bound<'_, Self>, charges: &Bound<'_, PyAny>) -> PyResult<()> {
        let total_atoms = slf.borrow().get_number_of_all_atoms();
        let charges_len = charges.len().map_err(|_| {
            let type_name = charges
                .get_type()
                .name()
                .map(|n| n.to_string())
                .unwrap_or_else(|_| "unknown".to_string());
            PyTypeError::new_err(format!("Expected sequence for charges, got {}", type_name))
        })?;
        if charges_len != total_atoms {
            return Err(PyValueError::new_err(format!(
                "charges length ({}) does not match atom count ({})",
                charges_len, total_atoms
            )));
        }

        let mut index = 0;
        assign_charges_rec(slf, charges, &mut index)?;
        if index != charges_len {
            return Err(PyValueError::new_err(format!(
                "charges length ({}) does not match atom count ({})",
                charges_len, index
            )));
        }
        Ok(())
    }

    /// Navigates the hierarchy tree and retrieves the group matching an absolute path.
    ///
    /// Searches descendants first; if not found within the subtree, traverses up
    /// through ancestor groups to search the tree. Does not support relative paths
    /// (such as ``"."`` or ``".."``).
    ///
    /// Args:
    ///     query_path (Union[str, int]): The absolute path to navigate to (e.g. ``"/protein/chainA"``).
    ///
    /// Returns:
    ///     Optional[AtomGroup]: The matched group, or None if not found.
    pub fn get_family(
        slf: &Bound<'_, Self>,
        query_path: &Bound<'_, PyAny>,
    ) -> PyResult<Option<Py<PyAtomGroup>>> {
        let mut query = extract_key(query_path)?;
        if !query.ends_with('/') {
            query.push('/');
        }

        let py = slf.py();
        let (self_path, parent_opt) = {
            let b = slf.borrow();
            (b.compute_path(py), b.get_parent(py))
        };

        if query == self_path {
            return Ok(Some(slf.clone().unbind()));
        }

        if query.starts_with(&self_path) {
            let mut curr = slf.clone().unbind();
            let mut curr_path = self_path;

            while curr_path != query {
                let curr_bound = curr.bind(py);
                let pairs: Vec<(String, Py<PyAtomGroup>)> = curr_bound
                    .borrow()
                    .groups
                    .iter()
                    .map(|(k, g)| (k.clone(), g.clone_ref(py)))
                    .collect();

                let mut matched = None;
                for (k, g) in pairs {
                    ensure_child_link(curr_bound, &k, &g)?;
                    let child_path = g.borrow(py).compute_path(py);
                    if query.starts_with(&child_path) {
                        matched = Some((g, child_path));
                        break;
                    }
                }

                if let Some((child_grp, child_path)) = matched {
                    curr = child_grp;
                    curr_path = child_path;
                } else {
                    return Ok(None);
                }
            }
            return Ok(Some(curr));
        }

        if let Some(parent_py) = parent_opt {
            return Self::get_family(parent_py.bind(py), query_path);
        }

        Ok(None)
    }

    /// Splits a hierarchical path string into individual component strings.
    ///
    /// Args:
    ///     path (Union[str, int]): The path string to split.
    ///
    /// Returns:
    ///     List[str]: List of path segment strings.
    #[staticmethod]
    pub fn divide_path(path: &Bound<'_, PyAny>) -> PyResult<Vec<String>> {
        let s = extract_key(path)?;
        let parts: Vec<String> = s
            .split('/')
            .filter(|p| !p.is_empty())
            .map(|p| p.to_string())
            .collect();
        Ok(parts)
    }

    /// Total electric charge of this group and all its descendant groups.
    ///
    /// Returns:
    ///     float: Sum of atomic partial charges.
    #[getter]
    pub fn charge(&self, py: Python<'_>) -> f64 {
        let mut sum = 0.0;
        for g in self.groups.values() {
            sum += g.borrow(py).charge(py);
        }
        for a in self.atoms.values() {
            sum += a.borrow(py).charge;
        }
        sum
    }

    /// Sum of atomic numbers of all atoms in this group and all descendant groups.
    ///
    /// Returns:
    ///     float: Sum of atomic numbers.
    #[getter]
    pub fn real_nuclei_charge(&self, py: Python<'_>) -> f64 {
        let mut sum = 0.0;
        for g in self.groups.values() {
            sum += g.borrow(py).real_nuclei_charge(py);
        }
        for a in self.atoms.values() {
            let atom = a.borrow(py);
            if atom.atomic_number > 0 {
                sum += atom.atomic_number as f64;
            }
        }
        sum
    }

    /// Nuclear charge of all atoms in this group and all descendant groups.
    ///
    /// For atoms where atomic number is zero, their partial charge is used instead.
    ///
    /// Returns:
    ///     float: Total nuclear charge.
    #[getter]
    pub fn nuclei_charge(&self, py: Python<'_>) -> f64 {
        let mut sum = 0.0;
        for g in self.groups.values() {
            sum += g.borrow(py).nuclei_charge(py);
        }
        for a in self.atoms.values() {
            let atom = a.borrow(py);
            if atom.atomic_number > 0 {
                sum += atom.atomic_number as f64;
            } else {
                sum += atom.charge;
            }
        }
        sum
    }

    /// Returns a set of unique chemical element symbols present in the hierarchy.
    ///
    /// Returns:
    ///     Set[str]: Set of element symbols (e.g. `{"C", "H", "N", "O"}`).
    pub fn get_atom_kinds<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PySet>> {
        let set = PySet::empty(py)?;
        self.collect_atom_kinds(py, &set)?;
        Ok(set)
    }

    /// Generates an XYZ-format coordinate block string for all atoms in the hierarchy.
    ///
    /// Returns:
    ///     str: Multi-line XYZ coordinate string.
    pub fn get_xyz(&self, py: Python<'_>) -> PyResult<String> {
        let count = self.get_number_of_all_atoms();
        let mut output = format!("{}\n# \n", count);
        self.collect_xyz(py, &mut output)?;
        Ok(output)
    }

    /// Returns a dictionary mapping element symbols to their occurrence counts.
    ///
    /// Returns:
    ///     Dict[str, int]: Counts of each chemical element.
    pub fn get_atom_kinds_count(&self) -> HashMap<String, usize> {
        Python::attach(|py| {
            let mut counts = HashMap::new();
            self.for_each_atom_recursive(py, &mut |a| {
                if let Ok(sym) = proteindf_bridge::periodic_table::PeriodicTable::get_symbol(
                    a.borrow().atomic_number,
                ) {
                    *counts.entry(sym.to_string()).or_insert(0) += 1;
                }
            });
            counts
        })
    }

    /// Computes the chemical formula formatted using the Hill system.
    ///
    /// If carbon is present, carbon is listed first, followed by hydrogen, and then
    /// all other elements in alphabetical order. Otherwise, all elements are listed
    /// in alphabetical order.
    ///
    /// Returns:
    ///     str: The Hill-system chemical formula (e.g. "C6H12O6").
    pub fn get_formula(&self) -> String {
        let kinds = self.get_atom_kinds_count();
        proteindf_bridge::atom_group::formula_from_kinds(&kinds)
    }

    /// Returns the chemical formula formatted using the Hill system.
    ///
    /// Alias for :meth:`get_formula`.
    ///
    /// Returns:
    ///     str: The Hill-system chemical formula.
    pub fn formula(&self) -> String {
        self.get_formula()
    }

    /// Total molecular weight of all atoms in this group and descendant groups.
    ///
    /// Returns:
    ///     float: Total molecular weight in atomic mass units (Da).
    #[getter]
    pub fn weight(&self) -> f64 {
        Python::attach(|py| {
            let mut sum = 0.0;
            self.for_each_atom_recursive(py, &mut |a| {
                sum += proteindf_bridge::periodic_table::PeriodicTable::atomic_weight(
                    a.borrow().atomic_number,
                )
                .unwrap_or(0.0);
            });
            sum
        })
    }

    /// Computes the geometric center (centroid) of all atoms in the hierarchy.
    ///
    /// Returns:
    ///     Position: The centroid coordinates.
    pub fn center(&self) -> PyPosition {
        Python::attach(|py| PyPosition::from_core(self.center_py(py)))
    }

    /// Computes the axis-aligned bounding box enclosing all atoms in the hierarchy.
    ///
    /// Returns:
    ///     Tuple[Position, Position]: A tuple of `(min_position, max_position)`.
    pub fn r#box(&self) -> (PyPosition, PyPosition) {
        Python::attach(|py| {
            let (bmin, bmax) = self.box_py(py);
            (PyPosition::from_core(bmin), PyPosition::from_core(bmax))
        })
    }

    /// Computes the axis-aligned bounding box enclosing all atoms in the hierarchy.
    ///
    /// Alias for :meth:`box`.
    ///
    /// Returns:
    ///     Tuple[Position, Position]: A tuple of `(min_position, max_position)`.
    pub fn get_box(&self) -> (PyPosition, PyPosition) {
        self.r#box()
    }

    /// Merges another `AtomGroup` hierarchy into this group in-place.
    ///
    /// Subgroups and atoms are merged recursively. Existing bonds and secondary
    /// structure information from `other` are also incorporated.
    ///
    /// Args:
    ///     other (AtomGroup): The atom group to merge into this group.
    pub fn merge(slf: &Bound<'_, Self>, other: &Bound<'_, PyAtomGroup>) -> PyResult<()> {
        let py = slf.py();
        let other_ref = other.borrow();

        // 1. groups
        for (k, child_py) in &other_ref.groups {
            let exists_key = slf.borrow().groups.contains_key(k);
            let matching_name_key = if !exists_key {
                let child_name = child_py.borrow(py).name.clone();
                if !child_name.is_empty() {
                    slf.borrow()
                        .groups
                        .iter()
                        .find(|(_, g)| g.borrow(py).name == child_name)
                        .map(|(existing_k, _)| existing_k.clone())
                } else {
                    None
                }
            } else {
                None
            };

            if exists_key {
                let existing_child = slf.borrow().groups.get(k).unwrap().clone_ref(py);
                Self::merge(existing_child.bind(py), child_py.bind(py))?;
            } else if let Some(matched_k) = matching_name_key {
                let existing_child = slf.borrow().groups.get(&matched_k).unwrap().clone_ref(py);
                Self::merge(existing_child.bind(py), child_py.bind(py))?;
            } else {
                let parent_py = slf.clone().unbind();

                let child_borrow = child_py.borrow(py);
                let copy_py =
                    Self::deep_clone_py(py, &child_borrow, Some(parent_py), Some(k.clone()))?;
                slf.borrow_mut().groups.insert(k.clone(), copy_py);
            }
        }

        // 2. atoms
        for (k, atom_py) in &other_ref.atoms {
            let mut atom_copy = atom_py.borrow(py).__copy__(py)?;
            let group_path = slf.borrow().compute_path(py);
            atom_copy.path = format!("{}{}", group_path, k);
            let new_atom_py = Py::new(py, atom_copy)?;
            slf.borrow_mut().atoms.insert(k.clone(), new_atom_py);
        }

        // 3. bonds
        for bond in &other_ref.bonds {
            if !slf.borrow().bonds.contains(bond) {
                slf.borrow_mut().bonds.push(bond.clone());
            }
        }

        // 4. secondary structure
        if other_ref.secondary_structure.is_some() {
            slf.borrow_mut().secondary_structure = other_ref.secondary_structure;
        }

        Ok(())
    }

    /// Extracts a subset of atoms/groups matching a selector into a new `AtomGroup`.
    ///
    /// Args:
    ///     selector (Union[Selector, Any]): A selector object (e.g. `SelectSymbol`,
    ///         `SelectName`, `SelectPath`, etc.) or any object providing an `is_match` method.
    ///
    /// Returns:
    ///     AtomGroup: A new atom group containing only the matched components.
    ///
    /// Raises:
    ///     TypeError: If `selector` is not a recognized selector type.
    pub fn select<'py>(slf: &Bound<'py, Self>, selector: &Bound<'py, PyAny>) -> PyResult<Self> {
        let py = slf.py();

        let sel = if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectSymbol>>() {
            ExtractedSelector::Symbol(s)
        } else if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectName>>() {
            ExtractedSelector::Name(s)
        } else if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectPathSimple>>() {
            ExtractedSelector::PathSimple(s)
        } else if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectPathWildcard>>() {
            ExtractedSelector::PathWildcard(s)
        } else if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectPathRegex>>() {
            ExtractedSelector::PathRegex(s)
        } else if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectPath>>() {
            ExtractedSelector::Path(s)
        } else if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectRange>>() {
            ExtractedSelector::Range(s)
        } else if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectAtom>>() {
            ExtractedSelector::Atom(s)
        } else if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectAtomGroup>>() {
            ExtractedSelector::AtomGroup(s)
        } else if selector.hasattr("is_match")? {
            ExtractedSelector::Custom(selector)
        } else {
            return Err(PyTypeError::new_err(format!(
                "Unsupported selector type: {}",
                selector.get_type()
            )));
        };

        let root_path = slf.borrow().compute_path(py);

        if sel.is_match_group(py, slf, &root_path)? {
            return slf.borrow().clone_for_select(py);
        }

        let child_groups: Vec<(String, Py<PyAtomGroup>)> = slf
            .borrow()
            .groups
            .iter()
            .map(|(k, g)| (k.clone(), g.clone_ref(py)))
            .collect();

        let direct_atoms: Vec<(String, Py<PyAtom>)> = slf
            .borrow()
            .atoms
            .iter()
            .map(|(k, a)| (k.clone(), a.clone_ref(py)))
            .collect();

        let mut answer_groups = IndexMap::new();
        for (k, child_py) in &child_groups {
            ensure_child_link(slf, k, child_py)?;
            let child_path = format!("{}{}/", root_path, k);
            if let Some(res_sub_py) = select_child_group(py, child_py, &child_path, k, None, &sel)?
            {
                answer_groups.insert(k.clone(), res_sub_py);
            }
        }

        let mut answer_atoms = IndexMap::new();
        for (k, atom_py) in &direct_atoms {
            let atom_path = format!("{}{}", root_path, k);
            if sel.is_match_atom(py, atom_py, &atom_path)? {
                let mut atom_copy = atom_py.borrow(py).__copy__(py)?;
                atom_copy.path = atom_path;
                answer_atoms.insert(k.clone(), Py::new(py, atom_copy)?);
            }
        }

        let self_borrow = slf.borrow();
        Ok(Self {
            name: self_borrow.name.clone(),
            atoms: answer_atoms,
            groups: answer_groups,
            bonds: Vec::new(),
            secondary_structure: None,
            parent: None,
            my_key: None,
            root_path: Some(root_path),
            sort_atoms: None,
            sort_groups: None,
        })
    }

    /// Registers a chemical bond between two atoms in the group hierarchy.
    ///
    /// The bond is stored at the lowest common ancestor group of the two atoms.
    ///
    /// Args:
    ///     atom1 (Atom): The first atom.
    ///     atom2 (Atom): The second atom.
    ///     order (int, optional): The bond order (1 for single, 2 for double, etc.). Defaults to 1.
    #[pyo3(signature = (atom1, atom2, order=1))]
    pub fn add_bond(
        slf: &Bound<'_, Self>,
        atom1: &PyAtom,
        atom2: &PyAtom,
        order: usize,
    ) -> PyResult<()> {
        let py = slf.py();
        let p1 = &atom1.path;
        let p2 = &atom2.path;

        let common_path = CoreAtomGroup::get_common_path(p1, p2);

        let query_str = PyString::new(py, &common_path);
        let target_group = if let Some(family) = Self::get_family(slf, query_str.as_any())? {
            family.into_bound(py)
        } else {
            slf.clone()
        };

        let target_path = target_group.borrow().compute_path(py);
        let common1 = CoreAtomGroup::get_common_path(&target_path, p1);
        let rel_p1 = if p1.starts_with(&common1) {
            &p1[common1.len()..]
        } else if p1.starts_with(&target_path) {
            &p1[target_path.len()..]
        } else {
            p1.as_str()
        }
        .trim_start_matches('/');

        let common2 = CoreAtomGroup::get_common_path(&target_path, p2);
        let rel_p2 = if p2.starts_with(&common2) {
            &p2[common2.len()..]
        } else if p2.starts_with(&target_path) {
            &p2[target_path.len()..]
        } else {
            p2.as_str()
        }
        .trim_start_matches('/');

        let b = BondRecord {
            atom1_path: rel_p1.to_string(),
            atom2_path: rel_p2.to_string(),
            order,
        };

        let mut target_borrow = target_group.borrow_mut();
        if !target_borrow.bonds.contains(&b) {
            target_borrow.bonds.push(b);
        }

        Ok(())
    }

    /// Returns a list of bond records stored directly in this group.
    ///
    /// Returns:
    ///     List[Tuple[str, str, int]]: List of `(atom1_path, atom2_path, order)` tuples.
    pub fn bonds(&self) -> Vec<(String, String, usize)> {
        self.bonds
            .iter()
            .map(|b| (b.atom1_path.clone(), b.atom2_path.clone(), b.order))
            .collect()
    }

    /// Returns a list of all bonds in this group and all its descendant groups.
    ///
    /// Returns:
    ///     List[Tuple[str, str, int]]: List of `(atom1_path, atom2_path, order)` tuples.
    pub fn get_bond_list(&self, py: Python<'_>) -> Vec<(String, String, usize)> {
        let base_path = self.compute_path(py);
        let mut bond_list = Vec::new();
        collect_bond_list_py(self, py, &base_path, &mut bond_list);
        bond_list
    }

    /// Assigns standard bonds, bond orders, and atom types using the built-in CCD template database.
    ///
    /// Raises:
    ///     BrValueError: If an atom has an unrecognized symbol or lacks covalent radius data.
    pub fn setup(&mut self, py: Python<'_>) -> PyResult<()> {
        let mut core = self.to_core(py);
        core.setup().map_err(to_py_err)?;
        self.update_from_core(py, core)?;
        Ok(())
    }

    /// Assigns bonds, bond orders, and atom types using a provided CCD template database.
    ///
    /// Args:
    ///     db (CcdTemplateDb): The Chemical Component Dictionary template database.
    ///
    /// Raises:
    ///     BrValueError: If an atom has an unrecognized symbol or lacks covalent radius data.
    pub fn setup_with_db(&mut self, py: Python<'_>, db: &PyCcdTemplateDb) -> PyResult<()> {
        let mut core = self.to_core(py);
        core.setup_with_db(&db.inner).map_err(to_py_err)?;
        self.update_from_core(py, core)?;
        Ok(())
    }

    /// Adds missing hydrogen atoms based on CCD residue templates.
    ///
    /// Args:
    ///     db (Optional[CcdTemplateDb], optional): The CCD template database to use.
    ///         If None, the global default database is used. Defaults to None.
    ///
    /// Returns:
    ///     OverallHydrogenationReport: Summary report of added hydrogens.
    #[pyo3(signature = (db=None))]
    pub fn add_missing_hydrogens(
        &mut self,
        py: Python<'_>,
        db: Option<&PyCcdTemplateDb>,
    ) -> PyResult<PyOverallHydrogenationReport> {
        let core_db = match db {
            Some(d) => &d.inner,
            None => CoreCcdTemplateDb::global(),
        };
        let mut core = self.to_core(py);
        let report = core.add_missing_hydrogens(core_db).map_err(to_py_err)?;
        self.update_from_core(py, core)?;
        PyOverallHydrogenationReport::new(py, report)
    }

    /// Checks for missing and unexpected atoms based on CCD residue templates.
    ///
    /// Args:
    ///     db (Optional[CcdTemplateDb], optional): The CCD template database to use.
    ///         If None, the global default database is used. Defaults to None.
    ///     include_hydrogens (bool, optional): Whether to include hydrogen atoms
    ///         in the check. Defaults to False.
    ///
    /// Returns:
    ///     MissingAtomsReport: Summary report of missing and unexpected atoms.
    #[pyo3(signature = (db=None, include_hydrogens=false))]
    pub fn check_missing_atoms(
        &self,
        py: Python<'_>,
        db: Option<&PyCcdTemplateDb>,
        include_hydrogens: bool,
    ) -> PyResult<PyMissingAtomsReport> {
        let core_db = match db {
            Some(d) => &d.inner,
            None => CoreCcdTemplateDb::global(),
        };
        let core = self.to_core(py);
        let options = MissingAtomsOptions { include_hydrogens };
        let report = core.check_missing_atoms(core_db, &options);
        PyMissingAtomsReport::new(py, report)
    }

    /// Reconstructs and adds missing heavy atoms based on CCD residue templates.
    ///
    /// Args:
    ///     db (Optional[CcdTemplateDb], optional): The CCD template database to use.
    ///         If None, the global default database is used. Defaults to None.
    ///     max_rmsd (float, optional): Maximum allowed anchor superposition RMSD
    ///         in Angstroms. Defaults to 1.0.
    ///
    /// Returns:
    ///     CompletionReport: Summary report of completed and skipped residues.
    #[pyo3(signature = (db=None, max_rmsd=1.0))]
    pub fn complete_missing_atoms(
        &mut self,
        py: Python<'_>,
        db: Option<&PyCcdTemplateDb>,
        max_rmsd: f64,
    ) -> PyResult<PyCompletionReport> {
        let core_db = match db {
            Some(d) => &d.inner,
            None => CoreCcdTemplateDb::global(),
        };
        let mut core = self.to_core(py);
        let options = CompleteOptions { max_rmsd };
        let report = core.complete_missing_atoms(core_db, &options);
        self.update_from_core(py, core)?;
        PyCompletionReport::new(py, report)
    }

    /// Detects and assigns secondary structure annotations ('H', 'E', '-') to residues in this chain.
    ///
    /// Modifications are applied in-place to residue groups.
    pub fn apply_secondary_structure(&mut self, py: Python<'_>) -> PyResult<()> {
        let mut core = self.to_core(py);
        core_apply_secondary_structure(&mut core);
        self.update_from_core(py, core)?;
        Ok(())
    }

    /// Translates all atomic coordinates in this hierarchy by a displacement vector.
    ///
    /// Args:
    ///     dir (Union[Position, str, Sequence[float]]): The translation displacement vector.
    ///
    /// Raises:
    ///     TypeError: If `dir` cannot be converted to a `Position`.
    pub fn shift_by(&mut self, py: Python<'_>, dir: &Bound<'_, PyAny>) -> PyResult<()> {
        let pos = extract_pos(dir)?;
        self.shift_by_position(py, pos);
        Ok(())
    }

    /// Rotates all atomic coordinates in this hierarchy using a 3x3 rotation matrix.
    ///
    /// Args:
    ///     rotmat (Union[Matrix, Any]): A 3x3 rotation matrix. Accepts a :class:`Matrix`
    ///         instance or any object providing `rows`, `cols`, and `get(row, col)`.
    ///
    /// Raises:
    ///     TypeError: If `rotmat` is not a :class:`Matrix` or lacks `rows`, `cols`, or `get`.
    ///     ValueError: If the dimensions of `rotmat` are not 3x3.
    pub fn rotate(&mut self, py: Python<'_>, rotmat: &Bound<'_, PyAny>) -> PyResult<()> {
        let core_mat = crate::matrix::extract_rotation_matrix(rotmat)?;
        self.rotate_core(py, &core_mat)
    }

    /// Returns the sum of atomic numbers of all atoms in this group and descendant groups.
    ///
    /// Returns:
    ///     float: Sum of atomic numbers.
    pub fn sum_of_atomic_number(&self) -> f64 {
        Python::attach(|py| {
            let mut sum = 0.0;
            self.for_each_atom_recursive(py, &mut |a| {
                sum += a.borrow().atomic_number as f64;
            });
            sum
        })
    }

    /// Reorganizes atoms to match the hierarchy of a reference `AtomGroup` based on coordinate proximity.
    ///
    /// Args:
    ///     reference (AtomGroup): Reference structure providing the target hierarchy.
    ///     range (float, optional): Proximity distance threshold in Angstroms. Defaults to 1.0e-5.
    ///
    /// Returns:
    ///     AtomGroup: A new restructured atom group.
    #[pyo3(signature = (reference, range=1.0e-5))]
    pub fn restructure(
        &self,
        py: Python<'_>,
        reference: &PyAtomGroup,
        range: f64,
    ) -> PyResult<PyAtomGroup> {
        let self_core = self.to_core(py);
        let ref_core = reference.to_core(py);
        let target_selector = proteindf_bridge::selector::SelectAtomGroup::new(&self_core, range);
        let mut restructured = ref_core.select(&target_selector);

        copy_attributes(&mut restructured, &self_core);

        let rest_of_target = get_rest_of_frame_molecule(&self_core, &restructured);
        assign_rest_molecule(&rest_of_target, &mut restructured, "model_1", "Z", "UNK")?;

        PyAtomGroup::from_core(py, restructured)
    }

    pub fn __getitem__<'py>(
        slf: &Bound<'py, Self>,
        key: &Bound<'_, PyAny>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let key_str = extract_key(key)?;
        let py = slf.py();

        // 1. Child group by key
        if slf.borrow().groups.contains_key(&key_str) {
            if let Some(grp) = Self::get_group(slf, key)? {
                return Ok(grp.into_bound(py).into_any());
            }
        }
        // 2. Child atom by key
        if slf.borrow().atoms.contains_key(&key_str) {
            if let Some(atom) = Self::get_atom(slf, key)? {
                return Ok(atom.into_bound(py).into_any());
            }
        }
        // 3. Search child groups by name
        let found_grp = slf
            .borrow()
            .groups
            .iter()
            .find(|(_, g)| g.borrow(py).name == key_str)
            .map(|(k, _)| k.clone());
        if let Some(k) = found_grp {
            let k_py = k.into_pyobject(py)?;
            if let Some(grp) = Self::get_group(slf, k_py.as_any())? {
                return Ok(grp.into_bound(py).into_any());
            }
        }
        // 4. Search child atoms by name
        let found_atom = slf
            .borrow()
            .atoms
            .iter()
            .find(|(_, a)| a.borrow(py).name.trim() == key_str.trim())
            .map(|(k, _)| k.clone());
        if let Some(k) = found_atom {
            let k_py = k.into_pyobject(py)?;
            if let Some(atom) = Self::get_atom(slf, k_py.as_any())? {
                return Ok(atom.into_bound(py).into_any());
            }
        }

        Err(PyKeyError::new_err(key_str))
    }

    pub fn __setitem__(
        slf: &Bound<'_, Self>,
        key: &Bound<'_, PyAny>,
        value: &Bound<'_, PyAny>,
    ) -> PyResult<()> {
        if let Ok(grp) = value.extract::<PyRef<PyAtomGroup>>() {
            let k = extract_key(key)?;
            let py = slf.py();
            let parent_py = slf.clone().unbind();

            let parent_path = slf.borrow().compute_path(py);
            let copy_py = Self::deep_clone_py(py, &grp, Some(parent_py), Some(k.clone()))?;
            let old = slf.borrow_mut().groups.insert(k.clone(), copy_py);
            if let Some(old_group) = old {
                let mut b = old_group.borrow_mut(py);
                let old_path = format!("{}{}/", parent_path, k);
                b.parent = None;
                b.my_key = None;
                b.root_path = Some(old_path);
            }
            Ok(())
        } else if let Ok(atom) = value.extract::<PyRef<PyAtom>>() {
            let k = extract_key(key)?;
            let py = slf.py();
            let mut atom_copy = atom.__copy__(py)?;
            let group_path = slf.borrow().compute_path(py);
            atom_copy.path = format!("{}{}", group_path, k);
            let atom_py = Py::new(py, atom_copy)?;
            slf.borrow_mut().atoms.insert(k, atom_py);
            Ok(())
        } else {
            Err(PyValueError::new_err(format!(
                "Expected AtomGroup or Atom, got {}",
                value.get_type()
            )))
        }
    }

    pub fn __and__(&self, py: Python<'_>, other: &PyAtomGroup) -> PyResult<PyAtomGroup> {
        let c1 = self.to_core(py);
        let c2 = other.to_core(py);
        PyAtomGroup::from_core(py, &c1 & &c2)
    }

    pub fn __or__(&self, py: Python<'_>, other: &PyAtomGroup) -> PyResult<PyAtomGroup> {
        let c1 = self.to_core(py);
        let c2 = other.to_core(py);
        PyAtomGroup::from_core(py, &c1 | &c2)
    }

    pub fn __xor__(&self, py: Python<'_>, other: &PyAtomGroup) -> PyResult<PyAtomGroup> {
        let c1 = self.to_core(py);
        let c2 = other.to_core(py);
        PyAtomGroup::from_core(py, &c1 ^ &c2)
    }

    /// Multiplies the coordinates of all atoms in the hierarchy by a scalar in-place.
    ///
    /// Args:
    ///     scalar (float): Scaling factor to multiply atomic coordinates by.
    ///
    /// Raises:
    ///     TypeError: If `scalar` is not a numeric type.
    pub fn __imul__(&mut self, py: Python<'_>, scalar: f64) {
        self.for_each_atom_recursive(py, &mut |a| {
            a.borrow().xyz.borrow_mut(py).inner *= scalar;
        });
    }

    pub fn __len__(&self) -> usize {
        self.atoms.len()
    }

    pub fn __repr__(&self) -> String {
        Python::attach(|py| {
            format!(
                "AtomGroup(name='{}', path='{}', atoms={}, groups={})",
                self.name,
                self.compute_path(py),
                self.atoms.len(),
                self.groups.len()
            )
        })
    }

    pub fn __str__(&self, py: Python<'_>) -> String {
        let mut out = String::new();
        let mut root_path = self.compute_path(py);
        if root_path.is_empty() || !root_path.ends_with('/') {
            root_path.push('/');
        }
        format_atom_group_indented(self, py, &mut out, "", 0, &root_path);
        out
    }

    /// Serializes the hierarchy tree into a nested dictionary representation.
    ///
    /// Returns:
    ///     dict: Nested dictionary containing group and atom data.
    pub fn get_raw_data<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        build_atomgroup_raw_data(self, py)
    }

    /// Replaces the group content with data loaded from a nested dictionary.
    ///
    /// Unrecognized keys at top level are ignored.
    ///
    /// Args:
    ///     data (dict): Nested dictionary representing tree data.
    ///
    /// Returns:
    ///     AtomGroup: Self, populated with the dictionary data.
    ///
    /// Raises:
    ///     TypeError: If `data` is not a dict.
    ///     BrValueError: If nested atom data within the dictionary is invalid.
    pub fn set_by_dict_data<'py>(
        slf: &Bound<'py, Self>,
        data: &Bound<'_, pyo3::types::PyDict>,
    ) -> PyResult<Bound<'py, Self>> {
        let py = slf.py();
        let val = crate::value::py_to_value(py, data.as_any())?;
        let mut core = CoreAtomGroup::new();
        proteindf_bridge::brd::atomgroup_set_by_dict_data(&mut core, &val).map_err(to_py_err)?;
        let mut borrow = slf.borrow_mut();
        borrow.update_from_core(py, core)?;
        drop(borrow);
        for (k, g) in &slf.borrow().groups {
            ensure_child_link(slf, k, g)?;
        }
        Ok(slf.clone())
    }

    pub fn __copy__(slf: &Bound<'_, Self>) -> PyResult<Py<Self>> {
        let py = slf.py();
        let clone = slf.borrow().deep_clone(py)?;
        let new_py = Py::new(py, clone)?;
        for (k, child_py) in &new_py.borrow(py).groups {
            let mut b = child_py.borrow_mut(py);
            if b.parent.is_none() {
                b.parent = Some(new_py.clone_ref(py));
            }
            b.my_key = Some(k.clone());
        }
        Ok(new_py)
    }

    pub fn __deepcopy__(slf: &Bound<'_, Self>, _memo: &Bound<'_, PyAny>) -> PyResult<Py<Self>> {
        Self::__copy__(slf)
    }

    pub fn __getstate__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let raw_py = build_atomgroup_raw_data(self, py)?;
        let path = self.compute_path(py);

        let mut ss_list = Vec::new();
        collect_secondary_structures_py(self, py, &path, &mut ss_list);

        let ss_list_py = pyo3::types::PyList::empty(py);
        for (p, ss) in ss_list {
            let tup = (p, ss).into_pyobject(py)?;
            ss_list_py.append(tup)?;
        }

        let mut sort_list = Vec::new();
        let mut key_path = Vec::new();
        collect_sort_settings(self, py, &mut key_path, &mut sort_list);
        let sort_list_py = pyo3::types::PyList::empty(py);
        for (kp, sa, sg) in sort_list {
            let tup = (kp, sa, sg).into_pyobject(py)?;
            sort_list_py.append(tup)?;
        }

        let state_tuple = pyo3::types::PyTuple::new(
            py,
            &[
                raw_py.into_any(),
                path.into_pyobject(py)?.into_any(),
                ss_list_py.into_any(),
                sort_list_py.into_any(),
            ],
        )?;
        Ok(state_tuple.into_any())
    }

    pub fn __setstate__(&mut self, py: Python<'_>, state: &Bound<'_, PyAny>) -> PyResult<()> {
        let tup = state.cast::<pyo3::types::PyTuple>()?;
        if tup.len() < 3 {
            return Err(PyValueError::new_err("Invalid state tuple for AtomGroup"));
        }
        let raw_obj = tup.get_item(0)?;
        let path_obj = tup.get_item(1)?;
        let ss_obj = tup.get_item(2)?;

        let val = crate::value::py_to_value(py, &raw_obj)?;
        let mut core = CoreAtomGroup::new();
        proteindf_bridge::brd::atomgroup_set_by_dict_data(&mut core, &val).map_err(to_py_err)?;

        let path = path_obj.extract::<String>()?;
        let depth = path
            .trim_matches('/')
            .split('/')
            .filter(|s| !s.is_empty())
            .count();
        core.set_path_with_depth(path.clone(), depth);

        let ss_list = ss_obj.extract::<Vec<(String, String)>>()?;
        for (p, ss_str) in ss_list {
            let ss_code = match ss_str.as_str() {
                "H" => Some(SsCode::Helix),
                "E" => Some(SsCode::Strand),
                "-" => Some(SsCode::Loop),
                _ => None,
            };
            if let Some(target) = find_group_mut_by_path(&mut core, &p) {
                target.set_secondary_structure(ss_code);
            }
        }

        let mut restored = PyAtomGroup::from_core(py, core)?;

        if tup.len() >= 4 {
            let sort_obj = tup.get_item(3)?;
            let sort_list =
                sort_obj.extract::<Vec<(Vec<String>, Option<String>, Option<String>)>>()?;
            for (kp, sa, sg) in sort_list {
                restore_sort_settings(&mut restored, py, &kp, sa, sg);
            }
        }

        *self = restored;
        Ok(())
    }

    pub fn __reduce__<'py>(
        slf: Bound<'py, Self>,
    ) -> PyResult<(
        Bound<'py, PyAny>,
        Bound<'py, pyo3::types::PyTuple>,
        Bound<'py, PyAny>,
    )> {
        let py = slf.py();
        let cls = slf.get_type().into_any();
        let empty_args = pyo3::types::PyTuple::empty(py);
        let state = slf.borrow().__getstate__(py)?;
        Ok((cls, empty_args, state))
    }

    /// Exports the hierarchical atom and group data to a CSV file.
    ///
    /// The CSV is written atomically via a temporary file.
    ///
    /// Args:
    ///     path (str): Destination file path for the CSV file.
    ///
    /// Raises:
    ///     BrInputError: If file creation, writing, or renaming fails.
    pub fn save_csv(&self, py: Python<'_>, path: &str) -> PyResult<()> {
        let path_ref = std::path::Path::new(path);
        let parent = path_ref
            .parent()
            .unwrap_or_else(|| std::path::Path::new("."));
        let file_name = path_ref
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("output.csv");
        let temp_name = format!(
            ".{file_name}.tmp_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap_or_default()
                .as_nanos()
        );
        let temp_path = parent.join(temp_name);

        struct TempFileGuard<'a>(&'a std::path::Path, bool);
        impl<'a> Drop for TempFileGuard<'a> {
            fn drop(&mut self) {
                if !self.1 {
                    let _ = std::fs::remove_file(self.0);
                }
            }
        }
        let mut guard = TempFileGuard(&temp_path, false);

        let file = std::fs::File::create(&temp_path).map_err(|e| {
            to_py_err(proteindf_bridge::error::BridgeError::input_error(
                temp_path.display().to_string(),
                e.to_string(),
            ))
        })?;
        let mut writer = std::io::BufWriter::new(file);

        let mut rows = Vec::new();
        collect_csv_rows(self, &[], &mut rows, py)?;

        use std::io::Write;
        for row in rows {
            for (i, field) in row.iter().enumerate() {
                if i > 0 {
                    writer.write_all(b",").map_err(|e| {
                        to_py_err(proteindf_bridge::error::BridgeError::input_error(
                            temp_path.display().to_string(),
                            e.to_string(),
                        ))
                    })?;
                }
                write_csv_field(&mut writer, field).map_err(|e| {
                    to_py_err(proteindf_bridge::error::BridgeError::input_error(
                        temp_path.display().to_string(),
                        e.to_string(),
                    ))
                })?;
            }
            writer.write_all(b"\r\n").map_err(|e| {
                to_py_err(proteindf_bridge::error::BridgeError::input_error(
                    temp_path.display().to_string(),
                    e.to_string(),
                ))
            })?;
        }

        writer.flush().map_err(|e| {
            to_py_err(proteindf_bridge::error::BridgeError::input_error(
                temp_path.display().to_string(),
                e.to_string(),
            ))
        })?;
        drop(writer);

        std::fs::rename(&temp_path, path_ref).map_err(|e| {
            to_py_err(proteindf_bridge::error::BridgeError::input_error(
                path_ref.display().to_string(),
                e.to_string(),
            ))
        })?;
        guard.1 = true;

        Ok(())
    }
}

fn format_atom_group_indented(
    group: &PyAtomGroup,
    py: Python<'_>,
    f: &mut String,
    key: &str,
    depth: usize,
    current_path: &str,
) {
    use std::fmt::Write;
    let indent = "  ".repeat(depth);
    let _ = writeln!(f, "{}# group key={} name={}", indent, key, group.name);
    for (sub_key, child_py) in &group.groups {
        let child = child_py.borrow(py);
        let child_path = format!("{}{}/", current_path, sub_key);
        format_atom_group_indented(&child, py, f, sub_key, depth + 1, &child_path);
    }
    for (atom_key, atom_py) in &group.atoms {
        let mut atom_core = atom_py.borrow(py).to_core(py);
        atom_core.path = format!("{}{}", current_path, atom_key);
        let _ = writeln!(f, "{}{}", indent, atom_core);
    }
    for bond in &group.bonds {
        let _ = writeln!(
            f,
            "{}bond {} <-{}-> {}",
            indent, bond.atom1_path, bond.order, bond.atom2_path
        );
    }
}

fn collect_csv_rows(
    ag: &PyAtomGroup,
    parents: &[String],
    rows: &mut Vec<Vec<String>>,
    py: Python<'_>,
) -> PyResult<()> {
    for (key, subgrp_py) in &ag.groups {
        let subgrp = subgrp_py.borrow(py);
        let mut new_parents = parents.to_vec();
        new_parents.push(key.clone());
        new_parents.push(subgrp.name.clone());
        collect_csv_rows(&subgrp, &new_parents, rows, py)?;
    }
    for (_key, atom_py) in &ag.atoms {
        let atom = atom_py.borrow(py);
        let mut row = parents.to_vec();
        let sym = proteindf_bridge::periodic_table::PeriodicTable::get_symbol(atom.atomic_number)
            .unwrap_or("X")
            .to_string();
        let xyz = atom.xyz.borrow(py).inner;
        row.push(sym);
        row.push(format_py_float(py, xyz.x)?);
        row.push(format_py_float(py, xyz.y)?);
        row.push(format_py_float(py, xyz.z)?);
        row.push(format_py_float(py, atom.charge)?);
        rows.push(row);
    }
    Ok(())
}

fn format_py_float(py: Python<'_>, v: f64) -> PyResult<String> {
    pyo3::types::PyFloat::new(py, v).repr()?.extract::<String>()
}

fn write_csv_field<W: std::io::Write>(w: &mut W, s: &str) -> std::io::Result<()> {
    let need_quote = s.contains(',') || s.contains('"') || s.contains('\n') || s.contains('\r');
    if need_quote {
        w.write_all(b"\"")?;
        for b in s.bytes() {
            if b == b'"' {
                w.write_all(b"\"\"")?;
            } else {
                w.write_all(&[b])?;
            }
        }
        w.write_all(b"\"")?;
    } else {
        w.write_all(s.as_bytes())?;
    }
    Ok(())
}

fn copy_attributes(target: &mut CoreAtomGroup, reference: &CoreAtomGroup) {
    for (_key, subgrp) in target.groups_mut() {
        copy_attributes(subgrp, reference);
    }
    for (_key, atom) in target.atoms_mut() {
        let selector = proteindf_bridge::selector::SelectAtom::new(atom.clone(), 0.1);
        let ref_selected = reference.select(&selector);
        let ref_atoms = ref_selected.get_atom_list();
        if !ref_atoms.is_empty() {
            atom.charge = ref_atoms[0].charge;
        }
    }
}

fn get_rest_of_frame_molecule(
    frame_molecule: &CoreAtomGroup,
    selected_molecule: &CoreAtomGroup,
) -> CoreAtomGroup {
    let selector = proteindf_bridge::selector::SelectAtomGroup::new(selected_molecule, 1.0e-5);
    let selected = frame_molecule.select(&selector);
    frame_molecule ^ &selected
}

fn assign_rest_molecule(
    rest_molecule: &CoreAtomGroup,
    output_atom_group: &mut CoreAtomGroup,
    model_id: &str,
    chain_id: &str,
    res_name: &str,
) -> PyResult<()> {
    let mut chain = CoreAtomGroup::new();
    let mut res = CoreAtomGroup::new();
    res.name = res_name.to_string();
    for (i, atom) in rest_molecule.get_atom_list().into_iter().enumerate() {
        let atom_id = i + 1;
        res.set_atom(&atom_id.to_string(), atom);
    }
    chain.set_group("1", res);

    if let Some(model) = output_atom_group.get_group_mut(model_id) {
        model.set_group(chain_id, chain);
        Ok(())
    } else {
        Err(PyKeyError::new_err(model_id.to_string()))
    }
}

fn collect_bond_list_py(
    group: &PyAtomGroup,
    py: Python<'_>,
    current_path: &str,
    bond_list: &mut Vec<(String, String, usize)>,
) {
    for (key, child) in &group.groups {
        let child_borrow = child.borrow(py);
        let child_path = format!("{}{}/", current_path, key);
        collect_bond_list_py(&child_borrow, py, &child_path, bond_list);
    }
    for b in &group.bonds {
        let path1 = if b.atom1_path.starts_with('/') {
            b.atom1_path.clone()
        } else {
            format!("{}{}", current_path, b.atom1_path)
        };
        let path2 = if b.atom2_path.starts_with('/') {
            b.atom2_path.clone()
        } else {
            format!("{}{}", current_path, b.atom2_path)
        };
        bond_list.push((path1, path2, b.order));
    }
}

fn collect_secondary_structures_py(
    group: &PyAtomGroup,
    py: Python<'_>,
    current_path: &str,
    out: &mut Vec<(String, String)>,
) {
    if let Some(ss) = group.secondary_structure {
        out.push((current_path.to_string(), ss.as_char().to_string()));
    }
    for (key, child) in &group.groups {
        let child_borrow = child.borrow(py);
        let child_path = format!("{}{}/", current_path, key);
        collect_secondary_structures_py(&child_borrow, py, &child_path, out);
    }
}

fn build_atomgroup_raw_data<'py>(
    group: &PyAtomGroup,
    py: Python<'py>,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    if !group.groups.is_empty() {
        let grp_dict = PyDict::new(py);
        for (key, grp) in &group.groups {
            grp_dict.set_item(key, build_atomgroup_raw_data(&grp.borrow(py), py)?)?;
        }
        dict.set_item("groups", grp_dict)?;
    }
    if !group.atoms.is_empty() {
        let atm_dict = PyDict::new(py);
        for (key, atm) in &group.atoms {
            atm_dict.set_item(key, atm.borrow(py).get_raw_data(py)?)?;
        }
        dict.set_item("atoms", atm_dict)?;
    }
    dict.set_item("name", &group.name)?;
    if !group.bonds.is_empty() {
        let bonds_list = PyList::empty(py);
        for b in &group.bonds {
            let entry = PyList::new(
                py,
                &[
                    b.atom1_path.as_str().into_pyobject(py)?.into_any(),
                    b.atom2_path.as_str().into_pyobject(py)?.into_any(),
                    b.order.into_pyobject(py)?.into_any(),
                ],
            )?;
            bonds_list.append(entry)?;
        }
        dict.set_item("bonds", bonds_list)?;
    }
    Ok(dict)
}

fn find_group_mut_by_path<'a>(
    group: &'a mut CoreAtomGroup,
    path: &str,
) -> Option<&'a mut CoreAtomGroup> {
    if group.path() == path {
        return Some(group);
    }
    for (_, child) in group.groups_mut() {
        if let Some(found) = find_group_mut_by_path(child, path) {
            return Some(found);
        }
    }
    None
}

fn collect_sort_settings(
    group: &PyAtomGroup,
    py: Python<'_>,
    current_key_path: &mut Vec<String>,
    out: &mut Vec<(Vec<String>, Option<String>, Option<String>)>,
) {
    if group.sort_atoms.is_some() || group.sort_groups.is_some() {
        out.push((
            current_key_path.clone(),
            group.sort_atoms.clone(),
            group.sort_groups.clone(),
        ));
    }
    for (k, child) in &group.groups {
        current_key_path.push(k.clone());
        collect_sort_settings(&child.borrow(py), py, current_key_path, out);
        current_key_path.pop();
    }
}

fn restore_sort_settings(
    group: &mut PyAtomGroup,
    py: Python<'_>,
    key_path: &[String],
    sort_atoms: Option<String>,
    sort_groups: Option<String>,
) {
    if key_path.is_empty() {
        group.sort_atoms = sort_atoms;
        group.sort_groups = sort_groups;
        return;
    }
    let next_key = &key_path[0];
    if let Some(child) = group.groups.get(next_key) {
        restore_sort_settings(
            &mut child.borrow_mut(py),
            py,
            &key_path[1..],
            sort_atoms,
            sort_groups,
        );
    }
}

enum ExtractedSelector<'a, 'py> {
    Symbol(PyRef<'a, crate::selector::PySelectSymbol>),
    Name(PyRef<'a, crate::selector::PySelectName>),
    PathSimple(PyRef<'a, crate::selector::PySelectPathSimple>),
    PathWildcard(PyRef<'a, crate::selector::PySelectPathWildcard>),
    PathRegex(PyRef<'a, crate::selector::PySelectPathRegex>),
    Path(PyRef<'a, crate::selector::PySelectPath>),
    Range(PyRef<'a, crate::selector::PySelectRange>),
    Atom(PyRef<'a, crate::selector::PySelectAtom>),
    AtomGroup(PyRef<'a, crate::selector::PySelectAtomGroup>),
    Custom(&'a Bound<'py, PyAny>),
}

impl<'a, 'py> ExtractedSelector<'a, 'py> {
    fn as_core_selector(&self) -> Option<&dyn Selector> {
        match self {
            Self::Symbol(s) => Some(&s.inner),
            Self::Name(s) => Some(&s.inner),
            Self::PathSimple(s) => Some(&s.inner),
            Self::PathWildcard(s) => Some(&s.inner),
            Self::PathRegex(s) => Some(&s.inner),
            Self::Path(s) => Some(&s.inner),
            Self::Range(s) => Some(&s.inner),
            Self::Atom(s) => Some(&s.inner),
            Self::AtomGroup(s) => Some(&s.inner),
            Self::Custom(_) => None,
        }
    }

    fn is_match_group(
        &self,
        _py: Python<'py>,
        group_bound: &Bound<'py, PyAtomGroup>,
        current_path: &str,
    ) -> PyResult<bool> {
        if let Some(core_sel) = self.as_core_selector() {
            let name = group_bound.borrow().name.clone();
            let mut dummy = CoreAtomGroup::with_name(&name);
            dummy.set_path(current_path.to_string());
            Ok(core_sel.is_match_group(&dummy))
        } else if let Self::Custom(obj) = self {
            let res = obj.call_method1("is_match", (group_bound,))?;
            res.is_truthy()
        } else {
            Ok(false)
        }
    }

    fn is_match_atom(
        &self,
        py: Python<'py>,
        atom_py: &Py<PyAtom>,
        atom_path: &str,
    ) -> PyResult<bool> {
        if let Some(core_sel) = self.as_core_selector() {
            let mut core_atom = atom_py.borrow(py).to_core(py);
            core_atom.path = atom_path.to_string();
            Ok(core_sel.is_match_atom(&core_atom))
        } else if let Self::Custom(obj) = self {
            atom_py.borrow_mut(py).path = atom_path.to_string();
            let res = obj.call_method1("is_match", (atom_py.bind(py),))?;
            res.is_truthy()
        } else {
            Ok(false)
        }
    }
}

fn select_child_group<'py>(
    py: Python<'py>,
    child_group_py: &Py<PyAtomGroup>,
    current_path: &str,
    my_key: &str,
    parent_py: Option<Py<PyAtomGroup>>,
    selector: &ExtractedSelector<'_, 'py>,
) -> PyResult<Option<Py<PyAtomGroup>>> {
    let child_bound = child_group_py.bind(py);

    if selector.is_match_group(py, child_bound, current_path)? {
        if child_group_py.borrow(py).get_number_of_all_atoms_py(py) == 0 {
            return Ok(None);
        }
        let cloned = PyAtomGroup::clone_subtree_for_select(
            py,
            &child_group_py.borrow(py),
            parent_py,
            Some(my_key.to_string()),
        )?;
        return Ok(Some(cloned));
    }

    let child_groups: Vec<(String, Py<PyAtomGroup>)> = child_group_py
        .borrow(py)
        .groups
        .iter()
        .map(|(k, g)| (k.clone(), g.clone_ref(py)))
        .collect();

    let direct_atoms: Vec<(String, Py<PyAtom>)> = child_group_py
        .borrow(py)
        .atoms
        .iter()
        .map(|(k, a)| (k.clone(), a.clone_ref(py)))
        .collect();

    let mut sub_groups = IndexMap::new();
    for (sub_k, sub_py) in &child_groups {
        ensure_child_link(child_bound, sub_k, sub_py)?;
        let sub_path = format!("{}{}/", current_path, sub_k);
        if let Some(res_sub_py) = select_child_group(py, sub_py, &sub_path, sub_k, None, selector)?
        {
            sub_groups.insert(sub_k.clone(), res_sub_py);
        }
    }

    let mut sub_atoms = IndexMap::new();
    for (atm_k, atm_py) in &direct_atoms {
        let atom_path = format!("{}{}", current_path, atm_k);
        if selector.is_match_atom(py, atm_py, &atom_path)? {
            let mut atom_copy = atm_py.borrow(py).__copy__(py)?;
            atom_copy.path = atom_path;
            sub_atoms.insert(atm_k.clone(), Py::new(py, atom_copy)?);
        }
    }

    let total_atoms = sub_atoms.len()
        + sub_groups
            .values()
            .map(|g| g.borrow(py).get_number_of_all_atoms_py(py))
            .sum::<usize>();

    if total_atoms == 0 {
        return Ok(None);
    }

    let child_borrow = child_group_py.borrow(py);
    let new_group = PyAtomGroup {
        name: child_borrow.name.clone(),
        atoms: sub_atoms,
        groups: IndexMap::new(),
        bonds: Vec::new(),
        secondary_structure: None,
        parent: parent_py,
        my_key: Some(my_key.to_string()),
        root_path: None,
        sort_atoms: None,
        sort_groups: None,
    };
    drop(child_borrow);

    let new_group_py = Py::new(py, new_group)?;

    for (sub_k, res_sub_py) in sub_groups {
        res_sub_py.borrow_mut(py).parent = Some(new_group_py.clone_ref(py));
        new_group_py.borrow_mut(py).groups.insert(sub_k, res_sub_py);
    }

    Ok(Some(new_group_py))
}
