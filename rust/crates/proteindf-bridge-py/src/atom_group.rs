// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom::PyAtom;
use crate::ccd_templates::PyCcdTemplateDb;
use crate::error::to_py_err;
use crate::hydrogenation::PyOverallHydrogenationReport;
use crate::position::PyPosition;
use crate::schema::PySchemaViolation;
use indexmap::IndexMap;
use proteindf_bridge::atom_group::{AtomGroup as CoreAtomGroup, BondRecord, Selector};
use proteindf_bridge::ccd_templates::CcdTemplateDb as CoreCcdTemplateDb;
use proteindf_bridge::position::Position;
use proteindf_bridge::secondary_structure::{
    apply_secondary_structure as core_apply_secondary_structure, SsCode,
};
use pyo3::exceptions::{PyDeprecationWarning, PyKeyError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PySet;
use std::collections::HashMap;
use std::str::FromStr;

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
        core.set_path_with_depth(path.clone(), depth);
        for (k, g) in &self.groups {
            let child_path = format!("{}{}/", path, k);
            let child_core = g.borrow(py).to_core_with_path(py, child_path, depth + 1);
            core.set_group(k, child_core);
        }
        for (k, a) in &self.atoms {
            core.set_atom(k, a.borrow(py).to_core(py));
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

    #[getter]
    pub fn name(&self) -> String {
        self.name.clone()
    }

    #[setter]
    pub fn set_name(&mut self, val: String) {
        self.name = val;
    }

    #[getter]
    pub fn parent(&self, py: Python<'_>) -> Option<Py<PyAtomGroup>> {
        self.get_parent(py)
    }

    #[getter]
    pub fn path(&self) -> String {
        Python::attach(|py| self.compute_path(py))
    }

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

    pub fn path_depth(&self) -> usize {
        Python::attach(|py| self.compute_depth(py))
    }

    pub fn is_model_level(&self) -> bool {
        self.path_depth() == 1
    }

    pub fn is_chain_level(&self) -> bool {
        self.path_depth() == 2
    }

    pub fn is_residue_level(&self) -> bool {
        self.path_depth() == 3
    }

    pub fn validate_schema(&self) -> Vec<PySchemaViolation> {
        Python::attach(|py| {
            self.to_core(py)
                .validate_schema()
                .into_iter()
                .map(PySchemaViolation::from)
                .collect()
        })
    }

    #[getter]
    pub fn secondary_structure(&self) -> Option<String> {
        self.secondary_structure.map(|ss| ss.to_string())
    }

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

    pub fn get_number_of_atoms(&self) -> usize {
        self.atoms.len()
    }

    pub fn get_number_of_all_atoms(&self) -> usize {
        Python::attach(|py| self.to_core(py).get_number_of_all_atoms())
    }

    pub fn get_number_of_groups(&self) -> usize {
        self.groups.len()
    }

    pub fn get_number_of_bonds(&self) -> usize {
        self.bonds.len()
    }

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

    pub fn has_atom(&self, py: Python<'_>, key: &Bound<'_, PyAny>) -> PyResult<bool> {
        let k = extract_key(key)?;
        Ok(self.atoms.contains_key(&k)
            || self
                .atoms
                .values()
                .any(|a| a.borrow(py).name.trim() == k.trim()))
    }

    pub fn has_atomkey(&self, key: &Bound<'_, PyAny>) -> PyResult<bool> {
        let k = extract_key(key)?;
        Ok(self.atoms.contains_key(&k))
    }

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

    pub fn del_atom(
        &mut self,
        _py: Python<'_>,
        key: &Bound<'_, PyAny>,
    ) -> PyResult<Option<Py<PyAtom>>> {
        let k = extract_key(key)?;
        Ok(self.atoms.shift_remove(&k))
    }

    pub fn remove_atom(
        &mut self,
        py: Python<'_>,
        key: &Bound<'_, PyAny>,
    ) -> PyResult<Option<Py<PyAtom>>> {
        self.del_atom(py, key)
    }

    pub fn erase_atom(&mut self, py: Python<'_>, key: &Bound<'_, PyAny>) -> PyResult<()> {
        warn_deprecation(py, "erase_atom() is deprecated, use remove_atom() instead.")?;
        self.remove_atom(py, key)?;
        Ok(())
    }

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

    pub fn has_group(&self, py: Python<'_>, key: &Bound<'_, PyAny>) -> PyResult<bool> {
        let k = extract_key(key)?;
        Ok(self.groups.contains_key(&k) || self.groups.values().any(|g| g.borrow(py).name == k))
    }

    pub fn has_groupkey(&self, key: &Bound<'_, PyAny>) -> PyResult<bool> {
        let k = extract_key(key)?;
        Ok(self.groups.contains_key(&k))
    }

    pub fn has_groupname(&self, py: Python<'_>, name: &Bound<'_, PyAny>) -> PyResult<bool> {
        let target = extract_key(name)?;
        for g in self.groups.values() {
            if g.borrow(py).name == target {
                return Ok(true);
            }
        }
        Ok(false)
    }

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

    pub fn remove_group(
        slf: &Bound<'_, Self>,
        key: &Bound<'_, PyAny>,
    ) -> PyResult<Option<Py<PyAtomGroup>>> {
        Self::del_group(slf, key)
    }

    pub fn erase_group(slf: &Bound<'_, Self>, key: &Bound<'_, PyAny>) -> PyResult<()> {
        let py = slf.py();
        warn_deprecation(
            py,
            "erase_group() is deprecated, use remove_group() instead.",
        )?;
        Self::remove_group(slf, key)?;
        Ok(())
    }

    #[getter]
    pub fn sort_atoms(&self) -> Option<String> {
        self.sort_atoms.clone()
    }

    #[setter]
    pub fn set_sort_atoms(&mut self, val: Option<String>) {
        self.sort_atoms = val;
    }

    #[getter]
    pub fn sort_groups(&self) -> Option<String> {
        self.sort_groups.clone()
    }

    #[setter]
    pub fn set_sort_groups(&mut self, val: Option<String>) {
        self.sort_groups = val;
    }

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

    pub fn get_group_list(&self) -> Vec<String> {
        let mut keys: Vec<String> = self.groups.keys().cloned().collect();
        if self.sort_groups.as_deref() == Some("nice") {
            proteindf_bridge::brd::sort_nicely(&mut keys);
        }
        keys
    }

    pub fn get_atom_keys(&self) -> Vec<String> {
        let mut keys: Vec<String> = self.atoms.keys().cloned().collect();
        if self.sort_atoms.as_deref() == Some("nice") {
            proteindf_bridge::brd::sort_nicely(&mut keys);
        }
        keys
    }

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

    pub fn get_path_list(slf: &Bound<'_, Self>) -> PyResult<Vec<String>> {
        let py = slf.py();
        let atoms = Self::get_atom_list(slf)?;
        Ok(atoms.iter().map(|a| a.borrow(py).path.clone()).collect())
    }

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

    pub fn get_atom_kinds<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PySet>> {
        let set = PySet::empty(py)?;
        self.collect_atom_kinds(py, &set)?;
        Ok(set)
    }

    pub fn get_xyz(&self, py: Python<'_>) -> PyResult<String> {
        let count = self.get_number_of_all_atoms();
        let mut output = format!("{}\n# \n", count);
        self.collect_xyz(py, &mut output)?;
        Ok(output)
    }

    pub fn get_atom_kinds_count(&self) -> HashMap<String, usize> {
        Python::attach(|py| self.to_core(py).get_atom_kinds_count())
    }

    pub fn get_formula(&self) -> String {
        Python::attach(|py| self.to_core(py).get_formula())
    }

    pub fn formula(&self) -> String {
        self.get_formula()
    }

    #[getter]
    pub fn weight(&self) -> f64 {
        Python::attach(|py| {
            self.to_core(py)
                .get_atom_list()
                .iter()
                .map(|a| a.weight().unwrap_or(0.0))
                .sum()
        })
    }

    pub fn center(&self) -> PyPosition {
        Python::attach(|py| PyPosition::from_core(self.to_core(py).center()))
    }

    pub fn r#box(&self) -> (PyPosition, PyPosition) {
        Python::attach(|py| {
            let (bmin, bmax) = self.to_core(py).r#box();
            (PyPosition::from_core(bmin), PyPosition::from_core(bmax))
        })
    }

    pub fn get_box(&self) -> (PyPosition, PyPosition) {
        self.r#box()
    }

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

    pub fn select(&self, py: Python<'_>, selector: &Bound<'_, PyAny>) -> PyResult<PyAtomGroup> {
        let core = self.to_core(py);
        if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectSymbol>>() {
            return PyAtomGroup::from_core(py, core.select(&s.inner));
        }
        if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectName>>() {
            return PyAtomGroup::from_core(py, core.select(&s.inner));
        }
        if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectPathSimple>>() {
            return PyAtomGroup::from_core(py, core.select(&s.inner));
        }
        if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectPathWildcard>>() {
            return PyAtomGroup::from_core(py, core.select(&s.inner));
        }
        if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectPathRegex>>() {
            return PyAtomGroup::from_core(py, core.select(&s.inner));
        }
        if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectPath>>() {
            return PyAtomGroup::from_core(py, core.select(&s.inner));
        }
        if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectRange>>() {
            return PyAtomGroup::from_core(py, core.select(&s.inner));
        }
        if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectAtom>>() {
            return PyAtomGroup::from_core(py, core.select(&s.inner));
        }
        if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectAtomGroup>>() {
            return PyAtomGroup::from_core(py, core.select(&s.inner));
        }

        if selector.hasattr("is_match")? {
            struct PySelector<'a> {
                obj: &'a Bound<'a, PyAny>,
                err: std::cell::RefCell<Option<PyErr>>,
            }

            impl<'a> Selector for PySelector<'a> {
                fn is_match_group(&self, group: &CoreAtomGroup) -> bool {
                    if self.err.borrow().is_some() {
                        return false;
                    }
                    let py = self.obj.py();
                    let py_group = match PyAtomGroup::from_core(py, group.clone()) {
                        Ok(g) => g,
                        Err(e) => {
                            *self.err.borrow_mut() = Some(e);
                            return false;
                        }
                    };
                    match self.obj.call_method1("is_match", (py_group,)) {
                        Ok(res) => res.is_truthy().unwrap_or(false),
                        Err(e) => {
                            *self.err.borrow_mut() = Some(e);
                            false
                        }
                    }
                }

                fn is_match_atom(&self, atom: &proteindf_bridge::atom::Atom) -> bool {
                    if self.err.borrow().is_some() {
                        return false;
                    }
                    let py = self.obj.py();
                    let py_atom = match PyAtom::from_core(py, atom.clone()) {
                        Ok(a) => a,
                        Err(e) => {
                            *self.err.borrow_mut() = Some(e);
                            return false;
                        }
                    };
                    match self.obj.call_method1("is_match", (py_atom,)) {
                        Ok(res) => res.is_truthy().unwrap_or(false),
                        Err(e) => {
                            *self.err.borrow_mut() = Some(e);
                            false
                        }
                    }
                }
            }

            let py_sel = PySelector {
                obj: selector,
                err: std::cell::RefCell::new(None),
            };
            let result = core.select(&py_sel);
            if let Some(err) = py_sel.err.into_inner() {
                return Err(err);
            }
            return PyAtomGroup::from_core(py, result);
        }

        Err(PyTypeError::new_err(format!(
            "Unsupported selector type: {}",
            selector.get_type()
        )))
    }

    #[pyo3(signature = (atom1, atom2, order=1))]
    pub fn add_bond(&mut self, _py: Python<'_>, atom1: &PyAtom, atom2: &PyAtom, order: usize) {
        let b = BondRecord {
            atom1_path: atom1.path.clone(),
            atom2_path: atom2.path.clone(),
            order,
        };
        if !self.bonds.contains(&b) {
            self.bonds.push(b);
        }
    }

    pub fn bonds(&self) -> Vec<(String, String, usize)> {
        self.bonds
            .iter()
            .map(|b| (b.atom1_path.clone(), b.atom2_path.clone(), b.order))
            .collect()
    }

    pub fn get_bond_list(&self, py: Python<'_>) -> Vec<(String, String, usize)> {
        let mut core = self.to_core(py);
        core.get_bond_list()
            .into_iter()
            .map(|b| (b.atom1_path, b.atom2_path, b.order))
            .collect()
    }

    pub fn setup(&mut self, py: Python<'_>) -> PyResult<()> {
        let mut core = self.to_core(py);
        core.setup().map_err(to_py_err)?;
        self.update_from_core(py, core)?;
        Ok(())
    }

    pub fn setup_with_db(&mut self, py: Python<'_>, db: &PyCcdTemplateDb) -> PyResult<()> {
        let mut core = self.to_core(py);
        core.setup_with_db(&db.inner).map_err(to_py_err)?;
        self.update_from_core(py, core)?;
        Ok(())
    }

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

    pub fn apply_secondary_structure(&mut self, py: Python<'_>) -> PyResult<()> {
        let mut core = self.to_core(py);
        core_apply_secondary_structure(&mut core);
        self.update_from_core(py, core)?;
        Ok(())
    }

    pub fn shift_by(&mut self, py: Python<'_>, dir: &Bound<'_, PyAny>) -> PyResult<()> {
        let pos = extract_pos(dir)?;
        self.shift_by_position(py, pos);
        Ok(())
    }

    pub fn rotate(&mut self, py: Python<'_>, rotmat: &Bound<'_, PyAny>) -> PyResult<()> {
        let core_mat = crate::matrix::extract_rotation_matrix(rotmat)?;
        self.rotate_core(py, &core_mat)
    }

    pub fn sum_of_atomic_number(&self) -> f64 {
        Python::attach(|py| self.to_core(py).sum_of_atomic_number())
    }

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

    pub fn __str__(&self) -> String {
        Python::attach(|py| format!("{}", self.to_core(py)))
    }

    pub fn get_raw_data<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let core = self.to_core(py);
        let val = proteindf_bridge::brd::atomgroup_get_raw_data(&core);
        crate::value::value_to_py(py, &val)
    }

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
        let core = self.to_core(py);
        let val = proteindf_bridge::brd::atomgroup_get_raw_data(&core);
        let raw_py = crate::value::value_to_py(py, &val)?;
        let path = self.compute_path(py);

        let mut ss_list = Vec::new();
        collect_secondary_structures(&core, &mut ss_list);

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
                raw_py,
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
        let core = self.to_core(py);
        collect_csv_rows(&core, &[], &mut rows, py)?;

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

fn collect_csv_rows(
    ag: &CoreAtomGroup,
    parents: &[String],
    rows: &mut Vec<Vec<String>>,
    py: Python<'_>,
) -> PyResult<()> {
    for (key, subgrp) in ag.groups() {
        let mut new_parents = parents.to_vec();
        new_parents.push(key.clone());
        new_parents.push(subgrp.name.clone());
        collect_csv_rows(subgrp, &new_parents, rows, py)?;
    }
    for (_key, atom) in ag.atoms() {
        let mut row = parents.to_vec();
        let sym = atom.symbol().unwrap_or("X").to_string();
        row.push(sym);
        row.push(format_py_float(py, atom.xyz.x)?);
        row.push(format_py_float(py, atom.xyz.y)?);
        row.push(format_py_float(py, atom.xyz.z)?);
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

fn collect_secondary_structures(group: &CoreAtomGroup, out: &mut Vec<(String, String)>) {
    if let Some(ss) = group.secondary_structure() {
        out.push((group.path().to_string(), ss.as_char().to_string()));
    }
    for (_, child) in group.groups() {
        collect_secondary_structures(child, out);
    }
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
