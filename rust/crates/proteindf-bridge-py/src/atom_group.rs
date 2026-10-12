// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom::PyAtom;
use crate::ccd_templates::PyCcdTemplateDb;
use crate::error::to_py_err;
use crate::hydrogenation::PyOverallHydrogenationReport;
use crate::matrix::PyMatrix;
use crate::position::PyPosition;
use crate::schema::PySchemaViolation;
use proteindf_bridge::atom_group::{AtomGroup as CoreAtomGroup, Selector};
use proteindf_bridge::ccd_templates::CcdTemplateDb as CoreCcdTemplateDb;
use proteindf_bridge::position::Position;
use proteindf_bridge::secondary_structure::{
    apply_secondary_structure as core_apply_secondary_structure, SsCode,
};
use pyo3::exceptions::{PyKeyError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use std::collections::HashMap;
use std::str::FromStr;

#[pyclass(name = "AtomGroup", module = "proteindf_bridge.rs", from_py_object)]
#[derive(Clone, Default)]
pub struct PyAtomGroup {
    pub(crate) inner: CoreAtomGroup,
}

impl PyAtomGroup {
    pub fn from_core(group: CoreAtomGroup) -> Self {
        Self { inner: group }
    }
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
    #[new]
    #[pyo3(signature = (*args, **kwargs))]
    pub fn new(
        args: &Bound<'_, pyo3::types::PyTuple>,
        kwargs: Option<&Bound<'_, pyo3::types::PyDict>>,
    ) -> PyResult<Self> {
        let mut group = CoreAtomGroup::new();

        if args.len() == 1 {
            let first = args.get_item(0)?;
            if let Ok(other) = first.extract::<PyRef<PyAtomGroup>>() {
                group = other.inner.clone();
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
        }

        Ok(Self { inner: group })
    }

    #[getter]
    pub fn name(&self) -> String {
        self.inner.name.clone()
    }

    #[setter]
    pub fn set_name(&mut self, val: String) {
        self.inner.name = val;
    }

    #[getter]
    pub fn path(&self) -> String {
        self.inner.path().to_string()
    }

    #[setter]
    pub fn set_path(&mut self, val: String) {
        self.inner.set_path(val);
    }

    /// Returns the path depth of this group in the hierarchy.
    pub fn path_depth(&self) -> usize {
        self.inner.path_depth()
    }

    /// Returns whether this group is at the model level (path depth 1).
    pub fn is_model_level(&self) -> bool {
        self.inner.is_model_level()
    }

    /// Returns whether this group is at the chain level (path depth 2).
    pub fn is_chain_level(&self) -> bool {
        self.inner.is_chain_level()
    }

    /// Returns whether this group is at the residue level (path depth 3).
    pub fn is_residue_level(&self) -> bool {
        self.inner.is_residue_level()
    }

    /// Validates this AtomGroup hierarchy against the standard protein schema (/model_N/chain_id/res_key/atom_key).
    pub fn validate_schema(&self) -> Vec<PySchemaViolation> {
        self.inner
            .validate_schema()
            .into_iter()
            .map(PySchemaViolation::from)
            .collect()
    }

    /// Returns the secondary structure code ('H', 'E', '-', or None) for this group.
    #[getter]
    pub fn secondary_structure(&self) -> Option<String> {
        self.inner.secondary_structure().map(|ss| ss.to_string())
    }

    /// Sets the secondary structure code ('H', 'E', '-', or None) for this group.
    #[setter]
    pub fn set_secondary_structure(&mut self, val: Option<&str>) -> PyResult<()> {
        match val {
            None => {
                self.inner.set_secondary_structure(None);
                Ok(())
            }
            Some("H") => {
                self.inner.set_secondary_structure(Some(SsCode::Helix));
                Ok(())
            }
            Some("E") => {
                self.inner.set_secondary_structure(Some(SsCode::Strand));
                Ok(())
            }
            Some("-") => {
                self.inner.set_secondary_structure(Some(SsCode::Loop));
                Ok(())
            }
            Some(other) => Err(PyValueError::new_err(format!(
                "Invalid secondary structure code '{other}', expected 'H', 'E', '-', or None"
            ))),
        }
    }

    pub fn get_number_of_atoms(&self) -> usize {
        self.inner.get_number_of_atoms()
    }

    pub fn get_number_of_all_atoms(&self) -> usize {
        self.inner.get_number_of_all_atoms()
    }

    pub fn get_number_of_groups(&self) -> usize {
        self.inner.get_number_of_groups()
    }

    pub fn get_number_of_bonds(&self) -> usize {
        self.inner.get_number_of_bonds()
    }

    pub fn get_atom(&self, key: &Bound<'_, PyAny>) -> PyResult<Option<PyAtom>> {
        let k = extract_key(key)?;
        Ok(self
            .inner
            .get_atom(&k)
            .map(|a| PyAtom::from_core(a.clone())))
    }

    pub fn set_atom(&mut self, key: &Bound<'_, PyAny>, atom: &PyAtom) -> PyResult<()> {
        let k = extract_key(key)?;
        self.inner.set_atom(&k, atom.inner.clone());
        Ok(())
    }

    pub fn has_atom(&self, key: &Bound<'_, PyAny>) -> PyResult<bool> {
        let k = extract_key(key)?;
        Ok(self.inner.has_atom(&k))
    }

    pub fn del_atom(&mut self, key: &Bound<'_, PyAny>) -> PyResult<Option<PyAtom>> {
        let k = extract_key(key)?;
        Ok(self.inner.remove_atom(&k).map(PyAtom::from_core))
    }

    pub fn remove_atom(&mut self, key: &Bound<'_, PyAny>) -> PyResult<Option<PyAtom>> {
        self.del_atom(key)
    }

    pub fn get_group(&self, key: &Bound<'_, PyAny>) -> PyResult<Option<PyAtomGroup>> {
        let k = extract_key(key)?;
        Ok(self
            .inner
            .get_group(&k)
            .map(|g| PyAtomGroup::from_core(g.clone())))
    }

    pub fn set_group(&mut self, key: &Bound<'_, PyAny>, group: &PyAtomGroup) -> PyResult<()> {
        let k = extract_key(key)?;
        self.inner.set_group(&k, group.inner.clone());
        Ok(())
    }

    pub fn has_group(&self, key: &Bound<'_, PyAny>) -> PyResult<bool> {
        let k = extract_key(key)?;
        Ok(self.inner.has_group(&k))
    }

    pub fn has_groupkey(&self, key: &Bound<'_, PyAny>) -> PyResult<bool> {
        let k = extract_key(key)?;
        Ok(self.inner.has_groupkey(&k))
    }

    pub fn del_group(&mut self, key: &Bound<'_, PyAny>) -> PyResult<Option<PyAtomGroup>> {
        let k = extract_key(key)?;
        Ok(self.inner.remove_group(&k).map(PyAtomGroup::from_core))
    }

    pub fn atoms(&self) -> Vec<(String, PyAtom)> {
        self.inner
            .atoms()
            .map(|(k, a)| (k.clone(), PyAtom::from_core(a.clone())))
            .collect()
    }

    pub fn groups(&self) -> Vec<(String, PyAtomGroup)> {
        self.inner
            .groups()
            .map(|(k, g)| (k.clone(), PyAtomGroup::from_core(g.clone())))
            .collect()
    }

    pub fn get_group_list(&self) -> Vec<String> {
        self.inner.get_group_list()
    }

    pub fn get_atom_list(&self) -> Vec<PyAtom> {
        self.inner
            .get_atom_list()
            .into_iter()
            .map(PyAtom::from_core)
            .collect()
    }

    pub fn get_path_list(&self) -> Vec<String> {
        self.inner.get_path_list()
    }

    pub fn get_atom_kinds_count(&self) -> HashMap<String, usize> {
        self.inner.get_atom_kinds_count()
    }

    pub fn get_formula(&self) -> String {
        self.inner.get_formula()
    }

    pub fn formula(&self) -> String {
        self.inner.get_formula()
    }

    #[getter]
    pub fn weight(&self) -> f64 {
        self.inner
            .get_atom_list()
            .iter()
            .map(|a| a.weight().unwrap_or(0.0))
            .sum()
    }

    pub fn center(&self) -> PyPosition {
        PyPosition::from_core(self.inner.center())
    }

    pub fn r#box(&self) -> (PyPosition, PyPosition) {
        let (bmin, bmax) = self.inner.r#box();
        (PyPosition::from_core(bmin), PyPosition::from_core(bmax))
    }

    pub fn get_box(&self) -> (PyPosition, PyPosition) {
        self.r#box()
    }

    pub fn merge(&mut self, other: &PyAtomGroup) {
        self.inner.merge(&other.inner);
    }

    pub fn select(&self, selector: &Bound<'_, PyAny>) -> PyResult<PyAtomGroup> {
        if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectSymbol>>() {
            return Ok(PyAtomGroup::from_core(self.inner.select(&s.inner)));
        }
        if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectName>>() {
            return Ok(PyAtomGroup::from_core(self.inner.select(&s.inner)));
        }
        if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectPathSimple>>() {
            return Ok(PyAtomGroup::from_core(self.inner.select(&s.inner)));
        }
        if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectPathWildcard>>() {
            return Ok(PyAtomGroup::from_core(self.inner.select(&s.inner)));
        }
        if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectPathRegex>>() {
            return Ok(PyAtomGroup::from_core(self.inner.select(&s.inner)));
        }
        if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectPath>>() {
            return Ok(PyAtomGroup::from_core(self.inner.select(&s.inner)));
        }
        if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectRange>>() {
            return Ok(PyAtomGroup::from_core(self.inner.select(&s.inner)));
        }
        if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectAtom>>() {
            return Ok(PyAtomGroup::from_core(self.inner.select(&s.inner)));
        }
        if let Ok(s) = selector.extract::<PyRef<crate::selector::PySelectAtomGroup>>() {
            return Ok(PyAtomGroup::from_core(self.inner.select(&s.inner)));
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
                    let py_group = PyAtomGroup::from_core(group.clone());
                    match py_group.into_pyobject(py) {
                        Ok(py_obj) => match self.obj.call_method1("is_match", (py_obj,)) {
                            Ok(res) => match res.extract::<bool>() {
                                Ok(b) => b,
                                Err(e) => {
                                    *self.err.borrow_mut() = Some(e);
                                    false
                                }
                            },
                            Err(e) => {
                                *self.err.borrow_mut() = Some(e);
                                false
                            }
                        },
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
                    let py_atom = PyAtom::from_core(atom.clone());
                    match py_atom.into_pyobject(py) {
                        Ok(py_obj) => match self.obj.call_method1("is_match", (py_obj,)) {
                            Ok(res) => match res.extract::<bool>() {
                                Ok(b) => b,
                                Err(e) => {
                                    *self.err.borrow_mut() = Some(e);
                                    false
                                }
                            },
                            Err(e) => {
                                *self.err.borrow_mut() = Some(e);
                                false
                            }
                        },
                        Err(e) => {
                            *self.err.borrow_mut() = Some(e);
                            false
                        }
                    }
                }
            }

            let adapter = PySelector {
                obj: selector,
                err: std::cell::RefCell::new(None),
            };
            let result = self.inner.select(&adapter);
            if let Some(err) = adapter.err.into_inner() {
                return Err(err);
            }
            return Ok(PyAtomGroup::from_core(result));
        }

        Err(PyTypeError::new_err(format!(
            "Expected Select, got {}",
            selector.get_type()
        )))
    }

    #[pyo3(signature = (atom1, atom2, order=1))]
    pub fn add_bond(&mut self, atom1: &PyAtom, atom2: &PyAtom, order: usize) {
        self.inner.add_bond(&atom1.inner, &atom2.inner, order);
    }

    pub fn bonds(&self) -> Vec<(String, String, usize)> {
        self.inner
            .bonds()
            .iter()
            .map(|b| (b.atom1_path.clone(), b.atom2_path.clone(), b.order))
            .collect()
    }

    pub fn get_bond_list(&mut self) -> Vec<(String, String, usize)> {
        self.inner
            .get_bond_list()
            .into_iter()
            .map(|b| (b.atom1_path, b.atom2_path, b.order))
            .collect()
    }

    /// Sets up chemical bonds using embedded CCD templates and covalent radius heuristics.
    ///
    /// Modifies this AtomGroup in place. Note that calling this on a subtree copy does not
    /// affect the original tree.
    pub fn setup(&mut self) -> PyResult<()> {
        self.inner.setup().map_err(to_py_err)
    }

    /// Sets up chemical bonds using the provided CCD template database and covalent radius heuristics.
    ///
    /// Modifies this AtomGroup in place. Note that calling this on a subtree copy does not
    /// affect the original tree.
    pub fn setup_with_db(&mut self, db: &PyCcdTemplateDb) -> PyResult<()> {
        self.inner.setup_with_db(&db.inner).map_err(to_py_err)
    }

    /// Adds missing hydrogens to this structure using standard polymer connectivity
    /// and CCD idealized templates.
    ///
    /// Modifies this AtomGroup in place. Note that calling this on a subtree copy does not
    /// affect the original tree.
    ///
    /// # Arguments
    /// * `db` - Optional CCD template database (`CcdTemplateDb`). If None, uses the standard embedded CCD database.
    ///
    /// # Returns
    /// An [`OverallHydrogenationReport`] summarizing added/removed hydrogens, modified residues,
    /// skipped components (e.g. water), and step errors.
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
        let report = self
            .inner
            .add_missing_hydrogens(core_db)
            .map_err(to_py_err)?;
        PyOverallHydrogenationReport::new(py, report)
    }

    /// Applies 3-state secondary structure assignments ('H', 'E', '-') directly to the residue groups of this chain in-place.
    ///
    /// Note: This modifies the AtomGroup in-place. If called on a copy or a sub-tree copy (such as one returned by indexing or get_group()), the original root tree will remain unchanged.
    pub fn apply_secondary_structure(&mut self) {
        core_apply_secondary_structure(&mut self.inner);
    }

    pub fn shift_by(&mut self, dir: &Bound<'_, PyAny>) -> PyResult<()> {
        let d = extract_pos(dir)?;
        self.inner.shift_by(d);
        Ok(())
    }

    pub fn rotate(&mut self, rotmat: &PyMatrix) -> PyResult<()> {
        self.inner.rotate(&rotmat.inner).map_err(to_py_err)
    }

    pub fn sum_of_atomic_number(&self) -> f64 {
        self.inner.sum_of_atomic_number()
    }

    /// Rebuilds the data structure by referring to the structure of reference.
    ///
    /// Matches atoms within `range` distance, copies charges from reference,
    /// and assigns any remaining atoms to model_1/Z/1 (res_name UNK).
    #[pyo3(signature = (reference, range=1.0e-5))]
    pub fn restructure(&self, reference: &PyAtomGroup, range: f64) -> PyResult<PyAtomGroup> {
        let target_selector = proteindf_bridge::selector::SelectAtomGroup::new(&self.inner, range);
        let mut restructured = reference.inner.select(&target_selector);

        copy_attributes(&mut restructured, &self.inner);

        let rest_of_target = get_rest_of_frame_molecule(&self.inner, &restructured);
        assign_rest_molecule(&rest_of_target, &mut restructured, "model_1", "Z", "UNK")?;

        Ok(PyAtomGroup::from_core(restructured))
    }

    pub fn __getitem__<'py>(
        &self,
        py: Python<'py>,
        key: &Bound<'_, PyAny>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let key_str = extract_key(key)?;
        let key_ref = key_str.as_str();
        // 1. Child group by key
        if self.inner.has_groupkey(key_ref) {
            if let Some(grp) = self.inner.get_group(key_ref) {
                let py_grp = PyAtomGroup::from_core(grp.clone());
                return py_grp.into_pyobject(py).map(|b| b.into_any());
            }
        }
        // 2. Child atom by key
        if self.inner.has_atomkey(key_ref) {
            if let Some(atom) = self.inner.get_atom(key_ref) {
                let py_atom = PyAtom::from_core(atom.clone());
                return py_atom.into_pyobject(py).map(|b| b.into_any());
            }
        }
        // 3. Search child groups by name
        for (_, grp) in self.inner.groups() {
            if grp.name == key_ref {
                let py_grp = PyAtomGroup::from_core(grp.clone());
                return py_grp.into_pyobject(py).map(|b| b.into_any());
            }
        }
        // 4. Search child atoms by name
        for (_, atom) in self.inner.atoms() {
            if atom.name == key_ref {
                let py_atom = PyAtom::from_core(atom.clone());
                return py_atom.into_pyobject(py).map(|b| b.into_any());
            }
        }
        Err(PyKeyError::new_err(key_str))
    }

    pub fn __setitem__(
        &mut self,
        key: &Bound<'_, PyAny>,
        value: &Bound<'_, PyAny>,
    ) -> PyResult<()> {
        let key_str = extract_key(key)?;
        if let Ok(grp) = value.extract::<PyRef<PyAtomGroup>>() {
            self.inner.set_group(&key_str, grp.inner.clone());
            Ok(())
        } else if let Ok(atom) = value.extract::<PyRef<PyAtom>>() {
            self.inner.set_atom(&key_str, atom.inner.clone());
            Ok(())
        } else {
            Err(pyo3::exceptions::PyValueError::new_err(format!(
                "Expected AtomGroup or Atom, got {}",
                value.get_type()
            )))
        }
    }

    pub fn __and__(&self, other: &PyAtomGroup) -> Self {
        Self {
            inner: &self.inner & &other.inner,
        }
    }

    pub fn __or__(&self, other: &PyAtomGroup) -> Self {
        Self {
            inner: &self.inner | &other.inner,
        }
    }

    pub fn __xor__(&self, other: &PyAtomGroup) -> Self {
        Self {
            inner: &self.inner ^ &other.inner,
        }
    }

    pub fn __len__(&self) -> usize {
        self.inner.get_number_of_atoms()
    }

    pub fn __repr__(&self) -> String {
        format!(
            "AtomGroup(name='{}', path='{}', atoms={}, groups={})",
            self.inner.name,
            self.inner.path(),
            self.inner.get_number_of_atoms(),
            self.inner.get_number_of_groups()
        )
    }

    pub fn __str__(&self) -> String {
        format!("{}", self.inner)
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
        collect_csv_rows(&self.inner, &[], &mut rows, py)?;

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
    if s.contains(',') || s.contains('"') || s.contains('\n') || s.contains('\r') {
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
