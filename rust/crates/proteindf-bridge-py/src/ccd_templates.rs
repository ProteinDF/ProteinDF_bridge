// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::error::to_py_err;
use proteindf_bridge::ccd_templates::{
    CcdAtom as CoreCcdAtom, CcdBondTemplate as CoreCcdBondTemplate,
    CcdTemplateDb as CoreCcdTemplateDb,
};
use proteindf_bridge::error::BridgeError;
use proteindf_bridge::format::mmcif::SimpleMmcif as CoreSimpleMmcif;
use pyo3::prelude::*;
use std::path::Path;

/// A single atom of a CCD component, with element and idealized geometry.
#[pyclass(name = "CcdAtom", module = "proteindf_bridge_rs", skip_from_py_object)]
#[derive(Clone)]
pub struct PyCcdAtom {
    pub(crate) inner: CoreCcdAtom,
}

#[pymethods]
impl PyCcdAtom {
    #[getter]
    pub fn name(&self) -> String {
        self.inner.name.clone()
    }

    #[getter]
    pub fn element(&self) -> String {
        self.inner.element.clone()
    }

    #[getter]
    pub fn ideal_xyz(&self) -> Option<(f64, f64, f64)> {
        self.inner.ideal_xyz
    }

    #[getter]
    pub fn is_hydrogen(&self) -> bool {
        self.inner.is_hydrogen()
    }

    pub fn __repr__(&self) -> String {
        format!(
            "CcdAtom(name='{}', element='{}', ideal_xyz={:?})",
            self.inner.name, self.inner.element, self.inner.ideal_xyz
        )
    }
}

/// A bond template for a chemical component in the CCD.
#[pyclass(
    name = "CcdBondTemplate",
    module = "proteindf_bridge_rs",
    skip_from_py_object
)]
#[derive(Clone)]
pub struct PyCcdBondTemplate {
    pub(crate) inner: CoreCcdBondTemplate,
}

#[pymethods]
impl PyCcdBondTemplate {
    #[getter]
    pub fn comp_id(&self) -> String {
        self.inner.comp_id.clone()
    }

    #[getter]
    pub fn atoms(&self) -> Vec<PyCcdAtom> {
        self.inner
            .atoms
            .iter()
            .cloned()
            .map(|a| PyCcdAtom { inner: a })
            .collect()
    }

    #[getter]
    pub fn bonds(&self) -> Vec<(String, String, usize)> {
        self.inner.bonds.clone()
    }

    pub fn get_atom(&self, name: &str) -> Option<PyCcdAtom> {
        self.inner
            .get_atom(name)
            .cloned()
            .map(|a| PyCcdAtom { inner: a })
    }

    pub fn __len__(&self) -> usize {
        self.inner.atoms.len()
    }

    pub fn __repr__(&self) -> String {
        format!(
            "CcdBondTemplate(comp_id='{}', atoms={}, bonds={})",
            self.inner.comp_id,
            self.inner.atoms.len(),
            self.inner.bonds.len()
        )
    }
}

/// In-memory lookup database of CCD bond templates.
#[pyclass(name = "CcdTemplateDb", module = "proteindf_bridge_rs", from_py_object)]
#[derive(Clone, Default)]
pub struct PyCcdTemplateDb {
    pub(crate) inner: CoreCcdTemplateDb,
}

#[pymethods]
impl PyCcdTemplateDb {
    #[new]
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns a new database initialized with standard 29 embedded residue templates.
    #[classmethod]
    pub fn builtin(_cls: &Bound<'_, pyo3::types::PyType>) -> Self {
        Self {
            inner: CoreCcdTemplateDb::default(),
        }
    }

    /// Returns a new, empty CCD template database.
    #[classmethod]
    pub fn empty(_cls: &Bound<'_, pyo3::types::PyType>) -> Self {
        Self {
            inner: CoreCcdTemplateDb::new(),
        }
    }

    /// Looks up a component template by its standard identifier (e.g. "ALA", "ARG").
    pub fn lookup(&self, comp_id: &str) -> Option<PyCcdBondTemplate> {
        self.inner
            .lookup(comp_id)
            .cloned()
            .map(|t| PyCcdBondTemplate { inner: t })
    }

    /// Returns the number of registered component templates.
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// Returns whether the database is empty.
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    pub fn __len__(&self) -> usize {
        self.inner.len()
    }

    pub fn __contains__(&self, comp_id: &str) -> bool {
        self.inner.lookup(comp_id).is_some()
    }

    /// Inserts a template into the database.
    pub fn insert(&mut self, template: &PyCcdBondTemplate) {
        self.inner.insert(template.inner.clone());
    }

    /// Merges another CcdTemplateDb into this one.
    ///
    /// Conflicting entries in `other` take precedence and overwrite existing entries.
    pub fn merge(&mut self, other: &PyCcdTemplateDb) {
        self.inner.merge(&other.inner);
    }

    /// Loads component templates from a user-supplied mmCIF file (which may contain multiple data blocks)
    /// and inserts them into this database.
    ///
    /// # Errors
    /// Returns an error if:
    /// - The file cannot be read or parsed.
    /// - No CCD data blocks are found in the file.
    /// - Any data block fails to parse (e.g. conflicting duplicate atoms, missing atom entries,
    ///   or macromolecular structure data blocks containing `_atom_site`).
    /// In case of error, this database is left completely unmodified.
    pub fn add_from_file(&mut self, file_path: &str) -> PyResult<usize> {
        let mut cif = CoreSimpleMmcif::new();
        cif.load(Path::new(file_path)).map_err(to_py_err)?;

        let names = cif.get_molecule_names();
        if names.is_empty() {
            return Err(to_py_err(BridgeError::input_error(
                "ccd_templates",
                format!("No data blocks found in file: {file_path}"),
            )));
        }

        let mut parsed_templates = Vec::new();
        let mut errors = Vec::new();

        for name in &names {
            match cif.get_data_block(name) {
                Some(block) => match CoreCcdBondTemplate::from_mmcif_block(block, name) {
                    Ok(template) => {
                        parsed_templates.push(template);
                    }
                    Err(e) => {
                        errors.push(format!("block '{name}': {e}"));
                    }
                },
                None => {
                    errors.push(format!("block '{name}': block not found in mmCIF data"));
                }
            }
        }

        if !errors.is_empty() {
            return Err(to_py_err(BridgeError::input_error(
                "ccd_templates",
                format!(
                    "Failed to parse CCD component blocks from {file_path}: {}",
                    errors.join("; ")
                ),
            )));
        }

        if parsed_templates.is_empty() {
            return Err(to_py_err(BridgeError::input_error(
                "ccd_templates",
                format!("No valid CCD component blocks found in file: {file_path}"),
            )));
        }

        let count = parsed_templates.len();
        for template in parsed_templates {
            self.inner.insert(template);
        }

        Ok(count)
    }

    pub fn __repr__(&self) -> String {
        format!("CcdTemplateDb(components={})", self.inner.len())
    }
}
