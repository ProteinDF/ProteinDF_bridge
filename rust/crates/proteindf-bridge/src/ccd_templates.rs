// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

//! wwPDB Chemical Component Dictionary (CCD) bond template database.
//!
//! Provides canonical atom names/elements/idealized geometry and bond topology
//! (with bond orders 1, 2, 3, 4) for standard amino acids, nucleic acids, and
//! water, embedded at compile time via MessagePack binary for zero-filesystem
//! portable access. The idealized geometry (including explicit hydrogens) is the
//! data foundation for CCD-reference-based hydrogen addition (`RUST_PORT_SPEC.md`
//! §3.16); this module itself only provides the template data and bond-order
//! resolution (`AtomGroup::apply_ccd_bond_templates`), not hydrogenation itself.
//!
//! # Runtime Extension with User-Supplied CCD Data
//! While the standard 29 residue templates are embedded at compile time via [`CcdTemplateDb::global`],
//! the database can be extended at runtime with user-supplied CCD files (e.g. `components.cif`
//! or individual ligand CIF files) using [`CcdBondTemplate::from_mmcif_block`] and [`CcdTemplateDb::insert`]
//! or [`CcdTemplateDb::merge`].
//!
//! ## Example
//! ```no_run
//! use std::path::Path;
//! use proteindf_bridge::ccd_templates::{CcdBondTemplate, CcdTemplateDb};
//! use proteindf_bridge::format::mmcif::SimpleMmcif;
//!
//! // 1. Load user-supplied CCD file via SimpleMmcif
//! let mut cif = SimpleMmcif::new();
//! cif.load(Path::new("path/to/custom_ligand.cif")).unwrap();
//!
//! // 2. Convert a CCD data block into CcdBondTemplate
//! let block = cif.get_data_block("comp_LIG").unwrap();
//! let template = CcdBondTemplate::from_mmcif_block(block, "LIG").unwrap();
//!
//! // 3. Insert into a mutable template DB (cloned from default or freshly created)
//! let mut db = CcdTemplateDb::default();
//! db.insert(template);
//!
//! // 4. Apply to AtomGroup
//! // atom_group.apply_ccd_bond_templates(&db);
//! ```

use std::collections::HashMap;
use std::sync::OnceLock;

use serde::{Deserialize, Serialize};

use crate::error::{BridgeError, Result};
use crate::format::mmcif::MmcifDataBlock;

/// Raw embedded MessagePack bytes for standard CCD bond templates.
pub const CCD_BOND_TEMPLATES_MSGPACK: &[u8] = include_bytes!("data/ccd_bond_templates.msgpack");

/// A single atom of a CCD component, with element and idealized geometry.
///
/// `ideal_xyz` is `None` when neither an idealized (`pdbx_model_Cartn_*_ideal`) nor a
/// model (`model_Cartn_*`) coordinate could be resolved for this atom (e.g. a
/// minimal user-supplied CCD file that only lists atom names/elements). Bond-order
/// resolution (`AtomGroup::apply_ccd_bond_templates`) does not need coordinates, so a
/// missing `ideal_xyz` does not prevent a template from being usable for that purpose.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CcdAtom {
    /// Canonical atom name (e.g. "CA", "HB2").
    pub name: String,
    /// Element symbol (e.g. "C", "N", "H"). Deuterium ("D") is normalized to "H".
    pub element: String,
    /// Idealized (or, failing that, model) Cartesian coordinates in Angstrom.
    pub ideal_xyz: Option<(f64, f64, f64)>,
}

impl CcdAtom {
    /// Returns whether this atom's element is hydrogen.
    pub fn is_hydrogen(&self) -> bool {
        self.element.eq_ignore_ascii_case("H")
    }
}

/// A bond template for a chemical component in the CCD.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CcdBondTemplate {
    /// 3-letter component identifier (e.g. "ALA", "ARG", "DA", "HOH").
    pub comp_id: String,
    /// List of canonical atoms (name, element, idealized geometry).
    pub atoms: Vec<CcdAtom>,
    /// List of intra-component bonds as (atom_id_1, atom_id_2, bond_order).
    pub bonds: Vec<(String, String, usize)>,
}

impl CcdBondTemplate {
    /// Looks up an atom by its canonical name.
    pub fn get_atom(&self, name: &str) -> Option<&CcdAtom> {
        self.atoms.iter().find(|a| a.name == name)
    }

    /// Extracts an idealized (falling back to model) coordinate axis from an mmCIF
    /// `_chem_comp_atom` row, mirroring `format::mmcif`'s coordinate-resolution
    /// convention (idealized coordinates take priority over model coordinates).
    fn get_coordinate(axis: &str, dict: &indexmap::IndexMap<String, String>) -> Option<f64> {
        let ideal_key = format!("_chem_comp_atom.pdbx_model_Cartn_{axis}_ideal");
        let model_key = format!("_chem_comp_atom.model_Cartn_{axis}");

        if let Some(val) = dict.get(&ideal_key) {
            if let Ok(v) = val.parse::<f64>() {
                return Some(v);
            }
        }
        if let Some(val) = dict.get(&model_key) {
            if let Ok(v) = val.parse::<f64>() {
                return Some(v);
            }
        }
        None
    }

    /// Builds a `CcdAtom` from an mmCIF `_chem_comp_atom` row, if it contains an `atom_id`.
    fn atom_from_row(row: &indexmap::IndexMap<String, String>) -> Option<CcdAtom> {
        let name = row.get("_chem_comp_atom.atom_id")?.clone();
        let mut element = row
            .get("_chem_comp_atom.type_symbol")
            .cloned()
            .unwrap_or_else(|| "X".to_string());
        if element == "D" {
            element = "H".to_string();
        }
        let ideal_xyz = match (
            Self::get_coordinate("x", row),
            Self::get_coordinate("y", row),
            Self::get_coordinate("z", row),
        ) {
            (Some(x), Some(y), Some(z)) => Some((x, y, z)),
            _ => None,
        };
        Some(CcdAtom {
            name,
            element,
            ideal_xyz,
        })
    }

    /// Builds a `CcdBondTemplate` from an mmCIF CCD data block.
    ///
    /// Extracts canonical atom names from `_chem_comp_atom.atom_id` and intra-component
    /// bonds and bond orders from `_chem_comp_bond`.
    ///
    /// # Errors
    /// Returns an error if:
    /// - The block contains `_atom_site` records (full macromolecular structure data rather than a CCD chemical component).
    /// - No `_chem_comp_atom` entries are found in the block.
    pub fn from_mmcif_block(block: &MmcifDataBlock, comp_id: &str) -> Result<Self> {
        if block.has_atom_site() {
            return Err(BridgeError::input_error(
                comp_id,
                "Cannot build CcdBondTemplate from an mmCIF data block containing _atom_site records (macromolecular structure)",
            ));
        }

        let actual_comp_id = block
            .key_values
            .get("_chem_comp.id")
            .cloned()
            .unwrap_or_else(|| {
                let trimmed = comp_id.trim();
                if trimmed.is_empty() {
                    "UNKNOWN".to_string()
                } else {
                    trimmed.to_string()
                }
            });

        // 1. Extract canonical atoms (name, element, idealized geometry; preserving appearance order)
        let mut atoms: Vec<CcdAtom> = Vec::new();
        if let Some(atom) = Self::atom_from_row(&block.key_values) {
            if !atoms.iter().any(|a| a.name == atom.name) {
                atoms.push(atom);
            }
        }
        for table in &block.tables {
            for row in table {
                if let Some(atom) = Self::atom_from_row(row) {
                    if !atoms.iter().any(|a| a.name == atom.name) {
                        atoms.push(atom);
                    }
                }
            }
        }

        if atoms.is_empty() {
            return Err(BridgeError::input_error(
                comp_id,
                format!("No _chem_comp_atom entries found in mmCIF data block for '{comp_id}'"),
            ));
        }

        // 2. Extract bonds and bond orders (reusing parse_chem_comp_bond_order)
        let mut bonds = Vec::new();
        if let (Some(a1), Some(a2)) = (
            block.key_values.get("_chem_comp_bond.atom_id_1"),
            block.key_values.get("_chem_comp_bond.atom_id_2"),
        ) {
            let order_str = block
                .key_values
                .get("_chem_comp_bond.value_order")
                .map(|s| s.as_str())
                .unwrap_or("");
            let order = crate::format::mmcif::parse_chem_comp_bond_order(order_str);
            if order > 0 {
                bonds.push((a1.clone(), a2.clone(), order));
            }
        }

        for table in &block.tables {
            for row in table {
                if row.contains_key("_chem_comp_bond.comp_id")
                    || row.contains_key("_chem_comp_bond.atom_id_1")
                {
                    if let (Some(a1), Some(a2)) = (
                        row.get("_chem_comp_bond.atom_id_1"),
                        row.get("_chem_comp_bond.atom_id_2"),
                    ) {
                        let order_str = row
                            .get("_chem_comp_bond.value_order")
                            .map(|s| s.as_str())
                            .unwrap_or("");
                        let order = crate::format::mmcif::parse_chem_comp_bond_order(order_str);
                        if order > 0 {
                            bonds.push((a1.clone(), a2.clone(), order));
                        }
                    }
                }
            }
        }

        Ok(Self {
            comp_id: actual_comp_id,
            atoms,
            bonds,
        })
    }
}

/// In-memory lookup database of CCD bond templates.
#[derive(Debug, Clone)]
pub struct CcdTemplateDb {
    templates: HashMap<String, CcdBondTemplate>,
}

impl CcdTemplateDb {
    /// Creates a new, empty CCD template database.
    pub fn new() -> Self {
        Self {
            templates: HashMap::new(),
        }
    }

    /// Returns the global singleton instance of the standard CCD template database.
    pub fn global() -> &'static Self {
        static INSTANCE: OnceLock<CcdTemplateDb> = OnceLock::new();
        INSTANCE.get_or_init(|| {
            Self::load_from_bytes(CCD_BOND_TEMPLATES_MSGPACK)
                .expect("failed to deserialize embedded CCD bond templates")
        })
    }

    /// Deserializes a CCD template database from raw MessagePack bytes.
    pub fn load_from_bytes(bytes: &[u8]) -> Result<Self> {
        let templates: HashMap<String, CcdBondTemplate> =
            rmp_serde::from_slice(bytes).map_err(|e| {
                BridgeError::input_error("ccd_templates", format!("MessagePack decode error: {e}"))
            })?;
        Ok(Self { templates })
    }

    /// Inserts a template into the database.
    ///
    /// If a template for the same `comp_id` already exists, it is replaced and the old
    /// template is returned.
    pub fn insert(&mut self, template: CcdBondTemplate) -> Option<CcdBondTemplate> {
        self.templates.insert(template.comp_id.clone(), template)
    }

    /// Merges another `CcdTemplateDb` into this one.
    ///
    /// # Conflict Resolution
    /// If a component with the same `comp_id` exists in both databases, the entry from
    /// `other` takes precedence and overwrites the existing entry ("last-write-wins").
    pub fn merge(&mut self, other: &CcdTemplateDb) {
        for (comp_id, template) in &other.templates {
            self.templates.insert(comp_id.clone(), template.clone());
        }
    }

    /// Looks up a component template by its standard identifier (e.g. "ALA", "ARG").
    pub fn lookup(&self, comp_id: &str) -> Option<&CcdBondTemplate> {
        self.templates.get(comp_id)
    }

    /// Returns the number of registered component templates.
    pub fn len(&self) -> usize {
        self.templates.len()
    }

    /// Returns whether the database is empty.
    pub fn is_empty(&self) -> bool {
        self.templates.is_empty()
    }
}

impl Default for CcdTemplateDb {
    fn default() -> Self {
        Self::global().clone()
    }
}
