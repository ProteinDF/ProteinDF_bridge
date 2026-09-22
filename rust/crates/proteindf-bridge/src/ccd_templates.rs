// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

//! wwPDB Chemical Component Dictionary (CCD) bond template database.
//!
//! Provides canonical atom names and bond topology (with bond orders 1, 2, 3, 4)
//! for standard amino acids, nucleic acids, and water, embedded at compile time
//! via MessagePack binary for zero-filesystem portable access.

use std::collections::HashMap;
use std::sync::OnceLock;

use serde::{Deserialize, Serialize};

use crate::error::{BridgeError, Result};

/// Raw embedded MessagePack bytes for standard CCD bond templates.
pub const CCD_BOND_TEMPLATES_MSGPACK: &[u8] = include_bytes!("data/ccd_bond_templates.msgpack");

/// A bond template for a chemical component in the CCD.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CcdBondTemplate {
    /// 3-letter component identifier (e.g. "ALA", "ARG", "DA", "HOH").
    pub comp_id: String,
    /// List of canonical atom names.
    pub atoms: Vec<String>,
    /// List of intra-component bonds as (atom_id_1, atom_id_2, bond_order).
    pub bonds: Vec<(String, String, usize)>,
}

/// In-memory lookup database of CCD bond templates.
#[derive(Debug, Clone)]
pub struct CcdTemplateDb {
    templates: HashMap<String, CcdBondTemplate>,
}

impl CcdTemplateDb {
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
