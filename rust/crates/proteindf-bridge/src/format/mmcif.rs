// Copyright (C) 2019 The ProteinDF development team.
// see also AUTHORS and README if provided.
//
// This file is a part of the ProteinDF software package.
//
// The ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

use std::fs;
use std::path::Path;
use std::str::FromStr;

use indexmap::IndexMap;

use crate::atom::Atom;
use crate::atom_group::AtomGroup;
use crate::error::{BridgeError, Result};
use crate::periodic_table::PeriodicTable;
use crate::position::Position;

/// Represents a CIF data block containing key-value pairs and loop tables.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MmcifDataBlock {
    pub key_values: IndexMap<String, String>,
    pub tables: Vec<Vec<IndexMap<String, String>>>,
}

/// Represents an atom record from the `_atom_site` category.
#[derive(Debug, Clone, PartialEq)]
pub struct AtomSiteRecord {
    pub group_pdb: String,
    pub id: usize,
    pub type_symbol: String,
    pub label_atom_id: String,
    pub label_alt_id: String,
    pub label_comp_id: String,
    pub label_asym_id: String,
    pub label_entity_id: String,
    pub label_seq_id: Option<i32>,
    pub pdbx_pdb_ins_code: String,
    pub cartn_x: f64,
    pub cartn_y: f64,
    pub cartn_z: f64,
    pub occupancy: f64,
    pub b_iso_or_equiv: f64,
    pub pdbx_formal_charge: f64,
    pub auth_seq_id: Option<i32>,
    pub auth_comp_id: String,
    pub auth_asym_id: String,
    pub auth_atom_id: String,
    pub pdbx_pdb_model_num: usize,
}

impl AtomSiteRecord {
    /// Parses an `AtomSiteRecord` from an `_atom_site` table row.
    pub fn from_row(row: &IndexMap<String, String>, default_id: usize) -> Result<Option<Self>> {
        let Some(group_pdb) = row.get("_atom_site.group_PDB") else {
            return Ok(None);
        };
        if group_pdb != "ATOM" && group_pdb != "HETATM" {
            return Ok(None);
        }
        let group_pdb = group_pdb.clone();

        let id = row
            .get("_atom_site.id")
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(default_id);

        let type_symbol = row
            .get("_atom_site.type_symbol")
            .cloned()
            .unwrap_or_else(|| "X".to_string());

        let label_atom_id = row
            .get("_atom_site.label_atom_id")
            .cloned()
            .unwrap_or_else(|| "X".to_string());

        let label_alt_id = row
            .get("_atom_site.label_alt_id")
            .map(|s| if s == "." || s == "?" { "" } else { s.as_str() })
            .unwrap_or("")
            .to_string();

        let label_comp_id = row
            .get("_atom_site.label_comp_id")
            .cloned()
            .unwrap_or_else(|| "UNK".to_string());

        let label_asym_id = row
            .get("_atom_site.label_asym_id")
            .filter(|s| *s != "." && *s != "?")
            .cloned()
            .unwrap_or_else(|| "_".to_string());

        let label_entity_id = row
            .get("_atom_site.label_entity_id")
            .cloned()
            .unwrap_or_default();

        let label_seq_id = row
            .get("_atom_site.label_seq_id")
            .and_then(|s| s.parse::<i32>().ok());

        let pdbx_pdb_ins_code = row
            .get("_atom_site.pdbx_PDB_ins_code")
            .filter(|s| *s != "." && *s != "?")
            .cloned()
            .unwrap_or_default();

        let cartn_x = row
            .get("_atom_site.Cartn_x")
            .ok_or_else(|| {
                BridgeError::input_error("_atom_site.Cartn_x", "missing Cartn_x coordinate")
            })?
            .parse::<f64>()
            .map_err(|e| {
                BridgeError::input_error("_atom_site.Cartn_x", format!("invalid float: {e}"))
            })?;
        let cartn_y = row
            .get("_atom_site.Cartn_y")
            .ok_or_else(|| {
                BridgeError::input_error("_atom_site.Cartn_y", "missing Cartn_y coordinate")
            })?
            .parse::<f64>()
            .map_err(|e| {
                BridgeError::input_error("_atom_site.Cartn_y", format!("invalid float: {e}"))
            })?;
        let cartn_z = row
            .get("_atom_site.Cartn_z")
            .ok_or_else(|| {
                BridgeError::input_error("_atom_site.Cartn_z", "missing Cartn_z coordinate")
            })?
            .parse::<f64>()
            .map_err(|e| {
                BridgeError::input_error("_atom_site.Cartn_z", format!("invalid float: {e}"))
            })?;

        let occupancy = row
            .get("_atom_site.occupancy")
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(1.0);

        let b_iso_or_equiv = row
            .get("_atom_site.B_iso_or_equiv")
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(0.0);

        let pdbx_formal_charge = row
            .get("_atom_site.pdbx_formal_charge")
            .filter(|s| *s != "." && *s != "?")
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(0.0);

        let auth_seq_id = row
            .get("_atom_site.auth_seq_id")
            .and_then(|s| s.parse::<i32>().ok());

        let auth_comp_id = row
            .get("_atom_site.auth_comp_id")
            .filter(|s| *s != "." && *s != "?")
            .cloned()
            .unwrap_or_else(|| label_comp_id.clone());

        let auth_asym_id = row
            .get("_atom_site.auth_asym_id")
            .filter(|s| *s != "." && *s != "?")
            .cloned()
            .unwrap_or_else(|| label_asym_id.clone());

        let auth_atom_id = row
            .get("_atom_site.auth_atom_id")
            .filter(|s| *s != "." && *s != "?")
            .cloned()
            .unwrap_or_else(|| label_atom_id.clone());

        let pdbx_pdb_model_num = row
            .get("_atom_site.pdbx_PDB_model_num")
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(1);

        Ok(Some(Self {
            group_pdb,
            id,
            type_symbol,
            label_atom_id,
            label_alt_id,
            label_comp_id,
            label_asym_id,
            label_entity_id,
            label_seq_id,
            pdbx_pdb_ins_code,
            cartn_x,
            cartn_y,
            cartn_z,
            occupancy,
            b_iso_or_equiv,
            pdbx_formal_charge,
            auth_seq_id,
            auth_comp_id,
            auth_asym_id,
            auth_atom_id,
            pdbx_pdb_model_num,
        }))
    }
}

/// Represents a connection record from `_struct_conn` (e.g. disulfide bonds).
#[derive(Debug, Clone, PartialEq)]
pub struct StructConnRecord {
    pub id: String,
    pub conn_type_id: String,
    pub ptnr1_asym_id: String,
    pub ptnr1_seq_id: Option<i32>,
    pub ptnr1_atom_id: String,
    pub ptnr2_asym_id: String,
    pub ptnr2_seq_id: Option<i32>,
    pub ptnr2_atom_id: String,
}

impl StructConnRecord {
    /// Parses a `StructConnRecord` from a `_struct_conn` table row.
    pub fn from_row(row: &IndexMap<String, String>) -> Option<Self> {
        let conn_type_id = row.get("_struct_conn.conn_type_id")?.clone();
        let id = row.get("_struct_conn.id").cloned().unwrap_or_default();

        let ptnr1_asym_id = row
            .get("_struct_conn.ptnr1_auth_asym_id")
            .filter(|s| *s != "." && *s != "?")
            .or_else(|| {
                row.get("_struct_conn.ptnr1_label_asym_id")
                    .filter(|s| *s != "." && *s != "?")
            })
            .cloned()
            .unwrap_or_default();
        let ptnr1_seq_id = row
            .get("_struct_conn.ptnr1_auth_seq_id")
            .and_then(|s| s.parse::<i32>().ok())
            .or_else(|| {
                row.get("_struct_conn.ptnr1_label_seq_id")
                    .and_then(|s| s.parse::<i32>().ok())
            });
        let ptnr1_atom_id = row
            .get("_struct_conn.ptnr1_label_atom_id")
            .cloned()
            .unwrap_or_default();

        let ptnr2_asym_id = row
            .get("_struct_conn.ptnr2_auth_asym_id")
            .filter(|s| *s != "." && *s != "?")
            .or_else(|| {
                row.get("_struct_conn.ptnr2_label_asym_id")
                    .filter(|s| *s != "." && *s != "?")
            })
            .cloned()
            .unwrap_or_default();
        let ptnr2_seq_id = row
            .get("_struct_conn.ptnr2_auth_seq_id")
            .and_then(|s| s.parse::<i32>().ok())
            .or_else(|| {
                row.get("_struct_conn.ptnr2_label_seq_id")
                    .and_then(|s| s.parse::<i32>().ok())
            });
        let ptnr2_atom_id = row
            .get("_struct_conn.ptnr2_label_atom_id")
            .cloned()
            .unwrap_or_default();

        Some(Self {
            id,
            conn_type_id,
            ptnr1_asym_id,
            ptnr1_seq_id,
            ptnr1_atom_id,
            ptnr2_asym_id,
            ptnr2_seq_id,
            ptnr2_atom_id,
        })
    }
}

impl MmcifDataBlock {
    /// Checks if this data block contains an `_atom_site` table or key-values.
    pub fn has_atom_site(&self) -> bool {
        self.key_values.keys().any(|k| k.starts_with("_atom_site."))
            || self.tables.iter().any(|table| {
                table
                    .first()
                    .is_some_and(|row| row.keys().any(|k| k.starts_with("_atom_site.")))
            })
    }

    /// Extracts all `AtomSiteRecord` entries from the `_atom_site` table or key-values.
    pub fn get_atom_site_records(&self) -> Result<Vec<AtomSiteRecord>> {
        let mut records = Vec::new();
        if self.key_values.keys().any(|k| k.starts_with("_atom_site.")) {
            if let Some(rec) = AtomSiteRecord::from_row(&self.key_values, 1)? {
                records.push(rec);
            }
        }
        for table in &self.tables {
            if let Some(first_row) = table.first() {
                if first_row.keys().any(|k| k.starts_with("_atom_site.")) {
                    for (idx, row) in table.iter().enumerate() {
                        if let Some(rec) = AtomSiteRecord::from_row(row, idx + 1)? {
                            records.push(rec);
                        }
                    }
                }
            }
        }
        Ok(records)
    }

    /// Extracts all `StructConnRecord` entries from the `_struct_conn` table or key-values.
    pub fn get_struct_conn_records(&self) -> Vec<StructConnRecord> {
        let mut records = Vec::new();
        if self
            .key_values
            .keys()
            .any(|k| k.starts_with("_struct_conn."))
        {
            if let Some(rec) = StructConnRecord::from_row(&self.key_values) {
                records.push(rec);
            }
        }
        for table in &self.tables {
            if let Some(first_row) = table.first() {
                if first_row.keys().any(|k| k.starts_with("_struct_conn.")) {
                    for row in table {
                        if let Some(rec) = StructConnRecord::from_row(row) {
                            records.push(rec);
                        }
                    }
                }
            }
        }
        records
    }
}

/// Parser for Crystallographic Information Files (CIF/mmCIF).
///
/// Ports `proteindf_bridge.mmcif.SimpleMmcif` for Chemical Component Dictionary (CCD) entries
/// and provides full-structure parsing via the `_atom_site` category.
#[derive(Debug, Clone, Default)]
pub struct SimpleMmcif {
    data: IndexMap<String, MmcifDataBlock>,
}

#[derive(Debug, PartialEq, Eq)]
enum Token {
    DataBlock(String),
    Loop,
    Tag(String),
    Value(String),
}

impl SimpleMmcif {
    /// Creates a new empty `SimpleMmcif`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Loads an mmCIF file from the given path.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        let mut mmcif = Self::new();
        mmcif.load(path)?;
        Ok(mmcif)
    }

    /// Parses mmCIF data from a string.
    #[allow(clippy::should_implement_trait)]
    pub fn from_str(s: &str) -> Result<Self> {
        <Self as FromStr>::from_str(s)
    }

    /// Loads an mmCIF file from the given path into this instance.
    pub fn load(&mut self, path: impl AsRef<Path>) -> Result<()> {
        let path_ref = path.as_ref();
        let content = fs::read_to_string(path_ref)
            .map_err(|e| BridgeError::input_error(path_ref.display().to_string(), e.to_string()))?;
        self.load_from_str(&content)
    }

    /// Parses mmCIF data from a string.
    pub fn load_from_str(&mut self, content: &str) -> Result<()> {
        let tokens = Self::tokenize(content);
        self.parse_tokens(&tokens)
    }

    /// Returns the list of molecule / data block names found in the file.
    pub fn get_molecule_names(&self) -> Vec<String> {
        self.data.keys().cloned().collect()
    }

    /// Returns a reference to the parsed data blocks map.
    pub fn data(&self) -> &IndexMap<String, MmcifDataBlock> {
        &self.data
    }

    /// Returns a specific data block by name.
    pub fn get_data_block(&self, name: &str) -> Option<&MmcifDataBlock> {
        self.data.get(name)
    }

    /// Constructs an `AtomGroup` from the given data block name.
    ///
    /// If the block contains `_atom_site` data, parses full structure coordinates.
    /// Otherwise, parses as Chemical Component Dictionary (CCD) data.
    pub fn get_atomgroup(&self, name: &str) -> Result<AtomGroup> {
        let block = self.data.get(name).ok_or_else(|| {
            BridgeError::input_error(name, format!("Invalid mmcif data: name={}", name))
        })?;

        if block.has_atom_site() {
            return self.get_structure_atomgroup_for_block(name, None, None);
        }

        let mut ag = AtomGroup::new();

        // 1. Determine residue / molecule name from _chem_comp.id
        if let Some(id) = block.key_values.get("_chem_comp.id") {
            ag.name = id.clone();
        }

        // Check key-values and tables for atoms and component name
        Self::extract_atoms_and_name(&block.key_values, &mut ag);

        for table in &block.tables {
            for row in table {
                Self::extract_atoms_and_name(row, &mut ag);
            }
        }

        // 2. Add bonds from _chem_comp_bond
        for table in &block.tables {
            for row in table {
                if row.contains_key("_chem_comp_bond.comp_id") {
                    if let (Some(atom1_name), Some(atom2_name)) = (
                        row.get("_chem_comp_bond.atom_id_1"),
                        row.get("_chem_comp_bond.atom_id_2"),
                    ) {
                        let bond_order_str = row
                            .get("_chem_comp_bond.value_order")
                            .map(|s| s.as_str())
                            .unwrap_or("");
                        let bond_order = match bond_order_str {
                            "SING" => 1,
                            "DOUB" => 2,
                            "TRIP" => 3,
                            _ => 0,
                        };

                        let atom1_opt = ag.get_atom(atom1_name).cloned();
                        let atom2_opt = ag.get_atom(atom2_name).cloned();
                        if let (Some(atom1), Some(atom2)) = (atom1_opt, atom2_opt) {
                            ag.add_bond(&atom1, &atom2, bond_order);
                        }
                    }
                }
            }
        }

        Ok(ag)
    }

    /// Builds an `AtomGroup` hierarchy representing the mmCIF structure for the specified data block.
    ///
    /// The resulting hierarchy matches `Pdb::get_atomgroup`:
    /// `root -> model_<serial> -> <chain_id> -> <res_seq> -> <serial>_<name>`
    ///
    /// If `select_model` is `None`, all models are included.
    /// Alternate location atoms matching `select_altloc` (default: "A") or blank are retained.
    pub fn get_structure_atomgroup_for_block(
        &self,
        block_name: &str,
        select_model: Option<usize>,
        select_altloc: Option<&str>,
    ) -> Result<AtomGroup> {
        let block = self.data.get(block_name).ok_or_else(|| {
            BridgeError::input_error(block_name, format!("Data block '{block_name}' not found"))
        })?;

        let records = block.get_atom_site_records()?;
        if records.is_empty() {
            return Err(BridgeError::input_error(
                block_name,
                format!("No _atom_site records found in block '{block_name}'"),
            ));
        }

        let conns = block.get_struct_conn_records();
        let altloc_filter = select_altloc.unwrap_or("A");

        // Group records by model_num preserving order
        let mut models_map: IndexMap<usize, Vec<&AtomSiteRecord>> = IndexMap::new();
        for rec in &records {
            if let Some(target) = select_model {
                if target != rec.pdbx_pdb_model_num {
                    continue;
                }
            }
            models_map
                .entry(rec.pdbx_pdb_model_num)
                .or_default()
                .push(rec);
        }

        let mut root = AtomGroup::new();

        for (&model_serial, model_records) in &models_map {
            let model_name = format!("model_{model_serial}");
            let mut model = AtomGroup::new();
            model.name = model_name.clone();

            for item in model_records {
                // Check altloc filter: retain matching altloc or empty
                let alt_loc = item.label_alt_id.trim();
                if !alt_loc.is_empty() && alt_loc != altloc_filter {
                    continue;
                }

                // Determine chain ID (auth_asym_id prioritized, fallback to label_asym_id, default to "_")
                let mut chain_id = item.auth_asym_id.trim().to_string();
                if chain_id.is_empty() || chain_id == "_" {
                    chain_id = item.label_asym_id.trim().to_string();
                }
                if chain_id.is_empty() || chain_id == " " {
                    chain_id = "_".to_string();
                }

                if !model.has_group(&chain_id) {
                    let mut chain = AtomGroup::new();
                    chain.name = chain_id.clone();
                    model.set_group(&chain_id, chain);
                }

                // Determine residue sequence key:
                // Prioritize auth_seq_id for both ATOM and HETATM, fallback to label_seq_id, default to 1
                let res_seq = item.auth_seq_id.or(item.label_seq_id).unwrap_or(1);
                let res_key = format!("{res_seq}");

                let mut res_name = item.auth_comp_id.clone();
                if matches!(res_name.as_str(), "HID" | "HIE" | "HIP") {
                    res_name = "HIS".to_string();
                }

                if let Some(chain) = model.get_group_mut(&chain_id) {
                    if !chain.has_group(&res_key) {
                        let mut residue = AtomGroup::new();
                        residue.name = res_name;
                        chain.set_group(&res_key, residue);
                    }
                }

                let mut element = item.type_symbol.clone();
                if element == "D" {
                    element = "H".to_string();
                }

                let mut atom = Atom::new();
                if let Ok(num) = PeriodicTable::get_atomic_number(&element) {
                    atom.set_atomic_number(num);
                }
                atom.xyz = Position::new(item.cartn_x, item.cartn_y, item.cartn_z);
                atom.name = item.auth_atom_id.clone();
                atom.charge = item.pdbx_formal_charge;

                let atom_key = format!("{}_{}", item.id, item.auth_atom_id);
                if let Some(chain) = model.get_group_mut(&chain_id) {
                    if let Some(residue) = chain.get_group_mut(&res_key) {
                        residue.set_atom(&atom_key, atom);
                    }
                }
            }

            // Link disulfide bonds from _struct_conn
            for conn in &conns {
                if conn.conn_type_id == "disulf" {
                    if let (Some(seq1), Some(seq2)) = (conn.ptnr1_seq_id, conn.ptnr2_seq_id) {
                        let res_key1 = format!("{seq1}");
                        let res_key2 = format!("{seq2}");

                        let sg1_opt = model
                            .get_group(&conn.ptnr1_asym_id)
                            .and_then(|c| c.get_group(&res_key1))
                            .and_then(|r| r.get_atom(&conn.ptnr1_atom_id))
                            .cloned();

                        let sg2_opt = model
                            .get_group(&conn.ptnr2_asym_id)
                            .and_then(|c| c.get_group(&res_key2))
                            .and_then(|r| r.get_atom(&conn.ptnr2_atom_id))
                            .cloned();

                        if let (Some(sg1), Some(sg2)) = (sg1_opt, sg2_opt) {
                            model.add_bond(&sg1, &sg2, 1);
                        }
                    }
                }
            }

            root.set_group(&model_name, model);
        }

        Ok(root)
    }

    /// Builds an `AtomGroup` hierarchy representing the mmCIF structure from the first data block.
    pub fn get_structure_atomgroup(
        &self,
        select_model: Option<usize>,
        select_altloc: Option<&str>,
    ) -> Result<AtomGroup> {
        let first_block = self
            .data
            .keys()
            .next()
            .ok_or_else(|| BridgeError::input_error("mmCIF", "No data blocks found in file"))?;
        self.get_structure_atomgroup_for_block(first_block, select_model, select_altloc)
    }

    fn extract_atoms_and_name(dict: &IndexMap<String, String>, ag: &mut AtomGroup) {
        if let Some(id) = dict.get("_chem_comp.id") {
            ag.name = id.clone();
        }

        if let Some(atom_id) = dict.get("_chem_comp_atom.atom_id") {
            let mut atom = Atom::new();
            atom.name = atom_id.clone();

            let mut symbol = dict
                .get("_chem_comp_atom.type_symbol")
                .cloned()
                .unwrap_or_else(|| "X".to_string());
            if symbol == "D" {
                symbol = "H".to_string();
            }
            if let Ok(num) = PeriodicTable::get_atomic_number(&symbol) {
                atom.set_atomic_number(num);
            }

            let x = Self::get_coordinate("x", dict);
            let y = Self::get_coordinate("y", dict);
            let z = Self::get_coordinate("z", dict);
            if let (Some(x), Some(y), Some(z)) = (x, y, z) {
                atom.xyz = Position::new(x, y, z);
            }

            if let Some(charge_str) = dict.get("_chem_comp_atom.charge") {
                if charge_str != "?" {
                    if let Ok(charge) = charge_str.parse::<f64>() {
                        atom.charge = charge;
                    }
                }
            }

            ag.set_atom(atom_id, atom);
        }
    }

    /// Extracts a coordinate axis (x, y, or z) prioritizing ideal coordinates over model coordinates.
    fn get_coordinate(axis: &str, dict: &IndexMap<String, String>) -> Option<f64> {
        let ideal_key = format!("_chem_comp_atom.pdbx_model_Cartn_{}_ideal", axis);
        let model_key = format!("_chem_comp_atom.model_Cartn_{}", axis);

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

    /// Parses tokens into data blocks, key-values, and loop tables.
    fn parse_tokens(&mut self, tokens: &[Token]) -> Result<()> {
        let mut idx = 0;
        let mut current_block_name: Option<String> = None;
        let mut current_block = MmcifDataBlock::default();

        while idx < tokens.len() {
            match &tokens[idx] {
                Token::DataBlock(name) => {
                    if let Some(prev_name) = current_block_name.take() {
                        self.data.insert(prev_name, current_block);
                        current_block = MmcifDataBlock::default();
                    }
                    current_block_name = Some(name.clone());
                    idx += 1;
                }
                Token::Tag(key) => {
                    idx += 1;
                    if idx < tokens.len() {
                        let value = match &tokens[idx] {
                            Token::Value(v) => v.clone(),
                            Token::Tag(v) => v.clone(),
                            _ => String::new(),
                        };
                        current_block.key_values.insert(key.clone(), value);
                        idx += 1;
                    }
                }
                Token::Loop => {
                    idx += 1;
                    // Read header tags
                    let mut headers = Vec::new();
                    while idx < tokens.len() {
                        if let Token::Tag(tag) = &tokens[idx] {
                            headers.push(tag.clone());
                            idx += 1;
                        } else {
                            break;
                        }
                    }

                    if headers.is_empty() {
                        continue;
                    }

                    let num_cols = headers.len();
                    let mut table = Vec::new();

                    // Read rows until next Tag, Loop, or DataBlock
                    while idx < tokens.len() {
                        match &tokens[idx] {
                            Token::Tag(_) | Token::Loop | Token::DataBlock(_) => {
                                break;
                            }
                            Token::Value(_) => {
                                // Read num_cols values
                                let mut row = IndexMap::new();
                                for header in &headers {
                                    if idx < tokens.len() {
                                        if let Token::Value(val) = &tokens[idx] {
                                            row.insert(header.clone(), val.clone());
                                            idx += 1;
                                        } else {
                                            break;
                                        }
                                    }
                                }
                                if row.len() == num_cols {
                                    table.push(row);
                                }
                            }
                        }
                    }

                    current_block.tables.push(table);
                }
                Token::Value(_) => {
                    idx += 1;
                }
            }
        }

        if let Some(name) = current_block_name {
            self.data.insert(name, current_block);
        }

        Ok(())
    }

    /// Tokenizes mmCIF file content into a token stream.
    fn tokenize(content: &str) -> Vec<Token> {
        let mut tokens = Vec::new();
        let lines: Vec<&str> = content.lines().collect();
        let mut line_idx = 0;

        while line_idx < lines.len() {
            let line = lines[line_idx];

            // 1. Semicolon block check: starts with ';' at line index 0
            if let Some(rest) = line.strip_prefix(';') {
                let mut block_text = String::new();
                if !rest.is_empty() {
                    block_text.push_str(rest);
                }
                line_idx += 1;
                while line_idx < lines.len() {
                    let cur_line = lines[line_idx];
                    if cur_line.starts_with(';') {
                        line_idx += 1;
                        break;
                    }
                    if !block_text.is_empty() {
                        block_text.push('\n');
                    }
                    block_text.push_str(cur_line);
                    line_idx += 1;
                }
                tokens.push(Token::Value(block_text));
                continue;
            }

            // 2. Line comment check
            let trimmed = line.trim();
            if trimmed.starts_with('#') || trimmed.is_empty() {
                line_idx += 1;
                continue;
            }

            // 3. Scan line character by character
            let chars: Vec<char> = line.chars().collect();
            let mut i = 0;
            while i < chars.len() {
                // Skip whitespace
                if chars[i].is_whitespace() {
                    i += 1;
                    continue;
                }

                // Comment outside quotes
                if chars[i] == '#' {
                    break;
                }

                // Quoted string: "..." or '...'
                if chars[i] == '"' || chars[i] == '\'' {
                    let quote_char = chars[i];
                    i += 1;
                    let mut s = String::new();
                    while i < chars.len() && chars[i] != quote_char {
                        s.push(chars[i]);
                        i += 1;
                    }
                    if i < chars.len() && chars[i] == quote_char {
                        i += 1;
                    }
                    tokens.push(Token::Value(s));
                    continue;
                }

                // Normal unquoted word
                let mut word = String::new();
                while i < chars.len() && !chars[i].is_whitespace() && chars[i] != '#' {
                    word.push(chars[i]);
                    i += 1;
                }

                if word.starts_with("data_") {
                    tokens.push(Token::DataBlock(word));
                } else if word == "loop_" {
                    tokens.push(Token::Loop);
                } else if word.starts_with('_') {
                    tokens.push(Token::Tag(word));
                } else {
                    tokens.push(Token::Value(word));
                }
            }

            line_idx += 1;
        }

        tokens
    }
}

impl FromStr for SimpleMmcif {
    type Err = BridgeError;

    fn from_str(s: &str) -> Result<Self> {
        let mut mmcif = Self::new();
        mmcif.load_from_str(s)?;
        Ok(mmcif)
    }
}
