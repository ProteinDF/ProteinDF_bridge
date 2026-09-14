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

/// Parser for Crystallographic Information Files (CIF/mmCIF).
///
/// Ports `proteindf_bridge.mmcif.SimpleMmcif` for Chemical Component Dictionary (CCD) entries.
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

    /// Constructs an `AtomGroup` from the Chemical Component Dictionary (CCD) data of the given molecule name.
    pub fn get_atomgroup(&self, name: &str) -> Result<AtomGroup> {
        let block = self.data.get(name).ok_or_else(|| {
            BridgeError::input_error(name, format!("Invalid mmcif data: name={}", name))
        })?;

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
