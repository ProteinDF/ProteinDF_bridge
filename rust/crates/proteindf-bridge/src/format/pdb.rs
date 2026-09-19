// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::collections::BTreeMap;
use std::fmt;
use std::fs;
use std::path::Path;

use crate::atom::Atom;
use crate::atom_group::AtomGroup;
use crate::error::{BridgeError, Result};
use crate::periodic_table::PeriodicTable;
use crate::position::Position;

/// Represents a single PDB ATOM, HETATM, or TER record.
#[derive(Debug, Clone, PartialEq)]
pub struct PdbRecord {
    pub record_name: String,
    pub serial: usize,
    pub name: String,
    pub alt_loc: String,
    pub res_name: String,
    pub chain_id: String,
    pub res_seq: i32,
    pub i_code: String,
    pub coord: [f64; 3],
    pub occupancy: f64,
    pub temp_factor: f64,
    pub element: String,
    pub charge: String,
}

impl Default for PdbRecord {
    fn default() -> Self {
        Self {
            record_name: "ATOM  ".to_string(),
            serial: 1,
            name: String::new(),
            alt_loc: " ".to_string(),
            res_name: String::new(),
            chain_id: " ".to_string(),
            res_seq: 1,
            i_code: " ".to_string(),
            coord: [0.0, 0.0, 0.0],
            occupancy: 1.0,
            temp_factor: 0.0,
            element: "  ".to_string(),
            charge: "  ".to_string(),
        }
    }
}

/// Represents a PDB SSBOND record (disulfide bond linkage).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SsBondRecord {
    pub chain_id1: String,
    pub seq_num1: i32,
    pub chain_id2: String,
    pub seq_num2: i32,
}

/// PDB file parser and serializer corresponding to `proteindf_bridge.biopdb.Pdb`.
#[derive(Debug, Clone)]
pub struct Pdb {
    data: BTreeMap<usize, Vec<PdbRecord>>,
    ssbonds: Vec<SsBondRecord>,
    conects: Vec<(usize, usize)>,
    mode: Option<String>,
}

impl Default for Pdb {
    fn default() -> Self {
        Self::new(None)
    }
}

impl Pdb {
    /// AmberTools compatibility version (default 22).
    pub const AMBER_TOOLS_VER: u32 = 22;

    /// Creates a new empty `Pdb` object.
    pub fn new(mode: Option<&str>) -> Self {
        let mut data = BTreeMap::new();
        data.insert(1, Vec::new());
        Self {
            data,
            ssbonds: Vec::new(),
            conects: Vec::new(),
            mode: mode.map(|m| m.to_ascii_uppercase()),
        }
    }

    /// Creates a `Pdb` object by loading a file from disk.
    pub fn from_file(path: impl AsRef<Path>, mode: Option<&str>) -> Result<Self> {
        let mut pdb = Self::new(mode);
        pdb.load(path)?;
        Ok(pdb)
    }

    /// Creates a `Pdb` object by parsing a PDB-formatted string.
    pub fn from_str(content: &str, mode: Option<&str>) -> Result<Self> {
        let mut pdb = Self::new(mode);
        pdb.parse_str(content)?;
        Ok(pdb)
    }

    /// Returns the active dialect mode (e.g. "AMBER", "FORMAL").
    pub fn mode(&self) -> Option<&str> {
        self.mode.as_deref()
    }

    /// Sets the active dialect mode.
    pub fn set_mode(&mut self, mode: Option<&str>) {
        self.mode = mode.map(|m| m.to_ascii_uppercase());
    }

    /// Returns a reference to the parsed records grouped by model serial number.
    pub fn data(&self) -> &BTreeMap<usize, Vec<PdbRecord>> {
        &self.data
    }

    /// Returns a reference to the parsed SSBOND records.
    pub fn ssbonds(&self) -> &[SsBondRecord] {
        &self.ssbonds
    }

    /// Returns a reference to the parsed CONECT records as deduplicated `(serial1, serial2)` pairs.
    pub fn conects(&self) -> &[(usize, usize)] {
        &self.conects
    }

    /// Renumbers atom serial numbers within each model starting from 1.
    pub fn renumber(&mut self) {
        for records in self.data.values_mut() {
            for (idx, record) in records.iter_mut().enumerate() {
                record.serial = idx + 1;
            }
        }
    }

    /// Loads and parses a PDB file from the given file path.
    pub fn load(&mut self, path: impl AsRef<Path>) -> Result<()> {
        let path_ref = path.as_ref();
        let content = fs::read_to_string(path_ref).map_err(|e| {
            BridgeError::input_error(
                path_ref.display().to_string(),
                format!("failed to read file: {e}"),
            )
        })?;
        self.parse_str(&content)
    }

    /// Parses PDB contents from a string.
    pub fn parse_str(&mut self, content: &str) -> Result<()> {
        self.data.clear();
        self.ssbonds.clear();
        self.conects.clear();

        let mut model_serial: usize = 1;
        let mut chain_serial: usize = 0;
        self.data.entry(model_serial).or_default();

        for line in content.lines() {
            let line = line.trim_end_matches(['\r', '\n']);
            let chars: Vec<char> = line.chars().collect();

            if chars.is_empty() {
                continue;
            }

            let record_name = slice_chars(&chars, 0, 6);

            if record_name == "SSBOND" {
                let chain_id1 = slice_chars(&chars, 15, 16);
                let seq_num1 = slice_chars(&chars, 17, 21)
                    .trim()
                    .parse::<i32>()
                    .map_err(|e| {
                        BridgeError::input_error("SSBOND seqNum1", format!("invalid integer: {e}"))
                    })?;
                let chain_id2 = slice_chars(&chars, 29, 30);
                let seq_num2 = slice_chars(&chars, 31, 35)
                    .trim()
                    .parse::<i32>()
                    .map_err(|e| {
                        BridgeError::input_error("SSBOND seqNum2", format!("invalid integer: {e}"))
                    })?;

                self.ssbonds.push(SsBondRecord {
                    chain_id1,
                    seq_num1,
                    chain_id2,
                    seq_num2,
                });
            } else if record_name == "ATOM  " || record_name == "HETATM" {
                // Pad line to 80 characters with spaces if needed
                let mut padded_chars = chars.clone();
                if padded_chars.len() < 80 {
                    padded_chars.resize(80, ' ');
                }

                let serial = slice_chars(&padded_chars, 6, 11)
                    .trim()
                    .parse::<usize>()
                    .map_err(|e| {
                        BridgeError::input_error("ATOM serial", format!("invalid integer: {e}"))
                    })?;
                let name4 = slice_chars(&padded_chars, 12, 16);
                let name = name4.trim().to_string();
                let alt_loc = slice_chars(&padded_chars, 16, 17);
                let mut res_name = slice_chars(&padded_chars, 17, 20);

                // Rename AMBER residue dialects
                if matches!(res_name.as_str(), "HID" | "HIE" | "HIP") {
                    res_name = "HIS".to_string();
                }

                let mut chain_id = slice_chars(&padded_chars, 21, 22);
                if chain_id == " " && res_name != "WAT" {
                    chain_id = ((b'A' + (chain_serial % 26) as u8) as char).to_string();
                }

                let res_seq = slice_chars(&padded_chars, 22, 26)
                    .trim()
                    .parse::<i32>()
                    .map_err(|e| {
                        BridgeError::input_error("ATOM res_seq", format!("invalid integer: {e}"))
                    })?;
                let i_code = slice_chars(&padded_chars, 26, 27);

                let coord_x = slice_chars(&padded_chars, 30, 38)
                    .trim()
                    .parse::<f64>()
                    .map_err(|e| {
                        BridgeError::input_error("ATOM coord_x", format!("invalid float: {e}"))
                    })?;
                let coord_y = slice_chars(&padded_chars, 38, 46)
                    .trim()
                    .parse::<f64>()
                    .map_err(|e| {
                        BridgeError::input_error("ATOM coord_y", format!("invalid float: {e}"))
                    })?;
                let coord_z = slice_chars(&padded_chars, 46, 54)
                    .trim()
                    .parse::<f64>()
                    .map_err(|e| {
                        BridgeError::input_error("ATOM coord_z", format!("invalid float: {e}"))
                    })?;

                let occ_str = slice_chars(&padded_chars, 54, 60);
                let occupancy = if occ_str.trim().is_empty() {
                    1.0
                } else {
                    occ_str.trim().parse::<f64>().map_err(|e| {
                        BridgeError::input_error("ATOM occupancy", format!("invalid float: {e}"))
                    })?
                };

                let temp_str = slice_chars(&padded_chars, 60, 66);
                let temp_factor = if temp_str.trim().is_empty() {
                    0.0
                } else {
                    temp_str.trim().parse::<f64>().map_err(|e| {
                        BridgeError::input_error("ATOM temp_factor", format!("invalid float: {e}"))
                    })?
                };

                let element_col = slice_chars(&padded_chars, 76, 78);
                let element = guess_element(&name4, &element_col);

                let charge_col = slice_chars(&padded_chars, 78, 80);
                let charge = format_charge(&charge_col);

                let record = PdbRecord {
                    record_name,
                    serial,
                    name,
                    alt_loc,
                    res_name,
                    chain_id,
                    res_seq,
                    i_code,
                    coord: [coord_x, coord_y, coord_z],
                    occupancy,
                    temp_factor,
                    element,
                    charge,
                };
                self.data.entry(model_serial).or_default().push(record);
            } else if record_name == "MODEL "
                || (record_name.starts_with("MODEL") && chars.len() >= 14)
            {
                let serial = slice_chars(&chars, 10, 14)
                    .trim()
                    .parse::<usize>()
                    .map_err(|e| {
                        BridgeError::input_error("MODEL serial", format!("invalid integer: {e}"))
                    })?;
                model_serial = serial;
                self.data.entry(model_serial).or_default();
                chain_serial = 0;
            } else if record_name == "TER   " || record_name.starts_with("TER") {
                let mut padded_chars = chars.clone();
                if padded_chars.len() < 27 {
                    padded_chars.resize(27, ' ');
                }
                let serial = slice_chars(&padded_chars, 6, 11)
                    .trim()
                    .parse::<usize>()
                    .unwrap_or(0);
                let res_name = slice_chars(&padded_chars, 17, 20);
                let mut chain_id = slice_chars(&padded_chars, 21, 22);
                if chain_id == " " && res_name != "WAT" {
                    chain_id = ((b'A' + (chain_serial % 26) as u8) as char).to_string();
                    chain_serial += 1;
                }
                let res_seq = slice_chars(&padded_chars, 22, 26)
                    .trim()
                    .parse::<i32>()
                    .unwrap_or(0);
                let i_code = slice_chars(&padded_chars, 26, 27);

                let record = PdbRecord {
                    record_name: "TER   ".to_string(),
                    serial,
                    name: String::new(),
                    alt_loc: " ".to_string(),
                    res_name,
                    chain_id,
                    res_seq,
                    i_code,
                    coord: [0.0, 0.0, 0.0],
                    occupancy: 0.0,
                    temp_factor: 0.0,
                    element: "  ".to_string(),
                    charge: "  ".to_string(),
                };
                self.data.entry(model_serial).or_default().push(record);
            } else if record_name == "CONECT" {
                let serial_str = slice_chars(&chars, 6, 11);
                let serial = serial_str.trim().parse::<usize>().map_err(|e| {
                    BridgeError::input_error("CONECT serial", format!("invalid integer: {e}"))
                })?;

                // Up to 4 bonded partners in columns 12-16, 17-21, 22-26, 27-31
                let partner_ranges = [(11, 16), (16, 21), (21, 26), (26, 31)];
                for (start, end) in partner_ranges {
                    let partner_str = slice_chars(&chars, start, end);
                    let trimmed = partner_str.trim();
                    if !trimmed.is_empty() {
                        let partner = trimmed.parse::<usize>().map_err(|e| {
                            BridgeError::input_error(
                                "CONECT partner serial",
                                format!("invalid integer: {e}"),
                            )
                        })?;
                        if partner != serial {
                            let bond_key = (serial.min(partner), serial.max(partner));
                            if !self.conects.contains(&bond_key) {
                                self.conects.push(bond_key);
                            }
                        }
                    }
                }
            }
        }

        Ok(())
    }

    /// Builds an `AtomGroup` hierarchy representing the PDB structure.
    ///
    /// The resulting hierarchy is structured as:
    /// `root -> model_<serial> -> <chain_id> -> <res_seq> -> <serial>_<name>`
    ///
    /// If the PDB file contains `SSBOND` (disulfide bonds) or `CONECT` records, the returned
    /// `AtomGroup` includes the explicit bond topology parsed from the file with deduplication.
    /// According to the bond priority policy (see `RUST_PORT_SPEC.md` §3.8), explicit
    /// file-derived bonds take precedence over heuristic estimation (`Bond::setup()`).
    ///
    /// If `select_model` is `None`, all models are included.
    /// Alternate location atoms matching `select_altloc` (default: "A") or blank are retained.
    pub fn get_atomgroup(
        &self,
        select_model: Option<usize>,
        select_altloc: Option<&str>,
    ) -> Result<AtomGroup> {
        let altloc_filter = select_altloc.unwrap_or("A");
        let mut root = AtomGroup::new();

        for (&model_serial, model_items) in &self.data {
            if let Some(target) = select_model {
                if target != model_serial {
                    continue;
                }
            }

            let model_name = format!("model_{model_serial}");
            let mut model = AtomGroup::new();
            model.name = model_name.clone();
            let mut serial_to_atom = std::collections::HashMap::new();

            for item in model_items {
                if item.record_name == "ATOM  " || item.record_name == "HETATM" {
                    let mut chain_id = item.chain_id.clone();
                    if chain_id == " " {
                        chain_id = "_".to_string();
                    }

                    if !model.has_group(&chain_id) {
                        let mut chain = AtomGroup::new();
                        chain.name = chain_id.clone();
                        model.set_group(&chain_id, chain);
                    }

                    let res_key = format!("{}", item.res_seq);
                    if let Some(chain) = model.get_group_mut(&chain_id) {
                        if !chain.has_group(&res_key) {
                            let mut residue = AtomGroup::new();
                            residue.name = item.res_name.clone();
                            chain.set_group(&res_key, residue);
                        }
                    }

                    let alt_loc_str = item.alt_loc.trim();
                    if alt_loc_str.is_empty() || alt_loc_str == altloc_filter {
                        let charge_val = item.charge.trim().parse::<f64>().unwrap_or(0.0);
                        let mut atom = Atom::new();
                        if let Ok(num) = PeriodicTable::get_atomic_number(&item.element) {
                            atom.set_atomic_number(num);
                        }
                        atom.xyz = Position::new(item.coord[0], item.coord[1], item.coord[2]);
                        atom.name = item.name.clone();
                        atom.charge = charge_val;

                        let atom_key = format!("{}_{}", item.serial, item.name);
                        if let Some(chain) = model.get_group_mut(&chain_id) {
                            if let Some(residue) = chain.get_group_mut(&res_key) {
                                residue.set_atom(&atom_key, atom);
                                if let Some(stored) = residue.get_atom(&atom_key) {
                                    serial_to_atom.insert(item.serial, stored.clone());
                                }
                            }
                        }
                    }
                }
            }

            // Track bonded atom path pairs to prevent duplicate bonds between SSBOND and CONECT
            let mut existing_bonds = std::collections::HashSet::new();

            // Link SSBOND disulfide bonds
            for ssbond in &self.ssbonds {
                let res_key1 = format!("{}", ssbond.seq_num1);
                let res_key2 = format!("{}", ssbond.seq_num2);

                let sg1_opt = model
                    .get_group(&ssbond.chain_id1)
                    .and_then(|c| c.get_group(&res_key1))
                    .and_then(|r| r.get_atom("SG"))
                    .cloned();

                let sg2_opt = model
                    .get_group(&ssbond.chain_id2)
                    .and_then(|c| c.get_group(&res_key2))
                    .and_then(|r| r.get_atom("SG"))
                    .cloned();

                if let (Some(sg1), Some(sg2)) = (sg1_opt, sg2_opt) {
                    let pair = if sg1.path < sg2.path {
                        (sg1.path.clone(), sg2.path.clone())
                    } else {
                        (sg2.path.clone(), sg1.path.clone())
                    };
                    if existing_bonds.insert(pair) {
                        model.add_bond(&sg1, &sg2, 1);
                    }
                }
            }

            // Link CONECT bonds
            for &(s1, s2) in &self.conects {
                if let (Some(a1), Some(a2)) = (serial_to_atom.get(&s1), serial_to_atom.get(&s2)) {
                    let pair = if a1.path < a2.path {
                        (a1.path.clone(), a2.path.clone())
                    } else {
                        (a2.path.clone(), a1.path.clone())
                    };
                    if existing_bonds.insert(pair) {
                        model.add_bond(a1, a2, 1);
                    }
                }
            }

            root.set_group(&model_name, model);
        }

        Ok(root)
    }

    /// Reconstructs the internal PDB records from an `AtomGroup`.
    pub fn set_by_atomgroup(
        &mut self,
        atomgroup: &AtomGroup,
        is_charge2tempfactor: bool,
    ) -> Result<()> {
        let modified_ag = self.get_modpdb_atomgroup(atomgroup);

        self.data.clear();
        let mut serial = 1;

        for (model_key, model) in modified_ag.groups() {
            let model_serial = model_key
                .strip_prefix("model_")
                .and_then(|s| {
                    let digits: String = s.chars().take_while(|c| c.is_ascii_digit()).collect();
                    digits.parse::<usize>().ok()
                })
                .unwrap_or(1);

            let records = self.data.entry(model_serial).or_default();

            for (chain_id, chain) in model.groups() {
                let rec_chain_id = if chain_id != "_" {
                    chain_id.clone()
                } else {
                    " ".to_string()
                };

                for (res_key, residue) in chain.groups() {
                    let digits: String =
                        res_key.chars().take_while(|c| c.is_ascii_digit()).collect();
                    let res_seq = digits.parse::<i32>().unwrap_or(0);

                    let res_name = residue.name.clone();
                    let mut has_oxt = false;

                    for (_atom_key, atom) in residue.atoms() {
                        let is_oxt = atom.name.trim() == "OXT";
                        if is_oxt {
                            has_oxt = true;
                        }

                        let temp_factor = if is_charge2tempfactor {
                            atom.charge
                        } else {
                            0.0
                        };

                        let charge_str = if atom.charge.abs() < 1e-6 {
                            "  ".to_string()
                        } else {
                            format!("{:+1.0}", atom.charge)
                        };

                        let element_str = atom.symbol().unwrap_or("X").to_string();

                        let record = PdbRecord {
                            record_name: "ATOM  ".to_string(),
                            serial,
                            name: atom.name.clone(),
                            alt_loc: " ".to_string(),
                            res_name: res_name.clone(),
                            chain_id: rec_chain_id.clone(),
                            res_seq,
                            i_code: " ".to_string(),
                            coord: [atom.xyz.x, atom.xyz.y, atom.xyz.z],
                            occupancy: 1.0,
                            temp_factor,
                            element: element_str,
                            charge: charge_str,
                        };
                        records.push(record);
                        serial += 1;
                    }

                    if has_oxt {
                        let ter_record = PdbRecord {
                            record_name: "TER   ".to_string(),
                            serial,
                            name: String::new(),
                            alt_loc: " ".to_string(),
                            res_name: res_name.clone(),
                            chain_id: rec_chain_id.clone(),
                            res_seq,
                            i_code: " ".to_string(),
                            coord: [0.0, 0.0, 0.0],
                            occupancy: 0.0,
                            temp_factor: 0.0,
                            element: "  ".to_string(),
                            charge: "  ".to_string(),
                        };
                        records.push(ter_record);
                        serial += 1;
                    }
                }

                // Add TER record at chain end if not already terminated
                if let Some(last) = records.last() {
                    if last.record_name != "TER   " {
                        let (res_name, res_seq, i_code) =
                            (last.res_name.clone(), last.res_seq, last.i_code.clone());
                        let ter_record = PdbRecord {
                            record_name: "TER   ".to_string(),
                            serial,
                            name: String::new(),
                            alt_loc: " ".to_string(),
                            res_name,
                            chain_id: rec_chain_id.clone(),
                            res_seq,
                            i_code,
                            coord: [0.0, 0.0, 0.0],
                            occupancy: 0.0,
                            temp_factor: 0.0,
                            element: "  ".to_string(),
                            charge: "  ".to_string(),
                        };
                        records.push(ter_record);
                        serial += 1;
                    }
                }
            }
        }

        self.sort_by_serial();
        Ok(())
    }

    /// Sorts records in each model by serial number.
    pub fn sort_by_serial(&mut self) {
        for records in self.data.values_mut() {
            records.sort_by_key(|r| r.serial);
        }
    }

    /// Applies dialect renaming rules (AMBER / FORMAL) to residues and atoms.
    pub fn get_modpdb_atomgroup(&self, ag_protein: &AtomGroup) -> AtomGroup {
        let mode = self.mode.as_deref();
        let mut modified = ag_protein.clone();

        for (_model_key, model) in modified.groups_mut() {
            for (_chain_key, chain) in model.groups_mut() {
                for (_res_key, res) in chain.groups_mut() {
                    self.modpdb_res(res, mode);
                    let res_name = res.name.clone();
                    for (_atom_key, atom) in res.atoms_mut() {
                        self.modpdb_resatom(&res_name, atom, mode);
                        self.modpdb_atom(atom, mode);
                    }
                }
            }
        }

        modified
    }

    fn modpdb_res(&self, res: &mut AtomGroup, mode: Option<&str>) {
        if mode == Some("AMBER") {
            Self::rename_to_amber_dialect(res);
        }

        let resname = res.name.trim().to_ascii_uppercase();
        if mode == Some("AMBER") {
            if resname == "NA" {
                res.name = "Na+".to_string();
            } else if resname == "CL" {
                res.name = "Cl-".to_string();
            }
        } else if resname == "NA" {
            res.name = "NA ".to_string();
        } else if resname == "CL" {
            res.name = "CL ".to_string();
        }
    }

    fn modpdb_atom(&self, atom: &mut Atom, mode: Option<&str>) {
        let atomname = atom.name.trim().to_ascii_uppercase();
        let symbol = atom.symbol().unwrap_or("X");

        if mode == Some("AMBER") {
            if atomname == "NA" && symbol == "Na" {
                atom.name = "Na+".to_string();
            } else if atomname == "CL" && symbol == "Cl" {
                atom.name = "Cl-".to_string();
            }
        } else if atomname == "NA" && symbol == "Na" {
            atom.name = "NA".to_string();
        } else if atomname == "CL" && symbol == "Cl" {
            atom.name = "CL".to_string();
        }
    }

    fn modpdb_resatom(&self, resname: &str, atom: &mut Atom, mode: Option<&str>) {
        let resname_upper = resname.trim().to_ascii_uppercase();
        let atomname_upper = atom.name.trim().to_ascii_uppercase();

        if resname_upper == "NME" {
            if mode == Some("AMBER") {
                if atomname_upper == "HN2" {
                    atom.name = "H".to_string();
                } else if atomname_upper == "CH3" {
                    atom.name = "C".to_string();
                }
            } else if atomname_upper == "H" {
                atom.name = "HN2".to_string();
            } else if atomname_upper == "C" {
                atom.name = "CH3".to_string();
            }
        }
    }

    fn rename_to_amber_dialect(res: &mut AtomGroup) {
        if res.name == "HIS" {
            let has_delta_h = res.has_atomname("HD1") && res.has_atomname("HD2");
            let has_epsilon_h = res.has_atomname("HE1") && res.has_atomname("HE2");

            if has_delta_h && has_epsilon_h {
                res.name = "HIP".to_string();
            } else if has_delta_h {
                res.name = "HID".to_string();
            } else if has_epsilon_h {
                res.name = "HIE".to_string();
            }
        }
    }

    /// Formats the PDB structure as an 80-column PDB string.
    pub fn get_text(&self) -> String {
        let mut output = String::new();

        for (&model_serial, records) in &self.data {
            output.push_str(&format!("MODEL     {:>4}\n", model_serial));

            for item in records {
                if item.record_name == "ATOM  " || item.record_name == "HETATM" {
                    let formatted_name = format_chimera_atom_name(&item.name, &item.element);
                    let charge_str = format_charge_display(&item.charge);
                    let element_upper = item.element.trim().to_ascii_uppercase();

                    output.push_str(&format!(
                        "{:<6}{:5} {:4}{:1}{:>3} {:1}{:4}{:1}   {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}          {:>2}{:2}\n",
                        item.record_name,
                        item.serial,
                        formatted_name,
                        item.alt_loc,
                        item.res_name,
                        item.chain_id,
                        item.res_seq,
                        item.i_code,
                        item.coord[0],
                        item.coord[1],
                        item.coord[2],
                        item.occupancy,
                        item.temp_factor,
                        element_upper,
                        charge_str,
                    ));
                } else if item.record_name == "TER   " {
                    output.push_str(&format!(
                        "TER   {:5}      {:3} {:1}{:4}{:1}\n",
                        item.serial, item.res_name, item.chain_id, item.res_seq, item.i_code,
                    ));
                }
            }
        }

        output
    }

    /// Writes the PDB structure to a file on disk.
    pub fn save(&self, path: impl AsRef<Path>) -> Result<()> {
        let path_ref = path.as_ref();
        fs::write(path_ref, self.get_text()).map_err(|e| {
            BridgeError::general(format!(
                "failed to write PDB to {}: {}",
                path_ref.display(),
                e
            ))
        })?;
        Ok(())
    }
}

impl fmt::Display for Pdb {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.get_text())
    }
}

// ---------------------------------------------------------------------------
// Helper functions for PDB fixed-column parsing and formatting
// ---------------------------------------------------------------------------

fn slice_chars(chars: &[char], start: usize, end: usize) -> String {
    if start >= chars.len() {
        return String::new();
    }
    let actual_end = end.min(chars.len());
    chars[start..actual_end].iter().collect()
}

fn guess_element(name4: &str, element_col: &str) -> String {
    let el = element_col.trim();
    if !el.is_empty() {
        if el == "D" {
            "H".to_string()
        } else {
            el.to_string()
        }
    } else {
        // Chimera element deduction rules
        let name4_upper = name4.to_ascii_uppercase();
        let chars: Vec<char> = name4_upper.chars().collect();
        let name_trimmed = name4_upper.trim();

        if name_trimmed.len() == 4 && chars.first().copied() == Some('H') {
            "H".to_string()
        } else if let Some(&c) = chars.first() {
            if c.is_ascii_digit() {
                chars.get(1).copied().unwrap_or(' ').to_string()
            } else {
                let name2s: String = chars.iter().take(2).collect();
                let name2s_trimmed = name2s.trim();
                if name2s_trimmed.len() == 2 {
                    let first = chars[0];
                    let second = chars[1].to_ascii_lowercase();
                    format!("{first}{second}")
                } else if chars.len() >= 2 && chars[1] != ' ' {
                    chars[1].to_string()
                } else {
                    chars.first().unwrap_or(&' ').to_string()
                }
            }
        } else {
            "X".to_string()
        }
    }
}

fn format_charge(charge_col: &str) -> String {
    let trimmed = charge_col.trim();
    if trimmed.is_empty() {
        "  ".to_string()
    } else {
        let last_char = trimmed.chars().last().unwrap();
        if last_char == '+' || last_char == '-' {
            let rest: String = trimmed.chars().take(trimmed.len() - 1).collect();
            format!("{last_char}{rest}")
        } else {
            trimmed.to_string()
        }
    }
}

fn format_chimera_atom_name(name: &str, element: &str) -> String {
    if name.len() < 4 {
        if element.trim().len() == 1 {
            format!(" {name:<3}")
        } else {
            format!("{name:<4}")
        }
    } else {
        name.to_string()
    }
}

fn format_charge_display(charge: &str) -> String {
    let trimmed = charge.trim();
    if trimmed.is_empty() {
        "  ".to_string()
    } else if let Ok(val) = trimmed.parse::<i32>() {
        if val == 0 {
            "  ".to_string()
        } else {
            format!("{val:+2}")
        }
    } else {
        trimmed.to_string()
    }
}
