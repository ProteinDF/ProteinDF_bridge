// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

//! Writer for PDBx/mmCIF structure files.
//!
//! Provides sequential stream-based writing (`std::io::Write`) and file saving
//! for [`AtomGroup`] hierarchies conforming to the protein schema
//! (`/model_N/chain_id/res_key/atom_key`).

use std::fs;
use std::io::Write;
use std::path::Path;

use crate::atom_group::AtomGroup;
use crate::error::{BridgeError, Result};
use crate::periodic_table::PeriodicTable;

/// Standard amino acids (20) defined in the wwPDB CCD.
pub const STANDARD_AMINO_ACIDS: &[&str] = &[
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
];

/// Standard nucleic acids (8) defined in the wwPDB CCD (4 DNA + 4 RNA).
pub const STANDARD_NUCLEIC_ACIDS: &[&str] = &["DA", "DC", "DG", "DT", "A", "C", "G", "U"];

/// Standard amino acids (20) and nucleic acids (8) defined in the wwPDB CCD.
/// Residues in this set are classified as `ATOM`; all other residues are `HETATM`.
pub const STANDARD_RESIDUES: &[&str] = &[
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "DA", "DC", "DG", "DT", "A", "C", "G", "U",
];

/// Checks if a residue name corresponds to a standard amino acid.
pub fn is_standard_amino_acid(res_name: &str) -> bool {
    STANDARD_AMINO_ACIDS.contains(&res_name)
}

/// Checks if a residue name corresponds to a standard nucleic acid.
pub fn is_standard_nucleic_acid(res_name: &str) -> bool {
    STANDARD_NUCLEIC_ACIDS.contains(&res_name)
}

/// Checks if a residue name corresponds to a standard amino acid or nucleic acid.
pub fn is_standard_residue(res_name: &str) -> bool {
    STANDARD_RESIDUES.contains(&res_name)
}

/// Options for writing an mmCIF structure.
#[derive(Debug, Clone, PartialEq)]
pub struct MmcifWriteOptions {
    /// Data block name (written as `data_<data_block_name>`).
    /// Defaults to `"structure"`. Must not be empty, and must not contain whitespace,
    /// control characters, quotes, semicolons, or hashes.
    pub data_block_name: String,
    /// If true, atom charges are written to `_atom_site.B_iso_or_equiv` (formatted
    /// to 4 decimal places, e.g. `-0.8341`, to preserve standard partial charge precision
    /// from RESP and force fields like AMBER/CHARMM, which mmCIF allows without fixed column limits)
    /// and `_atom_site.pdbx_formal_charge` is set to `?`.
    /// If false, formal charges are written to `_atom_site.pdbx_formal_charge`
    /// (which must be an integer, otherwise an error is returned) and
    /// `_atom_site.B_iso_or_equiv` is set to `0.00`.
    /// Defaults to false.
    pub charge_to_b_factor: bool,
}

impl Default for MmcifWriteOptions {
    fn default() -> Self {
        Self {
            data_block_name: "structure".to_string(),
            charge_to_b_factor: false,
        }
    }
}

/// CIF quotation style for a data value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CifQuoteKind {
    /// Value does not need quotes.
    None,
    /// Value is enclosed in single quotes `'...'`.
    Single,
    /// Value is enclosed in double quotes `"..."` (e.g. contains single quote like `O5'`).
    Double,
}

/// Determines the quotation style needed for a CIF data value according to mmCIF syntax rules.
///
/// Returns an error if the value contains both single and double quotes.
pub fn determine_cif_quote_kind(val: &str) -> Result<CifQuoteKind> {
    let has_single = val.contains('\'');
    let has_double = val.contains('"');
    if has_single && has_double {
        return Err(BridgeError::input_error(
            "cif_value",
            format!(
                "value contains both single and double quotes, which cannot be represented safely in mmCIF: {val:?}"
            ),
        ));
    }

    if has_single {
        return Ok(CifQuoteKind::Double);
    }

    let needs_quote = val.is_empty()
        || val.chars().any(|c| c.is_whitespace())
        || val.starts_with('_')
        || val.starts_with('#')
        || val.starts_with('$')
        || val.starts_with('\'')
        || val.starts_with('"')
        || val.starts_with('[')
        || val.starts_with(']')
        || val.starts_with(';')
        || val == "."
        || val == "?"
        || {
            let lower = val.to_ascii_lowercase();
            lower.starts_with("data_")
                || lower.starts_with("save_")
                || lower == "loop_"
                || lower == "global_"
                || lower == "stop_"
        };

    if needs_quote {
        Ok(CifQuoteKind::Single)
    } else {
        Ok(CifQuoteKind::None)
    }
}

/// Quotes a CIF string value according to mmCIF syntax rules.
///
/// Returns an error if the value contains both single and double quotes.
/// Encloses the string in double quotes if it contains a single quote.
/// Encloses the string in single quotes if it contains whitespace, is empty,
/// starts with a reserved CIF prefix character (`_#$'\"[];`), matches CIF special values (`.`, `?`),
/// or matches reserved keywords (`data_*`, `save_*`, `loop_`, `global_`, `stop_`).
pub fn quote_cif_value(val: &str) -> Result<String> {
    match determine_cif_quote_kind(val)? {
        CifQuoteKind::None => Ok(val.to_string()),
        CifQuoteKind::Single => Ok(format!("'{val}'")),
        CifQuoteKind::Double => Ok(format!("\"{val}\"")),
    }
}

pub(crate) fn write_cif_value(w: &mut impl Write, val: &str) -> Result<()> {
    match determine_cif_quote_kind(val)? {
        CifQuoteKind::None => {
            w.write_all(val.as_bytes())
                .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
        }
        CifQuoteKind::Single => {
            write!(w, "'{val}'").map_err(|e| BridgeError::input_error("write", e.to_string()))?;
        }
        CifQuoteKind::Double => {
            write!(w, "\"{val}\"").map_err(|e| BridgeError::input_error("write", e.to_string()))?;
        }
    }
    Ok(())
}

/// Validates a CIF data block name.
///
/// Data block names cannot be empty, cannot contain whitespace or control characters,
/// and cannot contain characters invalid in unquoted CIF identifiers (e.g. quotes, hash, semicolon).
pub fn validate_data_block_name(name: &str) -> Result<&str> {
    let stripped = name.strip_prefix("data_").unwrap_or(name);
    if stripped.is_empty() {
        return Err(BridgeError::input_error(
            "data_block_name",
            "data block name cannot be empty",
        ));
    }
    if stripped
        .chars()
        .any(|c| c.is_whitespace() || c.is_control())
    {
        return Err(BridgeError::input_error(
            "data_block_name",
            format!("data block name '{name}' cannot contain whitespace or control characters"),
        ));
    }
    if stripped
        .chars()
        .any(|c| c == '#' || c == '\'' || c == '"' || c == ';')
    {
        return Err(BridgeError::input_error(
            "data_block_name",
            format!("data block name '{name}' contains invalid characters ('#', ';', or quotes)"),
        ));
    }
    Ok(stripped)
}

/// Parses a residue key into an integer sequence number and an optional insertion code.
///
/// Accepts patterns such as `"1"`, `"-1"`, `"52A"`, `"-10B"`, `"+3C"`.
/// Returns an error if the key cannot be parsed into a signed integer prefix and an optional insertion code.
pub fn parse_residue_key(res_key: &str) -> Result<(i32, Option<&str>)> {
    if res_key.is_empty() {
        return Err(BridgeError::input_error(
            "residue_key",
            "residue key cannot be empty",
        ));
    }
    let bytes = res_key.as_bytes();
    let mut i = 0;
    if bytes[0] == b'+' || bytes[0] == b'-' {
        i += 1;
    }
    let num_start = i;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        i += 1;
    }
    if i == num_start {
        return Err(BridgeError::input_error(
            "residue_key",
            format!("residue key '{res_key}' does not start with a valid integer"),
        ));
    }
    let seq_str = &res_key[..i];
    let seq: i32 = seq_str.parse().map_err(|e| {
        BridgeError::input_error(
            "residue_key",
            format!("residue key '{res_key}' integer part '{seq_str}' is invalid: {e}"),
        )
    })?;
    let ins_code = if i < bytes.len() {
        Some(&res_key[i..])
    } else {
        None
    };
    Ok((seq, ins_code))
}

/// Parses a model key into a 1-based model number.
///
/// Expected format: `"model_N"`, where N is a positive integer >= 1.
pub fn parse_model_key(model_key: &str) -> Result<usize> {
    let num_str = model_key.strip_prefix("model_").ok_or_else(|| {
        BridgeError::input_error(
            "model_key",
            format!("model key '{model_key}' does not start with 'model_'"),
        )
    })?;
    let num: usize = num_str.parse().map_err(|e| {
        BridgeError::input_error(
            "model_key",
            format!("model key '{model_key}' does not contain a valid model integer: {e}"),
        )
    })?;
    if num == 0 {
        return Err(BridgeError::input_error(
            "model_key",
            format!("model key '{model_key}' has invalid model number 0 (must be >= 1)"),
        ));
    }
    Ok(num)
}

/// Represents an exported connection record in `_struct_conn`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StructConnExport {
    pub id: String,
    pub conn_type_id: String,
    pub pdbx_value_order: String,
    pub ptnr1_label_asym_id: String,
    pub ptnr1_label_comp_id: String,
    pub ptnr1_label_seq_id: String,
    pub ptnr1_label_atom_id: String,
    pub pdbx_ptnr1_ins_code: Option<String>,
    pub ptnr1_auth_asym_id: String,
    pub ptnr1_auth_seq_id: i32,
    pub ptnr2_label_asym_id: String,
    pub ptnr2_label_comp_id: String,
    pub ptnr2_label_seq_id: String,
    pub ptnr2_label_atom_id: String,
    pub pdbx_ptnr2_ins_code: Option<String>,
    pub ptnr2_auth_asym_id: String,
    pub ptnr2_auth_seq_id: i32,
}

/// Resolves chain key, residue key, atom name, and residue name from a hierarchical atom path.
fn resolve_atom_path_components(
    root: &AtomGroup,
    path: &str,
) -> Result<(String, String, String, String)> {
    let parts = AtomGroup::divide_path(path);
    if parts.len() < 4 {
        return Err(BridgeError::input_error(
            "struct_conn.atom_path",
            format!(
                "atom path '{path}' does not conform to protein schema (/model_N/chain/res/atom)"
            ),
        ));
    }
    let model_key = &parts[0];
    let chain_key = &parts[1];
    let res_key = &parts[2];
    let atom_key = &parts[3];

    let model = root.get_group(model_key).ok_or_else(|| {
        BridgeError::input_error(
            "struct_conn.atom_path",
            format!("model '{model_key}' not found for path '{path}'"),
        )
    })?;
    let chain = model.get_group(chain_key).ok_or_else(|| {
        BridgeError::input_error(
            "struct_conn.atom_path",
            format!("chain '{chain_key}' not found for path '{path}'"),
        )
    })?;
    let residue = chain.get_group(res_key).ok_or_else(|| {
        BridgeError::input_error(
            "struct_conn.atom_path",
            format!("residue '{res_key}' not found for path '{path}'"),
        )
    })?;
    let atom = residue.get_atom(atom_key).ok_or_else(|| {
        BridgeError::input_error(
            "struct_conn.atom_path",
            format!("atom '{atom_key}' not found for path '{path}'"),
        )
    })?;

    Ok((
        chain_key.clone(),
        res_key.clone(),
        atom.name.clone(),
        residue.name.clone(),
    ))
}

/// Collects and validates inter-residue bonds from the first model for mmCIF `_struct_conn` export.
///
/// **Known Limitations**:
/// - `_struct_conn` records do not support per-model connectivity in mmCIF. Only bonds
///   from the first model are exported (all models are assumed to share identical bond topology).
/// - Classification into metal coordination (`metalc`) is out of scope. Covalent heuristic bonds
///   between metal ions and coordinating residues are exported as `covale`.
pub fn collect_struct_conn_records(ag: &AtomGroup) -> Result<Vec<StructConnExport>> {
    let Some((_first_model_key, first_model)) = ag.groups().next() else {
        return Ok(Vec::new());
    };

    let all_bonds = first_model.get_bond_list_ref();
    let mut disulf_records = Vec::new();
    let mut covale_records = Vec::new();
    let mut seen_pairs = std::collections::HashSet::new();

    for bond in all_bonds {
        let (c1, r1_key, a1_name, r1_name) = resolve_atom_path_components(ag, &bond.atom1_path)?;
        let (c2, r2_key, a2_name, r2_name) = resolve_atom_path_components(ag, &bond.atom2_path)?;

        // Skip intra-residue bonds
        if c1 == c2 && r1_key == r2_key {
            continue;
        }

        // Avoid duplicate bonds in undirected graph
        let pair_key = if bond.atom1_path <= bond.atom2_path {
            (bond.atom1_path.clone(), bond.atom2_path.clone())
        } else {
            (bond.atom2_path.clone(), bond.atom1_path.clone())
        };
        if !seen_pairs.insert(pair_key) {
            continue;
        }

        // Exclude standard peptide bonds (C of standard amino acid to N of standard amino acid)
        let is_peptide = is_standard_amino_acid(&r1_name)
            && is_standard_amino_acid(&r2_name)
            && ((a1_name == "C" && a2_name == "N") || (a1_name == "N" && a2_name == "C"));
        if is_peptide {
            continue;
        }

        // Exclude standard nucleic acid backbone bonds (O3' of standard nucleic acid to P of standard nucleic acid)
        let is_nucleic_backbone = is_standard_nucleic_acid(&r1_name)
            && is_standard_nucleic_acid(&r2_name)
            && ((a1_name == "O3'" && a2_name == "P") || (a1_name == "P" && a2_name == "O3'"));
        if is_nucleic_backbone {
            continue;
        }

        let is_disulf = r1_name == "CYS" && a1_name == "SG" && r2_name == "CYS" && a2_name == "SG";
        let conn_type_id = if is_disulf { "disulf" } else { "covale" };

        let (auth_seq_1, ins_1) = parse_residue_key(&r1_key)?;
        let (auth_seq_2, ins_2) = parse_residue_key(&r2_key)?;

        let label_seq_1 = if is_standard_residue(&r1_name) {
            auth_seq_1.to_string()
        } else {
            ".".to_string()
        };
        let label_seq_2 = if is_standard_residue(&r2_name) {
            auth_seq_2.to_string()
        } else {
            ".".to_string()
        };

        let label_asym_1 = if c1 == "_" { "." } else { c1.as_str() };
        let label_asym_2 = if c2 == "_" { "." } else { c2.as_str() };
        let auth_asym_1 = label_asym_1;
        let auth_asym_2 = label_asym_2;

        let value_order = match bond.order {
            1 => "sing",
            2 => "doub",
            3 => "trip",
            4 => "quad",
            _ => "sing",
        };

        // Validate CIF quote kinds for all textual fields
        if label_asym_1 != "." {
            determine_cif_quote_kind(label_asym_1)?;
        }
        determine_cif_quote_kind(&r1_name)?;
        determine_cif_quote_kind(&a1_name)?;
        if let Some(ins) = ins_1 {
            determine_cif_quote_kind(ins)?;
        }

        if label_asym_2 != "." {
            determine_cif_quote_kind(label_asym_2)?;
        }
        determine_cif_quote_kind(&r2_name)?;
        determine_cif_quote_kind(&a2_name)?;
        if let Some(ins) = ins_2 {
            determine_cif_quote_kind(ins)?;
        }

        let rec = StructConnExport {
            id: String::new(),
            conn_type_id: conn_type_id.to_string(),
            pdbx_value_order: value_order.to_string(),
            ptnr1_label_asym_id: label_asym_1.to_string(),
            ptnr1_label_comp_id: r1_name.to_string(),
            ptnr1_label_seq_id: label_seq_1,
            ptnr1_label_atom_id: a1_name.to_string(),
            pdbx_ptnr1_ins_code: ins_1.map(|s| s.to_string()),
            ptnr1_auth_asym_id: auth_asym_1.to_string(),
            ptnr1_auth_seq_id: auth_seq_1,
            ptnr2_label_asym_id: label_asym_2.to_string(),
            ptnr2_label_comp_id: r2_name.to_string(),
            ptnr2_label_seq_id: label_seq_2,
            ptnr2_label_atom_id: a2_name.to_string(),
            pdbx_ptnr2_ins_code: ins_2.map(|s| s.to_string()),
            ptnr2_auth_asym_id: auth_asym_2.to_string(),
            ptnr2_auth_seq_id: auth_seq_2,
        };

        if is_disulf {
            disulf_records.push(rec);
        } else {
            covale_records.push(rec);
        }
    }

    // Assign IDs: disulf1, disulf2, ... then covale1, covale2, ...
    for (i, rec) in disulf_records.iter_mut().enumerate() {
        rec.id = format!("disulf{}", i + 1);
    }
    for (i, rec) in covale_records.iter_mut().enumerate() {
        rec.id = format!("covale{}", i + 1);
    }

    let mut result = disulf_records;
    result.extend(covale_records);
    Ok(result)
}

/// Validates the entire [`AtomGroup`] tree for mmCIF writing prior to emitting any bytes.
///
/// Guarantees that [`write_structure`] will not fail mid-stream due to invalid keys,
/// un-quotable values, unsupported elements, non-integer charges, or schema violations.
pub fn validate_for_mmcif_write(ag: &AtomGroup, opts: &MmcifWriteOptions) -> Result<()> {
    // 1. Schema check
    let violations = ag.validate_schema();
    if !violations.is_empty() {
        let msgs: Vec<String> = violations.iter().map(|v| v.to_string()).collect();
        return Err(BridgeError::input_error(
            "AtomGroup::validate_schema",
            format!("protein schema validation failed: {}", msgs.join("; ")),
        ));
    }

    // 2. Data block name
    validate_data_block_name(&opts.data_block_name)?;

    // 3. Hierarchy items validation
    for (model_key, model) in ag.groups() {
        parse_model_key(model_key)?;

        for (chain_key, chain) in model.groups() {
            if chain_key != "_" {
                determine_cif_quote_kind(chain_key)?;
            }

            for (res_key, residue) in chain.groups() {
                let (_, ins_code) = parse_residue_key(res_key)?;
                determine_cif_quote_kind(&residue.name)?;
                if let Some(ins) = ins_code {
                    determine_cif_quote_kind(ins)?;
                }

                for (_atom_key, atom) in residue.atoms() {
                    determine_cif_quote_kind(&atom.name)?;
                    PeriodicTable::get_symbol(atom.atomic_number())?;

                    if opts.charge_to_b_factor {
                        if !atom.charge.is_finite() {
                            return Err(BridgeError::input_error(
                                "atom.charge",
                                format!(
                                    "non-finite charge ({}) for atom '{}' at residue '{}' in chain '{}'",
                                    atom.charge, atom.name, res_key, chain_key
                                ),
                            ));
                        }
                    } else if !atom.charge.is_finite()
                        || (atom.charge.round() - atom.charge).abs() >= 1e-6
                    {
                        return Err(BridgeError::input_error(
                            "atom.charge",
                            format!(
                                "non-integer formal charge ({}) found for atom '{}' at residue '{}' in chain '{}'; formal charge in mmCIF must be an integer. Use MmcifWriteOptions {{ charge_to_b_factor: true, .. }} to write partial charges into B_iso_or_equiv instead",
                                atom.charge, atom.name, res_key, chain_key
                            ),
                        ));
                    }
                }
            }
        }
    }

    // 4. Inter-residue bond validation for _struct_conn
    collect_struct_conn_records(ag)?;

    Ok(())
}

/// Writes an [`AtomGroup`] hierarchy to an mmCIF stream (`std::io::Write`).
///
/// Performs full pre-validation before emitting any bytes. If any property violates
/// mmCIF specifications or protein schema, returns an error without writing anything to `w`.
pub fn write_structure(ag: &AtomGroup, w: &mut impl Write, opts: &MmcifWriteOptions) -> Result<()> {
    // 1. Full pre-validation (atomic error guarantee: no partial output written if error)
    validate_for_mmcif_write(ag, opts)?;
    let block_name = validate_data_block_name(&opts.data_block_name)?;
    let struct_conns = collect_struct_conn_records(ag)?;

    // 2. Data block header
    writeln!(w, "data_{block_name}")
        .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    writeln!(w, "#").map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    writeln!(w, "loop_").map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    writeln!(w, "_atom_site.group_PDB")
        .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    writeln!(w, "_atom_site.id").map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    writeln!(w, "_atom_site.type_symbol")
        .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    writeln!(w, "_atom_site.label_atom_id")
        .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    writeln!(w, "_atom_site.label_alt_id")
        .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    writeln!(w, "_atom_site.label_comp_id")
        .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    writeln!(w, "_atom_site.label_asym_id")
        .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    writeln!(w, "_atom_site.label_entity_id")
        .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    writeln!(w, "_atom_site.label_seq_id")
        .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    writeln!(w, "_atom_site.pdbx_PDB_ins_code")
        .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    writeln!(w, "_atom_site.Cartn_x")
        .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    writeln!(w, "_atom_site.Cartn_y")
        .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    writeln!(w, "_atom_site.Cartn_z")
        .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    writeln!(w, "_atom_site.occupancy")
        .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    writeln!(w, "_atom_site.B_iso_or_equiv")
        .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    writeln!(w, "_atom_site.pdbx_formal_charge")
        .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    writeln!(w, "_atom_site.auth_seq_id")
        .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    writeln!(w, "_atom_site.auth_comp_id")
        .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    writeln!(w, "_atom_site.auth_asym_id")
        .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    writeln!(w, "_atom_site.auth_atom_id")
        .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    writeln!(w, "_atom_site.pdbx_PDB_model_num")
        .map_err(|e| BridgeError::input_error("write", e.to_string()))?;

    let mut atom_serial: usize = 0;

    for (model_key, model) in ag.groups() {
        let model_num = parse_model_key(model_key)?;

        for (chain_key, chain) in model.groups() {
            let chain_is_empty = chain_key == "_";

            for (res_key, residue) in chain.groups() {
                let (auth_seq_id, ins_code) = parse_residue_key(res_key)?;
                let is_std = is_standard_residue(&residue.name);
                let group_pdb = if is_std { "ATOM" } else { "HETATM" };
                let label_seq_id = if is_std {
                    auth_seq_id.to_string()
                } else {
                    ".".to_string()
                };

                for (_atom_key, atom) in residue.atoms() {
                    atom_serial += 1;

                    let type_symbol = PeriodicTable::get_symbol(atom.atomic_number())?;

                    let (b_iso_str, charge_str) = if opts.charge_to_b_factor {
                        // 4 decimal places to preserve partial charge precision
                        (format!("{:.4}", atom.charge), "?".to_string())
                    } else {
                        let charge_int = atom.charge.round() as i32;
                        ("0.00".to_string(), charge_int.to_string())
                    };

                    let x = if atom.xyz.x.abs() < 1e-6 {
                        0.0
                    } else {
                        atom.xyz.x
                    };
                    let y = if atom.xyz.y.abs() < 1e-6 {
                        0.0
                    } else {
                        atom.xyz.y
                    };
                    let z = if atom.xyz.z.abs() < 1e-6 {
                        0.0
                    } else {
                        atom.xyz.z
                    };

                    write!(w, "{} {} {} ", group_pdb, atom_serial, type_symbol)
                        .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
                    write_cif_value(w, &atom.name)?;
                    write!(w, " . ")
                        .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
                    write_cif_value(w, &residue.name)?;
                    write!(w, " ").map_err(|e| BridgeError::input_error("write", e.to_string()))?;

                    if chain_is_empty {
                        write!(w, ". ? ")
                            .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
                    } else {
                        write_cif_value(w, chain_key)?;
                        write!(w, " ? ")
                            .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
                    }

                    write!(w, "{} ", label_seq_id)
                        .map_err(|e| BridgeError::input_error("write", e.to_string()))?;

                    if let Some(ins) = ins_code {
                        write_cif_value(w, ins)?;
                        write!(w, " ")
                            .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
                    } else {
                        write!(w, "? ")
                            .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
                    }

                    write!(
                        w,
                        "{:.3} {:.3} {:.3} 1.00 {} {} {} ",
                        x, y, z, b_iso_str, charge_str, auth_seq_id
                    )
                    .map_err(|e| BridgeError::input_error("write", e.to_string()))?;

                    write_cif_value(w, &residue.name)?;
                    write!(w, " ").map_err(|e| BridgeError::input_error("write", e.to_string()))?;

                    if chain_is_empty {
                        write!(w, ". ")
                            .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
                    } else {
                        write_cif_value(w, chain_key)?;
                        write!(w, " ")
                            .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
                    }

                    write_cif_value(w, &atom.name)?;
                    writeln!(w, " {}", model_num)
                        .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
                }
            }
        }
    }

    writeln!(w, "#").map_err(|e| BridgeError::input_error("write", e.to_string()))?;

    // 3. _struct_conn loop (if any inter-residue bonds exist)
    if !struct_conns.is_empty() {
        writeln!(w, "loop_").map_err(|e| BridgeError::input_error("write", e.to_string()))?;
        writeln!(w, "_struct_conn.id")
            .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
        writeln!(w, "_struct_conn.conn_type_id")
            .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
        writeln!(w, "_struct_conn.pdbx_value_order")
            .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
        writeln!(w, "_struct_conn.ptnr1_label_asym_id")
            .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
        writeln!(w, "_struct_conn.ptnr1_label_comp_id")
            .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
        writeln!(w, "_struct_conn.ptnr1_label_seq_id")
            .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
        writeln!(w, "_struct_conn.ptnr1_label_atom_id")
            .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
        writeln!(w, "_struct_conn.pdbx_ptnr1_PDB_ins_code")
            .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
        writeln!(w, "_struct_conn.ptnr1_auth_asym_id")
            .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
        writeln!(w, "_struct_conn.ptnr1_auth_seq_id")
            .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
        writeln!(w, "_struct_conn.ptnr2_label_asym_id")
            .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
        writeln!(w, "_struct_conn.ptnr2_label_comp_id")
            .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
        writeln!(w, "_struct_conn.ptnr2_label_seq_id")
            .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
        writeln!(w, "_struct_conn.ptnr2_label_atom_id")
            .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
        writeln!(w, "_struct_conn.pdbx_ptnr2_PDB_ins_code")
            .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
        writeln!(w, "_struct_conn.ptnr2_auth_asym_id")
            .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
        writeln!(w, "_struct_conn.ptnr2_auth_seq_id")
            .map_err(|e| BridgeError::input_error("write", e.to_string()))?;

        for conn in &struct_conns {
            write!(
                w,
                "{} {} {} ",
                conn.id, conn.conn_type_id, conn.pdbx_value_order
            )
            .map_err(|e| BridgeError::input_error("write", e.to_string()))?;

            if conn.ptnr1_label_asym_id == "." {
                write!(w, ". ").map_err(|e| BridgeError::input_error("write", e.to_string()))?;
            } else {
                write_cif_value(w, &conn.ptnr1_label_asym_id)?;
                write!(w, " ").map_err(|e| BridgeError::input_error("write", e.to_string()))?;
            }
            write_cif_value(w, &conn.ptnr1_label_comp_id)?;
            write!(w, " {} ", conn.ptnr1_label_seq_id)
                .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
            write_cif_value(w, &conn.ptnr1_label_atom_id)?;
            write!(w, " ").map_err(|e| BridgeError::input_error("write", e.to_string()))?;

            if let Some(ins) = &conn.pdbx_ptnr1_ins_code {
                write_cif_value(w, ins)?;
                write!(w, " ").map_err(|e| BridgeError::input_error("write", e.to_string()))?;
            } else {
                write!(w, "? ").map_err(|e| BridgeError::input_error("write", e.to_string()))?;
            }

            if conn.ptnr1_auth_asym_id == "." {
                write!(w, ". ").map_err(|e| BridgeError::input_error("write", e.to_string()))?;
            } else {
                write_cif_value(w, &conn.ptnr1_auth_asym_id)?;
                write!(w, " ").map_err(|e| BridgeError::input_error("write", e.to_string()))?;
            }
            write!(w, "{} ", conn.ptnr1_auth_seq_id)
                .map_err(|e| BridgeError::input_error("write", e.to_string()))?;

            if conn.ptnr2_label_asym_id == "." {
                write!(w, ". ").map_err(|e| BridgeError::input_error("write", e.to_string()))?;
            } else {
                write_cif_value(w, &conn.ptnr2_label_asym_id)?;
                write!(w, " ").map_err(|e| BridgeError::input_error("write", e.to_string()))?;
            }
            write_cif_value(w, &conn.ptnr2_label_comp_id)?;
            write!(w, " {} ", conn.ptnr2_label_seq_id)
                .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
            write_cif_value(w, &conn.ptnr2_label_atom_id)?;
            write!(w, " ").map_err(|e| BridgeError::input_error("write", e.to_string()))?;

            if let Some(ins) = &conn.pdbx_ptnr2_ins_code {
                write_cif_value(w, ins)?;
                write!(w, " ").map_err(|e| BridgeError::input_error("write", e.to_string()))?;
            } else {
                write!(w, "? ").map_err(|e| BridgeError::input_error("write", e.to_string()))?;
            }

            if conn.ptnr2_auth_asym_id == "." {
                write!(w, ". ").map_err(|e| BridgeError::input_error("write", e.to_string()))?;
            } else {
                write_cif_value(w, &conn.ptnr2_auth_asym_id)?;
                write!(w, " ").map_err(|e| BridgeError::input_error("write", e.to_string()))?;
            }
            writeln!(w, "{}", conn.ptnr2_auth_seq_id)
                .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
        }
        writeln!(w, "#").map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    }

    Ok(())
}

/// Saves an [`AtomGroup`] hierarchy to an mmCIF file at the given path.
///
/// Uses an atomic write strategy:
/// 1. Validates the structure and data block name up-front before touching disk.
/// 2. Writes to a temporary file in the same directory using buffered I/O.
/// 3. Atomically renames the temporary file to the destination upon success.
///
/// If an error occurs at any point, the destination file is untouched, and any
/// temporary file is cleaned up.
pub fn save_structure(
    ag: &AtomGroup,
    path: impl AsRef<Path>,
    opts: &MmcifWriteOptions,
) -> Result<()> {
    validate_for_mmcif_write(ag, opts)?;

    let path_ref = path.as_ref();
    let parent = path_ref.parent().unwrap_or_else(|| Path::new("."));

    let file_name = path_ref
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("structure.cif");
    let temp_name = format!(
        ".{file_name}.tmp_{}_{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos()
    );
    let temp_path = parent.join(temp_name);

    struct TempFileGuard<'a>(&'a Path, bool);
    impl<'a> Drop for TempFileGuard<'a> {
        fn drop(&mut self) {
            if !self.1 {
                let _ = fs::remove_file(self.0);
            }
        }
    }
    let mut guard = TempFileGuard(&temp_path, false);

    let file = fs::File::create(&temp_path)
        .map_err(|e| BridgeError::input_error(temp_path.display().to_string(), e.to_string()))?;
    let mut writer = std::io::BufWriter::new(file);
    write_structure(ag, &mut writer, opts)?;
    writer
        .flush()
        .map_err(|e| BridgeError::input_error(temp_path.display().to_string(), e.to_string()))?;

    fs::rename(&temp_path, path_ref)
        .map_err(|e| BridgeError::input_error(path_ref.display().to_string(), e.to_string()))?;
    guard.1 = true;

    Ok(())
}
