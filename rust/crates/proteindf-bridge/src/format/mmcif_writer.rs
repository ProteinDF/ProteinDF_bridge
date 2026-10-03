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

/// Standard amino acids (20) and nucleic acids (8) defined in the wwPDB CCD.
/// Residues in this set are classified as `ATOM`; all other residues are `HETATM`.
pub const STANDARD_RESIDUES: &[&str] = &[
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "DA", "DC", "DG", "DT", "A", "C", "G", "U",
];

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

/// Checks if a residue name corresponds to a standard amino acid or nucleic acid.
pub fn is_standard_residue(res_name: &str) -> bool {
    STANDARD_RESIDUES.contains(&res_name)
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
