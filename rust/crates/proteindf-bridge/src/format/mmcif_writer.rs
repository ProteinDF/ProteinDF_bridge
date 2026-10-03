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
    /// Defaults to `"structure"`.
    pub data_block_name: String,
    /// If true, atom charges are written to `_atom_site.B_iso_or_equiv` and
    /// `_atom_site.pdbx_formal_charge` is set to `?`.
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

/// Quotes a CIF string value according to mmCIF syntax rules.
///
/// Returns an error if the value contains both single and double quotes.
/// Encloses the string in double quotes if it contains a single quote.
/// Encloses the string in single quotes if it contains whitespace, is empty,
/// starts with a reserved CIF prefix character (`_#$'\"[];`), matches CIF special values (`.`, `?`),
/// or matches reserved keywords (`data_*`, `save_*`, `loop_`, `global_`, `stop_`).
pub fn quote_cif_value(val: &str) -> Result<String> {
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

    if has_single {
        Ok(format!("\"{val}\""))
    } else if needs_quote {
        Ok(format!("'{val}'"))
    } else {
        Ok(val.to_string())
    }
}

pub(crate) fn write_cif_value(w: &mut impl Write, val: &str) -> Result<()> {
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

    if has_single {
        write!(w, "\"{val}\"").map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    } else if needs_quote {
        write!(w, "'{val}'").map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    } else {
        w.write_all(val.as_bytes())
            .map_err(|e| BridgeError::input_error("write", e.to_string()))?;
    }
    Ok(())
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

/// Writes an [`AtomGroup`] hierarchy to an mmCIF stream (`std::io::Write`).
///
/// Validates that the hierarchy conforms to the protein schema (`/model_N/chain_id/res_key/atom_key`).
/// Returns an error if schema validation fails or any atom/residue/model properties violate CIF requirements.
pub fn write_structure(ag: &AtomGroup, w: &mut impl Write, opts: &MmcifWriteOptions) -> Result<()> {
    // 1. Validate schema
    let violations = ag.validate_schema();
    if !violations.is_empty() {
        let msgs: Vec<String> = violations.iter().map(|v| v.to_string()).collect();
        return Err(BridgeError::input_error(
            "AtomGroup::validate_schema",
            format!("protein schema validation failed: {}", msgs.join("; ")),
        ));
    }

    // 2. Data block header
    let block_name = if let Some(stripped) = opts.data_block_name.strip_prefix("data_") {
        stripped
    } else {
        opts.data_block_name.as_str()
    };
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
                        (format!("{:.2}", atom.charge), "?".to_string())
                    } else {
                        if !atom.charge.is_finite()
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
/// Uses buffered writing for high performance on large structures.
pub fn save_structure(
    ag: &AtomGroup,
    path: impl AsRef<Path>,
    opts: &MmcifWriteOptions,
) -> Result<()> {
    let path_ref = path.as_ref();
    let file = fs::File::create(path_ref)
        .map_err(|e| BridgeError::input_error(path_ref.display().to_string(), e.to_string()))?;
    let mut writer = std::io::BufWriter::new(file);
    write_structure(ag, &mut writer, opts)?;
    writer
        .flush()
        .map_err(|e| BridgeError::input_error(path_ref.display().to_string(), e.to_string()))?;
    Ok(())
}
