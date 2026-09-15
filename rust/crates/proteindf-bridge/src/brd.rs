// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

use rmpv::Value;

use crate::atom::Atom;
use crate::atom_group::{AtomGroup, BondRecord};
use crate::error::{BridgeError, Result};
use crate::position::Position;

/// Magic bytes identifying a YUI-compatible format.
pub const YUI_MAGIC: &[u8; 4] = b"YUI\0";
/// Current YUI format version.
pub const YUI_VERSION: u8 = 1;
/// Uncompressed payload flag.
pub const YUI_COMPRESSION_NONE: u8 = 0;
/// Zstandard compressed payload flag.
pub const YUI_COMPRESSION_ZSTD: u8 = 1;

/// Helper to convert a MessagePack Map key into a `&str`.
pub fn value_key_to_str(val: &Value) -> Option<&str> {
    match val {
        Value::String(s) => s.as_str(),
        Value::Binary(b) => std::str::from_utf8(b).ok(),
        _ => None,
    }
}

/// Helper to extract an `f64` from a Value (supports float and integer representations).
pub fn value_as_f64(val: &Value) -> Option<f64> {
    match val {
        Value::F64(v) => Some(*v),
        Value::F32(v) => Some(*v as f64),
        Value::Integer(i) => i.as_f64(),
        _ => None,
    }
}

/// Helper to extract a `usize` from a Value.
pub fn value_as_usize(val: &Value) -> Option<usize> {
    match val {
        Value::Integer(i) => i.as_u64().map(|v| v as usize),
        _ => None,
    }
}

/// Helper to extract a `Position` from a 3-element coordinate Value.
pub fn parse_position_from_value(val: &Value) -> Option<Position> {
    match val {
        Value::Array(arr) => {
            let x = arr.first().and_then(value_as_f64).unwrap_or(0.0);
            let y = arr.get(1).and_then(value_as_f64).unwrap_or(0.0);
            let z = arr.get(2).and_then(value_as_f64).unwrap_or(0.0);
            Some(Position::new(x, y, z))
        }
        _ => None,
    }
}

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
enum AlphanumChunk<'a> {
    Num(u64),
    Str(&'a str),
}

fn alphanum_chunks(s: &str) -> Vec<AlphanumChunk<'_>> {
    let mut chunks = Vec::new();
    let mut start = 0;
    let mut in_digit = false;
    let bytes = s.as_bytes();
    for (i, &b) in bytes.iter().enumerate() {
        let is_d = b.is_ascii_digit();
        if i == 0 {
            in_digit = is_d;
        } else if is_d != in_digit {
            let slice = &s[start..i];
            if in_digit {
                let n = slice.parse::<u64>().unwrap_or(0);
                chunks.push(AlphanumChunk::Num(n));
            } else {
                chunks.push(AlphanumChunk::Str(slice));
            }
            start = i;
            in_digit = is_d;
        }
    }
    if start < s.len() {
        let slice = &s[start..];
        if in_digit {
            let n = slice.parse::<u64>().unwrap_or(0);
            chunks.push(AlphanumChunk::Num(n));
        } else {
            chunks.push(AlphanumChunk::Str(slice));
        }
    }
    chunks
}

/// Sorts a slice of string keys using natural sort order (alphanumeric sorting),
/// reproducing Python's `StrUtils.sort_nicely`.
pub fn sort_nicely<T: AsRef<str>>(items: &mut [T]) {
    items.sort_by(|a, b| {
        let ca = alphanum_chunks(a.as_ref());
        let cb = alphanum_chunks(b.as_ref());
        ca.cmp(&cb)
    });
}

/// Returns the MessagePack raw data representation of an `Atom`.
/// Schema: `{"Z": int, "name": str, "Q": float, "xyz": [x,y,z], "force": [x,y,z]}`.
pub fn atom_get_raw_data(atom: &Atom) -> Value {
    let map = vec![
        (
            Value::String("Z".into()),
            Value::Integer((atom.atomic_number() as u64).into()),
        ),
        (
            Value::String("name".into()),
            Value::String(atom.name.clone().into()),
        ),
        (Value::String("Q".into()), Value::F64(atom.charge)),
        (
            Value::String("xyz".into()),
            Value::Array(vec![
                Value::F64(atom.xyz.x),
                Value::F64(atom.xyz.y),
                Value::F64(atom.xyz.z),
            ]),
        ),
        (
            Value::String("force".into()),
            Value::Array(vec![
                Value::F64(atom.force.x),
                Value::F64(atom.force.y),
                Value::F64(atom.force.z),
            ]),
        ),
    ];
    Value::Map(map)
}

/// Populates an `Atom` from a raw MessagePack dictionary Value.
/// Unknown keys are logged at debug level and ignored.
pub fn atom_set_by_raw_data(atom: &mut Atom, data: &Value) -> Result<()> {
    let entries = match data {
        Value::Map(entries) => entries,
        _ => {
            return Err(BridgeError::value_error(
                "data",
                "expected map for Atom raw data",
            ))
        }
    };
    for (k_val, v_val) in entries {
        let key = match value_key_to_str(k_val) {
            Some(k) => k,
            None => {
                log::debug!("bridge::Atom > unknown non-string key: {:?}", k_val);
                continue;
            }
        };
        match key {
            "Z" => {
                if let Some(z) = value_as_usize(v_val) {
                    atom.set_atomic_number(z);
                }
            }
            "name" => {
                if let Value::String(s) = v_val {
                    if let Some(str_val) = s.as_str() {
                        atom.name = str_val.to_string();
                    }
                }
            }
            "Q" => {
                if let Some(q) = value_as_f64(v_val) {
                    atom.charge = q;
                }
            }
            "xyz" => {
                if let Some(pos) = parse_position_from_value(v_val) {
                    atom.xyz = pos;
                }
            }
            "force" => {
                if let Some(pos) = parse_position_from_value(v_val) {
                    atom.force = pos;
                }
            }
            _ => {
                log::debug!("bridge::Atom > unknown key: {}={:?}", key, v_val);
            }
        }
    }
    Ok(())
}

/// Returns the MessagePack raw data representation of an `AtomGroup`.
/// Schema: `{"name": str, "groups": {...}, "atoms": {...}, "bonds": [...]}`.
/// Empty groups, atoms, or bonds are omitted matching Python behavior.
pub fn atomgroup_get_raw_data(group: &AtomGroup) -> Value {
    let mut map = Vec::new();

    if group.get_number_of_groups() > 0 {
        let mut grp_map = Vec::new();
        for (key, grp) in group.groups() {
            grp_map.push((
                Value::String(key.clone().into()),
                atomgroup_get_raw_data(grp),
            ));
        }
        map.push((Value::String("groups".into()), Value::Map(grp_map)));
    }

    if group.get_number_of_atoms() > 0 {
        let mut atm_map = Vec::new();
        for (key, atm) in group.atoms() {
            atm_map.push((Value::String(key.clone().into()), atom_get_raw_data(atm)));
        }
        map.push((Value::String("atoms".into()), Value::Map(atm_map)));
    }

    map.push((
        Value::String("name".into()),
        Value::String(group.name.clone().into()),
    ));

    if !group.bonds().is_empty() {
        let bonds_vec: Vec<Value> = group
            .bonds()
            .iter()
            .map(|b| {
                Value::Array(vec![
                    Value::String(b.atom1_path.clone().into()),
                    Value::String(b.atom2_path.clone().into()),
                    Value::Integer((b.order as u64).into()),
                ])
            })
            .collect();
        map.push((Value::String("bonds".into()), Value::Array(bonds_vec)));
    }

    Value::Map(map)
}

/// Populates an `AtomGroup` from a raw MessagePack dictionary Value.
/// Unknown keys (e.g. `Q=None`, `charge=0.0`) are logged as warnings and ignored.
/// Groups and atoms are inserted following `sort_nicely` order.
pub fn atomgroup_set_by_dict_data(group: &mut AtomGroup, data: &Value) -> Result<()> {
    let entries = match data {
        Value::Map(entries) => entries,
        _ => {
            return Err(BridgeError::value_error(
                "data",
                "expected map for AtomGroup dict data",
            ))
        }
    };

    let mut tmp_groups: Vec<(String, AtomGroup)> = Vec::new();
    let mut tmp_atoms: Vec<(String, Atom)> = Vec::new();
    let mut bonds: Vec<BondRecord> = Vec::new();

    for (k_val, v_val) in entries {
        let key = match value_key_to_str(k_val) {
            Some(k) => k,
            None => {
                log::warn!("AtomGroup::set_by_dict_data(): unknown key: {:?}", k_val);
                continue;
            }
        };

        match key {
            "name" => {
                if let Value::String(s) = v_val {
                    if let Some(str_val) = s.as_str() {
                        group.name = str_val.to_string();
                    }
                }
            }
            "groups" => {
                if let Value::Map(grp_entries) = v_val {
                    for (gk_val, gv_val) in grp_entries {
                        if let Some(gk) = value_key_to_str(gk_val) {
                            let mut sub_ag = AtomGroup::new();
                            atomgroup_set_by_dict_data(&mut sub_ag, gv_val)?;
                            tmp_groups.push((gk.to_string(), sub_ag));
                        }
                    }
                }
            }
            "atoms" => {
                if let Value::Map(atm_entries) = v_val {
                    for (ak_val, av_val) in atm_entries {
                        if let Some(ak) = value_key_to_str(ak_val) {
                            let mut atm = Atom::new();
                            atom_set_by_raw_data(&mut atm, av_val)?;
                            tmp_atoms.push((ak.to_string(), atm));
                        }
                    }
                }
            }
            "bonds" => {
                if let Value::Array(b_list) = v_val {
                    for b_item in b_list {
                        if let Value::Array(b_fields) = b_item {
                            let p1 = b_fields
                                .first()
                                .and_then(|v| {
                                    if let Value::String(s) = v {
                                        s.as_str()
                                    } else {
                                        None
                                    }
                                })
                                .unwrap_or("")
                                .to_string();
                            let p2 = b_fields
                                .get(1)
                                .and_then(|v| {
                                    if let Value::String(s) = v {
                                        s.as_str()
                                    } else {
                                        None
                                    }
                                })
                                .unwrap_or("")
                                .to_string();
                            let order = b_fields.get(2).and_then(value_as_usize).unwrap_or(1);
                            bonds.push(BondRecord {
                                atom1_path: p1,
                                atom2_path: p2,
                                order,
                            });
                        }
                    }
                }
            }
            _ => {
                log::warn!(
                    "AtomGroup::set_by_dict_data(): unknown key: {}={:?}",
                    key,
                    v_val
                );
            }
        }
    }

    // Sort and insert groups
    let mut grp_keys: Vec<String> = tmp_groups.iter().map(|(k, _)| k.clone()).collect();
    sort_nicely(&mut grp_keys);
    let mut grp_map: std::collections::HashMap<String, AtomGroup> =
        tmp_groups.into_iter().collect();
    for k in grp_keys {
        if let Some(g) = grp_map.remove(&k) {
            group.set_group(&k, g);
        }
    }

    // Sort and insert atoms
    let mut atm_keys: Vec<String> = tmp_atoms.iter().map(|(k, _)| k.clone()).collect();
    sort_nicely(&mut atm_keys);
    let mut atm_map: std::collections::HashMap<String, Atom> = tmp_atoms.into_iter().collect();
    for k in atm_keys {
        if let Some(a) = atm_map.remove(&k) {
            group.set_atom(&k, a);
        }
    }

    // Set bonds
    if !bonds.is_empty() {
        group.set_bonds(bonds);
    }

    group.set_path(group.path().to_string());
    Ok(())
}

// =============================================================================
// Plain MessagePack I/O (1:1 port of functions.py)
// =============================================================================

/// Loads a MessagePack file and returns the decoded Value.
pub fn load_msgpack<P: AsRef<Path>>(path: P) -> Result<Value> {
    let mut file = File::open(path)?;
    rmpv::decode::read_value(&mut file).map_err(|e| BridgeError::MsgPack(e.to_string()))
}

/// Encodes and saves a Value to a MessagePack file.
pub fn save_msgpack<P: AsRef<Path>>(val: &Value, path: P) -> Result<()> {
    let mut file = File::create(path)?;
    rmpv::encode::write_value(&mut file, val).map_err(|e| BridgeError::MsgPack(e.to_string()))
}

/// Loads a bridge (.brd) MessagePack file directly into an `AtomGroup`.
pub fn load_atomgroup<P: AsRef<Path>>(path: P) -> Result<AtomGroup> {
    let val = load_msgpack(path)?;
    let mut ag = AtomGroup::new();
    atomgroup_set_by_dict_data(&mut ag, &val)?;
    Ok(ag)
}

/// Loads a bridge (.brd) MessagePack byte slice directly into an `AtomGroup`.
pub fn load_atomgroup_from_bytes(bytes: &[u8]) -> Result<AtomGroup> {
    let val = rmpv::decode::read_value(&mut &bytes[..])
        .map_err(|e| BridgeError::MsgPack(e.to_string()))?;
    let mut ag = AtomGroup::new();
    atomgroup_set_by_dict_data(&mut ag, &val)?;
    Ok(ag)
}

/// Saves an `AtomGroup` into a plain bridge (.brd) MessagePack file.
pub fn save_atomgroup<P: AsRef<Path>>(group: &AtomGroup, path: P) -> Result<()> {
    let val = atomgroup_get_raw_data(group);
    save_msgpack(&val, path)
}

// =============================================================================
// YUI-compatible format I/O (RUST_PORT_SPEC.md §1.3)
// =============================================================================

/// Saves an `AtomGroup` with a YUI-compatible header (`[Magic: "YUI\0"(4B)] + [Version: 1(1B)] + [Compression Flag(1B)] + [Payload]`).
pub fn save_brd_yui<P: AsRef<Path>>(group: &AtomGroup, path: P, compress_zstd: bool) -> Result<()> {
    let val = atomgroup_get_raw_data(group);
    let mut raw_bytes = Vec::new();
    rmpv::encode::write_value(&mut raw_bytes, &val)
        .map_err(|e| BridgeError::MsgPack(e.to_string()))?;

    let (flag, payload) = if compress_zstd {
        let compressed =
            zstd::encode_all(&raw_bytes[..], 0).map_err(|e| BridgeError::Zstd(e.to_string()))?;
        (YUI_COMPRESSION_ZSTD, compressed)
    } else {
        (YUI_COMPRESSION_NONE, raw_bytes)
    };

    let mut file = File::create(path)?;
    file.write_all(YUI_MAGIC)?;
    file.write_all(&[YUI_VERSION, flag])?;
    file.write_all(&payload)?;
    Ok(())
}

/// Loads an `AtomGroup` from a file with a YUI-compatible header.
pub fn load_brd_yui<P: AsRef<Path>>(path: P) -> Result<AtomGroup> {
    let mut file = File::open(path)?;
    let mut header = [0u8; 6];
    file.read_exact(&mut header)?;

    if &header[0..4] != YUI_MAGIC {
        return Err(BridgeError::MsgPack(
            "Invalid YUI magic header in brd file".to_string(),
        ));
    }

    let version = header[4];
    if version != YUI_VERSION {
        return Err(BridgeError::MsgPack(format!(
            "Unsupported YUI version: {version}"
        )));
    }

    let comp_flag = header[5];
    let mut payload = Vec::new();
    file.read_to_end(&mut payload)?;

    let msgpack_bytes = match comp_flag {
        YUI_COMPRESSION_NONE => payload,
        YUI_COMPRESSION_ZSTD => {
            zstd::decode_all(&payload[..]).map_err(|e| BridgeError::Zstd(e.to_string()))?
        }
        other => {
            return Err(BridgeError::MsgPack(format!(
                "Unsupported YUI compression flag: {other}"
            )));
        }
    };

    let val = rmpv::decode::read_value(&mut &msgpack_bytes[..])
        .map_err(|e| BridgeError::MsgPack(e.to_string()))?;
    let mut group = AtomGroup::new();
    atomgroup_set_by_dict_data(&mut group, &val)?;
    Ok(group)
}
