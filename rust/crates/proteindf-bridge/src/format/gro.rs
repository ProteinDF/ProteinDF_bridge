// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::fmt;
use std::fs;
use std::path::Path;

use crate::atom::Atom;
use crate::atom_group::AtomGroup;
use crate::error::{BridgeError, Result};
use crate::periodic_table::PeriodicTable;
use crate::position::Position;

pub const MIN_RECORD_LEN: usize = 44;
pub const NUMBER_WRAP: usize = 10000;
pub const BOX_VECTORS_PAD_LEN: usize = 6;

/// A single atom record in GROMACS .gro format.
#[derive(Debug, Clone, PartialEq)]
pub struct GroAtom {
    pub residue_number: usize,
    pub residue_name: String,
    pub atom_name: String,
    pub atom_number: usize,
    /// Coordinates in nanometers (nm).
    pub position: Position,
    /// Velocities in nm/ps.
    pub velocity: [f64; 3],
}

/// GROMACS .gro format reader/writer corresponding to `proteindf_bridge.gro.SimpleGro`.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SimpleGro {
    pub title: String,
    pub num_of_atoms: usize,
    pub atoms: Vec<GroAtom>,
    pub box_vectors: Vec<f64>,
}

impl SimpleGro {
    /// Creates a new empty `SimpleGro` instance.
    pub fn new() -> Self {
        Self::default()
    }

    /// Loads a .gro file from the given path.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        let mut gro = Self::new();
        gro.load(path)?;
        Ok(gro)
    }

    /// Loads a .gro file.
    pub fn load(&mut self, file_path: impl AsRef<Path>) -> Result<()> {
        let path = file_path.as_ref();
        let content = fs::read_to_string(path).map_err(|e| {
            BridgeError::input_error(
                path.display().to_string(),
                format!("failed to read file: {}", e),
            )
        })?;
        self.parse_str(&content)
    }

    /// Parses .gro content from a string.
    pub fn parse_str(&mut self, content: &str) -> Result<()> {
        let mut lines = content.lines();

        // Line 1: title
        let title_line = lines
            .next()
            .ok_or_else(|| BridgeError::input_error("GRO", "empty content"))?;
        self.title = title_line.trim_end().to_string();

        // Line 2: number of atoms
        let count_line = lines
            .next()
            .ok_or_else(|| BridgeError::input_error("GRO", "missing atom count line"))?
            .trim();
        self.num_of_atoms = count_line
            .parse::<usize>()
            .map_err(|_| BridgeError::input_error(count_line, "invalid atom count"))?;

        self.atoms.clear();

        for i in 0..self.num_of_atoms {
            let line = lines.next().ok_or_else(|| {
                BridgeError::input_error("GRO", format!("unexpected EOF at atom {}", i + 1))
            })?;

            let chars: Vec<char> = line.chars().collect();
            if chars.len() < MIN_RECORD_LEN {
                return Err(BridgeError::input_error(
                    line,
                    format!("line {} too short for GRO record", i + 3),
                ));
            }

            let res_num_str = get_char_slice(&chars, 0, 5);
            let res_num = res_num_str.trim().parse::<usize>().map_err(|_| {
                BridgeError::input_error(
                    &res_num_str,
                    format!("invalid residue number at line {}", i + 3),
                )
            })?;
            let res_name = get_char_slice(&chars, 5, 10).trim().to_string();
            let atom_name = get_char_slice(&chars, 10, 15).trim().to_string();
            let atom_num_str = get_char_slice(&chars, 15, 20);
            let atom_num = atom_num_str.trim().parse::<usize>().map_err(|_| {
                BridgeError::input_error(
                    &atom_num_str,
                    format!("invalid atom number at line {}", i + 3),
                )
            })?;

            let px = parse_fixed_float(&chars, 20, 28)?;
            let py = parse_fixed_float(&chars, 28, 36)?;
            let pz = parse_fixed_float(&chars, 36, 44)?;

            let vx = parse_fixed_float(&chars, 44, 52).unwrap_or(0.0);
            let vy = parse_fixed_float(&chars, 52, 60).unwrap_or(0.0);
            let vz = parse_fixed_float(&chars, 60, 68).unwrap_or(0.0);

            self.atoms.push(GroAtom {
                residue_number: res_num,
                residue_name: res_name,
                atom_name,
                atom_number: atom_num,
                position: Position::new(px, py, pz),
                velocity: [vx, vy, vz],
            });
        }

        // Box vectors line
        if let Some(box_line) = lines.next() {
            let mut vectors: Vec<f64> = box_line
                .split_whitespace()
                .filter_map(|w| w.parse::<f64>().ok())
                .collect();
            while vectors.len() < BOX_VECTORS_PAD_LEN {
                vectors.push(0.0);
            }
            self.box_vectors = vectors;
        } else {
            self.box_vectors = vec![0.0; BOX_VECTORS_PAD_LEN];
        }

        Ok(())
    }

    /// Converts this GRO representation into a hierarchical `AtomGroup`.
    ///
    /// Hierarchical structure:
    /// `output` -> model "1" -> chain "_" -> residue `<residue_number>` -> atom `<atom_number>`
    /// Coordinates are converted from nanometers (nm) to Angstroms (Å) by multiplying by 10.0.
    pub fn get_atomgroup(&self) -> Result<AtomGroup> {
        let mut output = AtomGroup::with_name(&self.title);
        let mut model = AtomGroup::new();
        let mut chain = AtomGroup::new();

        let mut current_res_id = usize::MAX;
        let mut current_ag: Option<AtomGroup> = None;

        for atom_data in &self.atoms {
            if current_res_id != atom_data.residue_number {
                if let Some(ag) = current_ag.take() {
                    chain.set_group(&current_res_id.to_string(), ag);
                }
                current_res_id = atom_data.residue_number;
                let mut new_ag = AtomGroup::new();
                new_ag.name = atom_data.residue_name.clone();
                current_ag = Some(new_ag);
            }

            let mut atom = Atom::new();
            atom.name = atom_data.atom_name.clone();

            // Guess chemical element symbol from atom name
            let symbol = deduce_symbol(&atom_data.atom_name);
            atom.set_symbol(&symbol)?;

            // Convert nm to Angstroms (x 10.0)
            atom.xyz = atom_data.position * 10.0;

            if let Some(ref mut ag) = current_ag {
                ag.set_atom(&atom_data.atom_number.to_string(), atom);
            }
        }

        if let Some(ag) = current_ag {
            chain.set_group(&current_res_id.to_string(), ag);
            model.set_group("_", chain);
            output.set_group("1", model);
        }

        if output.get_bond_list().is_empty() {
            output.setup()?;
        }

        Ok(output)
    }

    /// Populates this GRO representation from an `AtomGroup`.
    /// Coordinates are converted from Angstroms (Å) to nanometers (nm) by multiplying by 0.1.
    pub fn set_by_atomgroup(&mut self, atomgroup: &AtomGroup) {
        self.title = atomgroup.name.clone();
        self.atoms.clear();

        let mut serial = 1;
        let mut residue_index = 1;

        for (_model_key, model) in atomgroup.groups() {
            for (_chain_key, chain) in model.groups() {
                for (_res_key, residue) in chain.groups() {
                    let residue_number = residue_index;
                    residue_index += 1;
                    let residue_name = residue.name.clone();

                    for (_atom_key, atom) in residue.atoms() {
                        let atom_name = atom.name.clone();
                        let atom_number = serial;
                        serial += 1;

                        // Angstroms -> nm (x 0.1)
                        let gro_pos = atom.xyz * 0.1;

                        self.atoms.push(GroAtom {
                            residue_number,
                            residue_name: residue_name.clone(),
                            atom_name,
                            atom_number,
                            position: gro_pos,
                            velocity: [0.0, 0.0, 0.0],
                        });
                    }
                }
            }
        }

        self.num_of_atoms = self.atoms.len();
        let (pos1, pos2) = atomgroup.r#box();
        self.box_vectors = vec![
            (pos2.x - pos1.x).abs(),
            (pos2.y - pos1.y).abs(),
            (pos2.z - pos1.z).abs(),
        ];
    }

    /// Generates the GRO format string representation.
    pub fn get_text(&self) -> String {
        let mut output = format!("{}\n{:>5}\n", self.title, self.atoms.len());
        for atom in &self.atoms {
            let res_num = atom.residue_number % NUMBER_WRAP;
            let atom_num = atom.atom_number % NUMBER_WRAP;
            let res_name = truncate_chars(&atom.residue_name, 5);
            let atom_name = truncate_chars(&atom.atom_name, 5);

            output.push_str(&format!(
                "{:>5}{:<5}{:>5}{:>5}{:8.3}{:8.3}{:8.3}{:8.4}{:8.4}{:8.4}\n",
                res_num,
                res_name,
                atom_name,
                atom_num,
                atom.position.x,
                atom.position.y,
                atom.position.z,
                atom.velocity[0],
                atom.velocity[1],
                atom.velocity[2]
            ));
        }

        if !self.box_vectors.is_empty() {
            let box_strs: Vec<String> = self
                .box_vectors
                .iter()
                .map(|v| format!("{:10.5}", v))
                .collect();
            output.push_str(&box_strs.join(" "));
            output.push('\n');
        }

        output
    }
}

impl fmt::Display for SimpleGro {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.get_text())
    }
}

fn get_char_slice(chars: &[char], start: usize, end: usize) -> String {
    if start >= chars.len() {
        return String::new();
    }
    let actual_end = end.min(chars.len());
    chars[start..actual_end].iter().collect()
}

fn truncate_chars(s: &str, max_len: usize) -> String {
    s.chars().take(max_len).collect()
}

fn parse_fixed_float(chars: &[char], start: usize, end: usize) -> Result<f64> {
    if start >= chars.len() {
        return Ok(0.0);
    }
    let actual_end = end.min(chars.len());
    let slice: String = chars[start..actual_end].iter().collect();
    let trimmed = slice.trim();
    if trimmed.is_empty() {
        return Ok(0.0);
    }
    trimmed.parse::<f64>().map_err(|_| {
        BridgeError::input_error(trimmed, "failed to parse coordinate/velocity in GRO line")
    })
}

fn deduce_symbol(name: &str) -> String {
    let trimmed = name.trim();
    if trimmed.len() <= 1 {
        return trimmed.to_uppercase();
    }
    let first_two = &trimmed[0..2];
    if PeriodicTable::contains(first_two) {
        first_two.to_string()
    } else {
        trimmed[0..1].to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    // Ported from tests/test_gro.py
    #[test]
    fn test_load_sample_gro() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let fixture_path = PathBuf::from(manifest_dir).join("tests/data/sample.gro");

        let gro = SimpleGro::from_file(&fixture_path).unwrap();
        assert_eq!(gro.num_of_atoms, 6);
        assert_eq!(gro.atoms.len(), 6);

        // Verify first atom
        let first = &gro.atoms[0];
        assert_eq!(first.residue_number, 1);
        assert_eq!(first.residue_name, "WATER");
        assert_eq!(first.atom_name, "OW1");
        assert_eq!(first.atom_number, 1);
        assert!((first.position.x - 0.126).abs() < 1e-5);
        assert!((first.position.y - 1.624).abs() < 1e-5);
        assert!((first.position.z - 1.679).abs() < 1e-5);
        assert!((first.velocity[0] - 0.1227).abs() < 1e-5);

        // Convert to AtomGroup
        let ag = gro.get_atomgroup().unwrap();
        assert_eq!(ag.get_number_of_all_atoms(), 6);

        // Verify nm -> Angstroms conversion (0.126 nm -> 1.26 Å)
        let model = ag.get_group("1").unwrap();
        let chain = model.get_group("_").unwrap();
        let res1 = chain.get_group("1").unwrap();
        let ow1 = res1.get_atom("1").unwrap();
        assert_eq!(ow1.name, "OW1");
        assert_eq!(ow1.symbol().unwrap(), "O");
        assert!((ow1.xyz.x - 1.26).abs() < 1e-4);
        assert!((ow1.xyz.y - 16.24).abs() < 1e-4);
        assert!((ow1.xyz.z - 16.79).abs() < 1e-4);

        // Verify box vectors padded to 6 elements
        assert_eq!(gro.box_vectors.len(), 6);
        assert!((gro.box_vectors[0] - 1.82060).abs() < 1e-5);
        assert!((gro.box_vectors[1] - 1.82060).abs() < 1e-5);
        assert!((gro.box_vectors[2] - 1.82060).abs() < 1e-5);
        assert_eq!(gro.box_vectors[3], 0.0);
        assert_eq!(gro.box_vectors[4], 0.0);
        assert_eq!(gro.box_vectors[5], 0.0);
    }

    #[test]
    fn test_gro_non_ascii() {
        // Non-ASCII characters in residue/atom names must not panic
        let gro_text = format!(
            "{}\n{:>5}\n{:>5}{:<5}{:>5}{:>5}{:8.3}{:8.3}{:8.3}\n 1.0 1.0 1.0\n",
            "Title with 日本語", 1, 1, "水分子", "酸素1", 1, 0.1, 0.2, 0.3
        );
        let mut gro = SimpleGro::new();
        let res = gro.parse_str(&gro_text);
        assert!(res.is_ok(), "parse failed: {:?}", res.err());
        assert_eq!(gro.atoms.len(), 1);
        assert_eq!(gro.atoms[0].residue_name, "水分子");
        assert_eq!(gro.atoms[0].atom_name, "酸素1");
    }

    #[test]
    fn test_gro_invalid_number_error() {
        let gro_text = "\
Test
    1
*****WATER  OW1    1   0.126   1.624   1.679
   1.0 1.0 1.0
";
        let mut gro = SimpleGro::new();
        let res = gro.parse_str(gro_text);
        assert!(res.is_err());
    }

    #[test]
    fn test_gro_number_wrap_and_box() {
        let mut ag = AtomGroup::with_name("test_wrap");
        let mut model = AtomGroup::new();
        let mut chain = AtomGroup::new();
        let mut res = AtomGroup::with_name("ALA");
        res.set_atom(
            "1",
            Atom::new_with_pos("C", Position::new(0.0, 0.0, 0.0)).unwrap(),
        );
        res.set_atom(
            "2",
            Atom::new_with_pos("C", Position::new(10.0, 20.0, 30.0)).unwrap(),
        );
        chain.set_group("1", res);
        model.set_group("A", chain);
        ag.set_group("1", model);

        let mut gro = SimpleGro::new();
        gro.set_by_atomgroup(&ag);

        assert_eq!(gro.box_vectors.len(), 3);
        assert!((gro.box_vectors[0] - 10.0).abs() < 1e-5);
        assert!((gro.box_vectors[1] - 20.0).abs() < 1e-5);
        assert!((gro.box_vectors[2] - 30.0).abs() < 1e-5);

        // Check number wrapping at 10000
        gro.atoms[0].residue_number = 10005;
        gro.atoms[0].atom_number = 20003;
        let text = gro.get_text();
        assert!(text.contains("    5"));
        assert!(text.contains("    3"));
    }
}
