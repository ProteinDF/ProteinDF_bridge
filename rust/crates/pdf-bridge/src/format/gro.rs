// Copyright (C) 2014 The ProteinDF development team.
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

use std::fmt;
use std::fs;
use std::path::Path;

use crate::atom::Atom;
use crate::atom_group::AtomGroup;
use crate::error::{BridgeError, Result};
use crate::periodic_table::PeriodicTable;
use crate::position::Position;

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

            if line.len() < 44 {
                return Err(BridgeError::input_error(
                    line,
                    format!("line {} too short for GRO record", i + 3),
                ));
            }

            let res_num = line[0..5.min(line.len())]
                .trim()
                .parse::<usize>()
                .unwrap_or(0);
            let res_name = line[5..10.min(line.len())].trim().to_string();
            let atom_name = line[10..15.min(line.len())].trim().to_string();
            let atom_num = line[15..20.min(line.len())]
                .trim()
                .parse::<usize>()
                .unwrap_or(0);

            let px = parse_fixed_float(line, 20, 28)?;
            let py = parse_fixed_float(line, 28, 36)?;
            let pz = parse_fixed_float(line, 36, 44)?;

            let vx = parse_fixed_float(line, 44, 52).unwrap_or(0.0);
            let vy = parse_fixed_float(line, 52, 60).unwrap_or(0.0);
            let vz = parse_fixed_float(line, 60, 68).unwrap_or(0.0);

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
            while vectors.len() < 3 {
                vectors.push(0.0);
            }
            self.box_vectors = vectors;
        } else {
            self.box_vectors = vec![0.0, 0.0, 0.0];
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
            atom.xyz = Position::new(
                atom_data.position.x * 10.0,
                atom_data.position.y * 10.0,
                atom_data.position.z * 10.0,
            );

            if let Some(ref mut ag) = current_ag {
                ag.set_atom(&atom_data.atom_number.to_string(), atom);
            }
        }

        if let Some(ag) = current_ag {
            chain.set_group(&current_res_id.to_string(), ag);
            model.set_group("_", chain);
            output.set_group("1", model);
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
                        let px = atom.xyz.x * 0.1;
                        let py = atom.xyz.y * 0.1;
                        let pz = atom.xyz.z * 0.1;

                        self.atoms.push(GroAtom {
                            residue_number,
                            residue_name: residue_name.clone(),
                            atom_name,
                            atom_number,
                            position: Position::new(px, py, pz),
                            velocity: [0.0, 0.0, 0.0],
                        });
                    }
                }
            }
        }

        self.num_of_atoms = self.atoms.len();
        self.box_vectors = vec![0.0, 0.0, 0.0];
    }

    /// Generates the GRO format string representation.
    pub fn get_text(&self) -> String {
        let mut output = format!("{}\n{:>5}\n", self.title, self.atoms.len());
        for atom in &self.atoms {
            let res_num = atom.residue_number % 100000;
            let atom_num = atom.atom_number % 100000;
            let res_name = if atom.residue_name.len() > 5 {
                &atom.residue_name[0..5]
            } else {
                &atom.residue_name
            };
            let atom_name = if atom.atom_name.len() > 5 {
                &atom.atom_name[0..5]
            } else {
                &atom.atom_name
            };

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

fn parse_fixed_float(line: &str, start: usize, end: usize) -> Result<f64> {
    if start >= line.len() {
        return Ok(0.0);
    }
    let actual_end = end.min(line.len());
    let slice = line[start..actual_end].trim();
    if slice.is_empty() {
        return Ok(0.0);
    }
    slice.parse::<f64>().map_err(|_| {
        BridgeError::input_error(slice, "failed to parse coordinate/velocity in GRO line")
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
    }
}
