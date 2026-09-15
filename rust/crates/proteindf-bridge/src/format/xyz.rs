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
use crate::position::Position;

/// A single atom entry in an XYZ file.
#[derive(Debug, Clone, PartialEq)]
pub struct XyzAtom {
    pub symbol: String,
    pub position: Position,
}

/// Reader and writer for XYZ coordinate files, corresponding to `proteindf_bridge.xyz.Xyz`.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Xyz {
    comment: String,
    atoms: Vec<XyzAtom>,
}

impl Xyz {
    /// Creates a new empty `Xyz` instance.
    pub fn new() -> Self {
        Self::default()
    }

    /// Loads an XYZ file from the given path.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        let mut xyz = Self::new();
        xyz.load(path)?;
        Ok(xyz)
    }

    /// Creates an `Xyz` instance from an `AtomGroup`.
    pub fn from_atomgroup(ag: &AtomGroup) -> Self {
        let mut xyz = Self::new();
        xyz.set_by_atomgroup(ag);
        xyz
    }

    /// Returns the comment / title line.
    pub fn comment(&self) -> &str {
        &self.comment
    }

    /// Sets the comment / title line.
    pub fn set_comment(&mut self, comment: impl Into<String>) {
        self.comment = comment.into();
    }

    /// Returns the slice of atoms.
    pub fn atoms(&self) -> &[XyzAtom] {
        &self.atoms
    }

    /// Loads XYZ data from a file.
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

    /// Parses XYZ format string.
    pub fn parse_str(&mut self, content: &str) -> Result<()> {
        let mut lines = content.lines();
        let num_str = lines
            .next()
            .ok_or_else(|| BridgeError::input_error("XYZ", "empty content"))?
            .trim();
        let num_of_atoms = num_str
            .parse::<usize>()
            .map_err(|_| BridgeError::input_error(num_str, "invalid atom count"))?;

        self.comment = lines.next().unwrap_or("").trim().to_string();
        self.atoms.clear();

        for i in 0..num_of_atoms {
            let line = lines.next().ok_or_else(|| {
                BridgeError::input_error("XYZ", format!("unexpected EOF at atom {}", i + 1))
            })?;
            let words: Vec<&str> = line.split_whitespace().collect();
            if words.len() < 4 {
                return Err(BridgeError::input_error(
                    format!("line {}: '{}'", i + 3, line),
                    "insufficient atom fields in XYZ",
                ));
            }
            let symbol = words[0].to_string();
            let x = words[1]
                .parse::<f64>()
                .map_err(|_| BridgeError::input_error(words[1], "invalid x coordinate"))?;
            let y = words[2]
                .parse::<f64>()
                .map_err(|_| BridgeError::input_error(words[2], "invalid y coordinate"))?;
            let z = words[3]
                .parse::<f64>()
                .map_err(|_| BridgeError::input_error(words[3], "invalid z coordinate"))?;

            self.atoms.push(XyzAtom {
                symbol,
                position: Position::new(x, y, z),
            });
        }

        Ok(())
    }

    /// Saves XYZ data to a file.
    pub fn save(&self, file_path: impl AsRef<Path>) -> Result<()> {
        let path = file_path.as_ref();
        fs::write(path, self.get_text()).map_err(|e| {
            BridgeError::general(format!("failed to write XYZ to {}: {}", path.display(), e))
        })
    }

    /// Converts this XYZ instance into an `AtomGroup`.
    pub fn get_atom_group(&self) -> Result<AtomGroup> {
        let mut root = AtomGroup::with_name(&self.comment);
        for (i, atom_data) in self.atoms.iter().enumerate() {
            let atom = Atom::new_with_pos(&atom_data.symbol, atom_data.position)?;
            root.set_atom(&i.to_string(), atom);
        }
        Ok(root)
    }

    /// Populates this XYZ instance from an `AtomGroup` recursively.
    pub fn set_by_atomgroup(&mut self, ag: &AtomGroup) {
        self.comment = ag.name.clone();
        self.atoms.clear();
        self.collect_atoms_recursive(ag);
    }

    fn collect_atoms_recursive(&mut self, ag: &AtomGroup) {
        for (_key, group) in ag.groups() {
            self.collect_atoms_recursive(group);
        }
        for (_key, atom) in ag.atoms() {
            let symbol = atom.symbol().unwrap_or("X").to_string();
            self.atoms.push(XyzAtom {
                symbol,
                position: atom.xyz,
            });
        }
    }

    /// Returns the XYZ format text representation.
    pub fn get_text(&self) -> String {
        let mut output = format!("{}\n{}\n", self.atoms.len(), self.comment);
        for atom in &self.atoms {
            output.push_str(&format!(
                "{} {:10.6} {:10.6} {:10.6}\n",
                atom.symbol, atom.position.x, atom.position.y, atom.position.z
            ));
        }
        output
    }
}

impl fmt::Display for Xyz {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.get_text())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    // Ported from tests/test_xyz.py
    #[test]
    fn test_init_with_atomgroup() {
        let mut ag = AtomGroup::with_name("water");
        ag.set_atom(
            "1",
            Atom::new_with_pos("O", Position::new(0.0, 0.0, 0.0)).unwrap(),
        );
        ag.set_atom(
            "2",
            Atom::new_with_pos("H", Position::new(0.0, 1.0, 0.0)).unwrap(),
        );
        ag.set_atom(
            "3",
            Atom::new_with_pos("H", Position::new(1.0, 0.0, 0.0)).unwrap(),
        );

        let xyz = Xyz::from_atomgroup(&ag);
        let text = xyz.get_text();
        assert!(text.starts_with("3\n"));
        assert!(text.contains("O"));
        assert!(text.contains("H"));
    }

    #[test]
    fn test_save_and_load() {
        let mut ag = AtomGroup::with_name("test_mol");
        ag.set_atom(
            "1",
            Atom::new_with_pos("C", Position::new(1.2, 3.4, 5.6)).unwrap(),
        );

        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("test_save_and_load.xyz");

        let xyz1 = Xyz::from_atomgroup(&ag);
        xyz1.save(&temp_path).unwrap();

        let xyz2 = Xyz::from_file(&temp_path).unwrap();
        let loaded_ag = xyz2.get_atom_group().unwrap();
        assert_eq!(loaded_ag.get_number_of_atoms(), 1);

        let atom = loaded_ag.get_atom("0").unwrap();
        assert_eq!(atom.symbol().unwrap(), "C");
        assert!((atom.xyz.x - 1.2).abs() < 1e-4);
        assert!((atom.xyz.y - 3.4).abs() < 1e-4);
        assert!((atom.xyz.z - 5.6).abs() < 1e-4);

        let _ = fs::remove_file(&temp_path);
    }

    #[test]
    fn test_load_fixture() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let fixture_path = PathBuf::from(manifest_dir).join("tests/data/ACE_ALA_NME.xyz");

        let xyz = Xyz::from_file(&fixture_path).unwrap();
        assert_eq!(xyz.atoms().len(), 22);

        let ag = xyz.get_atom_group().unwrap();
        assert_eq!(ag.get_number_of_all_atoms(), 22);
    }

    #[test]
    fn test_truncated_xyz_error() {
        let content = "3\ntest comment\nC 0.0 0.0 0.0\n";
        let mut xyz = Xyz::new();
        let result = xyz.parse_str(content);
        assert!(result.is_err());
    }
}
