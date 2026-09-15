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

use crate::atom_group::AtomGroup;
use crate::error::{BridgeError, Result};

/// Tripos Mol2 format writer corresponding to `proteindf_bridge.mol2.SimpleMol2`.
#[derive(Debug, Clone, Default)]
pub struct SimpleMol2 {
    atomgroup: AtomGroup,
    atom_index_table: Vec<String>,
}

impl SimpleMol2 {
    /// Creates a new empty `SimpleMol2` instance.
    pub fn new() -> Self {
        Self::default()
    }

    /// Creates a `SimpleMol2` instance populated from an `AtomGroup`.
    pub fn from_atomgroup(atomgroup: &AtomGroup) -> Self {
        let mut mol2 = Self::new();
        mol2.set_by_atomgroup(atomgroup);
        mol2
    }

    /// Sets the underlying `AtomGroup` and builds the atom index table.
    pub fn set_by_atomgroup(&mut self, atomgroup: &AtomGroup) {
        self.atomgroup = atomgroup.clone();
        self.atom_index_table = self
            .atomgroup
            .get_atom_list()
            .into_iter()
            .map(|a| a.name)
            .collect();
    }

    /// Saves the Mol2 representation to a file.
    pub fn save(&self, file_path: impl AsRef<Path>) -> Result<()> {
        let path = file_path.as_ref();
        let contents = self.get_text()?;
        fs::write(path, contents).map_err(|e| {
            BridgeError::input_error(
                path.display().to_string(),
                format!("failed to write file: {}", e),
            )
        })
    }

    /// Returns the `@<TRIPOS>MOLECULE` section string.
    fn get_contents_molecule(&self) -> String {
        let molecular_name = &self.atomgroup.name;
        let num_atoms = self.atomgroup.get_number_of_all_atoms();
        let num_bonds = self.atomgroup.get_number_of_bonds();
        let mol_type = "SMALL";
        let charge_type = "USER_CHARGES";

        format!(
            "@<TRIPOS>MOLECULE\n{}\n{} {}\n{}\n{}\n\n",
            molecular_name, num_atoms, num_bonds, mol_type, charge_type
        )
    }

    /// Returns the `@<TRIPOS>ATOM` section string.
    /// Note: if direct atoms are present, uses them; otherwise falls back to all atoms (get_atom_list),
    /// which allows Mol2 generation from hierarchical AtomGroups as an intentional enhancement over Python.
    fn get_contents_atom(&self) -> String {
        let mut output = String::from("@<TRIPOS>ATOM\n");
        let atoms = if self.atomgroup.get_number_of_atoms() > 0 {
            self.atomgroup
                .atoms()
                .map(|(_, a)| a.clone())
                .collect::<Vec<_>>()
        } else {
            self.atomgroup.get_atom_list()
        };

        for (i, atom) in atoms.iter().enumerate() {
            let atom_id = i + 1;
            let atom_name = &atom.name;
            let x = atom.xyz.x;
            let y = atom.xyz.y;
            let z = atom.xyz.z;
            let atom_type = atom.symbol().unwrap_or("");

            output.push_str(&format!(
                "{:<5} {:<3} {:8.3} {:8.3} {:8.3} {}\n",
                atom_id, atom_name, x, y, z, atom_type
            ));
        }

        output
    }

    /// Returns the `@<TRIPOS>BOND` section string.
    fn get_contents_bond(&self) -> Result<String> {
        let mut output = String::from("@<TRIPOS>BOND\n");
        let mut ag_clone = self.atomgroup.clone();
        let bond_list = ag_clone.get_bond_list();

        for (i, bond) in bond_list.iter().enumerate() {
            let bond_id = i + 1;

            let atom1 = ag_clone.get_atom_by_path(&bond.atom1_path);
            let atom2 = ag_clone.get_atom_by_path(&bond.atom2_path);

            let a1_name = atom1.map(|a| a.name.as_str()).unwrap_or_else(|| {
                bond.atom1_path
                    .rsplit('/')
                    .next()
                    .unwrap_or(&bond.atom1_path)
            });
            let a2_name = atom2.map(|a| a.name.as_str()).unwrap_or_else(|| {
                bond.atom2_path
                    .rsplit('/')
                    .next()
                    .unwrap_or(&bond.atom2_path)
            });

            let atom_id1 = self
                .atom_index_table
                .iter()
                .position(|n| n == a1_name)
                .ok_or_else(|| {
                    BridgeError::input_error(
                        a1_name,
                        format!("atom '{}' in bond not found in atom index table", a1_name),
                    )
                })?
                + 1;
            let atom_id2 = self
                .atom_index_table
                .iter()
                .position(|n| n == a2_name)
                .ok_or_else(|| {
                    BridgeError::input_error(
                        a2_name,
                        format!("atom '{}' in bond not found in atom index table", a2_name),
                    )
                })?
                + 1;
            let bond_type = bond.order;

            output.push_str(&format!(
                "{:<5} {:<5} {:<5} {}\n",
                bond_id, atom_id1, atom_id2, bond_type
            ));
        }
        output.push('\n');

        Ok(output)
    }

    /// Generates the complete Tripos Mol2 format string representation.
    pub fn get_text(&self) -> Result<String> {
        let mut output = String::new();
        output.push_str(&self.get_contents_molecule());
        output.push_str(&self.get_contents_atom());
        output.push_str(&self.get_contents_bond()?);
        Ok(output)
    }
}

impl fmt::Display for SimpleMol2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.get_text() {
            Ok(text) => write!(f, "{}", text),
            Err(e) => write!(f, "{}", e),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::position::Position;

    // Ported from tests/test_mol2.py
    #[test]
    fn test_save_and_format() {
        let mut ag = AtomGroup::new();
        ag.name = "water".to_string();

        let mut a1 = Atom::new_with_pos("O", Position::new(0.0, 0.0, 0.0)).unwrap();
        a1.name = "O1".to_string();
        a1.charge = -0.8;

        let mut a2 = Atom::new_with_pos("H", Position::new(0.0, 1.0, 0.0)).unwrap();
        a2.name = "H1".to_string();
        a2.charge = 0.4;

        let mut a3 = Atom::new_with_pos("H", Position::new(1.0, 0.0, 0.0)).unwrap();
        a3.name = "H2".to_string();
        a3.charge = 0.4;

        ag.set_atom("1", a1);
        ag.set_atom("2", a2);
        ag.set_atom("3", a3);

        let mol2 = SimpleMol2::from_atomgroup(&ag);
        let text = mol2.get_text().unwrap();

        assert!(text.contains("@<TRIPOS>MOLECULE"));
        assert!(text.contains("water"));
        assert!(text.contains("@<TRIPOS>ATOM"));
        assert!(text.contains("O1"));
        assert!(text.contains("H1"));
        assert!(text.contains("@<TRIPOS>BOND"));

        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("test_save_and_format.mol2");

        mol2.save(&temp_path).unwrap();
        assert!(temp_path.exists());

        let content = fs::read_to_string(&temp_path).unwrap();
        assert!(content.contains("@<TRIPOS>MOLECULE"));

        let _ = fs::remove_file(&temp_path);
    }

    #[test]
    fn test_mol2_with_bonds() {
        let mut ag = AtomGroup::with_name("ethane");
        let mut a1 = Atom::new_with_pos("C", Position::new(0.0, 0.0, 0.0)).unwrap();
        a1.name = "C1".to_string();
        let mut a2 = Atom::new_with_pos("C", Position::new(1.5, 0.0, 0.0)).unwrap();
        a2.name = "C2".to_string();

        ag.set_atom("1", a1.clone());
        ag.set_atom("2", a2.clone());
        ag.add_bond(&a1, &a2, 1);

        let mol2 = SimpleMol2::from_atomgroup(&ag);
        let text = mol2.get_text().unwrap();

        assert!(text.contains("ethane"));
        assert!(text.contains("2 1\n")); // 2 atoms, 1 bond
        assert!(text.contains("@<TRIPOS>BOND\n1     1     2     1\n"));
    }

    #[test]
    fn test_mol2_bond_atom_not_found_error() {
        let mut ag = AtomGroup::with_name("incomplete");
        let mut a1 = Atom::new_with_pos("C", Position::new(0.0, 0.0, 0.0)).unwrap();
        a1.name = "C1".to_string();
        let mut a2 = Atom::new_with_pos("C", Position::new(1.5, 0.0, 0.0)).unwrap();
        a2.name = "C2".to_string();

        // Only register a1, but add bond to non-existent a2
        ag.set_atom("1", a1.clone());
        ag.add_bond(&a1, &a2, 1);

        let mol2 = SimpleMol2::from_atomgroup(&ag);
        let res = mol2.get_text();
        assert!(res.is_err());
    }
}
