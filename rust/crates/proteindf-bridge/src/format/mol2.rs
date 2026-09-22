// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::collections::HashMap;
use std::fmt;
use std::fs;
use std::path::Path;

use crate::atom::Atom;
use crate::atom_group::AtomGroup;
use crate::error::{BridgeError, Result};
use crate::periodic_table::PeriodicTable;
use crate::position::Position;

/// Tripos Mol2 format reader and writer corresponding to `proteindf_bridge.mol2.SimpleMol2`.
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

    /// Loads a Mol2 file from the given path.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        let mut mol2 = Self::new();
        mol2.load(path)?;
        Ok(mol2)
    }

    /// Parses a Mol2 format string into a `SimpleMol2` instance.
    #[allow(clippy::should_implement_trait)]
    pub fn from_str(content: &str) -> Result<Self> {
        <Self as std::str::FromStr>::from_str(content)
    }

    /// Loads a Mol2 file from the specified path.
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

    /// Parses Mol2 content from a string.
    pub fn parse_str(&mut self, content: &str) -> Result<()> {
        enum Section {
            None,
            Molecule,
            Atom,
            Bond,
            Other,
        }

        let mut section = Section::None;
        let mut mol_name = String::new();
        let mut mol_lines = Vec::new();

        let mut parsed_atoms: Vec<(usize, Atom)> = Vec::new();
        let mut parsed_bonds: Vec<(usize, usize, usize)> = Vec::new();

        for line in content.lines() {
            let trimmed = line.trim();
            if trimmed.starts_with("@<TRIPOS>") {
                section = match trimmed {
                    "@<TRIPOS>MOLECULE" => {
                        mol_lines.clear();
                        Section::Molecule
                    }
                    "@<TRIPOS>ATOM" => Section::Atom,
                    "@<TRIPOS>BOND" => Section::Bond,
                    _ => Section::Other,
                };
                continue;
            }

            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }

            match section {
                Section::Molecule => {
                    mol_lines.push(trimmed);
                    if mol_lines.len() == 1 {
                        mol_name = trimmed.to_string();
                    }
                }
                Section::Atom => {
                    let tokens: Vec<&str> = trimmed.split_whitespace().collect();
                    if tokens.len() < 6 {
                        continue;
                    }
                    let atom_id = tokens[0].parse::<usize>().map_err(|_| {
                        BridgeError::input_error(tokens[0], "invalid atom_id in @<TRIPOS>ATOM")
                    })?;
                    let atom_name = tokens[1].to_string();
                    let x = tokens[2].parse::<f64>().map_err(|_| {
                        BridgeError::input_error(tokens[2], "invalid x coordinate in @<TRIPOS>ATOM")
                    })?;
                    let y = tokens[3].parse::<f64>().map_err(|_| {
                        BridgeError::input_error(tokens[3], "invalid y coordinate in @<TRIPOS>ATOM")
                    })?;
                    let z = tokens[4].parse::<f64>().map_err(|_| {
                        BridgeError::input_error(tokens[4], "invalid z coordinate in @<TRIPOS>ATOM")
                    })?;
                    let atom_type = tokens[5];

                    let charge = if tokens.len() >= 9 {
                        tokens[8].parse::<f64>().ok()
                    } else if tokens.len() == 7 || tokens.len() == 8 {
                        tokens.last().and_then(|tok| tok.parse::<f64>().ok())
                    } else {
                        None
                    };

                    let symbol = deduce_symbol_from_mol2(atom_type, &atom_name);
                    let mut atom = Atom::new_with_pos(&symbol, Position::new(x, y, z))?;
                    atom.name = atom_name;
                    if let Some(ch) = charge {
                        atom.charge = ch;
                    }

                    parsed_atoms.push((atom_id, atom));
                }
                Section::Bond => {
                    let tokens: Vec<&str> = trimmed.split_whitespace().collect();
                    if tokens.len() < 4 {
                        continue;
                    }
                    let atom1_id = tokens[1].parse::<usize>().map_err(|_| {
                        BridgeError::input_error(tokens[1], "invalid atom_id1 in @<TRIPOS>BOND")
                    })?;
                    let atom2_id = tokens[2].parse::<usize>().map_err(|_| {
                        BridgeError::input_error(tokens[2], "invalid atom_id2 in @<TRIPOS>BOND")
                    })?;
                    let order = match tokens[3].to_lowercase().as_str() {
                        "1" | "am" | "du" | "un" => 1,
                        "2" => 2,
                        "3" => 3,
                        "ar" => 1,
                        other => other.parse::<usize>().unwrap_or(1),
                    };
                    parsed_bonds.push((atom1_id, atom2_id, order));
                }
                _ => {}
            }
        }

        let mut ag = AtomGroup::with_name(&mol_name);
        let mut id_to_atom: HashMap<usize, Atom> = HashMap::new();

        for (id, atom) in &parsed_atoms {
            let key = id.to_string();
            ag.set_atom(&key, atom.clone());
            if let Some(stored) = ag.get_atom(&key) {
                id_to_atom.insert(*id, stored.clone());
            } else {
                let mut a = atom.clone();
                a.path = format!("{}{}", ag.path(), key);
                id_to_atom.insert(*id, a);
            }
        }

        for (id1, id2, order) in parsed_bonds {
            if let (Some(a1), Some(a2)) = (id_to_atom.get(&id1), id_to_atom.get(&id2)) {
                ag.add_bond(a1, a2, order);
            }
        }

        if ag.get_bond_list().is_empty() {
            ag.setup()?;
        }

        self.set_by_atomgroup(&ag);
        Ok(())
    }

    /// Returns a reference to the underlying `AtomGroup`.
    ///
    /// If the MOL2 source contains an `@<TRIPOS>BOND` section, the returned `AtomGroup`
    /// includes the explicit bond topology (pairs and orders) parsed from the file.
    /// According to the bond priority policy (see `RUST_PORT_SPEC.md` §3.8 & §3.14), explicit
    /// file-derived bonds take precedence. If no explicit bond records exist,
    /// [`AtomGroup::setup`] is automatically executed to resolve bonds.
    pub fn get_atomgroup(&self) -> &AtomGroup {
        &self.atomgroup
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

        let atoms = if self.atomgroup.get_number_of_atoms() > 0 {
            self.atomgroup
                .atoms()
                .map(|(_, a)| a.clone())
                .collect::<Vec<_>>()
        } else {
            self.atomgroup.get_atom_list()
        };

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

            let atom_id1 = if let Some(pos) = atoms.iter().position(|a| a.path == bond.atom1_path) {
                pos + 1
            } else {
                self.atom_index_table
                    .iter()
                    .position(|n| n == a1_name)
                    .ok_or_else(|| {
                        BridgeError::input_error(
                            a1_name,
                            format!("atom '{}' in bond not found in atom index table", a1_name),
                        )
                    })?
                    + 1
            };
            let atom_id2 = if let Some(pos) = atoms.iter().position(|a| a.path == bond.atom2_path) {
                pos + 1
            } else {
                self.atom_index_table
                    .iter()
                    .position(|n| n == a2_name)
                    .ok_or_else(|| {
                        BridgeError::input_error(
                            a2_name,
                            format!("atom '{}' in bond not found in atom index table", a2_name),
                        )
                    })?
                    + 1
            };
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

fn deduce_symbol_from_mol2(atom_type: &str, atom_name: &str) -> String {
    // 1. Try prefix of SYBYL atom type (e.g. "C" from "C.3", "Cl" from "Cl", "Fe" from "Fe")
    let raw_elem = atom_type.split('.').next().unwrap_or(atom_type);
    if !raw_elem.is_empty() {
        let mut chars = raw_elem.chars();
        let first = chars.next().unwrap().to_ascii_uppercase();
        let rest: String = chars.take(1).map(|c| c.to_ascii_lowercase()).collect();
        let candidate = format!("{}{}", first, rest);
        if PeriodicTable::contains(&candidate) {
            return candidate;
        }
        let single = first.to_string();
        if PeriodicTable::contains(&single) {
            return single;
        }
    }

    // 2. Fallback to atom_name deduction
    let trimmed = atom_name.trim();
    if trimmed.len() >= 2 {
        let mut chars = trimmed.chars();
        let first = chars.next().unwrap().to_ascii_uppercase();
        let second = chars.next().unwrap().to_ascii_lowercase();
        let two_char = format!("{}{}", first, second);
        if PeriodicTable::contains(&two_char) {
            return two_char;
        }
    }
    if !trimmed.is_empty() {
        let first = trimmed[..1].to_ascii_uppercase();
        if PeriodicTable::contains(&first) {
            return first;
        }
    }

    "C".to_string()
}

impl fmt::Display for SimpleMol2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.get_text() {
            Ok(text) => write!(f, "{}", text),
            Err(e) => write!(f, "{}", e),
        }
    }
}

impl std::str::FromStr for SimpleMol2 {
    type Err = BridgeError;

    fn from_str(s: &str) -> Result<Self> {
        let mut mol2 = Self::new();
        mol2.parse_str(s)?;
        Ok(mol2)
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

    #[test]
    fn test_mol2_roundtrip() {
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

        // Load back from string
        let loaded = SimpleMol2::from_str(&text).unwrap();
        let loaded_ag = loaded.get_atomgroup();

        assert_eq!(loaded_ag.name, "ethane");
        assert_eq!(loaded_ag.get_number_of_all_atoms(), 2);

        let at1 = loaded_ag.get_atom("C1").expect("C1 not found");
        let at2 = loaded_ag.get_atom("C2").expect("C2 not found");
        assert_eq!(at1.symbol().unwrap(), "C");
        assert_eq!(at2.symbol().unwrap(), "C");
        assert!((at1.xyz.x - 0.0).abs() < 1e-3);
        assert!((at2.xyz.x - 1.5).abs() < 1e-3);

        // Verify bond
        let mut ag_mut = loaded_ag.clone();
        let bonds = ag_mut.get_bond_list();
        assert_eq!(bonds.len(), 1);
        assert_eq!(bonds[0].order, 1);

        // Verify text identity
        let roundtrip_text = loaded.get_text().unwrap();
        assert_eq!(text, roundtrip_text);
    }

    #[test]
    fn test_mol2_load_synthetic_fixture() {
        let mol2_str = r#"@<TRIPOS>MOLECULE
benzene_derivative
3 2
SMALL
USER_CHARGES

@<TRIPOS>ATOM
      1 C1          0.0000    1.4000    0.0000 C.ar       1 BENZ   -0.1500
      2 C2          1.2124    0.7000    0.0000 C.ar       1 BENZ   -0.1500
      3 N1          2.4248    1.4000    0.0000 N.pl3      1 BENZ    0.3000
@<TRIPOS>BOND
     1     1     2 ar
     2     2     3 1
"#;

        let mol2 = SimpleMol2::from_str(mol2_str).unwrap();
        let ag = mol2.get_atomgroup();

        assert_eq!(ag.name, "benzene_derivative");
        assert_eq!(ag.get_number_of_all_atoms(), 3);

        let c1 = ag.get_atom("C1").unwrap();
        assert_eq!(c1.symbol().unwrap(), "C");
        assert!((c1.charge - (-0.15)).abs() < 1e-4);
        assert!((c1.xyz.y - 1.4).abs() < 1e-4);

        let n1 = ag.get_atom("N1").unwrap();
        assert_eq!(n1.symbol().unwrap(), "N");
        assert!((n1.charge - 0.30).abs() < 1e-4);

        let mut ag_mut = ag.clone();
        let bonds = ag_mut.get_bond_list();
        assert_eq!(bonds.len(), 2);
        // "ar" mapped to order 1
        assert_eq!(bonds[0].order, 1);
        assert_eq!(bonds[1].order, 1);
    }

    #[test]
    fn test_mol2_from_file_roundtrip() {
        let mut ag = AtomGroup::with_name("water");
        let mut a1 = Atom::new_with_pos("O", Position::new(0.0, 0.0, 0.0)).unwrap();
        a1.name = "O1".to_string();
        let mut a2 = Atom::new_with_pos("H", Position::new(0.0, 1.0, 0.0)).unwrap();
        a2.name = "H1".to_string();

        ag.set_atom("1", a1.clone());
        ag.set_atom("2", a2.clone());
        ag.add_bond(&a1, &a2, 1);

        let mol2 = SimpleMol2::from_atomgroup(&ag);

        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("test_roundtrip_file.mol2");
        mol2.save(&temp_path).unwrap();

        let loaded = SimpleMol2::from_file(&temp_path).unwrap();
        let _ = fs::remove_file(&temp_path);

        assert_eq!(loaded.get_atomgroup().name, "water");
        assert_eq!(loaded.get_atomgroup().get_number_of_all_atoms(), 2);
        assert_eq!(loaded.get_atomgroup().get_number_of_bonds(), 1);
    }

    #[test]
    fn test_mol2_duplicate_atom_names_bonds() {
        // Molecule with duplicate atom names (two hydrogen atoms both named "H")
        let mol2_str = r#"@<TRIPOS>MOLECULE
dup_test
4 2
SMALL
NO_CHARGES

@<TRIPOS>ATOM
      1 N           0.0000    0.0000    0.0000 N.3        1 RES      0.0000
      2 O           2.0000    0.0000    0.0000 O.3        1 RES      0.0000
      3 H           0.0000    1.0000    0.0000 H          1 RES      0.0000
      4 H           2.0000    1.0000    0.0000 H          1 RES      0.0000
@<TRIPOS>BOND
     1     1     3 1
     2     2     4 1
"#;

        let mol2 = SimpleMol2::from_str(mol2_str).unwrap();
        let ag = mol2.get_atomgroup();

        let mut ag_mut = ag.clone();
        let bonds = ag_mut.get_bond_list();
        assert_eq!(bonds.len(), 2);

        // Bond 1: connects /1 (N) and /3 (first H)
        let (b1_p1, b1_p2) = (&bonds[0].atom1_path, &bonds[0].atom2_path);
        assert!(
            (b1_p1 == "/1" && b1_p2 == "/3") || (b1_p1 == "/3" && b1_p2 == "/1"),
            "Bond 1 must connect /1 and /3, got: {} - {}",
            b1_p1,
            b1_p2
        );

        // Bond 2: connects /2 (O) and /4 (second H), MUST NOT erroneously connect to /3
        let (b2_p1, b2_p2) = (&bonds[1].atom1_path, &bonds[1].atom2_path);
        assert!(
            (b2_p1 == "/2" && b2_p2 == "/4") || (b2_p1 == "/4" && b2_p2 == "/2"),
            "Bond 2 must connect /2 and /4, got: {} - {}",
            b2_p1,
            b2_p2
        );

        // Ensure resolution points to correct distinct atoms
        let (a1_1, a1_2) = ag.resolve_bond(&bonds[0]).unwrap();
        let (a2_1, a2_2) = ag.resolve_bond(&bonds[1]).unwrap();

        assert_eq!(a1_1.name, "N");
        assert_eq!(a1_2.name, "H");
        assert_eq!(a2_1.name, "O");
        assert_eq!(a2_2.name, "H");

        // Coordinates of the bonded H atoms must differ (H3 at x=0, H4 at x=2)
        let h3_xyz = if a1_2.name == "H" { a1_2.xyz } else { a1_1.xyz };
        let h4_xyz = if a2_2.name == "H" { a2_2.xyz } else { a2_1.xyz };
        assert!((h3_xyz.x - 0.0).abs() < 1e-4);
        assert!((h4_xyz.x - 2.0).abs() < 1e-4);

        // Re-export text must preserve correct bond atom IDs (1 3 and 2 4)
        let reexported = mol2.get_text().unwrap();
        assert!(
            reexported.contains("1     1     3     1"),
            "Re-exported bond 1 must be '1 1 3 1'"
        );
        assert!(
            reexported.contains("2     2     4     1"),
            "Re-exported bond 2 must be '2 2 4 1' (not '2 2 3 1')"
        );
    }
}
