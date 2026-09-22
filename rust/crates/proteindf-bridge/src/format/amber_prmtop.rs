// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::fs;
use std::path::Path;

use crate::atom::Atom;
use crate::atom_group::AtomGroup;
use crate::error::{BridgeError, Result};
use crate::periodic_table::PeriodicTable;
use crate::position::Position;

/// Charge unit conversion factor for Amber PRMTOP files (18.2223).
pub const AMBER_CHARGE_FACTOR: f64 = 18.2223;

/// Amber PRMTOP topology and INPCRD coordinate reader corresponding to `proteindf_bridge.amber_prmtop.AmberPrmtop`.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct AmberPrmtop {
    atom_names: Vec<String>,
    charges: Vec<f64>,
    atomic_numbers: Vec<usize>,
    xyz: Vec<Position>,
    bonds: Vec<(usize, usize)>,
}

impl AmberPrmtop {
    /// Creates a new empty `AmberPrmtop` instance.
    pub fn new() -> Self {
        Self::default()
    }

    /// Loads and parses Amber PRMTOP and INPCRD files.
    pub fn from_files(
        prmtop_path: impl AsRef<Path>,
        inpcrd_path: impl AsRef<Path>,
    ) -> Result<Self> {
        let mut amber = Self::new();
        amber.load(prmtop_path, inpcrd_path)?;
        Ok(amber)
    }

    /// Parses Amber PRMTOP and INPCRD from content strings.
    pub fn from_strings(prmtop_content: &str, inpcrd_content: &str) -> Result<Self> {
        let mut amber = Self::new();
        amber.parse_prmtop_str(prmtop_content)?;
        amber.parse_inpcrd_str(inpcrd_content)?;
        amber.validate_data()?;
        Ok(amber)
    }

    /// Loads topology and coordinate files.
    pub fn load(
        &mut self,
        prmtop_path: impl AsRef<Path>,
        inpcrd_path: impl AsRef<Path>,
    ) -> Result<()> {
        let p_top = prmtop_path.as_ref();
        let top_content = fs::read_to_string(p_top).map_err(|e| {
            BridgeError::input_error(
                p_top.display().to_string(),
                format!("failed to read prmtop file: {}", e),
            )
        })?;

        let p_crd = inpcrd_path.as_ref();
        let crd_content = fs::read_to_string(p_crd).map_err(|e| {
            BridgeError::input_error(
                p_crd.display().to_string(),
                format!("failed to read inpcrd file: {}", e),
            )
        })?;

        self.parse_prmtop_str(&top_content)?;
        self.parse_inpcrd_str(&crd_content)?;
        self.validate_data()?;
        Ok(())
    }

    /// Parses the PRMTOP format string.
    pub fn parse_prmtop_str(&mut self, content: &str) -> Result<()> {
        self.bonds.clear();
        let mut lines = content.lines().peekable();

        while let Some(line) = lines.next() {
            let trimmed = line.trim_end();
            if trimmed == "%FLAG ATOM_NAME" {
                self.atom_names = Self::read_atom_names(&mut lines)?;
            } else if trimmed == "%FLAG CHARGE" {
                self.charges = Self::read_charges(&mut lines)?;
            } else if trimmed == "%FLAG ATOMIC_NUMBER" {
                self.atomic_numbers = Self::read_atomic_numbers(&mut lines)?;
            } else if trimmed == "%FLAG BONDS_INC_HYDROGEN"
                || trimmed == "%FLAG BONDS_WITHOUT_HYDROGEN"
            {
                let section_bonds = Self::read_bonds(&mut lines)?;
                self.bonds.extend(section_bonds);
            }
        }

        Ok(())
    }

    fn read_atom_names<'a>(
        lines: &mut std::iter::Peekable<impl Iterator<Item = &'a str>>,
    ) -> Result<Vec<String>> {
        let mut names = Vec::new();
        while let Some(&line) = lines.peek() {
            let trimmed = line.trim_end();
            if trimmed.starts_with("%FLAG") {
                break;
            }
            lines.next();
            if trimmed.starts_with("%FORMAT") {
                continue;
            }

            let chars: Vec<char> = trimmed.chars().collect();
            let mut start = 0;
            while start < chars.len() {
                let end = (start + 4).min(chars.len());
                let name: String = chars[start..end].iter().collect();
                names.push(name.trim_end().to_string());
                start += 4;
            }
        }
        Ok(names)
    }

    fn read_charges<'a>(
        lines: &mut std::iter::Peekable<impl Iterator<Item = &'a str>>,
    ) -> Result<Vec<f64>> {
        let mut charges = Vec::new();
        while let Some(&line) = lines.peek() {
            let trimmed = line.trim_end();
            if trimmed.starts_with("%FLAG") {
                break;
            }
            lines.next();
            if trimmed.starts_with("%FORMAT") {
                continue;
            }

            for val_str in trimmed.split_whitespace() {
                let val = val_str.parse::<f64>().map_err(|_| {
                    BridgeError::input_error(val_str, "failed to parse charge in PRMTOP")
                })?;
                charges.push(val / AMBER_CHARGE_FACTOR);
            }
        }
        Ok(charges)
    }

    fn read_atomic_numbers<'a>(
        lines: &mut std::iter::Peekable<impl Iterator<Item = &'a str>>,
    ) -> Result<Vec<usize>> {
        let mut numbers = Vec::new();
        while let Some(&line) = lines.peek() {
            let trimmed = line.trim_end();
            if trimmed.starts_with("%FLAG") {
                break;
            }
            lines.next();
            if trimmed.starts_with("%FORMAT") {
                continue;
            }

            for val_str in trimmed.split_whitespace() {
                let num = val_str.parse::<usize>().map_err(|_| {
                    BridgeError::input_error(val_str, "failed to parse atomic number in PRMTOP")
                })?;
                numbers.push(num);
            }
        }
        Ok(numbers)
    }

    fn read_bonds<'a>(
        lines: &mut std::iter::Peekable<impl Iterator<Item = &'a str>>,
    ) -> Result<Vec<(usize, usize)>> {
        let mut raw_ints = Vec::new();
        while let Some(&line) = lines.peek() {
            let trimmed = line.trim_end();
            if trimmed.starts_with("%FLAG") {
                break;
            }
            lines.next();
            if trimmed.starts_with("%FORMAT") {
                continue;
            }

            for val_str in trimmed.split_whitespace() {
                let val = val_str.parse::<usize>().map_err(|_| {
                    BridgeError::input_error(val_str, "failed to parse bond index in PRMTOP")
                })?;
                raw_ints.push(val);
            }
        }

        if raw_ints.len() % 3 != 0 {
            return Err(BridgeError::input_error(
                "BONDS",
                format!(
                    "invalid number of entries in bond section: expected multiple of 3, got {}",
                    raw_ints.len()
                ),
            ));
        }

        let mut bonds = Vec::with_capacity(raw_ints.len() / 3);
        for chunk in raw_ints.as_chunks::<3>().0 {
            let offset1 = chunk[0];
            let offset2 = chunk[1];
            // chunk[2] is bond_type_idx (force field parameter index)

            if offset1 % 3 != 0 || offset2 % 3 != 0 {
                return Err(BridgeError::input_error(
                    "BONDS",
                    format!(
                        "invalid coordinate offset in bond entry: ({}, {}) must be multiples of 3",
                        offset1, offset2
                    ),
                ));
            }

            let atom1 = offset1 / 3;
            let atom2 = offset2 / 3;
            bonds.push((atom1, atom2));
        }

        Ok(bonds)
    }

    /// Parses the INPCRD format string.
    pub fn parse_inpcrd_str(&mut self, content: &str) -> Result<()> {
        let mut lines = content.lines();

        // Line 1: Title
        let _title = lines
            .next()
            .ok_or_else(|| BridgeError::input_error("INPCRD", "empty content"))?;

        // Line 2: Number of atoms
        let count_line = lines
            .next()
            .ok_or_else(|| BridgeError::input_error("INPCRD", "missing atom count line"))?
            .trim();
        let first_token = count_line
            .split_whitespace()
            .next()
            .ok_or_else(|| BridgeError::input_error(count_line, "invalid atom count"))?;
        let num_of_atoms = first_token
            .parse::<usize>()
            .map_err(|_| BridgeError::input_error(first_token, "invalid atom count integer"))?;

        let mut positions = Vec::with_capacity(num_of_atoms);

        for line in lines {
            let trimmed = line.trim_end();
            if trimmed.is_empty() {
                continue;
            }
            let chars: Vec<char> = trimmed.chars().collect();
            let mut offset = 0;

            while offset + 36 <= chars.len() && positions.len() < num_of_atoms {
                let x_str: String = chars[offset..offset + 12].iter().collect();
                let y_str: String = chars[offset + 12..offset + 24].iter().collect();
                let z_str: String = chars[offset + 24..offset + 36].iter().collect();

                let x = x_str
                    .trim()
                    .parse::<f64>()
                    .map_err(|_| BridgeError::input_error(&x_str, "failed to parse x in INPCRD"))?;
                let y = y_str
                    .trim()
                    .parse::<f64>()
                    .map_err(|_| BridgeError::input_error(&y_str, "failed to parse y in INPCRD"))?;
                let z = z_str
                    .trim()
                    .parse::<f64>()
                    .map_err(|_| BridgeError::input_error(&z_str, "failed to parse z in INPCRD"))?;

                positions.push(Position::new(x, y, z));
                offset += 36;
            }

            if positions.len() >= num_of_atoms {
                break;
            }
        }

        self.xyz = positions;
        Ok(())
    }

    fn validate_data(&self) -> Result<()> {
        let n_atoms = self.xyz.len();
        if !self.atom_names.is_empty() && self.atom_names.len() != n_atoms {
            return Err(BridgeError::input_error(
                "AmberPrmtop",
                format!(
                    "mismatch between atom_names count ({}) and xyz count ({})",
                    self.atom_names.len(),
                    n_atoms
                ),
            ));
        }
        if !self.charges.is_empty() && self.charges.len() != n_atoms {
            return Err(BridgeError::input_error(
                "AmberPrmtop",
                format!(
                    "mismatch between charges count ({}) and xyz count ({})",
                    self.charges.len(),
                    n_atoms
                ),
            ));
        }
        if !self.atomic_numbers.is_empty() && self.atomic_numbers.len() != n_atoms {
            return Err(BridgeError::input_error(
                "AmberPrmtop",
                format!(
                    "mismatch between atomic_numbers count ({}) and xyz count ({})",
                    self.atomic_numbers.len(),
                    n_atoms
                ),
            ));
        }
        for &(a1, a2) in &self.bonds {
            if n_atoms > 0 && (a1 >= n_atoms || a2 >= n_atoms) {
                return Err(BridgeError::input_error(
                    "AmberPrmtop",
                    format!(
                        "bond index out of bounds: ({}, {}) for atom count {}",
                        a1, a2, n_atoms
                    ),
                ));
            }
        }
        Ok(())
    }

    /// Returns the parsed atom names.
    pub fn atom_names(&self) -> &[String] {
        &self.atom_names
    }

    /// Returns the parsed charges.
    pub fn charges(&self) -> &[f64] {
        &self.charges
    }

    /// Returns the parsed atomic numbers.
    pub fn atomic_numbers(&self) -> &[usize] {
        &self.atomic_numbers
    }

    /// Returns the parsed coordinates.
    pub fn xyz(&self) -> &[Position] {
        &self.xyz
    }

    /// Returns the parsed bonds as `(atom1_idx, atom2_idx)` pairs.
    pub fn bonds(&self) -> &[(usize, usize)] {
        &self.bonds
    }

    /// Converts the parsed Amber data into an `AtomGroup`.
    ///
    /// If the PRMTOP file contains `BONDS_WITHOUT_HYDROGEN` or `BONDS_INC_HYDROGEN`
    /// sections, the returned `AtomGroup` includes the explicit bond topology parsed from
    /// the file. According to the bond priority policy (see `RUST_PORT_SPEC.md` §3.8 & §3.14),
    /// explicit file-derived bonds take precedence. If no explicit bond records exist,
    /// [`AtomGroup::setup`] is automatically executed to resolve bonds.
    pub fn get_atomgroup(&self) -> Result<AtomGroup> {
        let mut atomgroup = AtomGroup::new();
        let num_of_atoms = self.xyz.len();

        for i in 0..num_of_atoms {
            let mut atom = Atom::new();
            if i < self.atomic_numbers.len() {
                let symbol = PeriodicTable::get_symbol(self.atomic_numbers[i])?;
                atom.set_symbol(symbol)?;
            }
            if i < self.xyz.len() {
                atom.xyz = self.xyz[i];
            }
            if i < self.charges.len() {
                atom.charge = self.charges[i];
            }
            if i < self.atom_names.len() {
                atom.name = self.atom_names[i].clone();
            }
            atomgroup.set_atom(&i.to_string(), atom);
        }

        for &(idx1, idx2) in &self.bonds {
            let key1 = idx1.to_string();
            let key2 = idx2.to_string();
            if let (Some(a1), Some(a2)) = (atomgroup.get_atom(&key1), atomgroup.get_atom(&key2)) {
                let a1_clone = a1.clone();
                let a2_clone = a2.clone();
                atomgroup.add_bond(&a1_clone, &a2_clone, 1);
            }
        }

        if atomgroup.get_bond_list().is_empty() {
            atomgroup.setup()?;
        }

        Ok(atomgroup)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const PRMTOP_CONTENT: &str = "\
%VERSION  VERSION_STAMP = V0001.000  DATE = 08/25/26  12:00:00
%FLAG ATOM_NAME
%FORMAT(20a4)
C1  H1  
%FLAG CHARGE
%FORMAT(5E16.8)
 0.00000000E+00 0.00000000E+00
%FLAG ATOMIC_NUMBER
%FORMAT(10I8)
       6       1
";

    const INPCRD_CONTENT: &str = "\
default_name
    2
   0.0000000   0.0000000   0.0000000   1.0900000   0.0000000   0.0000000
";

    // Ported from tests/test_amber_prmtop.py
    #[test]
    fn test_load_prmtop_and_inpcrd() {
        let amber = AmberPrmtop::from_strings(PRMTOP_CONTENT, INPCRD_CONTENT).unwrap();

        assert_eq!(amber.atom_names().len(), 2);
        assert_eq!(amber.atom_names()[0], "C1");
        assert_eq!(amber.atom_names()[1], "H1");
        assert_eq!(amber.atomic_numbers(), &[6, 1]);
        assert_eq!(amber.xyz().len(), 2);

        let ag = amber.get_atomgroup().unwrap();
        assert_eq!(ag.get_number_of_atoms(), 2);

        let c1 = ag.get_atom("0").unwrap();
        assert_eq!(c1.name, "C1");
        assert_eq!(c1.symbol().unwrap(), "C");
        assert_eq!(c1.xyz.x, 0.0);

        let h1 = ag.get_atom("1").unwrap();
        assert_eq!(h1.name, "H1");
        assert_eq!(h1.symbol().unwrap(), "H");
        assert!((h1.xyz.x - 1.09).abs() < 1e-5);
    }

    #[test]
    fn test_from_files() {
        let temp_dir = std::env::temp_dir();
        let top_path = temp_dir.join("test.prmtop");
        let crd_path = temp_dir.join("test.inpcrd");

        fs::write(&top_path, PRMTOP_CONTENT).unwrap();
        fs::write(&crd_path, INPCRD_CONTENT).unwrap();

        let amber = AmberPrmtop::from_files(&top_path, &crd_path).unwrap();
        assert_eq!(amber.atom_names().len(), 2);
        assert_eq!(amber.xyz().len(), 2);

        let _ = fs::remove_file(&top_path);
        let _ = fs::remove_file(&crd_path);
    }

    const PRMTOP_WITH_BONDS: &str = "\
%VERSION  VERSION_STAMP = V0001.000  DATE = 08/25/26  12:00:00
%FLAG ATOM_NAME
%FORMAT(20a4)
C0  C1  H2  H3  
%FLAG CHARGE
%FORMAT(5E16.8)
 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
%FLAG ATOMIC_NUMBER
%FORMAT(10I8)
       6       6       1       1
%FLAG BONDS_WITHOUT_HYDROGEN
%FORMAT(10I8)
       0       3       1
%FLAG BONDS_INC_HYDROGEN
%FORMAT(10I8)
       0       6       2       3       9       2
";

    const INPCRD_4ATOMS: &str = "\
default_name
    4
   0.0000000   0.0000000   0.0000000   1.5000000   0.0000000   0.0000000
  -0.5000000   0.9000000   0.0000000   2.0000000   0.9000000   0.0000000
";

    #[test]
    fn test_load_prmtop_with_bonds() {
        let amber = AmberPrmtop::from_strings(PRMTOP_WITH_BONDS, INPCRD_4ATOMS).unwrap();

        assert_eq!(amber.atom_names().len(), 4);
        assert_eq!(amber.bonds().len(), 3);
        assert_eq!(amber.bonds(), &[(0, 1), (0, 2), (1, 3)]);

        let mut ag = amber.get_atomgroup().unwrap();
        assert_eq!(ag.get_number_of_atoms(), 4);
        assert_eq!(ag.get_number_of_bonds(), 3);

        let bond_list = ag.get_bond_list();
        assert_eq!(bond_list.len(), 3);

        // Verify that bonds resolve correctly to the expected atoms
        for record in &bond_list {
            let (a1, a2) = ag.resolve_bond(record).expect("bond must resolve");
            match (a1.name.as_str(), a2.name.as_str()) {
                ("C0", "C1") | ("C1", "C0") => assert_eq!(record.order, 1),
                ("C0", "H2") | ("H2", "C0") => assert_eq!(record.order, 1),
                ("C1", "H3") | ("H3", "C1") => assert_eq!(record.order, 1),
                other => panic!("unexpected bond pair: {:?}", other),
            }
        }
    }

    #[test]
    fn test_bonds_boundary_first_and_last_atoms() {
        // System with 5 atoms: 0, 1, 2, 3, 4.
        // Bond connects first atom (0, offset 0) and last atom (4, offset 12).
        let prmtop = "\
%FLAG ATOM_NAME
%FORMAT(20a4)
A0  A1  A2  A3  A4  
%FLAG CHARGE
%FORMAT(5E16.8)
 0.0 0.0 0.0 0.0 0.0
%FLAG ATOMIC_NUMBER
%FORMAT(10I8)
       6       6       6       6       6
%FLAG BONDS_WITHOUT_HYDROGEN
%FORMAT(10I8)
       0      12       1
";
        let inpcrd = "\
boundary_test
    5
   0.0000000   0.0000000   0.0000000   1.0000000   0.0000000   0.0000000
   2.0000000   0.0000000   0.0000000   3.0000000   0.0000000   0.0000000
   4.0000000   0.0000000   0.0000000
";
        let amber = AmberPrmtop::from_strings(prmtop, inpcrd).unwrap();
        assert_eq!(amber.bonds(), &[(0, 4)]);

        let mut ag = amber.get_atomgroup().unwrap();
        assert_eq!(ag.get_number_of_bonds(), 1);
        let bond = &ag.get_bond_list()[0];
        let (a1, a2) = ag.resolve_bond(bond).unwrap();
        assert_eq!(a1.name, "A0");
        assert_eq!(a2.name, "A4");
    }

    #[test]
    fn test_bonds_invalid_offset_not_multiple_of_three() {
        let prmtop = "\
%FLAG ATOM_NAME
%FORMAT(20a4)
C1  H1  
%FLAG BONDS_WITHOUT_HYDROGEN
%FORMAT(10I8)
       0       2       1
";
        let mut amber = AmberPrmtop::new();
        let result = amber.parse_prmtop_str(prmtop);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("must be multiples of 3"));
    }

    #[test]
    fn test_bonds_invalid_entry_count() {
        let prmtop = "\
%FLAG ATOM_NAME
%FORMAT(20a4)
C1  H1  
%FLAG BONDS_WITHOUT_HYDROGEN
%FORMAT(10I8)
       0       3
";
        let mut amber = AmberPrmtop::new();
        let result = amber.parse_prmtop_str(prmtop);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("expected multiple of 3"));
    }

    #[test]
    fn test_bonds_out_of_bounds() {
        let prmtop = "\
%FLAG ATOM_NAME
%FORMAT(20a4)
C1  H1  
%FLAG ATOMIC_NUMBER
%FORMAT(10I8)
       6       1
%FLAG BONDS_WITHOUT_HYDROGEN
%FORMAT(10I8)
       0       6       1
";
        let inpcrd = "\
test
    2
   0.0000000   0.0000000   0.0000000   1.0000000   0.0000000   0.0000000
";
        let result = AmberPrmtop::from_strings(prmtop, inpcrd);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("bond index out of bounds"));
    }
}
