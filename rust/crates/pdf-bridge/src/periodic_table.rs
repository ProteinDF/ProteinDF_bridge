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

use crate::error::{BridgeError, Result};

pub struct PeriodicTable;

pub const TABLE: &[&str] = &[
    "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
    "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Fl", "Lv", "Uup", "Uuh", "Uus", "Uuo",
];

pub const ATOMIC_WEIGHTS: &[f64] = &[
    0.0, 1.008, 4.003, 6.941, 9.012, 10.81, 12.01, 14.01, 16.00, 19.00, 20.18, 22.00, 24.31, 26.98,
    28.09, 30.97, 32.07, 35.45, 39.95, 39.10, 40.08, 44.96, 47.87, 50.94, 52.00, 54.94, 55.85,
    58.93, 58.69, 63.55, 65.38, 69.72, 72.63, 74.92, 78.97, 79.90, 83.80, 85.47, 87.62, 88.91,
    91.22, 92.91, 95.95, 99.00, 101.1, 102.9, 106.4, 107.9, 112.4, 114.8, 118.7, 121.8, 127.6,
    169.9, 131.3, 132.9, 137.3, 138.9, 140.1, 140.9, 144.2, 145.0, 150.4, 152.0, 157.3, 158.9,
    162.5, 164.9, 167.3, 168.9, 173.1, 175.0, 178.5, 180.9, 183.8, 186.2, 190.2, 192.2, 195.1,
    197.0, 200.6, 204.4, 207.2, 209.0, 210.0, 210.0, 222.0, 223.0, 226.0, 227.0, 232.0, 231.0,
    238.0, 237.0, 239.0, 243.0, 247.0, 252.0, 252.0, 257.0, 258.0, 259.0, 262.0, 267.0, 268.0,
    271.0, 272.0, 277.0, 276.0, 281.0, 280.0, 285.0, 289.0, 293.0,
];

pub const VDW: &[f64] = &[
    0.0, 1.2, 1.4, 1.82, 0.0, 0.0, 1.70, 1.55, 1.52, 1.47, 1.54, 2.27, 1.73, 0.0, 2.10, 1.80, 1.80,
    1.75, 1.88, 2.75, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.63, 1.40, 1.39, 1.87, 0.0, 1.85,
    2.02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.63, 1.72, 1.58, 1.93, 2.17, 0.0, 2.06,
    1.98, 2.16,
];

/// Helper trait to accept either an atomic number (usize) or a symbol (&str) as an atom identifier.
pub trait IntoAtomId {
    fn to_atomic_number(self) -> Result<usize>;
}

impl IntoAtomId for usize {
    fn to_atomic_number(self) -> Result<usize> {
        Ok(self)
    }
}

impl IntoAtomId for &str {
    fn to_atomic_number(self) -> Result<usize> {
        PeriodicTable::get_atomic_number(self)
    }
}

impl IntoAtomId for &String {
    fn to_atomic_number(self) -> Result<usize> {
        PeriodicTable::get_atomic_number(self.as_str())
    }
}

impl IntoAtomId for String {
    fn to_atomic_number(self) -> Result<usize> {
        PeriodicTable::get_atomic_number(self.as_str())
    }
}

fn normalize_symbol(symbol: &str) -> String {
    let mut chars = symbol.chars();
    match chars.next() {
        None => String::new(),
        Some(first) => {
            let mut s = first.to_uppercase().to_string();
            for c in chars {
                s.extend(c.to_lowercase());
            }
            s
        }
    }
}

impl PeriodicTable {
    /// Returns the total number of entries in the periodic table (including 'X' at index 0).
    pub fn get_num_of_atoms() -> usize {
        TABLE.len()
    }

    /// Returns the element symbol for a given atomic number.
    pub fn get_symbol(atomic_number: usize) -> Result<&'static str> {
        TABLE
            .get(atomic_number)
            .copied()
            .ok_or(BridgeError::AtomicNumberNotFound(atomic_number))
    }

    /// Returns the atomic number for a given element symbol.
    /// Case-insensitive (e.g. "ca", "CA", "Ca" all match "Ca").
    pub fn get_atomic_number(symbol: &str) -> Result<usize> {
        let normalized = normalize_symbol(symbol);
        TABLE
            .iter()
            .position(|&s| s == normalized)
            .ok_or_else(|| BridgeError::SymbolNotFound(symbol.to_string()))
    }

    /// Checks whether the given symbol is present in the periodic table.
    pub fn contains(symbol: &str) -> bool {
        let normalized = normalize_symbol(symbol);
        TABLE.iter().any(|&s| s == normalized)
    }

    /// Returns the van der Waals radius for a given atom (atomic number or symbol).
    pub fn vdw(atom: impl IntoAtomId) -> Result<f64> {
        let num = atom.to_atomic_number()?;
        VDW.get(num)
            .copied()
            .ok_or(BridgeError::VdwRadiusNotFound(num))
    }

    /// Returns the atomic weight for a given atom (atomic number or symbol).
    pub fn atomic_weight(atom: impl IntoAtomId) -> Result<f64> {
        let num = atom.to_atomic_number()?;
        ATOMIC_WEIGHTS
            .get(num)
            .copied()
            .ok_or(BridgeError::AtomicWeightNotFound(num))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Ported from tests/test_periodictable.py
    #[test]
    fn test_get_symbol() {
        assert_eq!(PeriodicTable::get_symbol(1).unwrap(), "H");
        assert_eq!(PeriodicTable::get_symbol(6).unwrap(), "C");
        assert_eq!(PeriodicTable::get_symbol(8).unwrap(), "O");
        assert_eq!(PeriodicTable::get_symbol(20).unwrap(), "Ca");
    }

    #[test]
    fn test_get_atomic_number() {
        assert_eq!(PeriodicTable::get_atomic_number("H").unwrap(), 1);
        assert_eq!(PeriodicTable::get_atomic_number("C").unwrap(), 6);
        assert_eq!(PeriodicTable::get_atomic_number("N").unwrap(), 7);
        assert_eq!(PeriodicTable::get_atomic_number("Cu").unwrap(), 29);

        // Case insensitivity tests
        assert_eq!(PeriodicTable::get_atomic_number("h").unwrap(), 1);
        assert_eq!(PeriodicTable::get_atomic_number("cu").unwrap(), 29);
        assert_eq!(PeriodicTable::get_atomic_number("CU").unwrap(), 29);
    }

    #[test]
    fn test_get_weight_and_vdw() {
        let weight = PeriodicTable::atomic_weight(6).unwrap();
        assert!((weight - 12.01).abs() < 0.01);

        let vdw = PeriodicTable::vdw(6).unwrap();
        assert!((vdw - 1.70).abs() < 0.01);

        // Test with symbol input
        let weight_c = PeriodicTable::atomic_weight("C").unwrap();
        assert_eq!(weight, weight_c);

        let vdw_c = PeriodicTable::vdw("C").unwrap();
        assert_eq!(vdw, vdw_c);
    }

    #[test]
    fn test_get_num_of_atoms_and_contains() {
        assert_eq!(PeriodicTable::get_num_of_atoms(), 119);
        assert!(PeriodicTable::contains("H"));
        assert!(PeriodicTable::contains("he"));
        assert!(!PeriodicTable::contains("Invalid"));
    }

    #[test]
    fn test_error_handling() {
        assert!(matches!(
            PeriodicTable::get_symbol(999),
            Err(BridgeError::AtomicNumberNotFound(999))
        ));
        assert!(matches!(
            PeriodicTable::get_atomic_number("ZZ"),
            Err(BridgeError::SymbolNotFound(_))
        ));
        assert!(matches!(
            PeriodicTable::vdw(999),
            Err(BridgeError::VdwRadiusNotFound(999))
        ));
        assert!(matches!(
            PeriodicTable::atomic_weight(999),
            Err(BridgeError::AtomicWeightNotFound(999))
        ));
    }
}
