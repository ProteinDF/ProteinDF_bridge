// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

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

/// Covalent atomic radii in Angstroms from:
/// Cordero, B. et al. "Covalent radii revisited." Dalton Trans. 2008, 2832–2838.
/// DOI: 10.1039/B801115J
/// Indexed by atomic number (index 0 is dummy 'X' with value 0.0).
pub const COVALENT_RADIUS: &[f64] = &[
    0.0,  // X (dummy)
    0.31, // 1: H
    0.28, // 2: He
    1.28, // 3: Li
    0.96, // 4: Be
    0.84, // 5: B
    0.76, // 6: C
    0.71, // 7: N
    0.66, // 8: O
    0.57, // 9: F
    0.58, // 10: Ne
    1.66, // 11: Na
    1.41, // 12: Mg
    1.21, // 13: Al
    1.11, // 14: Si
    1.07, // 15: P
    1.05, // 16: S
    1.02, // 17: Cl
    1.06, // 18: Ar
    2.03, // 19: K
    1.76, // 20: Ca
    1.70, // 21: Sc
    1.60, // 22: Ti
    1.53, // 23: V
    1.39, // 24: Cr
    1.39, // 25: Mn
    1.32, // 26: Fe
    1.26, // 27: Co
    1.24, // 28: Ni
    1.32, // 29: Cu
    1.22, // 30: Zn
    1.22, // 31: Ga
    1.20, // 32: Ge
    1.19, // 33: As
    1.20, // 34: Se
    1.20, // 35: Br
    1.16, // 36: Kr
    2.20, // 37: Rb
    1.95, // 38: Sr
    1.90, // 39: Y
    1.75, // 40: Zr
    1.64, // 41: Nb
    1.54, // 42: Mo
    1.47, // 43: Tc
    1.46, // 44: Ru
    1.42, // 45: Rh
    1.39, // 46: Pd
    1.45, // 47: Ag
    1.44, // 48: Cd
    1.42, // 49: In
    1.39, // 50: Sn
    1.39, // 51: Sb
    1.38, // 52: Te
    1.39, // 53: I
    1.40, // 54: Xe
    2.44, // 55: Cs
    2.15, // 56: Ba
    2.07, // 57: La
    2.04, // 58: Ce
    2.03, // 59: Pr
    2.01, // 60: Nd
    1.99, // 61: Pm
    1.98, // 62: Sm
    1.98, // 63: Eu
    1.96, // 64: Gd
    1.94, // 65: Tb
    1.92, // 66: Dy
    1.92, // 67: Ho
    1.89, // 68: Er
    1.90, // 69: Tm
    1.87, // 70: Yb
    1.87, // 71: Lu
    1.75, // 72: Hf
    1.70, // 73: Ta
    1.62, // 74: W
    1.51, // 75: Re
    1.44, // 76: Os
    1.41, // 77: Ir
    1.36, // 78: Pt
    1.36, // 79: Au
    1.32, // 80: Hg
    1.45, // 81: Tl
    1.46, // 82: Pb
    1.48, // 83: Bi
    1.40, // 84: Po
    1.50, // 85: At
    1.50, // 86: Rn
    2.60, // 87: Fr
    2.21, // 88: Ra
    2.15, // 89: Ac
    2.06, // 90: Th
    2.00, // 91: Pa
    1.96, // 92: U
    1.90, // 93: Np
    1.87, // 94: Pu
    1.80, // 95: Am
    1.69, // 96: Cm
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

    /// Returns the covalent radius in Angstroms for a given atom (atomic number or symbol).
    ///
    /// Values are from Cordero et al., Dalton Trans. 2008, 2832–2838 (DOI: 10.1039/B801115J).
    pub fn covalent_radius(atom: impl IntoAtomId) -> Result<f64> {
        let num = atom.to_atomic_number()?;
        if num == 0 {
            return Err(BridgeError::CovalentRadiusNotFound(0));
        }
        COVALENT_RADIUS
            .get(num)
            .copied()
            .ok_or(BridgeError::CovalentRadiusNotFound(num))
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
