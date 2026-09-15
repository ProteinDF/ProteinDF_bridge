// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::fmt;

use crate::error::Result;
use crate::matrix::Matrix;
use crate::periodic_table::PeriodicTable;
use crate::position::Position;

/// Represents an atom corresponding to `proteindf_bridge.atom.Atom`.
#[derive(Debug, Clone)]
pub struct Atom {
    atomic_number: usize,
    pub xyz: Position,
    pub force: Position,
    pub name: String,
    pub label: String,
    pub charge: f64,
    pub path: String,
}

impl Default for Atom {
    fn default() -> Self {
        Self {
            atomic_number: PeriodicTable::get_atomic_number("X").unwrap_or(0),
            xyz: Position::default(),
            force: Position::default(),
            name: String::new(),
            label: String::new(),
            charge: 0.0,
            path: String::new(),
        }
    }
}

impl Atom {
    /// Creates a new default Atom (symbol "X", atomic number 0).
    pub fn new() -> Self {
        Self::default()
    }

    /// Creates an Atom with the specified symbol.
    pub fn from_symbol(symbol: &str) -> Result<Self> {
        let atomic_number = PeriodicTable::get_atomic_number(symbol)?;
        Ok(Self {
            atomic_number,
            ..Default::default()
        })
    }

    /// Creates an Atom with atomic number and position.
    pub fn new_with_pos(symbol: &str, xyz: Position) -> Result<Self> {
        let atomic_number = PeriodicTable::get_atomic_number(symbol)?;
        Ok(Self {
            atomic_number,
            xyz,
            ..Default::default()
        })
    }

    /// Returns the atomic number.
    pub fn atomic_number(&self) -> usize {
        self.atomic_number
    }

    /// Sets the atomic number.
    pub fn set_atomic_number(&mut self, num: usize) {
        self.atomic_number = num;
    }

    /// Returns the element symbol.
    pub fn symbol(&self) -> Result<&'static str> {
        PeriodicTable::get_symbol(self.atomic_number)
    }

    /// Sets the element symbol.
    pub fn set_symbol(&mut self, symbol: &str) -> Result<()> {
        self.atomic_number = PeriodicTable::get_atomic_number(symbol)?;
        Ok(())
    }

    /// Returns true if the atom is a real element (atomic number > 0).
    pub fn is_real(&self) -> bool {
        self.atomic_number > 0
    }

    /// Returns the atomic weight.
    pub fn weight(&self) -> Result<f64> {
        PeriodicTable::atomic_weight(self.atomic_number)
    }

    /// Returns the van der Waals radius.
    pub fn vdw(&self) -> Result<f64> {
        PeriodicTable::vdw(self.atomic_number)
    }

    /// Moves the atom to the given position.
    pub fn move_to(&mut self, position: Position) -> &mut Self {
        self.xyz = position;
        self
    }

    /// Shifts the atom's position by the given displacement.
    pub fn shift_by(&mut self, direction: Position) -> &mut Self {
        self.xyz += direction;
        self
    }

    /// Rotates the atom's position by a 3x3 rotation matrix.
    pub fn rotate(&mut self, rotmat: &Matrix) -> Result<()> {
        self.xyz.rotate(rotmat)
    }
}

// Equality matching Python's `Atom.__eq__`:
// self.atomic_number == rhs.atomic_number and self.xyz == rhs.xyz
impl PartialEq for Atom {
    fn eq(&self, other: &Self) -> bool {
        self.atomic_number == other.atomic_number && self.xyz == other.xyz
    }
}

impl fmt::Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let sym = self.symbol().unwrap_or("?");
        write!(
            f,
            "{:<2}({:<4}) {:8.3} {:8.3} {:8.3}, {:5.2}, {:8.3} {:8.3} {:8.3} {}",
            sym,
            self.name,
            self.xyz.x,
            self.xyz.y,
            self.xyz.z,
            self.charge,
            self.force.x,
            self.force.y,
            self.force.z,
            self.path
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Ported from tests/test_atom.py
    #[test]
    fn test_init() {
        let mut atom = Atom::new();
        atom.set_symbol("Fe").unwrap();
        assert_eq!(atom.atomic_number(), 26);
        assert_eq!(atom.symbol().unwrap(), "Fe");
    }

    #[test]
    fn test_init2() {
        let atom = Atom::new_with_pos("Ni", Position::new(0.0, 1.0, 2.0)).unwrap();
        assert_eq!(atom.atomic_number(), 28);
        assert!((atom.xyz.x - 0.0).abs() < 1e-10);
        assert!((atom.xyz.y - 1.0).abs() < 1e-10);
        assert!((atom.xyz.z - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_init3() {
        let pos: Position = "0.0 1.0 2.0".parse().unwrap();
        let atom = Atom::new_with_pos("C", pos).unwrap();
        assert_eq!(atom.atomic_number(), 6);
        assert!((atom.xyz.x - 0.0).abs() < 1e-10);
        assert!((atom.xyz.y - 1.0).abs() < 1e-10);
        assert!((atom.xyz.z - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_doctests() {
        let mut a = Atom::new();
        a.set_symbol("Fe").unwrap();
        assert_eq!(a.atomic_number(), 26);
        assert_eq!(a.symbol().unwrap(), "Fe");

        a.charge = -0.2;
        assert!((a.charge - (-0.2)).abs() < 1e-10);

        let mut b = a.clone();
        b.set_symbol("Na").unwrap();
        assert_eq!(a.symbol().unwrap(), "Fe");
        assert_eq!(b.symbol().unwrap(), "Na");
    }

    #[test]
    fn test_is_real_and_vdw() {
        let h = Atom::from_symbol("H").unwrap();
        assert!(h.is_real());
        assert!((h.vdw().unwrap() - 1.2).abs() < 1e-5);
        assert!((h.weight().unwrap() - 1.008).abs() < 1e-5);

        let x = Atom::new();
        assert!(!x.is_real());
    }
}
