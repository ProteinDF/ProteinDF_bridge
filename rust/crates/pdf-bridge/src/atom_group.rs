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

use std::collections::{HashMap, HashSet};
use std::fmt;
use std::ops::{BitAnd, BitAndAssign, BitOr, BitOrAssign, BitXor, BitXorAssign, Index};

use indexmap::IndexMap;

use crate::atom::Atom;
use crate::error::Result;
use crate::matrix::Matrix;
use crate::periodic_table::PeriodicTable;
use crate::position::Position;

/// A bond record storing two atom paths and the bond order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BondRecord {
    pub atom1_path: String,
    pub atom2_path: String,
    pub order: usize,
}

/// Trait for selecting atoms or groups matching a predicate.
pub trait Selector {
    fn is_match_group(&self, _group: &AtomGroup) -> bool {
        false
    }
    fn is_match_atom(&self, atom: &Atom) -> bool;
}

/// Hierarchical atom group corresponding to `proteindf_bridge.atomgroup.AtomGroup`.
#[derive(Debug, Clone, PartialEq)]
pub struct AtomGroup {
    pub name: String,
    path: String,
    atoms: IndexMap<String, Atom>,
    groups: IndexMap<String, AtomGroup>,
    bonds: Vec<BondRecord>,
}

impl Default for AtomGroup {
    fn default() -> Self {
        Self {
            name: String::new(),
            path: "/".to_string(),
            atoms: IndexMap::new(),
            groups: IndexMap::new(),
            bonds: Vec::new(),
        }
    }
}

impl AtomGroup {
    /// Creates a new empty `AtomGroup`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Creates a new `AtomGroup` with the given name.
    pub fn with_name(name: &str) -> Self {
        Self {
            name: name.to_string(),
            ..Default::default()
        }
    }

    /// Returns the group's current path (e.g. "/" or "/grp/").
    pub fn path(&self) -> &str {
        &self.path
    }

    /// Sets the path of this group and updates descendant paths.
    pub fn set_path(&mut self, mut new_path: String) {
        if new_path.is_empty() || !new_path.ends_with('/') {
            new_path.push('/');
        }
        self.path = new_path;
        self.update_paths();
    }

    fn update_paths(&mut self) {
        for (key, group) in self.groups.iter_mut() {
            group.set_path(format!("{}{}/", self.path, key));
        }
        for (key, atom) in self.atoms.iter_mut() {
            atom.path = format!("{}{}", self.path, key);
        }
    }

    /// Returns the number of direct child groups.
    pub fn get_number_of_groups(&self) -> usize {
        self.groups.len()
    }

    /// Returns the number of direct child atoms.
    pub fn get_number_of_atoms(&self) -> usize {
        self.atoms.len()
    }

    /// Returns the total number of atoms across all subgroups recursively.
    pub fn get_number_of_all_atoms(&self) -> usize {
        let mut count = self.atoms.len();
        for group in self.groups.values() {
            count += group.get_number_of_all_atoms();
        }
        count
    }

    /// Returns the sum of atomic numbers of all atoms recursively.
    pub fn sum_of_atomic_number(&self) -> f64 {
        let mut sum = 0.0;
        for atom in self.atoms.values() {
            sum += atom.atomic_number() as f64;
        }
        for group in self.groups.values() {
            sum += group.sum_of_atomic_number();
        }
        sum
    }

    /// Returns an iterator over direct child groups: `(key, &AtomGroup)`.
    pub fn groups(&self) -> impl Iterator<Item = (&String, &AtomGroup)> {
        self.groups.iter()
    }

    /// Returns an iterator over mutable direct child groups: `(key, &mut AtomGroup)`.
    pub fn groups_mut(&mut self) -> impl Iterator<Item = (&String, &mut AtomGroup)> {
        self.groups.iter_mut()
    }

    /// Returns an iterator over direct child atoms: `(key, &Atom)`.
    pub fn atoms(&self) -> impl Iterator<Item = (&String, &Atom)> {
        self.atoms.iter()
    }

    /// Returns an iterator over mutable direct child atoms: `(key, &mut Atom)`.
    pub fn atoms_mut(&mut self) -> impl Iterator<Item = (&String, &mut Atom)> {
        self.atoms.iter_mut()
    }

    /// Checks if a child group with the given key exists.
    pub fn has_groupkey(&self, key: &str) -> bool {
        self.groups.contains_key(key)
    }

    /// Checks if a child group with the given name exists.
    pub fn has_groupname(&self, name: &str) -> bool {
        self.groups.values().any(|g| g.name == name)
    }

    /// Checks if a child group exists matching key or name.
    pub fn has_group(&self, key_or_name: &str) -> bool {
        self.has_groupkey(key_or_name) || self.has_groupname(key_or_name)
    }

    /// Retrieves a child group by key or name.
    pub fn get_group(&self, key_or_name: &str) -> Option<&AtomGroup> {
        if let Some(g) = self.groups.get(key_or_name) {
            return Some(g);
        }
        self.groups.values().find(|g| g.name == key_or_name)
    }

    /// Retrieves a mutable child group by key or name.
    pub fn get_group_mut(&mut self, key_or_name: &str) -> Option<&mut AtomGroup> {
        if self.groups.contains_key(key_or_name) {
            return self.groups.get_mut(key_or_name);
        }
        self.groups.values_mut().find(|g| g.name == key_or_name)
    }

    /// Sets or adds a child group under `key`.
    pub fn set_group(&mut self, key: &str, mut group: AtomGroup) {
        group.set_path(format!("{}{}/", self.path, key));
        self.groups.insert(key.to_string(), group);
    }

    /// Removes a child group by key.
    pub fn remove_group(&mut self, key: &str) -> Option<AtomGroup> {
        self.groups.shift_remove(key)
    }

    /// Returns the list of child group keys.
    pub fn get_group_list(&self) -> Vec<String> {
        self.groups.keys().cloned().collect()
    }

    /// Checks if a child atom with the given key exists.
    pub fn has_atomkey(&self, key: &str) -> bool {
        self.atoms.contains_key(key)
    }

    /// Checks if a child atom with the given name exists.
    pub fn has_atomname(&self, name: &str) -> bool {
        self.atoms.values().any(|a| a.name.trim() == name.trim())
    }

    /// Checks if a child atom exists matching key or name.
    pub fn has_atom(&self, key_or_name: &str) -> bool {
        self.has_atomkey(key_or_name) || self.has_atomname(key_or_name)
    }

    /// Retrieves a child atom by key or name.
    pub fn get_atom(&self, key_or_name: &str) -> Option<&Atom> {
        if let Some(a) = self.atoms.get(key_or_name) {
            return Some(a);
        }
        self.atoms.values().find(|a| a.name == key_or_name)
    }

    /// Retrieves a mutable child atom by key or name.
    pub fn get_atom_mut(&mut self, key_or_name: &str) -> Option<&mut Atom> {
        if self.atoms.contains_key(key_or_name) {
            return self.atoms.get_mut(key_or_name);
        }
        self.atoms.values_mut().find(|a| a.name == key_or_name)
    }

    /// Retrieves an atom by hierarchical path (e.g. "/group_A/group_B/C3" or "C1").
    pub fn get_atom_by_path(&self, path: &str) -> Option<&Atom> {
        let trimmed = path.trim_start_matches('/');
        let parts: Vec<&str> = trimmed.splitn(2, '/').collect();
        if parts.len() == 1 {
            self.get_atom(parts[0])
        } else {
            let grp_key = parts[0];
            let rest = parts[1];
            self.groups
                .get(grp_key)
                .and_then(|g| g.get_atom_by_path(rest))
        }
    }

    /// Sets an atom directly under `key`.
    fn set_atom_direct(&mut self, key: &str, mut atom: Atom) {
        atom.path = format!("{}{}", self.path, key);
        self.atoms.insert(key.to_string(), atom);
    }

    /// Sets an atom using a key or path (e.g. "C1" or "/group_A/group_B/C3").
    /// Intermediate groups are automatically created if they do not exist.
    pub fn set_atom(&mut self, path_or_key: &str, atom: Atom) {
        let trimmed = path_or_key.trim_start_matches('/');
        let parts: Vec<&str> = trimmed.splitn(2, '/').collect();
        if parts.len() == 1 {
            self.set_atom_direct(parts[0], atom);
        } else {
            let grp_key = parts[0];
            let rest = parts[1];
            if !self.has_groupkey(grp_key) {
                self.set_group(grp_key, AtomGroup::new());
            }
            self.groups.get_mut(grp_key).unwrap().set_atom(rest, atom);
        }
    }

    /// Removes an atom by key.
    pub fn remove_atom(&mut self, key: &str) -> Option<Atom> {
        self.atoms.shift_remove(key)
    }

    /// Collects all atoms within this group and all subgroups into a flat `Vec<Atom>`.
    pub fn get_atom_list(&self) -> Vec<Atom> {
        let mut list = Vec::new();
        for group in self.groups.values() {
            list.extend(group.get_atom_list());
        }
        for atom in self.atoms.values() {
            list.push(atom.clone());
        }
        list
    }

    /// Returns a list of paths of all atoms in this group and subgroups.
    pub fn get_path_list(&self) -> Vec<String> {
        let mut paths = Vec::new();
        for group in self.groups.values() {
            paths.extend(group.get_path_list());
        }
        for atom in self.atoms.values() {
            paths.push(atom.path.clone());
        }
        paths
    }

    /// Returns the distinct atom symbols present in this group and subgroups.
    pub fn get_atom_kinds(&self) -> HashSet<String> {
        let mut kinds = HashSet::new();
        for group in self.groups.values() {
            kinds.extend(group.get_atom_kinds());
        }
        for atom in self.atoms.values() {
            if let Ok(sym) = atom.symbol() {
                kinds.insert(sym.to_string());
            }
        }
        kinds
    }

    /// Returns a map from atom symbol to atom count.
    pub fn get_atom_kinds_count(&self) -> HashMap<String, usize> {
        let mut counts = HashMap::new();
        for group in self.groups.values() {
            for (sym, cnt) in group.get_atom_kinds_count() {
                *counts.entry(sym).or_insert(0) += cnt;
            }
        }
        for atom in self.atoms.values() {
            if let Ok(sym) = atom.symbol() {
                *counts.entry(sym.to_string()).or_insert(0) += 1;
            }
        }
        counts
    }

    /// Returns the molecular formula (Hill system ordering: H, C, then others by atomic number).
    pub fn get_formula(&self) -> String {
        let kinds = self.get_atom_kinds_count();
        let mut formula = String::new();

        let num_elements = PeriodicTable::get_num_of_atoms();
        for atomic_num in 1..num_elements {
            if let Ok(symbol) = PeriodicTable::get_symbol(atomic_num) {
                if let Some(&count) = kinds.get(symbol) {
                    formula.push_str(&format!("{}{}", symbol, count));
                }
            }
        }
        if let Some(&count) = kinds.get("X") {
            formula.push_str(&format!("X{}", count));
        }

        formula
    }

    /// Merges another `AtomGroup` into this group.
    pub fn merge(&mut self, other: &AtomGroup) {
        for (key, group) in other.groups() {
            if let Some(existing) = self.groups.get_mut(key) {
                existing.merge(group);
            } else if let Some(existing) = self
                .groups
                .values_mut()
                .find(|g| !g.name.is_empty() && g.name == group.name)
            {
                existing.merge(group);
            } else {
                self.set_group(key, group.clone());
            }
        }
        for (key, atom) in other.atoms() {
            self.set_atom_direct(key, atom.clone());
        }
        for bond in &other.bonds {
            if !self.bonds.contains(bond) {
                self.bonds.push(bond.clone());
            }
        }
    }

    /// Shifts all atoms in this group and subgroups by `direction`.
    pub fn shift_by(&mut self, direction: Position) {
        for group in self.groups.values_mut() {
            group.shift_by(direction);
        }
        for atom in self.atoms.values_mut() {
            atom.shift_by(direction);
        }
    }

    /// Rotates all atoms in this group and subgroups by a 3x3 rotation matrix.
    pub fn rotate(&mut self, rotmat: &Matrix) -> Result<()> {
        for group in self.groups.values_mut() {
            group.rotate(rotmat)?;
        }
        for atom in self.atoms.values_mut() {
            atom.rotate(rotmat)?;
        }
        Ok(())
    }

    /// Returns the geometric center of all atoms in this group.
    pub fn center(&self) -> Position {
        let mut sum = Position::default();
        let mut total_atoms = 0;
        for atom in self.atoms.values() {
            sum += atom.xyz;
            total_atoms += 1;
        }
        for group in self.groups.values() {
            let n = group.get_number_of_all_atoms();
            if n > 0 {
                sum += group.center() * (n as f64);
                total_atoms += n;
            }
        }
        if total_atoms > 0 {
            sum /= total_atoms as f64;
        }
        sum
    }

    /// Returns the bounding box `(box_min, box_max)` of the atom group.
    pub fn r#box(&self) -> (Position, Position) {
        let mut box_min = self.center();
        let mut box_max = box_min;

        for group in self.groups.values() {
            let (grp_min, grp_max) = group.r#box();
            box_min.x = box_min.x.min(grp_min.x);
            box_min.y = box_min.y.min(grp_min.y);
            box_min.z = box_min.z.min(grp_min.z);
            box_max.x = box_max.x.max(grp_max.x);
            box_max.y = box_max.y.max(grp_max.y);
            box_max.z = box_max.z.max(grp_max.z);
        }

        for atom in self.atoms.values() {
            box_min.x = box_min.x.min(atom.xyz.x);
            box_min.y = box_min.y.min(atom.xyz.y);
            box_min.z = box_min.z.min(atom.xyz.z);
            box_max.x = box_max.x.max(atom.xyz.x);
            box_max.y = box_max.y.max(atom.xyz.y);
            box_max.z = box_max.z.max(atom.xyz.z);
        }

        (box_min, box_max)
    }

    /// Alias for `r#box()` for Rust idiom.
    pub fn get_box(&self) -> (Position, Position) {
        self.r#box()
    }

    /// Selects matching atoms and groups into a new `AtomGroup`.
    pub fn select<S: Selector>(&self, selector: &S) -> Self {
        if selector.is_match_group(self) {
            return self.clone();
        }

        let mut answer = AtomGroup::with_name(&self.name);
        answer.path = self.path.clone();

        for (key, group) in self.groups() {
            let sub = group.select(selector);
            if sub.get_number_of_all_atoms() > 0 {
                answer.set_group(key, sub);
            }
        }
        for (key, atom) in self.atoms() {
            if selector.is_match_atom(atom) {
                answer.set_atom_direct(key, atom.clone());
            }
        }

        answer
    }

    /// Adds a bond between two atoms.
    pub fn add_bond(&mut self, atom1: &Atom, atom2: &Atom, order: usize) {
        let p1 = if atom1.path.is_empty() {
            self.atoms
                .iter()
                .find(|(_, a)| a.name == atom1.name)
                .map(|(k, _)| format!("{}{}", self.path, k))
                .unwrap_or_else(|| atom1.name.clone())
        } else {
            atom1.path.clone()
        };
        let p2 = if atom2.path.is_empty() {
            self.atoms
                .iter()
                .find(|(_, a)| a.name == atom2.name)
                .map(|(k, _)| format!("{}{}", self.path, k))
                .unwrap_or_else(|| atom2.name.clone())
        } else {
            atom2.path.clone()
        };
        self.bonds.push(BondRecord {
            atom1_path: p1,
            atom2_path: p2,
            order,
        });
    }

    /// Returns the number of bonds in this group.
    pub fn get_number_of_bonds(&self) -> usize {
        self.bonds.len()
    }

    /// Returns the list of bonds directly defined in this group.
    pub fn bonds(&self) -> &[BondRecord] {
        &self.bonds
    }

    /// Recursively returns the list of all bonds in this group and its subgroups.
    pub fn get_bond_list(&mut self) -> Vec<BondRecord> {
        self.update_paths();
        let mut bond_list = Vec::new();
        self.collect_bond_list(&mut bond_list);
        bond_list
    }

    fn collect_bond_list(&self, bond_list: &mut Vec<BondRecord>) {
        for group in self.groups.values() {
            group.collect_bond_list(bond_list);
        }
        for b in &self.bonds {
            let path1 = if b.atom1_path.starts_with('/') {
                b.atom1_path.clone()
            } else {
                format!("{}{}", self.path, b.atom1_path)
            };
            let path2 = if b.atom2_path.starts_with('/') {
                b.atom2_path.clone()
            } else {
                format!("{}{}", self.path, b.atom2_path)
            };
            bond_list.push(BondRecord {
                atom1_path: path1,
                atom2_path: path2,
                order: b.order,
            });
        }
    }

    /// Splits a path string into components, removing empty parts.
    pub fn divide_path(path: &str) -> Vec<String> {
        path.split('/')
            .filter(|s| !s.is_empty())
            .map(|s| s.to_string())
            .collect()
    }
}

// Indexing: group["key"] returns child AtomGroup
impl Index<&str> for AtomGroup {
    type Output = AtomGroup;
    fn index(&self, key: &str) -> &Self::Output {
        self.get_group(key)
            .unwrap_or_else(|| panic!("Group key not found: {}", key))
    }
}

// BitAnd: Intersection of two AtomGroups
impl BitAnd for &AtomGroup {
    type Output = AtomGroup;
    fn bitand(self, rhs: Self) -> Self::Output {
        let mut result = AtomGroup::new();
        result.path = self.path.clone();

        for (key, group) in self.groups() {
            if let Some(rhs_group) = rhs.get_group(key) {
                let inter = group & rhs_group;
                if inter.get_number_of_all_atoms() > 0 {
                    result.set_group(key, inter);
                }
            }
        }

        for (key, atom) in self.atoms() {
            if rhs.has_atom(key) {
                result.set_atom_direct(key, atom.clone());
            }
        }

        for bond in &self.bonds {
            if rhs.bonds.contains(bond) && !result.bonds.contains(bond) {
                result.bonds.push(bond.clone());
            }
        }

        result
    }
}

impl BitAnd for AtomGroup {
    type Output = AtomGroup;
    fn bitand(self, rhs: Self) -> Self::Output {
        &self & &rhs
    }
}

impl BitAndAssign<&AtomGroup> for AtomGroup {
    fn bitand_assign(&mut self, rhs: &AtomGroup) {
        *self = &*self & rhs;
    }
}

impl BitAndAssign for AtomGroup {
    fn bitand_assign(&mut self, rhs: Self) {
        *self = &*self & &rhs;
    }
}

// BitOr: Union of two AtomGroups
impl BitOr for &AtomGroup {
    type Output = AtomGroup;
    fn bitor(self, rhs: Self) -> Self::Output {
        let mut result = self.clone();
        result.merge(rhs);
        result
    }
}

impl BitOr for AtomGroup {
    type Output = AtomGroup;
    fn bitor(self, rhs: Self) -> Self::Output {
        &self | &rhs
    }
}

impl BitOrAssign<&AtomGroup> for AtomGroup {
    fn bitor_assign(&mut self, rhs: &AtomGroup) {
        self.merge(rhs);
    }
}

impl BitOrAssign for AtomGroup {
    fn bitor_assign(&mut self, rhs: Self) {
        self.merge(&rhs);
    }
}

// BitXor: Symmetric difference of two AtomGroups
impl BitXor for &AtomGroup {
    type Output = AtomGroup;
    fn bitxor(self, rhs: Self) -> Self::Output {
        let mut result = AtomGroup::new();
        result.path = self.path.clone();

        // Subgroups (preserve order: self first, then rhs)
        let mut all_group_keys: Vec<String> = self.groups.keys().cloned().collect();
        for k in rhs.groups.keys() {
            if !all_group_keys.contains(k) {
                all_group_keys.push(k.clone());
            }
        }

        for key in all_group_keys {
            let in_self = self.get_group(&key);
            let in_rhs = rhs.get_group(&key);

            match (in_self, in_rhs) {
                (Some(g1), Some(g2)) => {
                    let diff = g1 ^ g2;
                    if diff.get_number_of_all_atoms() > 0 {
                        result.set_group(&key, diff);
                    }
                }
                (Some(g1), None) => {
                    result.set_group(&key, g1.clone());
                }
                (None, Some(g2)) => {
                    result.set_group(&key, g2.clone());
                }
                (None, None) => {}
            }
        }

        // Atoms (preserve order: self first, then rhs)
        let mut all_atom_keys: Vec<String> = self.atoms.keys().cloned().collect();
        for k in rhs.atoms.keys() {
            if !all_atom_keys.contains(k) {
                all_atom_keys.push(k.clone());
            }
        }

        for key in all_atom_keys {
            let in_self = self.get_atom(&key);
            let in_rhs = rhs.get_atom(&key);

            match (in_self, in_rhs) {
                (Some(a1), None) => {
                    result.set_atom_direct(&key, a1.clone());
                }
                (None, Some(a2)) => {
                    result.set_atom_direct(&key, a2.clone());
                }
                _ => {} // present in both -> omitted
            }
        }

        // Bonds
        for bond in &self.bonds {
            if !rhs.bonds.contains(bond) && !result.bonds.contains(bond) {
                result.bonds.push(bond.clone());
            }
        }
        for bond in &rhs.bonds {
            if !self.bonds.contains(bond) && !result.bonds.contains(bond) {
                result.bonds.push(bond.clone());
            }
        }

        result
    }
}

impl BitXor for AtomGroup {
    type Output = AtomGroup;
    fn bitxor(self, rhs: Self) -> Self::Output {
        &self ^ &rhs
    }
}

impl BitXorAssign<&AtomGroup> for AtomGroup {
    fn bitxor_assign(&mut self, rhs: &AtomGroup) {
        *self = &*self ^ rhs;
    }
}

impl BitXorAssign for AtomGroup {
    fn bitxor_assign(&mut self, rhs: Self) {
        *self = &*self ^ &rhs;
    }
}

impl fmt::Display for AtomGroup {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "# group name={}", self.name)?;
        for (key, group) in self.groups() {
            write!(f, "  group key={}: {}", key, group)?;
        }
        for (key, atom) in self.atoms() {
            writeln!(f, "  {} {}", key, atom)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Test-only helper selector matching atoms within a Euclidean sphere.
    struct SelectRange {
        center: Position,
        radius: f64,
    }

    impl SelectRange {
        fn new(center: Position, radius: f64) -> Self {
            Self { center, radius }
        }
    }

    impl Selector for SelectRange {
        fn is_match_atom(&self, atom: &Atom) -> bool {
            atom.xyz.distance_from(&self.center) <= self.radius
        }
    }

    // Ported from tests/test_atomgroup.py
    #[test]
    fn test_init() {
        let mut group1 = AtomGroup::new();
        let atom1 = Atom::from_symbol("C").unwrap();
        let atom2 = Atom::from_symbol("H").unwrap();
        let atom3 = Atom::from_symbol("N").unwrap();

        let mut subgrp = AtomGroup::new();
        subgrp.set_atom("C1", atom1);
        subgrp.set_atom("H1", atom2);
        subgrp.set_atom("N1", atom3);
        group1.set_group("grp", subgrp);

        assert_eq!(group1["grp"].get_atom("C1").unwrap().symbol().unwrap(), "C");
        assert_eq!(group1["grp"].get_atom("C1").unwrap().path, "/grp/C1");
        assert_eq!(group1["grp"].get_atom("H1").unwrap().symbol().unwrap(), "H");
        assert_eq!(group1["grp"].get_atom("H1").unwrap().path, "/grp/H1");

        assert_eq!(group1.get_number_of_atoms(), 0);
        assert_eq!(group1.get_number_of_groups(), 1);
        assert_eq!(group1.get_number_of_all_atoms(), 3);
        assert!((group1.sum_of_atomic_number() - 14.0).abs() < 1e-10);
    }

    #[test]
    fn test_op_and() {
        let c1 = Atom::from_symbol("C").unwrap();
        let h1 = Atom::from_symbol("H").unwrap();
        let n1 = Atom::from_symbol("N").unwrap();

        let mut group1 = AtomGroup::new();
        let mut subgrp1 = AtomGroup::new();
        subgrp1.set_atom("C1", c1.clone());
        subgrp1.set_atom("H1", h1.clone());
        subgrp1.set_atom("N1", n1);
        group1.set_group("grp", subgrp1);

        let mut group2 = AtomGroup::new();
        let mut subgrp2 = AtomGroup::new();
        subgrp2.set_atom("C1", c1);
        subgrp2.set_atom("H1", h1);
        group2.set_group("grp", subgrp2);

        let group3 = &group1 & &group2;

        assert_eq!(group3.get_number_of_all_atoms(), 2);
        assert_eq!(group3.get_number_of_atoms(), 0);
        assert!(group3.has_groupkey("grp"));
        assert_eq!(group3["grp"].get_number_of_all_atoms(), 2);
        assert_eq!(group3["grp"].get_number_of_atoms(), 2);
        assert!(group3["grp"].has_atom("C1"));
        assert!(group3["grp"].has_atom("H1"));
    }

    #[test]
    fn test_op_iand() {
        let c1 = Atom::from_symbol("C").unwrap();
        let h1 = Atom::from_symbol("H").unwrap();
        let n1 = Atom::from_symbol("N").unwrap();

        let mut group1 = AtomGroup::new();
        let mut subgrp1 = AtomGroup::new();
        subgrp1.set_atom("C1", c1.clone());
        subgrp1.set_atom("H1", h1.clone());
        subgrp1.set_atom("N1", n1);
        group1.set_group("grp", subgrp1);

        let mut group2 = AtomGroup::new();
        let mut subgrp2 = AtomGroup::new();
        subgrp2.set_atom("C1", c1);
        subgrp2.set_atom("H1", h1);
        group2.set_group("grp", subgrp2);

        group1 &= &group2;

        assert_eq!(group1.get_number_of_all_atoms(), 2);
        assert_eq!(group1.get_number_of_atoms(), 0);
        assert!(group1.has_groupkey("grp"));
        assert_eq!(group1["grp"].get_number_of_all_atoms(), 2);
        assert_eq!(group1["grp"].get_number_of_atoms(), 2);
        assert!(group1["grp"].has_atom("C1"));
        assert!(group1["grp"].has_atom("H1"));
    }

    #[test]
    fn test_op_or() {
        let c1 = Atom::from_symbol("C").unwrap();
        let h1 = Atom::from_symbol("H").unwrap();
        let n1 = Atom::from_symbol("N").unwrap();

        let mut group1 = AtomGroup::new();
        let mut subgrp1 = AtomGroup::new();
        subgrp1.set_atom("C1", c1.clone());
        subgrp1.set_atom("H1", h1.clone());
        subgrp1.set_atom("N1", n1);
        group1.set_group("grp", subgrp1);

        let mut group2 = AtomGroup::new();
        let mut subgrp2 = AtomGroup::new();
        subgrp2.set_atom("C1", c1);
        subgrp2.set_atom("H1", h1);
        group2.set_group("grp", subgrp2);

        let group3 = &group1 | &group2;

        assert_eq!(group3.get_number_of_all_atoms(), 3);
        assert_eq!(group3.get_number_of_atoms(), 0);
        assert!(group3.has_groupkey("grp"));
        assert_eq!(group3["grp"].get_number_of_all_atoms(), 3);
        assert_eq!(group3["grp"].get_number_of_atoms(), 3);
        assert!(group3["grp"].has_atom("C1"));
        assert!(group3["grp"].has_atom("H1"));
        assert!(group3["grp"].has_atom("N1"));
    }

    #[test]
    fn test_op_ior() {
        let c1 = Atom::from_symbol("C").unwrap();
        let h1 = Atom::from_symbol("H").unwrap();
        let n1 = Atom::from_symbol("N").unwrap();

        let mut group1 = AtomGroup::new();
        let mut subgrp1 = AtomGroup::new();
        subgrp1.set_atom("C1", c1.clone());
        subgrp1.set_atom("H1", h1.clone());
        subgrp1.set_atom("N1", n1);
        group1.set_group("grp", subgrp1);

        let mut group2 = AtomGroup::new();
        let mut subgrp2 = AtomGroup::new();
        subgrp2.set_atom("C1", c1);
        subgrp2.set_atom("H1", h1);
        group2.set_group("grp", subgrp2);

        group1 |= &group2;

        assert_eq!(group1.get_number_of_all_atoms(), 3);
        assert_eq!(group1.get_number_of_atoms(), 0);
        assert!(group1.has_groupkey("grp"));
        assert_eq!(group1["grp"].get_number_of_all_atoms(), 3);
        assert_eq!(group1["grp"].get_number_of_atoms(), 3);
        assert!(group1["grp"].has_atom("C1"));
        assert!(group1["grp"].has_atom("H1"));
        assert!(group1["grp"].has_atom("N1"));
    }

    #[test]
    fn test_op_xor() {
        let c1 = Atom::from_symbol("C").unwrap();
        let h1 = Atom::from_symbol("H").unwrap();
        let n1 = Atom::from_symbol("N").unwrap();
        let o1 = Atom::from_symbol("O").unwrap();

        let mut group1 = AtomGroup::new();
        let mut subgrp1 = AtomGroup::new();
        subgrp1.set_atom("C1", c1.clone());
        subgrp1.set_atom("H1", h1.clone());
        subgrp1.set_atom("N1", n1);
        group1.set_group("grp", subgrp1);

        let mut group2 = AtomGroup::new();
        let mut subgrp2 = AtomGroup::new();
        subgrp2.set_atom("C1", c1);
        subgrp2.set_atom("H1", h1);
        subgrp2.set_atom("O1", o1);
        group2.set_group("grp", subgrp2);

        let group3 = &group1 ^ &group2;

        assert_eq!(group3.get_number_of_all_atoms(), 2);
        assert_eq!(group3.get_number_of_atoms(), 0);
        assert!(group3.has_groupkey("grp"));
        assert_eq!(group3["grp"].get_number_of_all_atoms(), 2);
        assert_eq!(group3["grp"].get_number_of_atoms(), 2);
        assert!(group3["grp"].has_atom("N1"));
        assert!(group3["grp"].has_atom("O1"));
    }

    #[test]
    fn test_op_ixor() {
        let c1 = Atom::from_symbol("C").unwrap();
        let h1 = Atom::from_symbol("H").unwrap();
        let n1 = Atom::from_symbol("N").unwrap();
        let o1 = Atom::from_symbol("O").unwrap();

        let mut group1 = AtomGroup::new();
        let mut subgrp1 = AtomGroup::new();
        subgrp1.set_atom("C1", c1.clone());
        subgrp1.set_atom("H1", h1.clone());
        subgrp1.set_atom("N1", n1);
        group1.set_group("grp", subgrp1);

        let mut group2 = AtomGroup::new();
        let mut subgrp2 = AtomGroup::new();
        subgrp2.set_atom("C1", c1);
        subgrp2.set_atom("H1", h1);
        subgrp2.set_atom("O1", o1);
        group2.set_group("grp", subgrp2);

        group1 ^= &group2;

        assert_eq!(group1.get_number_of_all_atoms(), 2);
        assert_eq!(group1.get_number_of_atoms(), 0);
        assert!(group1.has_groupkey("grp"));
        assert_eq!(group1["grp"].get_number_of_all_atoms(), 2);
        assert_eq!(group1["grp"].get_number_of_atoms(), 2);
        assert!(group1["grp"].has_atom("N1"));
        assert!(group1["grp"].has_atom("O1"));
    }

    #[test]
    fn test_atom_list() {
        let mut group1 = AtomGroup::new();
        let mut subgrp = AtomGroup::new();
        subgrp.set_atom("C1", Atom::from_symbol("C").unwrap());
        subgrp.set_atom("H1", Atom::from_symbol("H").unwrap());
        subgrp.set_atom("N1", Atom::from_symbol("N").unwrap());
        group1.set_group("grp", subgrp);

        let atom_list = group1.get_atom_list();
        assert_eq!(atom_list.len(), 3);
    }

    #[test]
    fn test_set_atom_by_path() {
        let atom1 = Atom::new_with_pos("C", Position::new(1.1, 2.1, 3.1)).unwrap();
        let atom2 = Atom::new_with_pos("C", Position::new(1.2, 2.2, 3.2)).unwrap();
        let atom3 = Atom::new_with_pos("C", Position::new(1.3, 2.3, 3.3)).unwrap();

        let mut atomgroup = AtomGroup::new();
        atomgroup.set_atom("/C1", atom1);
        atomgroup.set_atom("/group_A/C2", atom2);
        atomgroup.set_atom("/group_A/group_B/C3", atom3);

        assert_eq!(atomgroup.get_number_of_all_atoms(), 3);
        assert_eq!(atomgroup.get_number_of_atoms(), 1);
        assert!(atomgroup.has_groupkey("group_A"));
        assert_eq!(atomgroup["group_A"].get_number_of_all_atoms(), 2);
        assert_eq!(atomgroup["group_A"].get_number_of_atoms(), 1);
        assert!(atomgroup["group_A"].has_groupkey("group_B"));
        assert_eq!(atomgroup["group_A"]["group_B"].get_number_of_atoms(), 1);
    }

    #[test]
    fn test_path() {
        let mut atomgroup1 = AtomGroup::new();
        atomgroup1.set_atom("C", Atom::from_symbol("C").unwrap());
        atomgroup1.set_atom("H1", Atom::from_symbol("H").unwrap());
        atomgroup1.set_atom("H2", Atom::from_symbol("H").unwrap());
        atomgroup1.set_atom("H3", Atom::from_symbol("H").unwrap());

        assert_eq!(atomgroup1.path(), "/");
        assert_eq!(atomgroup1.get_atom("C").unwrap().path, "/C");
        assert_eq!(atomgroup1.get_atom("H1").unwrap().path, "/H1");
        assert_eq!(atomgroup1.get_atom("H2").unwrap().path, "/H2");
        assert_eq!(atomgroup1.get_atom("H3").unwrap().path, "/H3");

        let mut atomgroup2 = AtomGroup::new();
        atomgroup2.set_group("Me", atomgroup1);
        assert_eq!(atomgroup2.path(), "/");
        assert_eq!(atomgroup2["Me"].get_atom("C").unwrap().path, "/Me/C");

        let mut atomgroup3 = AtomGroup::new();
        atomgroup3.set_group("grp3", atomgroup2);
        assert_eq!(
            atomgroup3["grp3"]["Me"].get_atom("H3").unwrap().path,
            "/grp3/Me/H3"
        );
    }

    #[test]
    fn test_path_copy() {
        let atom10 = Atom::from_symbol("C").unwrap();
        let mut atomgroup1 = AtomGroup::new();
        atomgroup1.set_atom("C", atom10);

        assert_eq!(atomgroup1.path(), "/");
        assert_eq!(atomgroup1.get_atom("C").unwrap().path, "/C");

        let grp_cp = atomgroup1.clone();
        assert_eq!(grp_cp.path(), "/");
        assert_eq!(grp_cp.get_atom("C").unwrap().path, "/C");

        let mut atomgroup2 = AtomGroup::new();
        atomgroup2.set_group("Me", atomgroup1);
        let grp_cp2 = atomgroup2.clone();
        assert_eq!(grp_cp2.path(), "/");
        assert_eq!(grp_cp2["Me"].get_atom("C").unwrap().path, "/Me/C");
    }

    #[test]
    fn test_select_range() {
        let mut grp1 = AtomGroup::new();
        grp1.set_atom(
            "H0",
            Atom::new_with_pos("H", Position::new(1.0, 0.0, 0.0)).unwrap(),
        );
        grp1.set_atom(
            "H1",
            Atom::new_with_pos("H", Position::new(1.1, 0.0, 0.0)).unwrap(),
        );
        grp1.set_atom(
            "H2",
            Atom::new_with_pos("H", Position::new(1.2, 0.0, 0.0)).unwrap(),
        );
        grp1.set_atom(
            "H3",
            Atom::new_with_pos("H", Position::new(1.3, 0.0, 0.0)).unwrap(),
        );

        let mut grp10 = AtomGroup::new();
        grp10.set_group("g1", grp1);
        let mut grp100 = AtomGroup::new();
        grp100.set_group("g10", grp10);
        assert_eq!(
            grp100["g10"]["g1"].get_atom("H0").unwrap().path,
            "/g10/g1/H0"
        );

        let selector1 = SelectRange::new(Position::new(0.0, 0.0, 0.0), 1.01);
        let part1 = grp100.select(&selector1);
        assert_eq!(part1.get_number_of_all_atoms(), 1);
        assert_eq!(
            part1["g10"]["g1"].get_atom("H0").unwrap().path,
            "/g10/g1/H0"
        );

        let selector2 = SelectRange::new(Position::new(0.0, 0.0, 0.0), 2.00);
        let part2 = grp100.select(&selector2);
        assert_eq!(part2.get_number_of_all_atoms(), 4);
    }

    #[test]
    fn test_get_path_list() {
        let mut group1 = AtomGroup::new();
        let mut subgrp = AtomGroup::new();
        subgrp.set_atom("C1", Atom::from_symbol("C").unwrap());
        subgrp.set_atom("H1", Atom::from_symbol("H").unwrap());
        subgrp.set_atom("N1", Atom::from_symbol("N").unwrap());
        group1.set_group("grp", subgrp);

        let path_list = group1.get_path_list();
        assert_eq!(path_list.len(), 3);
        assert_eq!(path_list[0], "/grp/C1");
        assert_eq!(path_list[1], "/grp/H1");
        assert_eq!(path_list[2], "/grp/N1");
    }

    #[test]
    fn test_divide_path() {
        let parts = AtomGroup::divide_path("atom0");
        assert_eq!(parts.len(), 1);
        assert_eq!(parts[0], "atom0");

        let parts = AtomGroup::divide_path("/res1/atom2");
        assert_eq!(parts[0], "res1");
        assert_eq!(parts[1], "atom2");
    }

    #[test]
    fn test_get_formula() {
        let mut group1 = AtomGroup::new();
        let mut subgrp = AtomGroup::new();
        subgrp.set_atom("C1", Atom::from_symbol("C").unwrap());
        subgrp.set_atom("H1", Atom::from_symbol("H").unwrap());
        subgrp.set_atom("H2", Atom::from_symbol("H").unwrap());
        group1.set_group("grp", subgrp);

        let formula = group1.get_formula();
        assert_eq!(formula, "H2C1");
    }

    #[test]
    fn test_ixor_operator() {
        let mut ag1 = AtomGroup::new();
        ag1.set_atom("C1", Atom::from_symbol("C").unwrap());
        let mut ag2 = AtomGroup::new();
        ag2.set_atom("N1", Atom::from_symbol("N").unwrap());

        ag1 ^= &ag2;
        assert_eq!(ag1.get_number_of_all_atoms(), 2);
        assert!(ag1.has_atom("C1"));
        assert!(ag1.has_atom("N1"));
    }

    #[test]
    fn test_box() {
        let mut ag = AtomGroup::new();
        let a1 = Atom::new_with_pos("C", Position::new(-1.0, 2.0, 0.0)).unwrap();
        let a2 = Atom::new_with_pos("C", Position::new(3.0, -4.0, 5.0)).unwrap();
        ag.set_atom("1", a1);
        ag.set_atom("2", a2);

        let (bmin, bmax) = ag.r#box();
        assert_eq!(bmin.x, -1.0);
        assert_eq!(bmin.y, -4.0);
        assert_eq!(bmin.z, 0.0);
        assert_eq!(bmax.x, 3.0);
        assert_eq!(bmax.y, 2.0);
        assert_eq!(bmax.z, 5.0);
    }

    #[test]
    fn test_bonds_and_bond_list() {
        let mut ag = AtomGroup::with_name("mol");
        let mut sub = AtomGroup::with_name("sub");
        let mut a1 = Atom::from_symbol("C").unwrap();
        a1.name = "C1".to_string();
        let mut a2 = Atom::from_symbol("C").unwrap();
        a2.name = "C2".to_string();
        sub.set_atom("1", a1.clone());
        sub.set_atom("2", a2.clone());
        sub.add_bond(&a1, &a2, 2);
        assert_eq!(sub.get_number_of_bonds(), 1);

        ag.set_group("grp", sub);
        let bond_list = ag.get_bond_list();
        assert_eq!(bond_list.len(), 1);
        assert_eq!(bond_list[0].order, 2);
    }
}
