// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

//! Hierarchical molecular structure representation for `proteindf-bridge`.
//!
//! # Protein Path Schema Convention
//!
//! Standard protein structures in `proteindf-bridge` follow a 4-level hierarchical path convention:
//!
//! ```text
//! /model_N/chain_id/res_key/atom_key
//! ```
//!
//! - **Level 0 (Root / Models)**: Path depth 0 (e.g. `"/"`). Contains model subgroups. No direct atoms.
//! - **Level 1 (Model)**: Path depth 1 (e.g. `"/model_1/"`). Checked by [`AtomGroup::is_model_level`].
//!   Contains chain subgroups. No direct atoms.
//! - **Level 2 (Chain)**: Path depth 2 (e.g. `"/model_1/A/"`). Checked by [`AtomGroup::is_chain_level`].
//!   Contains residue subgroups. No direct atoms.
//! - **Level 3 (Residue)**: Path depth 3 (e.g. `"/model_1/A/6/"`). Checked by [`AtomGroup::is_residue_level`].
//!   Contains atoms. Must not contain subgroups.
//! - **Level 4 (Atom)**: Path (e.g. `"/model_1/A/6/CA"`). Leaf entity in the tree.
//!
//! ## Positional vs. Structural Checks
//!
//! Note the distinction between positional level checks and structural checks:
//! - **Positional checks** ([`AtomGroup::is_model_level`], [`AtomGroup::is_chain_level`], [`AtomGroup::is_residue_level`])
//!   inspect only the **path depth** of the group within the hierarchy.
//! - **Structural checks** ([`crate::Format::is_protein`], [`crate::Format::is_chain`], [`crate::Format::is_residue`])
//!   inspect the **content** (e.g. ensuring no direct atoms, all children are residues, etc.).
//!
//! In valid structures, positional and structural checks coincide. In malformed data (such as HETATM
//! or water molecules placed directly in a chain without a residue wrapper), a group's path depth remains
//! at the chain level (depth 2) while its structural validity ([`crate::Format::is_chain`]) fails.
//! Use [`AtomGroup::validate_schema`] to detect such violations across the entire hierarchy.

use std::collections::{HashMap, HashSet};
use std::fmt;
use std::ops::{BitAnd, BitAndAssign, BitOr, BitOrAssign, BitXor, BitXorAssign, Index};

use indexmap::IndexMap;

use crate::atom::Atom;
use crate::error::Result;
use crate::matrix::Matrix;
use crate::periodic_table::PeriodicTable;
use crate::position::Position;
use crate::secondary_structure::SsCode;

/// A violation of the standard protein schema (`/model_N/chain_id/res_key/atom_key`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SchemaViolation {
    /// Direct atoms found at a non-residue level (e.g., under root, model, or chain).
    DirectAtomsAtNonResidueLevel {
        path: String,
        depth: usize,
        atom_keys: Vec<String>,
    },
    /// Subgroups found inside a residue-level group (residues must be leaf groups).
    SubgroupsInResidue {
        path: String,
        group_keys: Vec<String>,
    },
    /// Group nested deeper than the residue level (depth > 3).
    ExcessiveDepth { path: String, depth: usize },
}

impl fmt::Display for SchemaViolation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DirectAtomsAtNonResidueLevel {
                path,
                depth,
                atom_keys,
            } => {
                write!(
                    f,
                    "Group '{}' at depth {} has {} direct atom(s) ({:?}), but atoms are only allowed at residue level (depth 3)",
                    path,
                    depth,
                    atom_keys.len(),
                    atom_keys
                )
            }
            Self::SubgroupsInResidue { path, group_keys } => {
                write!(
                    f,
                    "Residue-level group '{}' contains {} subgroup(s) ({:?}), but residues must not contain subgroups",
                    path,
                    group_keys.len(),
                    group_keys
                )
            }
            Self::ExcessiveDepth { path, depth } => {
                write!(
                    f,
                    "Group '{}' has path depth {}, exceeding maximum standard depth 3",
                    path, depth
                )
            }
        }
    }
}

impl std::error::Error for SchemaViolation {}

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
    depth: usize,
    atoms: IndexMap<String, Atom>,
    groups: IndexMap<String, AtomGroup>,
    bonds: Vec<BondRecord>,
    secondary_structure: Option<SsCode>,
}

impl Default for AtomGroup {
    fn default() -> Self {
        Self {
            name: String::new(),
            path: "/".to_string(),
            depth: 0,
            atoms: IndexMap::new(),
            groups: IndexMap::new(),
            bonds: Vec::new(),
            secondary_structure: None,
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

    /// Returns the depth of this group in the hierarchy tree.
    ///
    /// Root (`"/"`) has depth 0, model (`"/model_1/"`) has depth 1,
    /// chain (`"/model_1/A/"`) has depth 2, and residue (`"/model_1/A/6/"`) has depth 3.
    ///
    /// This reflects the true tree nesting depth rather than counting slashes
    /// in the path string, ensuring robustness against keys containing slashes
    /// or empty segments.
    pub fn path_depth(&self) -> usize {
        self.depth
    }

    /// Checks whether this group is at the model level based on path depth (depth == 1).
    ///
    /// # Positional vs. Structural Check
    /// This is a positional check based on the group's `path()` (e.g. `"/model_1/"`).
    /// In contrast, [`crate::Format::is_protein`] is a structural check that verifies
    /// there are no direct atoms and all subgroups satisfy chain conditions.
    pub fn is_model_level(&self) -> bool {
        self.path_depth() == 1
    }

    /// Checks whether this group is at the chain level based on path depth (depth == 2).
    ///
    /// # Positional vs. Structural Check
    /// This is a positional check based on the group's `path()` (e.g. `"/model_1/A/"`).
    /// In contrast, [`crate::Format::is_chain`] is a structural check that verifies
    /// there are no direct atoms and all subgroups satisfy residue conditions.
    pub fn is_chain_level(&self) -> bool {
        self.path_depth() == 2
    }

    /// Checks whether this group is at the residue level based on path depth (depth == 3).
    ///
    /// # Positional vs. Structural Check
    /// This is a positional check based on the group's `path()` (e.g. `"/model_1/A/6/"`).
    /// In contrast, [`crate::Format::is_residue`] is a structural check that verifies
    /// there are no child groups (subgroups).
    pub fn is_residue_level(&self) -> bool {
        self.path_depth() == 3
    }

    /// Validates this `AtomGroup` tree against the standard protein schema
    /// (`/model_N/chain_id/res_key/atom_key`).
    ///
    /// Traverses the entire tree recursively and collects all detected
    /// [`SchemaViolation`]s, including:
    /// - Direct atoms found at non-residue levels (depth != 3, e.g. HETATM placed directly under a chain).
    /// - Subgroups found inside a residue-level group (depth == 3).
    /// - Groups nested deeper than standard residue level (depth > 3).
    pub fn validate_schema(&self) -> Vec<SchemaViolation> {
        let mut violations = Vec::new();
        self.collect_schema_violations(self.depth, &mut violations);
        violations
    }

    fn collect_schema_violations(
        &self,
        current_depth: usize,
        violations: &mut Vec<SchemaViolation>,
    ) {
        // Atoms are only allowed at residue level (depth 3).
        if !self.atoms.is_empty() && current_depth != 3 {
            violations.push(SchemaViolation::DirectAtomsAtNonResidueLevel {
                path: self.path.clone(),
                depth: current_depth,
                atom_keys: self.atoms.keys().cloned().collect(),
            });
        }

        // Residues (depth 3) must not have subgroups.
        if current_depth == 3 && !self.groups.is_empty() {
            violations.push(SchemaViolation::SubgroupsInResidue {
                path: self.path.clone(),
                group_keys: self.groups.keys().cloned().collect(),
            });
        }

        // Nesting depth must not exceed 3.
        if current_depth > 3 {
            violations.push(SchemaViolation::ExcessiveDepth {
                path: self.path.clone(),
                depth: current_depth,
            });
        }

        for group in self.groups.values() {
            group.collect_schema_violations(current_depth + 1, violations);
        }
    }

    /// Sets the path of this group and updates descendant paths.
    ///
    /// Note: This method only updates `self.path` and propagates paths down the tree;
    /// it preserves `self.depth`. If repositioning or detaching a subtree where
    /// the root depth changes, use [`set_path_with_depth`](Self::set_path_with_depth)
    /// or re-attach the group using [`set_group`](Self::set_group).
    pub fn set_path(&mut self, mut new_path: String) {
        if new_path.is_empty() || !new_path.ends_with('/') {
            new_path.push('/');
        }
        self.path = new_path;
        self.update_paths();
    }

    /// Sets the path and depth of this group, updating descendant paths and depths.
    ///
    /// Use this when detaching a subtree or creating a standalone group at an explicit depth.
    pub fn set_path_with_depth(&mut self, mut new_path: String, depth: usize) {
        if new_path.is_empty() || !new_path.ends_with('/') {
            new_path.push('/');
        }
        self.path = new_path;
        self.depth = depth;
        self.update_paths();
    }

    fn update_paths(&mut self) {
        for (key, group) in self.groups.iter_mut() {
            group.depth = self.depth + 1;
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
        group.depth = self.depth + 1;
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

    /// Retrieves an atom by hierarchical path (e.g. "/model_1/A/1/CA" or "C1").
    ///
    /// The lookup traverses `IndexMap` groups level by level, taking O(depth) time
    /// (typically 4 levels in standard protein schema) independent of the total
    /// atom count in the structure.
    pub fn get_atom_by_path(&self, path: &str) -> Option<&Atom> {
        let trimmed = path.trim_start_matches('/');
        if let Some((grp_key, rest)) = trimmed.split_once('/') {
            self.groups
                .get(grp_key)
                .and_then(|g| g.get_atom_by_path(rest))
        } else {
            self.get_atom(trimmed)
        }
    }

    /// Resolves both endpoint atoms of a [`BondRecord`] by their hierarchical paths.
    ///
    /// Performs an O(depth) lookup for each endpoint. Returns `None` if either
    /// atom cannot be found.
    pub fn resolve_bond<'a>(&'a self, record: &BondRecord) -> Option<(&'a Atom, &'a Atom)> {
        let a1 = self.get_atom_by_path(&record.atom1_path)?;
        let a2 = self.get_atom_by_path(&record.atom2_path)?;
        Some((a1, a2))
    }

    /// Collects all atoms within this group and subgroups whose key or name matches.
    pub fn pickup_atoms(&self, key_or_name: &str) -> Vec<Atom> {
        let mut result = Vec::new();
        for subgrp in self.groups.values() {
            result.extend(subgrp.pickup_atoms(key_or_name));
        }
        for (atm_key, atm) in &self.atoms {
            if atm_key == key_or_name || atm.name == key_or_name {
                result.push(atm.clone());
            }
        }
        result
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
        if other.secondary_structure.is_some() {
            self.secondary_structure = other.secondary_structure;
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

    /// Computes the longest common directory path ending with '/' between two paths.
    pub fn get_common_path(path1: &str, path2: &str) -> String {
        let common_prefix: String = path1
            .chars()
            .zip(path2.chars())
            .take_while(|(c1, c2)| c1 == c2)
            .map(|(c, _)| c)
            .collect();

        let mut common = if common_prefix.ends_with('/') {
            common_prefix
        } else if let Some(last_slash) = common_prefix.rfind('/') {
            common_prefix[..=last_slash].to_string()
        } else {
            String::from("/")
        };

        if !common.starts_with('/') {
            common.insert(0, '/');
        }
        if !common.ends_with('/') {
            common.push('/');
        }
        common
    }

    /// Returns a reference to the descendant group matching `query_path`.
    ///
    /// Note: Unlike Python's `get_family` which traverses upwards via `self.parent`,
    /// this Rust implementation only searches downwards within `self`'s subtree
    /// because Rust's tree is ownership-based without parent back-references.
    pub fn get_family(&self, query_path: &str) -> Option<&AtomGroup> {
        let q = if !query_path.ends_with('/') {
            format!("{}/", query_path)
        } else {
            query_path.to_string()
        };
        if self.path == q {
            return Some(self);
        }
        for group in self.groups.values() {
            if let Some(found) = group.get_family(&q) {
                return Some(found);
            }
        }
        None
    }

    /// Returns a mutable reference to the descendant group matching `query_path`.
    ///
    /// Note: Downward search only, matching `get_family`.
    pub fn get_family_mut(&mut self, query_path: &str) -> Option<&mut AtomGroup> {
        let q = if !query_path.ends_with('/') {
            format!("{}/", query_path)
        } else {
            query_path.to_string()
        };
        if self.path == q {
            return Some(self);
        }
        for group in self.groups.values_mut() {
            if let Some(found) = group.get_family_mut(&q) {
                return Some(found);
            }
        }
        None
    }

    /// Adds a bond between two atoms, delegating storage to the nearest common ancestor group.
    /// Bonds are stored with relative paths with respect to the common ancestor group.
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

        let common_path = Self::get_common_path(&p1, &p2);

        if self.path == common_path {
            self.add_bond_direct(&p1, &p2, order);
            return;
        }

        let mut found = false;
        for group in self.groups.values_mut() {
            if let Some(family) = group.get_family_mut(&common_path) {
                family.add_bond_direct(&p1, &p2, order);
                found = true;
                break;
            }
        }

        if !found {
            self.add_bond_direct(&p1, &p2, order);
        }
    }

    /// Directly records a bond with relative paths stripped of `self.path` prefix.
    fn add_bond_direct(&mut self, p1: &str, p2: &str, order: usize) {
        let common1 = Self::get_common_path(&self.path, p1);
        let rel_p1 = if p1.starts_with(&common1) {
            &p1[common1.len()..]
        } else if p1.starts_with(&self.path) {
            &p1[self.path.len()..]
        } else {
            p1.trim_start_matches('/')
        };

        let common2 = Self::get_common_path(&self.path, p2);
        let rel_p2 = if p2.starts_with(&common2) {
            &p2[common2.len()..]
        } else if p2.starts_with(&self.path) {
            &p2[self.path.len()..]
        } else {
            p2.trim_start_matches('/')
        };

        self.bonds.push(BondRecord {
            atom1_path: rel_p1.to_string(),
            atom2_path: rel_p2.to_string(),
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

    /// Sets the list of bonds directly defined in this group.
    pub fn set_bonds(&mut self, bonds: Vec<BondRecord>) {
        self.bonds = bonds;
    }

    /// Returns the secondary structure code assigned to this group (typically at residue level).
    pub fn secondary_structure(&self) -> Option<SsCode> {
        self.secondary_structure
    }

    /// Sets the secondary structure code for this group.
    pub fn set_secondary_structure(&mut self, ss: Option<SsCode>) {
        self.secondary_structure = ss;
    }

    /// Applies 3-state secondary structure assignments to each residue in this chain group.
    pub fn apply_secondary_structure(&mut self) {
        crate::secondary_structure::apply_secondary_structure(self);
    }

    /// Returns the raw MessagePack Value representation of this AtomGroup.
    pub fn get_raw_data(&self) -> rmpv::Value {
        crate::brd::atomgroup_get_raw_data(self)
    }

    /// Populates this AtomGroup from a raw MessagePack dictionary Value.
    pub fn set_by_dict_data(&mut self, data: &rmpv::Value) -> Result<&mut Self> {
        crate::brd::atomgroup_set_by_dict_data(self, data)?;
        Ok(self)
    }

    /// Creates an AtomGroup from a raw MessagePack dictionary Value.
    pub fn from_dict_data(data: &rmpv::Value) -> Result<Self> {
        let mut group = AtomGroup::new();
        group.set_by_dict_data(data)?;
        Ok(group)
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

        result.secondary_structure = match (self.secondary_structure, rhs.secondary_structure) {
            (Some(s), Some(r)) if s == r => Some(s),
            _ => None,
        };

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

        result.secondary_structure = match (self.secondary_structure, rhs.secondary_structure) {
            (Some(s), None) => Some(s),
            (None, Some(r)) => Some(r),
            _ => None,
        };

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

    use crate::selector::SelectRange;

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

    #[test]
    fn test_get_common_path() {
        assert_eq!(AtomGroup::get_common_path("/A/1/1_SG", "/A/6/6_SG"), "/A/");
        assert_eq!(AtomGroup::get_common_path("/A/1/1_SG", "/B/2/2_SG"), "/");
        assert_eq!(
            AtomGroup::get_common_path("/model_1/A/1/1_SG", "/model_1/A/6/6_SG"),
            "/model_1/A/"
        );
        assert_eq!(
            AtomGroup::get_common_path("/model_1/A/1/1_SG", "/model_1/B/2/2_SG"),
            "/model_1/"
        );
    }

    #[test]
    fn test_hierarchical_bond_routing_and_reparent() {
        // Setup model with chain A and two residues (1 and 6)
        let mut model = AtomGroup::new();
        let mut chain_a = AtomGroup::new();
        let mut res1 = AtomGroup::new();
        let mut res6 = AtomGroup::new();

        let mut sg1 = Atom::from_symbol("S").unwrap();
        sg1.name = "SG".to_string();
        res1.set_atom("1_SG", sg1);

        let mut sg2 = Atom::from_symbol("S").unwrap();
        sg2.name = "SG".to_string();
        res6.set_atom("6_SG", sg2);

        chain_a.set_group("1", res1);
        chain_a.set_group("6", res6);
        model.set_group("A", chain_a);

        // Retrieve atoms through model to get current paths (/A/1/1_SG and /A/6/6_SG)
        let sg1_ref = model.get_atom_by_path("/A/1/1_SG").unwrap().clone();
        let sg2_ref = model.get_atom_by_path("/A/6/6_SG").unwrap().clone();

        // Add bond at model level (should route to chain A as nearest common ancestor)
        model.add_bond(&sg1_ref, &sg2_ref, 1);

        // Check that bond was routed into chain A
        let chain_a_grp = model.get_group("A").unwrap();
        assert_eq!(chain_a_grp.get_number_of_bonds(), 1);
        assert_eq!(model.get_number_of_bonds(), 0);

        // Now reparent model into root (root.set_group("model_1", model))
        let mut root = AtomGroup::new();
        root.set_group("model_1", model);

        // Collect bond list from root
        let bonds = root.get_bond_list();
        assert_eq!(bonds.len(), 1);
        assert_eq!(bonds[0].atom1_path, "/model_1/A/1/1_SG");
        assert_eq!(bonds[0].atom2_path, "/model_1/A/6/6_SG");
        assert_eq!(bonds[0].order, 1);

        // Verify that root can resolve both bonded atoms by their paths
        let resolved_sg1 = root.get_atom_by_path(&bonds[0].atom1_path);
        let resolved_sg2 = root.get_atom_by_path(&bonds[0].atom2_path);
        assert!(resolved_sg1.is_some());
        assert!(resolved_sg2.is_some());
        assert_eq!(resolved_sg1.unwrap().name, "SG");
        assert_eq!(resolved_sg2.unwrap().name, "SG");
    }

    #[test]
    fn test_inter_chain_bond_routing_and_reparent() {
        // Setup model with chain A and chain B
        let mut model = AtomGroup::new();
        let mut chain_a = AtomGroup::new();
        let mut chain_b = AtomGroup::new();
        let mut res_a = AtomGroup::new();
        let mut res_b = AtomGroup::new();

        let mut sg_a = Atom::from_symbol("S").unwrap();
        sg_a.name = "SG".to_string();
        res_a.set_atom("1_SG", sg_a);

        let mut sg_b = Atom::from_symbol("S").unwrap();
        sg_b.name = "SG".to_string();
        res_b.set_atom("2_SG", sg_b);

        chain_a.set_group("1", res_a);
        chain_b.set_group("2", res_b);
        model.set_group("A", chain_a);
        model.set_group("B", chain_b);

        let sg_a_ref = model.get_atom_by_path("/A/1/1_SG").unwrap().clone();
        let sg_b_ref = model.get_atom_by_path("/B/2/2_SG").unwrap().clone();

        // Add inter-chain bond at model level (nearest common ancestor is model itself)
        model.add_bond(&sg_a_ref, &sg_b_ref, 1);
        assert_eq!(model.get_number_of_bonds(), 1);

        // Reparent into root
        let mut root = AtomGroup::new();
        root.set_group("model_1", model);

        let bonds = root.get_bond_list();
        assert_eq!(bonds.len(), 1);
        assert_eq!(bonds[0].atom1_path, "/model_1/A/1/1_SG");
        assert_eq!(bonds[0].atom2_path, "/model_1/B/2/2_SG");

        let resolved_a = root.get_atom_by_path(&bonds[0].atom1_path);
        let resolved_b = root.get_atom_by_path(&bonds[0].atom2_path);
        assert!(resolved_a.is_some());
        assert!(resolved_b.is_some());
    }

    #[test]
    fn test_schema_level_checks_normal_hierarchy() {
        use crate::format::Format;

        let mut root = AtomGroup::new();
        let mut model = AtomGroup::with_name("model_1");
        let mut chain = AtomGroup::with_name("A");
        let mut residue = AtomGroup::with_name("6");
        let atom = Atom::from_symbol("C").unwrap();

        residue.set_atom("CA", atom);
        chain.set_group("6", residue);
        model.set_group("A", chain);
        root.set_group("model_1", model);

        // Root (depth 0)
        assert_eq!(root.path_depth(), 0);
        assert!(!root.is_model_level());
        assert!(!root.is_chain_level());
        assert!(!root.is_residue_level());

        // Model (depth 1)
        let m = root.get_group("model_1").unwrap();
        assert_eq!(m.path_depth(), 1);
        assert!(m.is_model_level());
        assert!(!m.is_chain_level());
        assert!(!m.is_residue_level());
        assert!(Format::is_protein(m));

        // Chain (depth 2)
        let c = m.get_group("A").unwrap();
        assert_eq!(c.path_depth(), 2);
        assert!(!c.is_model_level());
        assert!(c.is_chain_level());
        assert!(!c.is_residue_level());
        assert!(Format::is_chain(c));

        // Residue (depth 3)
        let r = c.get_group("6").unwrap();
        assert_eq!(r.path_depth(), 3);
        assert!(!r.is_model_level());
        assert!(!r.is_chain_level());
        assert!(r.is_residue_level());
        assert!(Format::is_residue(r));

        // Validation on clean hierarchy produces no violations
        let violations = root.validate_schema();
        assert!(
            violations.is_empty(),
            "Expected no violations, got {:?}",
            violations
        );
    }

    #[test]
    fn test_schema_violations_direct_atoms_in_chain() {
        use crate::format::Format;

        let mut root = AtomGroup::new();
        let mut model = AtomGroup::with_name("model_1");
        let mut chain = AtomGroup::with_name("A");
        let mut residue = AtomGroup::with_name("6");

        let ca = Atom::from_symbol("C").unwrap();
        residue.set_atom("CA", ca);
        chain.set_group("6", residue);

        // Intentionally violate schema: attach HETATM / water directly under chain without a residue
        let mut water = Atom::from_symbol("O").unwrap();
        water.name = "O".to_string();
        chain.set_atom("HOH_1", water);

        model.set_group("A", chain);
        root.set_group("model_1", model);

        let c = root.get_group("model_1").unwrap().get_group("A").unwrap();
        // Positional check still sees depth 2 (chain level)
        assert!(c.is_chain_level());
        // But structural check fails due to direct atoms
        assert!(!Format::is_chain(c));

        // validate_schema() detects the violation
        let violations = root.validate_schema();
        assert_eq!(violations.len(), 1);
        match &violations[0] {
            SchemaViolation::DirectAtomsAtNonResidueLevel {
                path,
                depth,
                atom_keys,
            } => {
                assert_eq!(path, "/model_1/A/");
                assert_eq!(*depth, 2);
                assert_eq!(atom_keys, &vec!["HOH_1".to_string()]);
            }
            other => panic!("Unexpected violation: {:?}", other),
        }
    }

    #[test]
    fn test_schema_violations_subgroups_in_residue_and_excessive_depth() {
        let mut root = AtomGroup::new();
        let mut model = AtomGroup::with_name("model_1");
        let mut chain = AtomGroup::with_name("A");
        let mut residue = AtomGroup::with_name("6");
        let mut sub_residue = AtomGroup::with_name("sub");

        let atom = Atom::from_symbol("C").unwrap();
        sub_residue.set_atom("C1", atom);
        residue.set_group("sub", sub_residue);
        chain.set_group("6", residue);
        model.set_group("A", chain);
        root.set_group("model_1", model);

        let violations = root.validate_schema();
        // Should detect:
        // 1. SubgroupsInResidue at "/model_1/A/6/"
        // 2. ExcessiveDepth at "/model_1/A/6/sub/" (depth 4)
        // 3. DirectAtomsAtNonResidueLevel at "/model_1/A/6/sub/" (depth 4 has direct atoms)
        assert!(violations.iter().any(|v| matches!(
            v,
            SchemaViolation::SubgroupsInResidue { path, group_keys }
                if path == "/model_1/A/6/" && group_keys == &vec!["sub".to_string()]
        )));
        assert!(violations.iter().any(|v| matches!(
            v,
            SchemaViolation::ExcessiveDepth { path, depth }
                if path == "/model_1/A/6/sub/" && *depth == 4
        )));
        assert!(violations.iter().any(|v| matches!(
            v,
            SchemaViolation::DirectAtomsAtNonResidueLevel { path, depth, atom_keys }
                if path == "/model_1/A/6/sub/" && *depth == 4 && atom_keys == &vec!["C1".to_string()]
        )));
    }

    #[test]
    fn test_schema_regression_key_with_slash() {
        use crate::format::Format;

        let mut root = AtomGroup::new();
        let mut model = AtomGroup::with_name("model_1");
        // Key with slash: "A/B"
        let mut chain = AtomGroup::with_name("A/B");

        // Directly attach atom under chain (schema violation)
        let mut atom = Atom::from_symbol("O").unwrap();
        atom.name = "O".to_string();
        chain.set_atom("HOH_1", atom);

        model.set_group("A/B", chain);
        root.set_group("model_1", model);

        let c = root.get_group("model_1").unwrap().get_group("A/B").unwrap();
        // The tree depth of this chain is 2 (root=0 -> model=1 -> chain=2)
        assert_eq!(c.path_depth(), 2);
        assert!(c.is_chain_level());
        assert!(!c.is_residue_level());
        assert!(!Format::is_chain(c));

        // validate_schema() must NOT report 0 violations due to slash miscount;
        // it must detect DirectAtomsAtNonResidueLevel at depth 2
        let violations = root.validate_schema();
        assert_eq!(violations.len(), 1);
        match &violations[0] {
            SchemaViolation::DirectAtomsAtNonResidueLevel {
                path,
                depth,
                atom_keys,
            } => {
                assert_eq!(path, "/model_1/A/B/");
                assert_eq!(*depth, 2);
                assert_eq!(atom_keys, &vec!["HOH_1".to_string()]);
            }
            other => panic!("Unexpected violation: {:?}", other),
        }
    }

    #[test]
    fn test_schema_regression_empty_string_keys() {
        let mut root = AtomGroup::new();
        let mut g1 = AtomGroup::new(); // depth 1
        let mut g2 = AtomGroup::new(); // depth 2
        let mut g3 = AtomGroup::new(); // depth 3
        let mut g4 = AtomGroup::new(); // depth 4 (excessive depth)

        let atom = Atom::from_symbol("C").unwrap();
        g4.set_atom("C1", atom);

        g3.set_group("", g4);
        g2.set_group("", g3);
        g1.set_group("", g2);
        root.set_group("", g1);

        let violations = root.validate_schema();
        // g4 is at depth 4: must detect ExcessiveDepth with depth == 4
        // and DirectAtomsAtNonResidueLevel with depth == 4
        assert!(violations.iter().any(|v| matches!(
            v,
            SchemaViolation::ExcessiveDepth { depth, .. } if *depth == 4
        )));
        assert!(violations.iter().any(|v| matches!(
            v,
            SchemaViolation::DirectAtomsAtNonResidueLevel { depth, atom_keys, .. }
                if *depth == 4 && atom_keys == &vec!["C1".to_string()]
        )));
    }
}
