// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::str::FromStr;

use regex::Regex;

use crate::atom::Atom;
use crate::atom_group::{AtomGroup, Selector};
use crate::error::{BridgeError, Result};
use crate::position::Position;

/// Selects atoms by their atomic symbol (case-insensitive).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SelectSymbol {
    atom_symbol: String,
}

impl SelectSymbol {
    pub fn new(symbol: &str) -> Self {
        Self {
            atom_symbol: symbol.trim().to_ascii_uppercase(),
        }
    }
}

impl Selector for SelectSymbol {
    fn is_match_atom(&self, atom: &Atom) -> bool {
        atom.symbol()
            .map(|s| s.to_ascii_uppercase() == self.atom_symbol)
            .unwrap_or(false)
    }
}

/// Selects atoms or groups by name (exact trimmed match).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SelectName {
    query: String,
}

impl SelectName {
    pub fn new(query: &str) -> Self {
        Self {
            query: query.trim().to_string(),
        }
    }
}

impl Selector for SelectName {
    fn is_match_atom(&self, atom: &Atom) -> bool {
        atom.name.trim() == self.query
    }

    fn is_match_group(&self, group: &AtomGroup) -> bool {
        group.name.trim() == self.query
    }
}

/// Selects atoms or groups by exact path match.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SelectPathSimple {
    query: String,
}

impl SelectPathSimple {
    pub fn new(query: &str) -> Self {
        Self {
            query: query.to_string(),
        }
    }
}

impl Selector for SelectPathSimple {
    fn is_match_atom(&self, atom: &Atom) -> bool {
        atom.path == self.query
    }

    fn is_match_group(&self, group: &AtomGroup) -> bool {
        group.path() == self.query
    }
}

/// Selects atoms or groups by path matching a regular expression.
#[derive(Debug, Clone)]
pub struct SelectPathRegex {
    query: String,
    regex: Regex,
}

impl SelectPathRegex {
    pub fn new(query: &str) -> Result<Self> {
        let regex = Regex::new(query).map_err(|e| {
            BridgeError::input_error(query, format!("invalid regular expression: {e}"))
        })?;
        Ok(Self {
            query: query.to_string(),
            regex,
        })
    }

    pub fn query(&self) -> &str {
        &self.query
    }
}

impl Selector for SelectPathRegex {
    fn is_match_atom(&self, atom: &Atom) -> bool {
        self.regex.is_match(&atom.path)
    }

    fn is_match_group(&self, group: &AtomGroup) -> bool {
        self.regex.is_match(group.path())
    }
}

/// Selects atoms or groups by path with wildcard support (`*` and `?`).
#[derive(Debug, Clone)]
pub struct SelectPathWildcard {
    query: String,
    regex_selector: SelectPathRegex,
}

impl SelectPathWildcard {
    pub fn new(query: &str) -> Result<Self> {
        let pattern = prepare_wildcard(query);
        let regex_selector = SelectPathRegex::new(&pattern)?;
        Ok(Self {
            query: query.to_string(),
            regex_selector,
        })
    }

    pub fn query(&self) -> &str {
        &self.query
    }
}

impl Selector for SelectPathWildcard {
    fn is_match_atom(&self, atom: &Atom) -> bool {
        self.regex_selector.is_match_atom(atom)
    }

    fn is_match_group(&self, group: &AtomGroup) -> bool {
        self.regex_selector.is_match_group(group)
    }
}

fn prepare_wildcard(query: &str) -> String {
    let mut result = String::from("^");
    let mut chars = query.chars().peekable();
    while let Some(c) = chars.next() {
        if c == '\\' {
            if let Some(&next) = chars.peek() {
                if next == '*' || next == '?' {
                    result.push('\\');
                    result.push(next);
                    chars.next();
                    continue;
                }
            }
            result.push('\\');
        } else if c == '*' {
            result.push_str(".*");
        } else if c == '?' {
            result.push('?');
        } else {
            result.push(c);
        }
    }
    result.push('$');
    result
}

/// Deprecated path selector corresponding to Python's `Select_Path`.
#[derive(Debug, Clone)]
pub struct SelectPath {
    wildcard_selector: Option<SelectPathWildcard>,
    simple_selector: Option<SelectPathSimple>,
}

impl SelectPath {
    pub fn new(query: &str, use_wildcard: bool) -> Result<Self> {
        if use_wildcard {
            Ok(Self {
                wildcard_selector: Some(SelectPathWildcard::new(query)?),
                simple_selector: None,
            })
        } else {
            Ok(Self {
                wildcard_selector: None,
                simple_selector: Some(SelectPathSimple::new(query)),
            })
        }
    }
}

impl Selector for SelectPath {
    fn is_match_atom(&self, atom: &Atom) -> bool {
        if let Some(ref w) = self.wildcard_selector {
            w.is_match_atom(atom)
        } else if let Some(ref s) = self.simple_selector {
            s.is_match_atom(atom)
        } else {
            false
        }
    }

    fn is_match_group(&self, group: &AtomGroup) -> bool {
        if let Some(ref w) = self.wildcard_selector {
            w.is_match_group(group)
        } else if let Some(ref s) = self.simple_selector {
            s.is_match_group(group)
        } else {
            false
        }
    }
}

/// Selects atoms within a spherical distance from a reference position.
#[derive(Debug, Clone, PartialEq)]
pub struct SelectRange {
    pos: Position,
    d: f64,
    d2: f64,
}

impl SelectRange {
    pub fn new(pos: Position, d: f64) -> Self {
        Self { pos, d, d2: d * d }
    }

    pub fn from_str(pos_str: &str, d: f64) -> Result<Self> {
        let pos = Position::from_str(pos_str)?;
        Ok(Self::new(pos, d))
    }

    pub fn center(&self) -> &Position {
        &self.pos
    }

    pub fn radius(&self) -> f64 {
        self.d
    }
}

impl Selector for SelectRange {
    fn is_match_atom(&self, atom: &Atom) -> bool {
        self.pos.square_distance_from(&atom.xyz) < self.d2
    }
}

/// Selects atoms having the same atomic number and within a distance from a reference atom.
#[derive(Debug, Clone, PartialEq)]
pub struct SelectAtom {
    atom: Atom,
    distance2: f64,
}

impl SelectAtom {
    pub fn new(atom: Atom, distance: f64) -> Self {
        Self {
            atom,
            distance2: distance * distance,
        }
    }
}

impl Selector for SelectAtom {
    fn is_match_atom(&self, atom: &Atom) -> bool {
        self.atom.atomic_number() == atom.atomic_number()
            && self.atom.xyz.square_distance_from(&atom.xyz) < self.distance2
    }
}

/// Selects atoms that also exist in a reference `AtomGroup` (same atomic number within distance range).
#[derive(Debug, Clone, PartialEq)]
pub struct SelectAtomGroup {
    ref_atoms: Vec<Atom>,
    range: f64,
}

impl SelectAtomGroup {
    pub fn new(ref_atomgroup: &AtomGroup, range: f64) -> Self {
        Self {
            ref_atoms: ref_atomgroup.get_atom_list(),
            range,
        }
    }

    pub fn with_default_range(ref_atomgroup: &AtomGroup) -> Self {
        Self::new(ref_atomgroup, 1.0e-5)
    }
}

impl Selector for SelectAtomGroup {
    fn is_match_atom(&self, atom: &Atom) -> bool {
        self.ref_atoms.iter().any(|ref_atom| {
            ref_atom.atomic_number() == atom.atomic_number()
                && ref_atom.xyz.distance_from(&atom.xyz) < self.range
        })
    }
}

// ---------------------------------------------------------------------------
// Python 1:1 naming compatibility aliases
// ---------------------------------------------------------------------------
#[allow(non_camel_case_types)]
pub type Select_Symbol = SelectSymbol;
#[allow(non_camel_case_types)]
pub type Select_Name = SelectName;
#[allow(non_camel_case_types)]
pub type Select_Path_simple = SelectPathSimple;
#[allow(non_camel_case_types)]
pub type Select_Path_wildcard = SelectPathWildcard;
#[allow(non_camel_case_types)]
pub type Select_Path = SelectPath;
#[allow(non_camel_case_types)]
pub type Select_PathRegex = SelectPathRegex;
#[allow(non_camel_case_types)]
pub type Select_Range = SelectRange;
#[allow(non_camel_case_types)]
pub type Select_Atom = SelectAtom;
#[allow(non_camel_case_types)]
pub type Select_AtomGroup = SelectAtomGroup;
