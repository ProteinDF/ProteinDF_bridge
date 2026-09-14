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

pub mod gro;
pub mod xyz;

pub use gro::SimpleGro;
pub use xyz::Xyz;

use crate::atom_group::AtomGroup;

/// Hierarchy validation helper corresponding to `proteindf_bridge.format.Format`.
pub struct Format;

impl Format {
    /// Checks whether an AtomGroup represents a residue (no subgroups; atom count does not affect validity in Python).
    pub fn is_residue(res: &AtomGroup) -> bool {
        res.get_number_of_groups() == 0
    }

    /// Checks whether an AtomGroup represents a chain (no direct atoms, all subgroups are residues).
    pub fn is_chain(chain: &AtomGroup) -> bool {
        if chain.get_number_of_atoms() != 0 {
            return false;
        }
        chain.groups().all(|(_, res)| Self::is_residue(res))
    }

    /// Checks whether an AtomGroup represents a protein/model (no direct atoms, all subgroups are chains).
    pub fn is_protein(model: &AtomGroup) -> bool {
        if model.get_number_of_atoms() != 0 {
            return false;
        }
        model.groups().all(|(_, chain)| Self::is_chain(chain))
    }

    /// Checks whether an AtomGroup represents a multi-model collection (no direct atoms, all subgroups are proteins).
    pub fn is_models(models: &AtomGroup) -> bool {
        if models.get_number_of_atoms() != 0 {
            return false;
        }
        models.groups().all(|(_, model)| Self::is_protein(model))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;

    #[test]
    fn test_format_empty() {
        let empty = AtomGroup::new();
        // Empty group has no subgroups and no atoms -> matches residue, chain, protein, and models in Python
        assert!(Format::is_residue(&empty));
        assert!(Format::is_chain(&empty));
        assert!(Format::is_protein(&empty));
        assert!(Format::is_models(&empty));
    }

    #[test]
    fn test_format_hierarchy() {
        // Build hierarchy: models -> model_1 -> A -> 3 -> atom
        let mut models = AtomGroup::with_name("models");
        let mut model = AtomGroup::with_name("model_1");
        let mut chain = AtomGroup::with_name("A");
        let mut res = AtomGroup::with_name("3");

        let atom = Atom::from_symbol("C").unwrap();
        res.set_atom("CA", atom);

        chain.set_group("3", res.clone());
        model.set_group("A", chain.clone());
        models.set_group("model_1", model.clone());

        // Test models
        assert!(Format::is_models(&models));
        assert!(!Format::is_protein(&models));
        assert!(!Format::is_chain(&models));
        assert!(!Format::is_residue(&models));

        // Test model
        assert!(!Format::is_models(&model));
        assert!(Format::is_protein(&model));
        assert!(!Format::is_chain(&model));
        assert!(!Format::is_residue(&model));

        // Test chain
        assert!(!Format::is_models(&chain));
        assert!(!Format::is_protein(&chain));
        assert!(Format::is_chain(&chain));
        assert!(!Format::is_residue(&chain));

        // Test residue
        assert!(!Format::is_models(&res));
        assert!(!Format::is_protein(&res));
        assert!(!Format::is_chain(&res));
        assert!(Format::is_residue(&res));
    }
}
