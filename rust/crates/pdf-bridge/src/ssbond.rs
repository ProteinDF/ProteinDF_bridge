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

use crate::atom::Atom;
use crate::atom_group::AtomGroup;

/// Finds disulfide bonds in protein models corresponding to `proteindf_bridge.ssbond.SSBond`.
#[derive(Debug, Clone)]
pub struct SSBond {
    model: AtomGroup,
    ssbonds: Vec<(String, String)>,
    is_checked: bool,
}

impl SSBond {
    /// Default maximum distance between sulfur atoms (SG-SG) to qualify as a disulfide bond: 2.1 * 1.1 = 2.31 Å.
    pub const SS_BOND_MAX_LENGTH: f64 = 2.1 * 1.1;

    /// Creates a new `SSBond` detector for the given protein model.
    pub fn new(model: &AtomGroup) -> Self {
        Self {
            model: model.clone(),
            ssbonds: Vec::new(),
            is_checked: false,
        }
    }

    /// Returns the list of detected disulfide bond residue path pairs `(path1, path2)`.
    pub fn get_bonds(&mut self) -> &[(String, String)] {
        if !self.is_checked {
            self.check();
        }
        &self.ssbonds
    }

    /// Static convenience function to find disulfide bonds in a model.
    pub fn find_bonds(model: &AtomGroup) -> Vec<(String, String)> {
        let mut ssbond = Self::new(model);
        ssbond.get_bonds().to_vec()
    }

    fn check(&mut self) {
        let mut sgs: Vec<(String, Atom)> = Vec::new();

        for (_chain_id, chain) in self.model.groups() {
            for (_res_key, res) in chain.groups() {
                if res.name == "CYS" || res.name == "CYX" {
                    if let Some(sg) = res.get_atom("SG") {
                        sgs.push((res.path().to_string(), sg.clone()));
                    }
                }
            }
        }

        self.check_sgs(&sgs);
        self.is_checked = true;
    }

    fn check_sgs(&mut self, sgs: &[(String, Atom)]) {
        for (i, (path1, sg1)) in sgs.iter().enumerate() {
            for (path2, sg2) in sgs.iter().skip(i + 1) {
                let distance = sg1.xyz.distance_from(&sg2.xyz);
                if distance < Self::SS_BOND_MAX_LENGTH {
                    self.ssbonds.push((path1.clone(), path2.clone()));
                }
            }
        }
    }
}
