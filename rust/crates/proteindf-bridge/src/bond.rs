// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom::Atom;
use crate::atom_group::AtomGroup;
use crate::error::Result;
use crate::matrix::SymmetricMatrix;

/// Bond detector based on VDW radii, corresponding to `proteindf_bridge.bond.Bond`.
#[derive(Debug, Clone, Default)]
pub struct Bond {
    pub atoms: Vec<Atom>,
    pub distmat: Option<SymmetricMatrix>,
    pub bondmat: Option<SymmetricMatrix>,
}

impl Bond {
    /// Creates a new `Bond` instance.
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets up bonds for the given `AtomGroup`.
    pub fn setup(&mut self, mol: &mut AtomGroup) -> Result<()> {
        self.atoms = mol.get_atom_list();
        self.make_distance_matrix();
        self.make_bond_matrix()?;

        let num_of_atoms = self.atoms.len();
        if let Some(ref bondmat) = self.bondmat {
            for p in 0..num_of_atoms {
                for q in 0..p {
                    let b = bondmat.get(p, q).unwrap_or(0.0);
                    if b > 0.0 {
                        mol.add_bond(&self.atoms[p], &self.atoms[q], b as usize);
                    }
                }
            }
        }
        Ok(())
    }

    fn make_distance_matrix(&mut self) {
        let n = self.atoms.len();
        let mut distmat = SymmetricMatrix::new(n);
        for p in 0..n {
            for q in 0..p {
                let d = self.atoms[p].xyz.distance_from(&self.atoms[q].xyz);
                distmat.set(p, q, d);
            }
        }
        self.distmat = Some(distmat);
    }

    fn make_bond_matrix(&mut self) -> Result<()> {
        let n = self.atoms.len();
        let mut bondmat = SymmetricMatrix::new(n);
        let distmat = self
            .distmat
            .as_ref()
            .expect("Distance matrix must be computed first");

        for p in 0..n {
            let vdw_p = self.atoms[p].vdw()?;
            for q in 0..p {
                let vdw_q = self.atoms[q].vdw()?;
                let r = distmat.get(p, q).unwrap_or(f64::INFINITY);
                if r <= (vdw_p + vdw_q) + 0.4 {
                    bondmat.set(p, q, 1.0);
                } else {
                    bondmat.set(p, q, 0.0);
                }
            }
        }
        self.bondmat = Some(bondmat);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::position::Position;

    // Ported from tests/test_bond.py
    #[test]
    fn test_bond_setup() {
        let mut ag = AtomGroup::with_name("mol");
        // Two carbon atoms at covalent bonding distance (~1.5 Å)
        let c1 = Atom::new_with_pos("C", Position::new(0.0, 0.0, 0.0)).unwrap();
        let c2 = Atom::new_with_pos("C", Position::new(1.5, 0.0, 0.0)).unwrap();
        // A distant carbon atom (10.0 Å)
        let c3 = Atom::new_with_pos("C", Position::new(10.0, 0.0, 0.0)).unwrap();

        ag.set_atom("1", c1);
        ag.set_atom("2", c2);
        ag.set_atom("3", c3);

        let mut bond = Bond::new();
        bond.setup(&mut ag).unwrap();

        assert!(bond.distmat.is_some());
        assert!(bond.bondmat.is_some());

        let distmat = bond.distmat.as_ref().unwrap();
        let bondmat = bond.bondmat.as_ref().unwrap();

        // Distance between c1 and c2 is 1.5
        assert!((distmat.get(1, 0).unwrap() - 1.5).abs() < 1e-10);

        // c1 and c2 are bonded (1), c1/c2 and c3 are not (0)
        assert_eq!(bondmat.get(1, 0).unwrap() as usize, 1);
        assert_eq!(bondmat.get(2, 0).unwrap() as usize, 0);
        assert_eq!(bondmat.get(2, 1).unwrap() as usize, 0);

        // Verify that bonds were added to ag
        assert_eq!(ag.bonds().len(), 1);
    }

    #[test]
    fn test_bond_setup_empty() {
        let mut ag = AtomGroup::with_name("empty");
        let mut bond = Bond::new();
        bond.setup(&mut ag).unwrap();

        assert!(bond.distmat.is_some());
        assert!(bond.bondmat.is_some());
        assert_eq!(bond.atoms.len(), 0);
        assert_eq!(ag.bonds().len(), 0);
    }
}
