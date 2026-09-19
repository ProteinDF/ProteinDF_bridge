// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom::Atom;
use crate::atom_group::AtomGroup;
use crate::error::Result;
use crate::matrix::SymmetricMatrix;
use crate::position::Position;
use crate::spatial::CellList;

/// Maximum number of atoms for which dense $O(N^2)$ matrices (`distmat` and `bondmat`) are constructed.
///
/// For structures with $N > \text{MAX_DENSE_MATRIX_ATOMS}$, `distmat` and `bondmat` are set to `None`
/// to protect against excessive memory consumption ($O(N^2)$ allocation). Bond detection continues
/// to operate in $O(N)$ time and memory via spatial cell lists ([`CellList`]).
pub const MAX_DENSE_MATRIX_ATOMS: usize = 2000;

/// Bond detector based on VDW radii, corresponding to `proteindf_bridge.bond.Bond`.
///
/// Bonds are established between pairs of atoms $(p, q)$ satisfying:
/// $$r_{pq} \le \text{vdw}_p + \text{vdw}_q + 0.4$$
///
/// # Performance & Scalability
/// Bond detection uses an $O(N)$ uniform spatial cell list ([`CellList`]) with dynamically
/// determined cell size based on the maximum VDW radius in the atom set.
///
/// Dense matrices (`distmat` and `bondmat`) are populated only when $N \le \text{MAX_DENSE_MATRIX_ATOMS}$
/// (up to 2,000 atoms, ~16 MB). For larger structures, these fields remain `None` to prevent out-of-memory
/// conditions, while the bond topology is directly added to the [`AtomGroup`].
#[derive(Debug, Clone, Default)]
pub struct Bond {
    pub atoms: Vec<Atom>,
    /// Dense distance matrix between atom pairs.
    ///
    /// Set to `Some` if $N \le \text{MAX_DENSE_MATRIX_ATOMS}$, or `None` for larger structures
    /// to avoid $O(N^2)$ memory consumption.
    pub distmat: Option<SymmetricMatrix>,
    /// Dense bond connectivity matrix (1.0 for bonded, 0.0 otherwise).
    ///
    /// Set to `Some` if $N \le \text{MAX_DENSE_MATRIX_ATOMS}$, or `None` for larger structures
    /// to avoid $O(N^2)$ memory consumption.
    pub bondmat: Option<SymmetricMatrix>,
}

impl Bond {
    /// Creates a new `Bond` instance.
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets up bonds for the given `AtomGroup` based on VDW radii.
    ///
    /// Uses an $O(N)$ spatial cell list with dynamically calculated cell size.
    /// Dense matrices (`distmat`/`bondmat`) are allocated only when $N \le \text{MAX_DENSE_MATRIX_ATOMS}$.
    pub fn setup(&mut self, mol: &mut AtomGroup) -> Result<()> {
        self.atoms = mol.get_atom_list();
        let n = self.atoms.len();
        if n == 0 {
            self.distmat = Some(SymmetricMatrix::new(0));
            self.bondmat = Some(SymmetricMatrix::new(0));
            return Ok(());
        }

        // Collect VDW radii and find maximum VDW radius to dynamically size the cell list
        let mut vdws = Vec::with_capacity(n);
        let mut max_vdw = 0.0_f64;
        for atom in &self.atoms {
            let v = atom.vdw()?;
            if v > max_vdw {
                max_vdw = v;
            }
            vdws.push(v);
        }

        // Dynamically determine cell size: must cover max possible cutoff (2 * max_vdw + 0.4)
        let max_cutoff = 2.0 * max_vdw + 0.4;
        let cell_size = max_cutoff.max(3.0);

        // Build dense matrices only for small structures (backward compatibility)
        if n <= MAX_DENSE_MATRIX_ATOMS {
            let mut distmat = SymmetricMatrix::new(n);
            let mut bondmat = SymmetricMatrix::new(n);
            for p in 0..n {
                let vdw_p = vdws[p];
                for (q, &vdw_q) in vdws.iter().enumerate().take(p) {
                    let d = self.atoms[p].xyz.distance_from(&self.atoms[q].xyz);
                    distmat.set(p, q, d);
                    if d <= (vdw_p + vdw_q) + 0.4 {
                        bondmat.set(p, q, 1.0);
                    } else {
                        bondmat.set(p, q, 0.0);
                    }
                }
            }
            self.distmat = Some(distmat);
            self.bondmat = Some(bondmat);
        } else {
            self.distmat = None;
            self.bondmat = None;
        }

        // Build spatial CellList and detect bonds in O(N) time
        let positions: Vec<Position> = self.atoms.iter().map(|a| a.xyz).collect();
        let cell_list = CellList::new(&positions, cell_size);

        let mut bonds_to_add = Vec::new();
        cell_list.for_each_neighbor_pair(max_cutoff, |p, q, dist| {
            // p < q is guaranteed by CellList
            if dist <= (vdws[p] + vdws[q]) + 0.4 {
                bonds_to_add.push((p, q));
            }
        });

        for (p, q) in bonds_to_add {
            mol.add_bond(&self.atoms[p], &self.atoms[q], 1);
        }

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

    #[test]
    fn test_bond_setup_large_threshold() {
        // N > MAX_DENSE_MATRIX_ATOMS (2001 atoms)
        let mut ag = AtomGroup::with_name("large");
        let n = MAX_DENSE_MATRIX_ATOMS + 1;
        for i in 0..n {
            // Line of carbon atoms spaced by 2.5 Å (cutoff is 3.8 Å): only adjacent pairs are bonded
            let atom = Atom::new_with_pos("C", Position::new(i as f64 * 2.5, 0.0, 0.0)).unwrap();
            ag.set_atom(&i.to_string(), atom);
        }

        let mut bond = Bond::new();
        bond.setup(&mut ag).unwrap();

        // Dense matrices should be None to protect memory
        assert!(bond.distmat.is_none());
        assert!(bond.bondmat.is_none());

        // Only adjacent atoms (distance 2.5 Å <= 3.8 Å, next is 5.0 Å > 3.8 Å) are bonded: exactly n - 1 bonds
        assert_eq!(ag.bonds().len(), n - 1);
    }

    #[test]
    fn test_bond_setup_equivalence_with_brute_force() {
        // Equivalence test with ~3,000 atoms (15 x 15 x 14 grid)
        // Verify that CellList detection produces the exact same bond pair set as brute force O(N^2)
        let mut ag = AtomGroup::with_name("equiv_test");
        let mut atoms = Vec::new();

        let nx = 15;
        let ny = 15;
        let nz = 14; // total = 3150 atoms
        let mut idx = 0;
        for x in 0..nx {
            for y in 0..ny {
                for z in 0..nz {
                    let symbol = if (x + y + z) % 3 == 0 {
                        "C"
                    } else if (x + y + z) % 3 == 1 {
                        "N"
                    } else {
                        "O"
                    };
                    let pos = Position::new(x as f64 * 1.8, y as f64 * 1.8, z as f64 * 1.8);
                    let atom = Atom::new_with_pos(symbol, pos).unwrap();
                    ag.set_atom(&idx.to_string(), atom.clone());
                    atoms.push(atom);
                    idx += 1;
                }
            }
        }

        // 1. Brute-force O(N^2) bond collection
        let n = atoms.len();
        let mut brute_bonds = Vec::new();
        for i in 0..n {
            let vdw_i = atoms[i].vdw().unwrap();
            for j in i + 1..n {
                let vdw_j = atoms[j].vdw().unwrap();
                let dist = atoms[i].xyz.distance_from(&atoms[j].xyz);
                if dist <= (vdw_i + vdw_j) + 0.4 {
                    brute_bonds.push((i, j));
                }
            }
        }
        brute_bonds.sort();

        // 2. CellList-based Bond::setup()
        let mut bond = Bond::new();
        bond.setup(&mut ag).unwrap();

        let mut cell_bonds: Vec<(usize, usize)> = ag
            .bonds()
            .iter()
            .map(|b| {
                let i: usize = b.atom1_path.trim_start_matches('/').parse().unwrap();
                let j: usize = b.atom2_path.trim_start_matches('/').parse().unwrap();
                if i < j {
                    (i, j)
                } else {
                    (j, i)
                }
            })
            .collect();
        cell_bonds.sort();

        assert_eq!(
            cell_bonds.len(),
            brute_bonds.len(),
            "bond count mismatch between CellList and brute force"
        );
        assert_eq!(
            cell_bonds, brute_bonds,
            "bond pair set mismatch between CellList and brute force"
        );
    }

    /// Benchmark test for 1,000,000 atoms (run with `cargo test -- --ignored --nocapture`)
    #[test]
    #[ignore]
    fn test_benchmark_1m_atoms() {
        println!("\n--- 1,000,000 Atoms Benchmark ---");
        let start_gen = std::time::Instant::now();

        // 100 x 100 x 100 grid = 1,000,000 atoms
        let n_side = 100;
        let mut ag = AtomGroup::with_name("bench_1m");
        let mut idx = 0;
        for x in 0..n_side {
            for y in 0..n_side {
                for z in 0..n_side {
                    let pos = Position::new(x as f64 * 2.0, y as f64 * 2.0, z as f64 * 2.0);
                    let atom = Atom::new_with_pos("C", pos).unwrap();
                    ag.set_atom(&idx.to_string(), atom);
                    idx += 1;
                }
            }
        }
        let gen_time = start_gen.elapsed();
        println!("Generation time for 1,000,000 atoms: {:.2?}", gen_time);

        let start_bond = std::time::Instant::now();
        let mut bond = Bond::new();
        bond.setup(&mut ag).unwrap();
        let bond_time = start_bond.elapsed();

        let num_bonds = ag.bonds().len();
        println!("Bond::setup() time for 1,000,000 atoms: {:.2?}", bond_time);
        println!("Detected bonds count: {}", num_bonds);
        assert!(num_bonds > 0);
    }
}
