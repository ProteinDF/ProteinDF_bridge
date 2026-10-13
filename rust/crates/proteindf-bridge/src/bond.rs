// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::collections::HashSet;

use crate::atom::Atom;
use crate::atom_group::AtomGroup;
use crate::error::Result;
use crate::format::Format;
use crate::matrix::SymmetricMatrix;
use crate::position::Position;
use crate::spatial::CellList;

/// Maximum number of atoms for which dense $O(N^2)$ matrices (`distmat` and `bondmat`) are constructed.
///
/// For structures with $N > \text{MAX_DENSE_MATRIX_ATOMS}$, `distmat` and `bondmat` are set to `None`
/// to protect against excessive memory consumption ($O(N^2)$ allocation). Bond detection continues
/// to operate in $O(N)$ time and memory via spatial cell lists ([`CellList`]).
pub const MAX_DENSE_MATRIX_ATOMS: usize = 2000;

/// Default tolerance added to sum of covalent radii for bond detection (in Å).
///
/// Based on OpenBabel's standard convention (`OBAtom::ConnectsTo`, which uses
/// $r_A + r_B + 0.45$ Å). This tolerance accommodates thermal vibration and experimental
/// uncertainty in macromolecular structures while strictly preventing non-bonded van der Waals
/// contacts and hydrogen bonds (~2.6–3.5 Å) from being falsely detected as covalent bonds.
pub const COVALENT_BOND_TOLERANCE: f64 = 0.45;

/// Bond detector based on covalent radii, corresponding to `proteindf_bridge.bond.Bond`.
///
/// Bonds are established between pairs of atoms $(p, q)$ satisfying:
/// $$r_{pq} \le \text{cov}_p + \text{cov}_q + \text{COVALENT_BOND_TOLERANCE}$$
///
/// # Intentional Deviation from Python Version
/// In Python `proteindf_bridge.bond.Bond`, bond detection was based on van der Waals radii
/// ($r_{pq} \le \text{vdw}_p + \text{vdw}_q + 0.4$). Because VDW radii represent non-bonded contact
/// distances (e.g. C-C cutoff was $1.70 + 1.70 + 0.4 = 3.8$ Å), that heuristic frequently falsely
/// identified hydrogen bonds (~2.6–3.5 Å) and steric VDW packing as covalent bonds.
/// In this Rust port, following modern cheminformatics standards (OpenBabel, ASE `natural_cutoffs`,
/// pymatgen), detection is based on Cordero et al. (2008) covalent radii plus a tolerance of 0.45 Å.
///
/// # Performance & Scalability
/// Bond detection uses an $O(N)$ uniform spatial cell list ([`CellList`]) with dynamically
/// determined cell size based on the maximum covalent radius in the atom set.
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

    /// Sets up bonds for the given `AtomGroup` based on covalent radii heuristics.
    ///
    /// Uses an $O(N)$ spatial cell list with dynamically calculated cell size.
    /// Dense matrices (`distmat`/`bondmat`) are allocated only when $N \le \text{MAX_DENSE_MATRIX_ATOMS}$.
    ///
    /// Pre-existing bonds in `mol` (e.g., from CCD templates or file-derived CONECT/bonds) are
    /// preserved without duplicate registration or overwriting existing bond orders.
    pub fn setup_heuristic(&mut self, mol: &mut AtomGroup) -> Result<()> {
        mol.update_paths();
        self.atoms = mol.get_atom_list();
        let n = self.atoms.len();
        if n == 0 {
            self.distmat = Some(SymmetricMatrix::new(0));
            self.bondmat = Some(SymmetricMatrix::new(0));
            return Ok(());
        }

        // Collect existing bonds in mol to avoid duplicate registrations
        let mut existing_bonds: HashSet<(String, String)> = HashSet::new();
        for b in mol.get_bond_list() {
            let p1 = b.atom1_path.clone();
            let p2 = b.atom2_path.clone();
            if p1 <= p2 {
                existing_bonds.insert((p1, p2));
            } else {
                existing_bonds.insert((p2, p1));
            }
        }

        // Collect covalent radii and find maximum covalent radius to dynamically size the cell list
        let mut cov_radii = Vec::with_capacity(n);
        let mut max_cov = 0.0_f64;
        for atom in &self.atoms {
            let r = atom.covalent_radius()?;
            if r > max_cov {
                max_cov = r;
            }
            cov_radii.push(r);
        }

        // Dynamically determine cell size: must cover max possible cutoff (2 * max_cov + COVALENT_BOND_TOLERANCE)
        let max_cutoff = 2.0 * max_cov + COVALENT_BOND_TOLERANCE;
        let cell_size = max_cutoff.max(3.0);

        // Determine whether mol represents a multi-model collection
        let model_ids: Option<Vec<usize>> = if Format::is_models(mol) {
            let mut ids = Vec::with_capacity(n);
            for (model_idx, (_, child)) in mol.groups().enumerate() {
                ids.resize(ids.len() + child.get_number_of_all_atoms(), model_idx);
            }
            Some(ids)
        } else {
            None
        };

        // Build dense matrices only for small structures (backward compatibility)
        if n <= MAX_DENSE_MATRIX_ATOMS {
            let mut distmat = SymmetricMatrix::new(n);
            let mut bondmat = SymmetricMatrix::new(n);
            for p in 0..n {
                let cov_p = cov_radii[p];
                for (q, &cov_q) in cov_radii.iter().enumerate().take(p) {
                    let d = self.atoms[p].xyz.distance_from(&self.atoms[q].xyz);
                    distmat.set(p, q, d);
                    let cross_model = match &model_ids {
                        Some(ids) => ids[p] != ids[q],
                        None => false,
                    };
                    if !cross_model && d <= (cov_p + cov_q) + COVALENT_BOND_TOLERANCE {
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
            let cross_model = match &model_ids {
                Some(ids) => ids[p] != ids[q],
                None => false,
            };
            if !cross_model && dist <= (cov_radii[p] + cov_radii[q]) + COVALENT_BOND_TOLERANCE {
                let p1 = &self.atoms[p].path;
                let p2 = &self.atoms[q].path;
                let key = if p1 <= p2 {
                    (p1.clone(), p2.clone())
                } else {
                    (p2.clone(), p1.clone())
                };
                if !existing_bonds.contains(&key) {
                    existing_bonds.insert(key);
                    bonds_to_add.push((p, q));
                }
            }
        });

        bonds_to_add.sort_unstable();

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
        bond.setup_heuristic(&mut ag).unwrap();

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
        bond.setup_heuristic(&mut ag).unwrap();

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
            // Line of carbon atoms spaced by 1.5 Å (covalent cutoff = 0.76 + 0.76 + 0.45 = 1.97 Å):
            // Only adjacent pairs (1.5 Å <= 1.97 Å) are bonded; next-nearest (3.0 Å > 1.97 Å) are not.
            // (Updated from 2.5 Å in legacy VDW test, as 2.5 Å is non-bonded under covalent radius)
            let atom = Atom::new_with_pos("C", Position::new(i as f64 * 1.5, 0.0, 0.0)).unwrap();
            ag.set_atom(&i.to_string(), atom);
        }

        let mut bond = Bond::new();
        bond.setup_heuristic(&mut ag).unwrap();

        // Dense matrices should be None to protect memory
        assert!(bond.distmat.is_none());
        assert!(bond.bondmat.is_none());

        // Only adjacent atoms (distance 1.5 Å <= 1.97 Å, next is 3.0 Å > 1.97 Å) are bonded: exactly n - 1 bonds
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
                    let pos = Position::new(x as f64 * 1.5, y as f64 * 1.5, z as f64 * 1.5);
                    let atom = Atom::new_with_pos(symbol, pos).unwrap();
                    ag.set_atom(&idx.to_string(), atom.clone());
                    atoms.push(atom);
                    idx += 1;
                }
            }
        }

        // 1. Brute-force O(N^2) bond collection using covalent radii
        let n = atoms.len();
        let mut brute_bonds = Vec::new();
        for i in 0..n {
            let cov_i = atoms[i].covalent_radius().unwrap();
            for j in i + 1..n {
                let cov_j = atoms[j].covalent_radius().unwrap();
                let dist = atoms[i].xyz.distance_from(&atoms[j].xyz);
                if dist <= (cov_i + cov_j) + COVALENT_BOND_TOLERANCE {
                    brute_bonds.push((i, j));
                }
            }
        }
        brute_bonds.sort();

        // 2. CellList-based Bond::setup_heuristic()
        let mut bond = Bond::new();
        bond.setup_heuristic(&mut ag).unwrap();

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

        // 100 x 100 x 100 grid = 1,000,000 atoms spaced by 1.5 Å (within C-C covalent cutoff)
        let n_side = 100;
        let mut ag = AtomGroup::with_name("bench_1m");
        let mut idx = 0;
        for x in 0..n_side {
            for y in 0..n_side {
                for z in 0..n_side {
                    let pos = Position::new(x as f64 * 1.5, y as f64 * 1.5, z as f64 * 1.5);
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
        bond.setup_heuristic(&mut ag).unwrap();
        let bond_time = start_bond.elapsed();

        let num_bonds = ag.bonds().len();
        println!(
            "Bond::setup_heuristic() time for 1,000,000 atoms: {:.2?}",
            bond_time
        );
        println!("Detected bonds count: {}", num_bonds);
        assert!(num_bonds > 0);
    }

    fn make_test_protein(name: &str, pos1: Position, pos2: Position) -> AtomGroup {
        let mut model = AtomGroup::with_name(name);
        let mut chain = AtomGroup::with_name("A");
        let mut res = AtomGroup::with_name("ALA_1");
        let c1 = Atom::new_with_pos("C", pos1).unwrap();
        let c2 = Atom::new_with_pos("C", pos2).unwrap();
        res.set_atom("C1", c1);
        res.set_atom("C2", c2);
        chain.set_group("ALA_1", res);
        model.set_group("A", chain);
        model
    }

    #[test]
    fn test_bond_cross_model_avoided_with_custom_names() {
        // (a) 子グループの名前が model_* でなくても、is_models を満たす2モデルを重ねると別モデル間に結合しない
        let single_m = make_test_protein(
            "prot_A",
            Position::new(0.0, 0.0, 0.0),
            Position::new(1.5, 0.0, 0.0),
        );
        let mut single_root = AtomGroup::with_name("single_root");
        single_root.set_group("prot_A", single_m);
        assert!(Format::is_models(&single_root));
        let mut bond_single = Bond::new();
        bond_single.setup_heuristic(&mut single_root).unwrap();
        let single_bonds = single_root.get_bond_list();
        assert_eq!(single_bonds.len(), 1);

        // 2つのモデル (名前は "prot_A", "prot_B") を同じ座標に重ねる
        let m1 = make_test_protein(
            "prot_A",
            Position::new(0.0, 0.0, 0.0),
            Position::new(1.5, 0.0, 0.0),
        );
        let m2 = make_test_protein(
            "prot_B",
            Position::new(0.0, 0.0, 0.0),
            Position::new(1.5, 0.0, 0.0),
        );
        let mut root = AtomGroup::with_name("root");
        root.set_group("prot_A", m1);
        root.set_group("prot_B", m2);
        assert!(Format::is_models(&root));

        let mut bond = Bond::new();
        bond.setup_heuristic(&mut root).unwrap();

        let bonds = root.get_bond_list();
        assert_eq!(bonds.len(), single_bonds.len() * 2);

        // 各結合が同一モデル内の原子どうしであることを確認
        for b in &bonds {
            let m_a = b
                .atom1_path
                .split('/')
                .find(|s| !s.is_empty() && *s != "root");
            let m_b = b
                .atom2_path
                .split('/')
                .find(|s| !s.is_empty() && *s != "root");
            assert_eq!(m_a, m_b);
        }

        // 密行列 (bondmat) でも別モデル間が 0 であることを確認
        assert!(bond.bondmat.is_some());
        let bondmat = bond.bondmat.as_ref().unwrap();
        for p in 0..bond.atoms.len() {
            let m_p = bond.atoms[p]
                .path
                .split('/')
                .find(|s| !s.is_empty() && *s != "root");
            for q in 0..p {
                let m_q = bond.atoms[q]
                    .path
                    .split('/')
                    .find(|s| !s.is_empty() && *s != "root");
                if m_p != m_q {
                    assert_eq!(bondmat.get(p, q).unwrap() as usize, 0);
                }
            }
        }
    }

    #[test]
    fn test_bond_cross_group_allowed_when_not_models() {
        // (b) is_models を満たさない木では、別の最上位グループの原子どうしでも距離が近ければ結合する
        let mut root = AtomGroup::with_name("root");
        let mut g1 = AtomGroup::with_name("g1");
        let mut g2 = AtomGroup::with_name("g2");
        let c1 = Atom::new_with_pos("C", Position::new(0.0, 0.0, 0.0)).unwrap();
        let c2 = Atom::new_with_pos("C", Position::new(1.5, 0.0, 0.0)).unwrap();
        g1.set_atom("C1", c1);
        g2.set_atom("C2", c2);
        root.set_group("g1", g1);
        root.set_group("g2", g2);
        // ルート直下に原子を追加して is_models を確実に false にする
        let c_extra = Atom::new_with_pos("C", Position::new(100.0, 100.0, 100.0)).unwrap();
        root.set_atom("C_root", c_extra);

        assert!(!Format::is_models(&root));

        let mut bond = Bond::new();
        bond.setup_heuristic(&mut root).unwrap();

        let bonds = root.get_bond_list();
        let cross_bond = bonds.iter().find(|b| {
            let p1 = &b.atom1_path;
            let p2 = &b.atom2_path;
            (p1.contains("g1") && p2.contains("g2")) || (p1.contains("g2") && p2.contains("g1"))
        });
        assert!(
            cross_bond.is_some(),
            "Bonds across groups must form when is_models is false"
        );
    }

    #[test]
    fn test_bond_2mgo_cross_model_avoided_rust() {
        use crate::format::Pdb;
        use std::path::Path;

        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("data")
            .join("2MGO.pdb");
        if !path.exists() {
            return;
        }

        let pdb = Pdb::from_file(&path, None).expect("failed to load 2MGO.pdb");
        let mut ag = pdb
            .get_atomgroup(None, None)
            .expect("failed to get atomgroup");
        assert!(Format::is_models(&ag));

        let mut bond = Bond::new();
        bond.setup_heuristic(&mut ag).unwrap();

        let bonds = ag.get_bond_list();
        // (c) 2MGO で 2720本・別モデル間0本
        assert_eq!(bonds.len(), 2720);

        for b in &bonds {
            let m_a = b.atom1_path.split('/').find(|s| s.starts_with("model_"));
            let m_b = b.atom2_path.split('/').find(|s| s.starts_with("model_"));
            assert_eq!(m_a, m_b);
        }

        // (d) ag["model_3"] だけを渡すと 136本
        let mut model_3 = ag.get_group("model_3").unwrap().clone();
        assert!(!Format::is_models(&model_3));
        let mut bond3 = Bond::new();
        bond3.setup_heuristic(&mut model_3).unwrap();
        assert_eq!(model_3.get_bond_list().len(), 136);
    }
}
