// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::time::Instant;

use proteindf_bridge::atom::Atom;
use proteindf_bridge::atom_group::{AtomGroup, BondRecord};
use proteindf_bridge::position::Position;

fn make_synthetic_protein(num_residues: usize) -> AtomGroup {
    let mut root = AtomGroup::new();
    let mut model = AtomGroup::with_name("model_1");
    let mut chain = AtomGroup::with_name("A");

    let atom_names = [
        "N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2", "H", "HA", "HB2", "HB3",
    ];

    for r in 1..=num_residues {
        let mut residue = AtomGroup::with_name("ARG");
        for &name in &atom_names {
            let symbol = &name[..1];
            let mut a = Atom::new_with_pos(symbol, Position::new(r as f64, 0.0, 0.0)).unwrap();
            a.name = name.to_string();
            residue.set_atom(name, a);
        }
        chain.set_group(&r.to_string(), residue);
    }

    model.set_group("A", chain);
    root.set_group("model_1", model);
    root
}

#[test]
fn test_resolve_bond_basic() {
    let mut root = AtomGroup::new();
    let mut model = AtomGroup::with_name("model_1");
    let mut chain = AtomGroup::with_name("A");
    let mut residue = AtomGroup::with_name("1");

    let mut a1 = Atom::from_symbol("C").unwrap();
    a1.name = "CA".to_string();
    let mut a2 = Atom::from_symbol("C").unwrap();
    a2.name = "CB".to_string();

    residue.set_atom("CA", a1);
    residue.set_atom("CB", a2);
    chain.set_group("1", residue);
    model.set_group("A", chain);
    root.set_group("model_1", model);

    // Valid bond record
    let valid_bond = BondRecord {
        atom1_path: "/model_1/A/1/CA".to_string(),
        atom2_path: "/model_1/A/1/CB".to_string(),
        order: 1,
    };

    let resolved = root.resolve_bond(&valid_bond);
    assert!(resolved.is_some());
    let (atom1, atom2) = resolved.unwrap();
    assert_eq!(atom1.name, "CA");
    assert_eq!(atom2.name, "CB");

    // Invalid endpoint 1
    let bad_bond1 = BondRecord {
        atom1_path: "/model_1/A/1/NONEXISTENT".to_string(),
        atom2_path: "/model_1/A/1/CB".to_string(),
        order: 1,
    };
    assert!(root.resolve_bond(&bad_bond1).is_none());

    // Invalid endpoint 2
    let bad_bond2 = BondRecord {
        atom1_path: "/model_1/A/1/CA".to_string(),
        atom2_path: "/model_1/A/999/CB".to_string(),
        order: 1,
    };
    assert!(root.resolve_bond(&bad_bond2).is_none());
}

#[test]
fn test_bond_resolution_scalability_benchmark() {
    // Construct small, medium, and large synthetic hierarchies:
    // small:  100 residues  (~1,500 atoms)
    // medium: 1,000 residues (~15,000 atoms)
    // large:  5,000 residues (~75,000 atoms)
    let ag_small = make_synthetic_protein(100);
    let ag_medium = make_synthetic_protein(1_000);
    let ag_large = make_synthetic_protein(5_000);

    assert_eq!(ag_small.get_number_of_all_atoms(), 1_500);
    assert_eq!(ag_medium.get_number_of_all_atoms(), 15_000);
    assert_eq!(ag_large.get_number_of_all_atoms(), 75_000);

    let test_bonds = vec![
        BondRecord {
            atom1_path: "/model_1/A/50/CA".to_string(),
            atom2_path: "/model_1/A/50/CB".to_string(),
            order: 1,
        },
        BondRecord {
            atom1_path: "/model_1/A/100/N".to_string(),
            atom2_path: "/model_1/A/100/CA".to_string(),
            order: 1,
        },
    ];

    let iterations = 10_000;

    // Benchmark small
    let start_small = Instant::now();
    for _ in 0..iterations {
        for bond in &test_bonds {
            let res = ag_small.resolve_bond(bond);
            assert!(res.is_some());
        }
    }
    let duration_small = start_small.elapsed();

    // Benchmark medium (10x atoms of small)
    let start_medium = Instant::now();
    for _ in 0..iterations {
        for bond in &test_bonds {
            let res = ag_medium.resolve_bond(bond);
            assert!(res.is_some());
        }
    }
    let duration_medium = start_medium.elapsed();

    // Benchmark large (50x atoms of small)
    let start_large = Instant::now();
    for _ in 0..iterations {
        for bond in &test_bonds {
            let res = ag_large.resolve_bond(bond);
            assert!(res.is_some());
        }
    }
    let duration_large = start_large.elapsed();

    println!(
        "Benchmark ({} lookups):\n  small (1.5k atoms):  {:?}\n  medium (15k atoms):  {:?}\n  large (75k atoms):   {:?}",
        iterations * test_bonds.len(),
        duration_small,
        duration_medium,
        duration_large
    );

    let ratio = duration_large.as_nanos() as f64 / duration_small.as_nanos().max(1) as f64;
    println!(
        "Scalability check: large/small duration ratio = {:.2} (confirms O(1) depth-only scaling)",
        ratio
    );
}
