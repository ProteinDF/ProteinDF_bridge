// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::path::PathBuf;

use proteindf_bridge::atom::Atom;
use proteindf_bridge::atom_group::AtomGroup;
use proteindf_bridge::format::pdb::Pdb;
use proteindf_bridge::hydrogen_bond::{
    calc_backbone_hbonds, calc_kabsch_sander_energy, calc_pseudo_hydrogen,
};
use proteindf_bridge::position::Position;

fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

fn make_atom(name: &str, symbol: &str, pos: Position) -> Atom {
    let mut a = Atom::new_with_pos(symbol, pos).unwrap();
    a.name = name.to_string();
    a
}

#[test]
fn test_kabsch_sander_energy_benchmark() {
    let path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&path, None).expect("failed to load 1hls.pdb");
    let ag = pdb.get_atomgroup(None, None).expect("failed to get ag");

    let model_1 = ag.get_group("model_1").expect("model_1 not found");
    let chain_a = model_1.get_group("A").expect("chain A not found");

    // Extract residues 1..=8 for benchmark testing
    let get_res_atoms = |key: &str| {
        let res = chain_a.get_group(key).unwrap();
        (
            res.get_atom("N").unwrap().xyz,
            res.get_atom("CA").unwrap().xyz,
            res.get_atom("C").unwrap().xyz,
            res.get_atom("O").unwrap().xyz,
        )
    };

    // Calculate pseudo-H for residues 2..=8
    let mut pseudo_h_map = std::collections::HashMap::new();
    for r in 2..=8 {
        let (_, _, c_prev, _) = get_res_atoms(&(r - 1).to_string());
        let (n_curr, ca_curr, _, _) = get_res_atoms(&r.to_string());
        let h = calc_pseudo_hydrogen(&c_prev, &n_curr, &ca_curr).expect("failed to calc pseudo-H");
        pseudo_h_map.insert(r, h);
    }

    // Benchmark energy values from docs/rust-port-handoff.md:
    // (donor, acceptor, expected_energy)
    let benchmark_cases: Vec<(usize, usize, f64)> = vec![
        (3, 2, -3.8074),
        (4, 3, -3.7912),
        (5, 4, -3.7696),
        (6, 5, -3.7979),
        (7, 6, -3.7884),
        (8, 7, -3.7836),
    ];

    for (d, a, expected_e) in benchmark_cases {
        let (n_d, _, _, _) = get_res_atoms(&d.to_string());
        let h_d = pseudo_h_map.get(&d).unwrap();
        let (_, _, c_a, o_a) = get_res_atoms(&a.to_string());

        let e = calc_kabsch_sander_energy(&n_d, h_d, &c_a, &o_a);
        assert!(
            (e - expected_e).abs() < 1e-3,
            "Energy mismatch for ({}, {}): expected {}, got {}",
            d,
            a,
            expected_e,
            e
        );
    }
}

#[test]
fn test_adjacent_exclusion_and_first_residue() {
    let path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&path, None).expect("failed to load 1hls.pdb");
    let ag = pdb.get_atomgroup(None, None).expect("failed to get ag");

    let model_1 = ag.get_group("model_1").expect("model_1 not found");
    let chain_a = model_1.get_group("A").expect("chain A not found");

    let hbonds = calc_backbone_hbonds(chain_a);

    // 1. First residue (1) cannot be a donor because pseudo-H cannot be computed
    assert!(
        !hbonds.iter().any(|hb| hb.donor_residue_key == "1"),
        "Residue 1 must not be a donor"
    );

    // 2. Adjacent pairs (|d - a| <= 2) must be excluded even if energy < -0.5
    // Specifically, check the benchmark pairs with strong energy (e.g. 3->2, 4->3, etc.)
    let excluded_pairs = [
        ("2", "1"),
        ("3", "2"),
        ("4", "3"),
        ("5", "4"),
        ("6", "5"),
        ("7", "6"),
        ("8", "7"),
        ("3", "1"),
        ("4", "2"),
        ("5", "3"),
    ];
    for (d, a) in excluded_pairs {
        assert!(
            !hbonds
                .iter()
                .any(|hb| hb.donor_residue_key == d && hb.acceptor_residue_key == a),
            "Pair ({} -> {}) has |d - a| <= 2 and must be excluded from hbonds",
            d,
            a
        );
    }

    // Every detected hydrogen bond must satisfy |d - a| > 2
    for hb in &hbonds {
        let d_num: isize = hb.donor_residue_key.parse().unwrap();
        let a_num: isize = hb.acceptor_residue_key.parse().unwrap();
        assert!(
            (d_num - a_num).abs() > 2,
            "Detected hbond ({} -> {}) violates |d - a| > 2 constraint",
            hb.donor_residue_key,
            hb.acceptor_residue_key
        );
        assert!(
            hb.energy < -0.5,
            "Detected hbond must have energy < -0.5, got {}",
            hb.energy
        );
    }
}

#[test]
fn test_calc_backbone_hbonds_1hls() {
    let path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&path, None).expect("failed to load 1hls.pdb");
    let ag = pdb.get_atomgroup(None, None).expect("failed to get ag");

    let model_1 = ag.get_group("model_1").expect("model_1 not found");
    let chain_a = model_1.get_group("A").expect("chain A not found");

    let hbonds = calc_backbone_hbonds(chain_a);
    assert_eq!(hbonds.len(), 9);

    // Expected 9 hbonds in chain A with their expected energies
    let expected: Vec<(&str, &str, f64)> = vec![
        ("6", "2", -1.6882),
        ("7", "3", -0.6522),
        ("8", "3", -1.6268),
        ("15", "12", -0.5502),
        ("16", "12", -0.5110),
        ("17", "14", -0.7758),
        ("18", "14", -0.5546),
        ("19", "16", -0.8396),
        ("20", "17", -1.4639),
    ];

    for (exp_d, exp_a, exp_e) in expected {
        let found = hbonds
            .iter()
            .find(|hb| hb.donor_residue_key == exp_d && hb.acceptor_residue_key == exp_a);
        assert!(
            found.is_some(),
            "Expected hbond ({} -> {}) not found",
            exp_d,
            exp_a
        );
        let hb = found.unwrap();
        assert!(
            (hb.energy - exp_e).abs() < 1e-3,
            "Energy mismatch for ({} -> {}): expected {}, got {}",
            exp_d,
            exp_a,
            exp_e,
            hb.energy
        );
    }
}

#[test]
fn test_unsorted_chain_keys_hbonds() {
    let path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&path, None).expect("failed to load 1hls.pdb");
    let ag = pdb.get_atomgroup(None, None).expect("failed to get ag");

    let model_1 = ag.get_group("model_1").expect("model_1 not found");
    let chain_a = model_1.get_group("A").expect("chain A not found");

    let hbonds_sorted = calc_backbone_hbonds(chain_a);

    // Build scrambled chain with deliberately shuffled insertion order
    let mut chain_scrambled = AtomGroup::with_name("A");
    let mut keys = chain_a.get_group_list();
    // Reverse or scramble order
    keys.reverse();
    for k in keys {
        let res = chain_a.get_group(&k).unwrap();
        chain_scrambled.set_group(&k, res.clone());
    }

    let hbonds_scrambled = calc_backbone_hbonds(&chain_scrambled);
    assert_eq!(hbonds_sorted.len(), hbonds_scrambled.len());
    assert_eq!(hbonds_sorted, hbonds_scrambled);
}

#[test]
fn test_missing_atoms_safe_skip() {
    let mut chain = AtomGroup::with_name("A");

    // Residue 1: complete
    let mut res1 = AtomGroup::with_name("1");
    res1.set_atom("N", make_atom("N", "N", Position::new(0.0, 0.0, 0.0)));
    res1.set_atom("CA", make_atom("CA", "C", Position::new(1.0, 0.0, 0.0)));
    res1.set_atom("C", make_atom("C", "C", Position::new(1.5, 1.0, 0.0)));
    res1.set_atom("O", make_atom("O", "O", Position::new(1.5, 2.0, 0.0)));
    chain.set_group("1", res1);

    // Residue 2: missing O
    let mut res2 = AtomGroup::with_name("2");
    res2.set_atom("N", make_atom("N", "N", Position::new(2.5, 1.0, 0.0)));
    res2.set_atom("CA", make_atom("CA", "C", Position::new(3.0, 2.0, 0.0)));
    res2.set_atom("C", make_atom("C", "C", Position::new(4.0, 2.0, 0.0)));
    // O is missing
    chain.set_group("2", res2);

    // Residue 3: complete
    let mut res3 = AtomGroup::with_name("3");
    res3.set_atom("N", make_atom("N", "N", Position::new(5.0, 2.0, 0.0)));
    res3.set_atom("CA", make_atom("CA", "C", Position::new(6.0, 2.0, 0.0)));
    res3.set_atom("C", make_atom("C", "C", Position::new(6.5, 3.0, 0.0)));
    res3.set_atom("O", make_atom("O", "O", Position::new(6.5, 4.0, 0.0)));
    chain.set_group("3", res3);

    // Residue 4: complete
    let mut res4 = AtomGroup::with_name("4");
    res4.set_atom("N", make_atom("N", "N", Position::new(7.5, 3.0, 0.0)));
    res4.set_atom("CA", make_atom("CA", "C", Position::new(8.5, 3.0, 0.0)));
    res4.set_atom("C", make_atom("C", "C", Position::new(9.0, 4.0, 0.0)));
    res4.set_atom("O", make_atom("O", "O", Position::new(9.0, 5.0, 0.0)));
    chain.set_group("4", res4);

    // Should not panic, skips residue 2 safely
    let hbonds = calc_backbone_hbonds(&chain);
    assert!(!hbonds
        .iter()
        .any(|hb| hb.donor_residue_key == "2" || hb.acceptor_residue_key == "2"));
}
