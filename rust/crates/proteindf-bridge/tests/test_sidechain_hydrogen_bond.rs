// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::path::PathBuf;

use proteindf_bridge::atom::Atom;
use proteindf_bridge::atom_group::AtomGroup;
use proteindf_bridge::format::pdb::Pdb;
use proteindf_bridge::hydrogen_bond::{calc_sidechain_hbonds, calc_sidechain_hbonds_with_options};
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
fn test_sidechain_hbonds_1hls_benchmark() {
    let path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&path, None).expect("failed to load 1hls.pdb");
    let ag = pdb.get_atomgroup(None, None).expect("failed to get ag");

    let model_1 = ag.get_group("model_1").expect("model_1 not found");
    let hbonds = calc_sidechain_hbonds(model_1);

    // Benchmark cases from docs/rust-port-handoff.md Phase 9:
    // (donor_subpath, donor_atom, acceptor_subpath, acceptor_atom, expected_distance)
    let benchmark_cases = [
        ("A/8", "OG1", "A/4", "O", 2.5853),
        ("B/16", "NE2", "B/13", "O", 2.6420),
        ("A/2", "N", "A/19", "OH", 2.7983),
        ("B/10", "ND1", "B/7", "O", 2.9566),
        ("B/5", "NE2", "A/8", "O", 3.1916),
    ];

    for (d_sub, d_atom, a_sub, a_atom, exp_dist) in benchmark_cases {
        let found = hbonds.iter().find(|hb| {
            hb.donor_path.trim_end_matches('/').ends_with(d_sub)
                && hb.donor_atom == d_atom
                && hb.acceptor_path.trim_end_matches('/').ends_with(a_sub)
                && hb.acceptor_atom == a_atom
        });

        assert!(
            found.is_some(),
            "Expected sidechain hbond {}.{} -> {}.{} not found in results",
            d_sub,
            d_atom,
            a_sub,
            a_atom
        );

        let hb = found.unwrap();
        assert!(
            (hb.distance - exp_dist).abs() < 1e-3,
            "Distance mismatch for {}.{} -> {}.{}: expected {}, got {}",
            d_sub,
            d_atom,
            a_sub,
            a_atom,
            exp_dist,
            hb.distance
        );
        // Heavy-atom-only mode (no explicit H in 1hls.pdb) -> angle must be None
        assert!(
            hb.angle.is_none(),
            "Expected angle to be None in heavy-atom-only mode, got {:?}",
            hb.angle
        );
    }
}

#[test]
fn test_sidechain_hbonds_explicit_hydrogen_angle() {
    // Test explicit hydrogen mode with synthetic data
    // Acceptor residue 1 (GLU): O at (0.0, 0.0, 0.0)
    // Donor residue 2 (SER): OG at (2.5, 0.0, 0.0)

    // Case 1: Ideal hydrogen angle (180 deg > 120 deg)
    // H placed between OG and O at (1.5, 0.0, 0.0) -> angle = 180°
    {
        let mut model = AtomGroup::with_name("model");
        let mut chain = AtomGroup::with_name("A");

        let mut res1 = AtomGroup::with_name("GLU");
        res1.set_atom("O", make_atom("O", "O", Position::new(0.0, 0.0, 0.0)));
        chain.set_group("1", res1);

        let mut res2 = AtomGroup::with_name("SER");
        res2.set_atom("OG", make_atom("OG", "O", Position::new(2.5, 0.0, 0.0)));
        res2.set_atom("HG", make_atom("HG", "H", Position::new(1.5, 0.0, 0.0)));
        chain.set_group("2", res2);

        model.set_group("A", chain);

        // In explicit-H mode, should detect 1 hbond when angle > 120°
        let hbonds = calc_sidechain_hbonds_with_options(&model, true);
        assert_eq!(hbonds.len(), 1, "Should detect 1 hbond when angle > 120°");
        let hb = &hbonds[0];
        assert_eq!(hb.donor_atom, "OG");
        assert_eq!(hb.acceptor_atom, "O");
        assert!((hb.distance - 2.5).abs() < 1e-6);
        assert!(hb.angle.is_some());
        assert!((hb.angle.unwrap() - 180.0).abs() < 1e-3);
    }

    // Case 2: Unfavorable hydrogen angle (~68 deg <= 120 deg)
    // H placed at (2.5, 1.0, 0.0) -> angle ~68.2°
    {
        let mut model = AtomGroup::with_name("model");
        let mut chain = AtomGroup::with_name("A");

        let mut res1 = AtomGroup::with_name("GLU");
        res1.set_atom("O", make_atom("O", "O", Position::new(0.0, 0.0, 0.0)));
        chain.set_group("1", res1);

        let mut res2 = AtomGroup::with_name("SER");
        res2.set_atom("OG", make_atom("OG", "O", Position::new(2.5, 0.0, 0.0)));
        res2.set_atom("HG", make_atom("HG", "H", Position::new(2.5, 1.0, 0.0)));
        chain.set_group("2", res2);

        model.set_group("A", chain);

        // In explicit-H mode, angle <= 120° must be excluded
        let hbonds_explicit = calc_sidechain_hbonds_with_options(&model, true);
        assert_eq!(
            hbonds_explicit.len(),
            0,
            "Should exclude hbond when explicit H angle <= 120°"
        );

        // In heavy-atom-only mode, distance < 3.5 Å is sufficient so it is detected
        let hbonds_heavy = calc_sidechain_hbonds(&model);
        assert_eq!(
            hbonds_heavy.len(),
            1,
            "Heavy-atom-only mode should detect pair purely by distance < 3.5 Å"
        );
    }
}

#[test]
fn test_same_residue_exclusion() {
    // Residue containing both a donor and an acceptor (e.g. SER with backbone N and sidechain OG)
    // Even if distance < 3.5 Å, it must be excluded.
    let mut model = AtomGroup::with_name("model");
    let mut chain = AtomGroup::with_name("A");

    let mut res1 = AtomGroup::with_name("SER");
    res1.set_atom("N", make_atom("N", "N", Position::new(0.0, 0.0, 0.0)));
    res1.set_atom("OG", make_atom("OG", "O", Position::new(2.0, 0.0, 0.0)));
    chain.set_group("1", res1);

    model.set_group("A", chain);

    let hbonds = calc_sidechain_hbonds(&model);
    assert!(
        hbonds.is_empty(),
        "Intra-residue donor-acceptor pairs must be excluded"
    );
}
