// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::path::PathBuf;

use proteindf_bridge::atom::Atom;
use proteindf_bridge::atom_group::AtomGroup;
use proteindf_bridge::ch_pi::{
    calc_ch_pi_interactions, calc_ch_pi_interactions_with_thresholds, calc_ring_geometry,
    AROMATIC_RINGS,
};
use proteindf_bridge::format::pdb::Pdb;
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
fn test_ring_centroids_benchmark() {
    let path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&path, None).expect("failed to load 1hls.pdb");
    let ag = pdb.get_atomgroup(None, None).expect("failed to get ag");

    let model_1 = ag.get_group("model_1").expect("model_1 not found");
    let chain_b = model_1.get_group("B").expect("chain B not found");

    // Benchmark ring centroids from docs/rust-port-handoff.md Phase 9:
    let benchmark_rings = [
        ("5", "HIS", Position::new(11.6482, 3.7998, 0.8784)),
        ("16", "HIS", Position::new(-2.8782, -4.1082, -7.6430)),
        ("24", "PHE", Position::new(-6.9067, -0.9103, -3.5685)),
        ("26", "TYR", Position::new(-3.3625, 4.6848, -0.3693)),
    ];

    for (res_key, res_name, exp_center) in benchmark_rings {
        let res = chain_b.get_group(res_key).expect("residue not found");
        assert_eq!(res.name.trim(), res_name);

        let def = AROMATIC_RINGS
            .iter()
            .find(|d| d.residue_name == res_name)
            .unwrap();

        let (center, normal) =
            calc_ring_geometry(res, def.atom_names).expect("failed to calc ring geometry");

        assert!(
            (center.x - exp_center.x).abs() < 1e-3,
            "Center X mismatch for {}: expected {}, got {}",
            res_key,
            exp_center.x,
            center.x
        );
        assert!(
            (center.y - exp_center.y).abs() < 1e-3,
            "Center Y mismatch for {}: expected {}, got {}",
            res_key,
            exp_center.y,
            center.y
        );
        assert!(
            (center.z - exp_center.z).abs() < 1e-3,
            "Center Z mismatch for {}: expected {}, got {}",
            res_key,
            exp_center.z,
            center.z
        );

        // Normal vector must be normalized
        assert!((normal.length() - 1.0).abs() < 1e-6);
    }
}

#[test]
fn test_ch_pi_candidate_pairs_benchmark() {
    let path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&path, None).expect("failed to load 1hls.pdb");
    let ag = pdb.get_atomgroup(None, None).expect("failed to get ag");

    let model_1 = ag.get_group("model_1").expect("model_1 not found");
    let interactions = calc_ch_pi_interactions(model_1);

    // Benchmark candidate pairs from docs/rust-port-handoff.md Phase 9:
    // (carbon_subpath, carbon_atom, ring_subpath, exp_dist, exp_angle, should_detect)
    let candidates = [
        ("A/10", "CG1", "B/5", 3.400, 23.79, true),
        ("B/15", "CB", "B/24", 3.793, 23.69, true),
        ("B/17", "CA", "B/16", 4.059, 15.62, true),
        ("A/10", "CD1", "B/5", 4.074, 41.11, false), // angle > 40 deg
        ("B/27", "C", "B/26", 4.578, 11.75, false),  // distance > 4.5 A
    ];

    for (c_sub, c_atom, r_sub, exp_dist, exp_angle, should_detect) in candidates {
        let found = interactions.iter().find(|i| {
            i.carbon_path.trim_end_matches('/').ends_with(c_sub)
                && i.carbon_atom == c_atom
                && i.ring_path.trim_end_matches('/').ends_with(r_sub)
        });

        if should_detect {
            assert!(
                found.is_some(),
                "Expected CH-pi interaction {}.{} -> {} not detected",
                c_sub,
                c_atom,
                r_sub
            );
            let item = found.unwrap();
            assert!(
                (item.distance - exp_dist).abs() < 1e-3,
                "Distance mismatch for {}.{} -> {}: expected {}, got {}",
                c_sub,
                c_atom,
                r_sub,
                exp_dist,
                item.distance
            );
            assert!(
                (item.angle - exp_angle).abs() < 1e-2,
                "Angle mismatch for {}.{} -> {}: expected {}, got {}",
                c_sub,
                c_atom,
                r_sub,
                exp_angle,
                item.angle
            );
        } else {
            assert!(
                found.is_none(),
                "Candidate {}.{} -> {} should not be detected at default thresholds (dist={}, angle={})",
                c_sub,
                c_atom,
                r_sub,
                exp_dist,
                exp_angle
            );
        }
    }
}

#[test]
fn test_ch_pi_tunable_thresholds() {
    let path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&path, None).expect("failed to load 1hls.pdb");
    let ag = pdb.get_atomgroup(None, None).expect("failed to get ag");

    let model_1 = ag.get_group("model_1").expect("model_1 not found");

    // 1. Expand angle threshold to 42.0°: A:10.CD1 -> B:5 (angle 41.11°) should now be detected
    let interactions_angle_42 = calc_ch_pi_interactions_with_thresholds(model_1, 4.5, 42.0);
    let found_angle_exceeded = interactions_angle_42.iter().find(|i| {
        i.carbon_path.trim_end_matches('/').ends_with("A/10")
            && i.carbon_atom == "CD1"
            && i.ring_path.trim_end_matches('/').ends_with("B/5")
    });
    assert!(
        found_angle_exceeded.is_some(),
        "A:10.CD1 -> B:5 must be detected when angle threshold is relaxed to 42.0°"
    );

    // 2. Expand distance threshold to 4.6 Å: B:27.C -> B:26 (dist 4.578 Å) should now be detected
    let interactions_dist_46 = calc_ch_pi_interactions_with_thresholds(model_1, 4.6, 40.0);
    let found_dist_exceeded = interactions_dist_46.iter().find(|i| {
        i.carbon_path.trim_end_matches('/').ends_with("B/27")
            && i.carbon_atom == "C"
            && i.ring_path.trim_end_matches('/').ends_with("B/26")
    });
    assert!(
        found_dist_exceeded.is_some(),
        "B:27.C -> B:26 must be detected when distance threshold is relaxed to 4.6 Å"
    );
}

#[test]
fn test_same_residue_exclusion() {
    // Synthetic PHE residue: its own CA, C, CB or ring atoms should not interact with its own ring
    let mut model = AtomGroup::with_name("model");
    let mut chain = AtomGroup::with_name("A");

    let mut phe = AtomGroup::with_name("PHE");
    phe.set_atom("N", make_atom("N", "N", Position::new(0.0, 0.0, 0.0)));
    phe.set_atom("CA", make_atom("CA", "C", Position::new(1.0, 0.0, 0.0)));
    phe.set_atom("C", make_atom("C", "C", Position::new(2.0, 0.0, 0.0)));
    phe.set_atom("O", make_atom("O", "O", Position::new(2.5, 1.0, 0.0)));
    phe.set_atom("CB", make_atom("CB", "C", Position::new(1.0, 1.5, 0.0)));
    phe.set_atom("CG", make_atom("CG", "C", Position::new(1.0, 2.5, 0.0)));
    phe.set_atom("CD1", make_atom("CD1", "C", Position::new(0.0, 3.0, 0.0)));
    phe.set_atom("CD2", make_atom("CD2", "C", Position::new(2.0, 3.0, 0.0)));
    phe.set_atom("CE1", make_atom("CE1", "C", Position::new(0.0, 4.0, 0.0)));
    phe.set_atom("CE2", make_atom("CE2", "C", Position::new(2.0, 4.0, 0.0)));
    phe.set_atom("CZ", make_atom("CZ", "C", Position::new(1.0, 4.5, 0.0)));
    chain.set_group("1", phe);

    model.set_group("A", chain);

    let interactions = calc_ch_pi_interactions(&model);
    assert!(
        interactions.is_empty(),
        "Intra-residue carbons must not form CH-pi interactions with own ring"
    );
}
