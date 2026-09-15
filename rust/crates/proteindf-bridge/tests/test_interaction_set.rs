// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::path::PathBuf;

use proteindf_bridge::format::Pdb;
use proteindf_bridge::interaction_set::{Interaction, InteractionKind, InteractionSet};

fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

#[test]
fn test_interaction_set_detect_all_1hls() {
    let path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&path, None).expect("failed to load 1hls.pdb");
    let ag = pdb.get_atomgroup(None, None).expect("failed to get ag");
    let model_1 = ag.get_group("model_1").expect("model_1 not found");

    let set = InteractionSet::detect_all(model_1, None, None);

    // 1. Disulfide bonds: 3 bonds in insulin (A:6-A:11, A:7-B:7, A:20-B:19)
    let disulfides = set.filter_by_kind(InteractionKind::Disulfide);
    assert_eq!(disulfides.len(), 3);
    assert_eq!(set.count_by_kind(InteractionKind::Disulfide), 3);
    for d in &disulfides {
        assert_eq!(d.atoms.len(), 2);
        assert!(d.distance.is_some());
        let dist = d.distance.unwrap();
        // SG-SG distance should be around 2.0 - 2.1 Å
        assert!(dist > 1.9 && dist < 2.3);
    }

    // 2. Salt bridges: 0 in 1hls
    assert_eq!(set.count_by_kind(InteractionKind::SaltBridge), 0);

    // 3. Hydrogen bonds:
    // Backbone: chain A has 9 hbonds, chain B has 13 hbonds (total 22 backbone hbonds)
    // Sidechain: 6 hbonds detected (including the 5 benchmark cases + A:9 SER.N -> A:8 THR.OG1)
    let hbonds = set.filter_by_kind(InteractionKind::HydrogenBond);
    let bb_count = hbonds
        .iter()
        .filter(|h| {
            h.donor_acceptor_role
                .as_ref()
                .is_some_and(|r| r.starts_with("backbone"))
        })
        .count();
    let sc_count = hbonds
        .iter()
        .filter(|h| {
            h.donor_acceptor_role
                .as_ref()
                .is_some_and(|r| r == "sidechain")
        })
        .count();

    assert_eq!(bb_count, 22);
    assert_eq!(sc_count, 6);
    assert_eq!(hbonds.len(), 28);
    assert_eq!(set.count_by_kind(InteractionKind::HydrogenBond), 28);

    // Verify sidechain benchmark cases are present in InteractionSet
    let sc_benchmark_cases = [
        ("A/8/OG1", "A/4/O"),
        ("B/16/NE2", "B/13/O"),
        ("A/2/N", "A/19/OH"),
        ("B/10/ND1", "B/7/O"),
        ("B/5/NE2", "A/8/O"),
    ];
    for (d_suffix, a_suffix) in sc_benchmark_cases {
        let found = hbonds.iter().find(|h| {
            h.atoms.len() == 2 && h.atoms[0].ends_with(d_suffix) && h.atoms[1].ends_with(a_suffix)
        });
        assert!(
            found.is_some(),
            "Expected sidechain hbond ({} -> {}) not found in InteractionSet",
            d_suffix,
            a_suffix
        );
    }

    // 4. CH-pi interactions: 7 detected with default thresholds (4.5 Å, 40.0°)
    let ch_pi = set.filter_by_kind(InteractionKind::ChPi);
    assert_eq!(ch_pi.len(), 7);
    assert_eq!(set.count_by_kind(InteractionKind::ChPi), 7);

    let ch_pi_cases = [
        ("A/10/CG1", "B/5/"),
        ("B/15/CB", "B/24/"),
        ("B/17/CA", "B/16/"),
    ];
    for (c_suffix, ring_suffix) in ch_pi_cases {
        let found = ch_pi.iter().find(|ch| {
            ch.atoms.len() == 2
                && ch.atoms[0].ends_with(c_suffix)
                && ch.atoms[1].ends_with(ring_suffix)
        });
        assert!(
            found.is_some(),
            "Expected CH-pi ({} -> {}) not found in InteractionSet",
            c_suffix,
            ring_suffix
        );
    }

    // Total interactions: 3 (SS) + 28 (HB) + 7 (CH-pi) = 38
    assert_eq!(set.len(), 3 + 28 + 7);
}

#[test]
fn test_interaction_set_ch_pi_tuning() {
    let path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&path, None).expect("failed to load 1hls.pdb");
    let ag = pdb.get_atomgroup(None, None).expect("failed to get ag");
    let model_1 = ag.get_group("model_1").expect("model_1 not found");

    // Default: 7
    let set_default = InteractionSet::detect_all(model_1, None, None);
    assert_eq!(set_default.count_by_kind(InteractionKind::ChPi), 7);

    // Relaxed thresholds (5.0 Å, 50.0°): should detect additional candidates
    let set_relaxed = InteractionSet::detect_all(model_1, Some(5.0), Some(50.0));
    assert!(set_relaxed.count_by_kind(InteractionKind::ChPi) > 7);

    // Strict thresholds (3.5 Å, 20.0°): fewer candidates
    let set_strict = InteractionSet::detect_all(model_1, Some(3.5), Some(20.0));
    assert!(set_strict.count_by_kind(InteractionKind::ChPi) < 7);
}

#[test]
fn test_interaction_set_msgpack_roundtrip() {
    let path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&path, None).expect("failed to load 1hls.pdb");
    let ag = pdb.get_atomgroup(None, None).expect("failed to get ag");
    let model_1 = ag.get_group("model_1").expect("model_1 not found");

    let original = InteractionSet::detect_all(model_1, None, None);
    assert!(!original.is_empty());

    // In-memory roundtrip
    let bytes = original.to_msgpack().expect("to_msgpack failed");
    let restored = InteractionSet::from_msgpack(&bytes).expect("from_msgpack failed");
    assert_eq!(original, restored);

    // File roundtrip
    let tmp_dir = std::env::temp_dir();
    let file_path = tmp_dir.join("interactions_test_roundtrip.msgpack");
    original
        .save_msgpack(&file_path)
        .expect("save_msgpack failed");
    let loaded = InteractionSet::load_msgpack(&file_path).expect("load_msgpack failed");
    let _ = std::fs::remove_file(&file_path);
    assert_eq!(original, loaded);
}

#[test]
fn test_interaction_set_yaml_roundtrip() {
    let path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&path, None).expect("failed to load 1hls.pdb");
    let ag = pdb.get_atomgroup(None, None).expect("failed to get ag");
    let model_1 = ag.get_group("model_1").expect("model_1 not found");

    let original = InteractionSet::detect_all(model_1, None, None);
    assert!(!original.is_empty());

    // In-memory roundtrip
    let yaml_str = original.to_yaml().expect("to_yaml failed");
    assert!(yaml_str.contains("disulfide"));
    assert!(yaml_str.contains("hydrogen_bond"));
    assert!(yaml_str.contains("ch_pi"));

    let restored = InteractionSet::from_yaml(&yaml_str).expect("from_yaml failed");
    assert_eq!(original, restored);

    // File roundtrip
    let tmp_dir = std::env::temp_dir();
    let file_path = tmp_dir.join("interactions_test_roundtrip.yaml");
    original.save_yaml(&file_path).expect("save_yaml failed");
    let loaded = InteractionSet::load_yaml(&file_path).expect("load_yaml failed");
    let _ = std::fs::remove_file(&file_path);
    assert_eq!(original, loaded);
}

#[test]
fn test_interaction_set_salt_bridge_inclusion() {
    // Verify that SaltBridge interactions are correctly serialized and filtered
    let mut set = InteractionSet::new();
    set.push(Interaction::with_metrics(
        InteractionKind::SaltBridge,
        vec!["/model_1/A/4/".to_string(), "/model_1/A/29/".to_string()],
        Some(3.2),
        None,
        Some("GLU - LYS".to_string()),
    ));

    assert_eq!(set.len(), 1);
    assert_eq!(set.count_by_kind(InteractionKind::SaltBridge), 1);

    // YAML check
    let yaml = set.to_yaml().unwrap();
    assert!(yaml.contains("salt_bridge"));
    assert!(yaml.contains("GLU - LYS"));

    let deserialized = InteractionSet::from_yaml(&yaml).unwrap();
    assert_eq!(set, deserialized);
}
