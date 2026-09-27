// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::path::PathBuf;

use proteindf_bridge::atom::Atom;
use proteindf_bridge::atom_group::AtomGroup;
use proteindf_bridge::ccd_templates::CcdTemplateDb;
use proteindf_bridge::format::Pdb;
use proteindf_bridge::position::Position;

fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

/// Recursively removes all hydrogens from an `AtomGroup`.
fn strip_all_hydrogens(group: &mut AtomGroup) {
    let h_keys: Vec<String> = group
        .atoms()
        .filter(|(_, a)| a.atomic_number() == 1 || a.symbol() == Ok("H"))
        .map(|(k, _)| k.clone())
        .collect();
    for k in h_keys {
        group.remove_atom(&k);
    }
    for (_, sub) in group.groups_mut() {
        strip_all_hydrogens(sub);
    }
}

/// Counts all hydrogen atoms in an `AtomGroup` recursively.
fn count_all_hydrogens(group: &AtomGroup) -> usize {
    let direct_h = group
        .atoms()
        .filter(|(_, a)| a.atomic_number() == 1 || a.symbol() == Ok("H"))
        .count();
    let sub_h: usize = group
        .groups()
        .map(|(_, sub)| count_all_hydrogens(sub))
        .sum();
    direct_h + sub_h
}

#[test]
fn test_hydrogenate_1hls_full_structure() {
    let pdb_path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&pdb_path, None).expect("Failed to read 1hls.pdb");
    let mut structure = pdb
        .get_atomgroup(None, None)
        .expect("Failed to convert to AtomGroup");
    let db = CcdTemplateDb::global();

    // 1hls in tests/data is an X-ray structure without hydrogens (or with initial hydrogens)
    let initial_h_count = count_all_hydrogens(&structure);

    // Strip all hydrogens to ensure uniform starting point
    if initial_h_count > 0 {
        strip_all_hydrogens(&mut structure);
    }
    assert_eq!(count_all_hydrogens(&structure), 0);

    // Run hydrogenation orchestrator
    let report = structure
        .add_missing_hydrogens(db)
        .expect("Hydrogenation failed on 1hls");

    assert_eq!(report.skipped_residues.len(), 0);
    assert_eq!(report.total_removed_hydrogens, 0);
    assert!(report.total_added_hydrogens > 0);
    assert_eq!(report.hydrogenated_residues, 51); // 21 (chain A) + 30 (chain B)

    // Verify chain A N-terminus (GIVEQCCTSICSLYQLENYCN, 21 res)
    let model = structure
        .get_group("model_1")
        .or_else(|| structure.get_group("1"))
        .expect("Model group must exist");

    let chain_a = model.get_group("A").expect("Chain A must exist");
    // Residue 1 (GLY) is N-terminus -> should have H1, H2, H3
    let res_a1 = chain_a.get_group("1").expect("Residue A1 must exist");
    assert!(res_a1.has_atom("H1"));
    assert!(res_a1.has_atom("H2"));
    assert!(res_a1.has_atom("H3"));
    assert!(!res_a1.has_atom("H"));

    // Residue 2 (ILE) is internal -> should have H and HA
    let res_a2 = chain_a.get_group("2").expect("Residue A2 must exist");
    assert!(res_a2.has_atom("H"));
    assert!(res_a2.has_atom("HA"));
    assert!(!res_a2.has_atom("H1"));

    // Verify chain B (30 res)
    let chain_b = model.get_group("B").expect("Chain B must exist");
    // Residue 1 (PHE) is N-terminus -> should have H1, H2, H3
    let res_b1 = chain_b.get_group("1").expect("Residue B1 must exist");
    assert!(res_b1.has_atom("H1"));
    assert!(res_b1.has_atom("H2"));
    assert!(res_b1.has_atom("H3"));

    // Residue 28 (PRO) is internal -> should NOT have H
    let res_b28 = chain_b.get_group("28").expect("Residue B28 must exist");
    assert_eq!(res_b28.name, "PRO");
    assert!(!res_b28.has_atom("H"));
    assert!(!res_b28.has_atom("H1"));
    assert!(res_b28.has_atom("HA"));
    assert!(res_b28.has_atom("HB2") || res_b28.has_atom("HB3"));

    // Total hydrogens added should match report.total_added_hydrogens
    // 1hls has 51 residues (21 in chain A, 30 in chain B, 1 Proline B28).
    // Standard total hydrogens = 394:
    // - 2 N-termini (GLY A1: H1,H2,H3, PHE B1: H1,H2,H3) = 6 backbone NH
    // - 48 internal non-Proline residues (H) = 48 backbone NH (B28 is PRO: 0 backbone H)
    // - Total backbone NH = 54
    // - Sidechain & CA hydrogens across 51 residues = 340
    // Total = 394 hydrogens.
    let final_h_count = count_all_hydrogens(&structure);
    assert_eq!(final_h_count, report.total_added_hydrogens);
    assert_eq!(final_h_count, 394);
}

#[test]
fn test_hydrogenate_partial_existing_hydrogens() {
    let pdb_path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&pdb_path, None).expect("Failed to read 1hls.pdb");
    let mut structure = pdb
        .get_atomgroup(None, None)
        .expect("Failed to convert to AtomGroup");
    let db = CcdTemplateDb::global();

    // 1hls.pdb has 379 hydrogens natively (lacking N-terminal NH3+ and some terminal H)
    let initial_h = count_all_hydrogens(&structure);
    assert_eq!(initial_h, 379);

    // Run hydrogenation on structure that already has partial hydrogens
    let report = structure
        .add_missing_hydrogens(db)
        .expect("Hydrogenation on partially hydrogenated structure should succeed");

    // Missing hydrogens (394 - 379 = 15) should be added
    assert_eq!(report.total_added_hydrogens, 15);
    let final_h = count_all_hydrogens(&structure);
    assert_eq!(final_h, 394);
}

#[test]
fn test_skip_unknown_residues() {
    let db = CcdTemplateDb::global();

    let mut chain = AtomGroup::with_name("A");
    chain.set_path("/model_1/A/".to_string());

    // Residue 1: ALA (known)
    let mut res_ala = AtomGroup::with_name("ALA");
    res_ala.set_atom(
        "N",
        Atom::new_with_pos("N", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    res_ala.set_atom(
        "CA",
        Atom::new_with_pos("C", Position::new(1.46, 0.0, 0.0)).unwrap(),
    );
    res_ala.set_atom(
        "C",
        Atom::new_with_pos("C", Position::new(2.0, 1.4, 0.0)).unwrap(),
    );
    res_ala.set_atom(
        "O",
        Atom::new_with_pos("O", Position::new(3.2, 1.5, 0.0)).unwrap(),
    );
    res_ala.set_atom(
        "CB",
        Atom::new_with_pos("C", Position::new(2.0, -0.7, 1.2)).unwrap(),
    );
    chain.set_group("1", res_ala);

    // Residue 2: UNK (unknown)
    let mut res_unk = AtomGroup::with_name("UNK");
    res_unk.set_atom(
        "N",
        Atom::new_with_pos("N", Position::new(1.3, 2.5, 0.0)).unwrap(),
    );
    res_unk.set_atom(
        "CA",
        Atom::new_with_pos("C", Position::new(1.8, 3.8, 0.0)).unwrap(),
    );
    res_unk.set_atom(
        "C",
        Atom::new_with_pos("C", Position::new(3.3, 3.9, 0.0)).unwrap(),
    );
    chain.set_group("2", res_unk);

    // Residue 3: GLY (known)
    let mut res_gly = AtomGroup::with_name("GLY");
    res_gly.set_atom(
        "N",
        Atom::new_with_pos("N", Position::new(4.0, 5.0, 0.0)).unwrap(),
    );
    res_gly.set_atom(
        "CA",
        Atom::new_with_pos("C", Position::new(5.4, 5.2, 0.0)).unwrap(),
    );
    res_gly.set_atom(
        "C",
        Atom::new_with_pos("C", Position::new(6.0, 6.6, 0.0)).unwrap(),
    );
    res_gly.set_atom(
        "O",
        Atom::new_with_pos("O", Position::new(7.2, 6.7, 0.0)).unwrap(),
    );
    chain.set_group("3", res_gly);

    let report = chain
        .add_missing_hydrogens(db)
        .expect("Hydrogenation should proceed even with unknown residue");

    // UNK should be recorded in skipped_residues
    assert_eq!(report.skipped_residues.len(), 1);
    let (skipped_path, reason) = &report.skipped_residues[0];
    assert!(skipped_path.contains("2"));
    assert!(reason.contains("UNK"));

    // Residue 1 (ALA) and Residue 3 (GLY) should still be hydrogenated
    let res1 = chain.get_group("1").unwrap();
    assert!(res1.has_atom("H1")); // N-terminus
    assert!(res1.has_atom("HA")); // Sidechain

    let res3 = chain.get_group("3").unwrap();
    assert!(res3.has_atom("H")); // Internal backbone
    assert!(res3.has_atom("HA2") || res3.has_atom("HA3")); // GLY alpha hydrogens
}

#[test]
fn test_chain_break_handling() {
    let db = CcdTemplateDb::global();

    let mut chain = AtomGroup::with_name("A");
    chain.set_path("/model_1/A/".to_string());

    // Residue 1: ALA
    let mut res1 = AtomGroup::with_name("ALA");
    res1.set_atom(
        "N",
        Atom::new_with_pos("N", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    res1.set_atom(
        "CA",
        Atom::new_with_pos("C", Position::new(1.46, 0.0, 0.0)).unwrap(),
    );
    res1.set_atom(
        "C",
        Atom::new_with_pos("C", Position::new(2.0, 1.4, 0.0)).unwrap(),
    );
    res1.set_atom(
        "O",
        Atom::new_with_pos("O", Position::new(3.2, 1.5, 0.0)).unwrap(),
    );
    res1.set_atom(
        "CB",
        Atom::new_with_pos("C", Position::new(2.0, -0.7, 1.2)).unwrap(),
    );
    chain.set_group("1", res1);

    // Residue 2: ALA with a gap (distance between res1.C and res2.N > 2.5 A, e.g. 10.0 A away)
    let mut res2 = AtomGroup::with_name("ALA");
    res2.set_atom(
        "N",
        Atom::new_with_pos("N", Position::new(12.0, 1.4, 0.0)).unwrap(),
    );
    res2.set_atom(
        "CA",
        Atom::new_with_pos("C", Position::new(13.46, 1.4, 0.0)).unwrap(),
    );
    res2.set_atom(
        "C",
        Atom::new_with_pos("C", Position::new(14.0, 2.8, 0.0)).unwrap(),
    );
    res2.set_atom(
        "O",
        Atom::new_with_pos("O", Position::new(15.2, 2.9, 0.0)).unwrap(),
    );
    res2.set_atom(
        "CB",
        Atom::new_with_pos("C", Position::new(14.0, 0.7, 1.2)).unwrap(),
    );
    chain.set_group("2", res2);

    let report = chain
        .add_missing_hydrogens(db)
        .expect("Hydrogenation across chain break should succeed");

    assert_eq!(report.skipped_residues.len(), 0);

    // Both residues must be treated as N-terminal fragments because of the chain break
    let r1 = chain.get_group("1").unwrap();
    assert!(r1.has_atom("H1"));
    assert!(r1.has_atom("H2"));
    assert!(r1.has_atom("H3"));
    assert!(!r1.has_atom("H"));

    let r2 = chain.get_group("2").unwrap();
    assert!(r2.has_atom("H1"));
    assert!(r2.has_atom("H2"));
    assert!(r2.has_atom("H3"));
    assert!(!r2.has_atom("H"));
}
