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

/// Regression test for Bug 1: Water molecules (HOH) have only 1 heavy atom ('O') in CCD templates,
/// which is less than MIN_SUPERPOSE_HEAVY_ATOMS (3).
/// The orchestrator must not propagate this error and abort the whole structure; instead, it must
/// record HOH into `skipped_residues` and successfully hydrogenate other residues.
#[test]
fn test_hydrogenate_water_hoh_resilient() {
    let db = CcdTemplateDb::global();

    let mut chain = AtomGroup::with_name("A");
    chain.set_path("/model_1/A/".to_string());

    // Residue 1: ALA (standard amino acid)
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

    // Residue 2: HOH (water molecule with single 'O' heavy atom)
    let mut res_hoh = AtomGroup::with_name("HOH");
    res_hoh.set_atom(
        "O",
        Atom::new_with_pos("O", Position::new(10.0, 10.0, 10.0)).unwrap(),
    );
    chain.set_group("2", res_hoh);

    let report = chain
        .add_missing_hydrogens(db)
        .expect("Hydrogenation must not abort when encountering HOH water molecules");

    // ALA should be hydrogenated
    let res1 = chain.get_group("1").unwrap();
    assert!(res1.has_atom("H1"));
    assert!(res1.has_atom("HA"));
    assert_eq!(report.hydrogenated_residues, 1);
    assert!(report.residue_reports.contains_key("/model_1/A/1/"));

    // HOH must be recorded in skipped_residues (not in residue_reports)
    assert_eq!(report.skipped_residues.len(), 1);
    let (skipped_path, reason) = &report.skipped_residues[0];
    assert_eq!(skipped_path, "/model_1/A/2/");
    assert!(
        reason.contains("common heavy atoms") || reason.contains("Sidechain/general hydrogenation"),
        "Reason should explain heavy atom deficiency: {reason}"
    );
    assert!(!report.residue_reports.contains_key("/model_1/A/2/"));
}

/// Regression test for Bug 1 on real fixture: adding crystal waters to 1hls.pdb
/// ensures real protein structures with water molecules complete without error.
#[test]
fn test_1hls_with_crystal_waters() {
    let pdb_path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&pdb_path, None).expect("Failed to read 1hls.pdb");
    let mut structure = pdb
        .get_atomgroup(None, None)
        .expect("Failed to convert to AtomGroup");
    let db = CcdTemplateDb::global();

    // Add simulated HOH water molecule under model_1
    let model = if structure.has_group("model_1") {
        structure.get_group_mut("model_1").unwrap()
    } else {
        structure.get_group_mut("1").unwrap()
    };

    let mut water_chain = AtomGroup::with_name("W");
    water_chain.set_path("/model_1/W/".to_string());
    let mut hoh1 = AtomGroup::with_name("HOH");
    hoh1.set_atom(
        "O",
        Atom::new_with_pos("O", Position::new(50.0, 50.0, 50.0)).unwrap(),
    );
    water_chain.set_group("1", hoh1);
    model.set_group("W", water_chain);

    // Strip all hydrogens
    strip_all_hydrogens(&mut structure);

    let report = structure
        .add_missing_hydrogens(db)
        .expect("Hydrogenation must succeed on 1hls with crystal waters");

    // All 51 protein residues must be hydrogenated (394 hydrogens)
    assert_eq!(report.hydrogenated_residues, 51);
    assert_eq!(report.total_added_hydrogens, 394);

    // The water molecule must be in skipped_residues and not in residue_reports
    assert_eq!(report.skipped_residues.len(), 1);
    assert_eq!(report.skipped_residues[0].0, "/model_1/W/1/");
    assert!(!report.residue_reports.contains_key("/model_1/W/1/"));
}

/// Regression test for Bug 2: Test that a completely unmodified residue (e.g. unknown ligand with no CCD template)
/// is recorded in `skipped_residues` and NEVER appears in `residue_reports`.
#[test]
fn test_unmodified_unknown_ligand_recorded_in_skipped_only() {
    let db = CcdTemplateDb::global();

    let mut chain = AtomGroup::with_name("A");
    chain.set_path("/model_1/A/".to_string());

    // Residue 1: ALA (known amino acid)
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

    // Residue 2: LIG (unknown ligand with no N/CA backbone atoms and no CCD template)
    let mut res_lig = AtomGroup::with_name("LIG");
    res_lig.set_atom(
        "C1",
        Atom::new_with_pos("C", Position::new(20.0, 20.0, 20.0)).unwrap(),
    );
    res_lig.set_atom(
        "C2",
        Atom::new_with_pos("C", Position::new(21.4, 20.0, 20.0)).unwrap(),
    );
    chain.set_group("2", res_lig);

    let report = chain
        .add_missing_hydrogens(db)
        .expect("Hydrogenation should proceed with unknown ligand");

    // LIG was not modified: must be in skipped_residues, must NOT be in residue_reports
    assert_eq!(report.skipped_residues.len(), 1);
    let (skipped_path, reason) = &report.skipped_residues[0];
    assert_eq!(skipped_path, "/model_1/A/2/");
    assert!(reason.contains("No CCD template found for residue 'LIG'"));
    assert!(!report.residue_reports.contains_key("/model_1/A/2/"));

    // ALA was modified: must be in residue_reports, must NOT be in skipped_residues
    assert_eq!(report.hydrogenated_residues, 1);
    assert!(report.residue_reports.contains_key("/model_1/A/1/"));
}

/// Regression test for Bug 2: Test that an unknown amino acid residue with backbone N/CA
/// that has backbone hydrogens added is recorded in `residue_reports` and NEVER in `skipped_residues`.
#[test]
fn test_partially_modified_unknown_residue_recorded_in_reports_only() {
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

    // Residue 2: UNK (unknown amino acid: has N, CA, C so backbone H is added, but no CCD template for sidechain)
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

    let report = chain
        .add_missing_hydrogens(db)
        .expect("Hydrogenation should proceed with unknown residue");

    // UNK had backbone H added (modified): must be in residue_reports, must NOT be in skipped_residues
    let res_unk_after = chain.get_group("2").unwrap();
    assert!(res_unk_after.has_atom("H"));
    assert!(
        report.residue_reports.contains_key("/model_1/A/2/"),
        "Modified UNK residue must be recorded in residue_reports"
    );
    assert!(
        !report
            .skipped_residues
            .iter()
            .any(|(p, _)| p == "/model_1/A/2/"),
        "Modified UNK residue must NEVER appear in skipped_residues (no double-counting)"
    );
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

/// Regression test for Round 2 Bug 1: Running hydrogenation on an already fully hydrogenated structure
/// must NOT record residues as skipped/failed; they must be treated as successful (with 0 additions)
/// in `residue_reports`, and `skipped_residues` must remain empty.
#[test]
fn test_idempotent_hydrogenation_already_complete() {
    let pdb_path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&pdb_path, None).expect("Failed to read 1hls.pdb");
    let mut structure = pdb
        .get_atomgroup(None, None)
        .expect("Failed to convert to AtomGroup");
    let db = CcdTemplateDb::global();

    // 1st run: adds missing hydrogens
    let report1 = structure
        .add_missing_hydrogens(db)
        .expect("1st hydrogenation should succeed");
    assert_eq!(report1.skipped_residues.len(), 0);
    assert_eq!(report1.step_errors.len(), 0);

    // 2nd run: already fully hydrogenated
    let report2 = structure
        .add_missing_hydrogens(db)
        .expect("2nd hydrogenation should succeed");

    // Must have 0 additions and 0 removals
    assert_eq!(report2.total_added_hydrogens, 0);
    assert_eq!(report2.total_removed_hydrogens, 0);

    // Critically: must NOT treat residues as skipped or failed!
    assert_eq!(
        report2.skipped_residues.len(),
        0,
        "Already hydrogenated residues must not be marked as skipped: {:?}",
        report2.skipped_residues
    );
    assert_eq!(report2.step_errors.len(), 0);

    // All 51 residues must be present in residue_reports as completed
    assert_eq!(report2.residue_reports.len(), 51);
}

/// Regression test for Round 2 Bug 2: When Step 1 (backbone) succeeds and modifies a residue,
/// but Step 2 (sidechain/general) fails, the error information must NOT be swallowed;
/// it must be preserved in `step_errors` while the backbone addition is preserved in `residue_reports`.
#[test]
fn test_step2_error_recorded_when_step1_succeeds() {
    let db = CcdTemplateDb::global();

    let mut chain = AtomGroup::with_name("A");
    chain.set_path("/model_1/A/".to_string());

    // Residue 1: ALA (standard preceding residue)
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

    // Residue 2: ALA with only N and CA (missing C, O, CB).
    // - Step 1: has N and CA, distance C_prev(2.0, 1.4, 0.0) to N(2.0, 2.7, 0.0) is 1.3 A (valid peptide bond).
    //   -> Backbone amide 'H' is successfully added!
    // - Step 2: ALA CCD template requires at least MIN_SUPERPOSE_HEAVY_ATOMS (3) heavy atoms, but only 2 (N, CA) exist.
    //   -> Step 2 fails with common_heavy_atoms error.
    let mut res2 = AtomGroup::with_name("ALA");
    res2.set_atom(
        "N",
        Atom::new_with_pos("N", Position::new(2.0, 2.7, 0.0)).unwrap(),
    );
    res2.set_atom(
        "CA",
        Atom::new_with_pos("C", Position::new(3.0, 3.5, 0.0)).unwrap(),
    );
    chain.set_group("2", res2);

    let report = chain
        .add_missing_hydrogens(db)
        .expect("Overall traversal must succeed even with partial residue failure");

    // 1. Backbone amide 'H' was added to residue 2
    let res2_after = chain.get_group("2").unwrap();
    assert!(
        res2_after.has_atom("H"),
        "Backbone amide H must have been added to residue 2"
    );

    // 2. Residue 2 must be recorded in residue_reports because it was modified
    assert!(
        report.residue_reports.contains_key("/model_1/A/2/"),
        "Residue 2 must be recorded in residue_reports"
    );

    // 3. Critically: Step 2 failure must NOT be silenced! It must be recorded in step_errors.
    let step2_err = report
        .step_errors
        .iter()
        .find(|(path, _)| path == "/model_1/A/2/");
    assert!(
        step2_err.is_some(),
        "Step 2 error for residue 2 must be recorded in step_errors: {:?}",
        report.step_errors
    );
    let err_msg = &step2_err.unwrap().1;
    assert!(
        err_msg.contains("common heavy atoms")
            || err_msg.contains("Sidechain/general hydrogenation"),
        "Error message must describe the superposition failure: {err_msg}"
    );
}
