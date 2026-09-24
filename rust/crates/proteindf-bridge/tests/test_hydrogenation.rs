// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::path::PathBuf;

use proteindf_bridge::atom::Atom;
use proteindf_bridge::atom_group::AtomGroup;
use proteindf_bridge::ccd_templates::CcdTemplateDb;
use proteindf_bridge::format::Pdb;
use proteindf_bridge::hydrogenation::{
    add_hydrogens_to_component, add_hydrogens_to_component_in_place,
    add_hydrogens_to_component_with_options, HydrogenationOptions, MIN_SUPERPOSE_HEAVY_ATOMS,
};
use proteindf_bridge::position::Position;
use proteindf_bridge::superposer::Superposer;

fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

/// Helper: calculates Euclidean distance between two positions.
fn distance(p1: &Position, p2: &Position) -> f64 {
    let diff = *p1 - *p2;
    diff.length()
}

/// Helper: creates a copy of a residue group containing only heavy atoms (strips hydrogens).
fn strip_hydrogens(res: &AtomGroup) -> AtomGroup {
    let mut res_heavy = AtomGroup::with_name(&res.name);
    for (k, atom) in res.atoms() {
        if atom.atomic_number() != 1 {
            res_heavy.set_atom(k, atom.clone());
        }
    }
    res_heavy
}

/// Helper: retrieves a residue group by chain and residue key from a loaded structure.
fn get_residue<'a>(protein: &'a AtomGroup, chain_id: &str, res_key: &str) -> &'a AtomGroup {
    let model = protein
        .get_group("model_1")
        .or_else(|| protein.get_group("1"))
        .expect("model group must exist");
    let chain = model
        .get_group(chain_id)
        .unwrap_or_else(|| panic!("chain {chain_id} must exist"));
    chain
        .get_group(res_key)
        .unwrap_or_else(|| panic!("residue {res_key} in chain {chain_id} must exist"))
}

// 1. Identity superposition test:
// When the component's heavy atoms have the exact idealized coordinates from the CCD template,
// the superposed hydrogen coordinates must match the template's ideal coordinates with
// numerical precision (tolerance 1e-6).
#[test]
fn test_hydrogenation_identity_superposition() {
    let db = CcdTemplateDb::global();
    let template = db.lookup("ALA").expect("ALA template must exist in DB");

    // Construct component with only heavy atoms at ideal coordinates
    let mut ala_heavy = AtomGroup::with_name("ALA");
    let mut expected_h_count = 0;

    for atom in &template.atoms {
        if atom.is_hydrogen() {
            expected_h_count += 1;
        } else {
            let (x, y, z) = atom.ideal_xyz.expect("ALA heavy atoms must have ideal_xyz");
            let mut a = Atom::new_with_pos(&atom.element, Position::new(x, y, z))
                .expect("valid atom creation");
            a.name = atom.name.clone();
            ala_heavy.set_atom(&atom.name, a);
        }
    }

    assert!(expected_h_count > 0);
    assert_eq!(ala_heavy.pickup_atoms("H").len(), 0);

    // Run hydrogenation
    let hydrogenated = add_hydrogens_to_component(&ala_heavy, template)
        .expect("hydrogenation should succeed on identical coordinates");

    // In PR#38 single-component hydrogenation, backbone N 'H2' is unconditionally skipped to avoid
    // over-protonating internal/C-terminal peptide residues.
    assert!(
        !hydrogenated.has_atom("H2"),
        "Backbone N 'H2' should be skipped unconditionally"
    );

    // Verify all applicable hydrogens were added and coordinates match exactly
    for atom in &template.atoms {
        if atom.is_hydrogen() && atom.name != "H2" && atom.name != "H3" {
            let (ix, iy, iz) = atom.ideal_xyz.unwrap();
            let ideal_pos = Position::new(ix, iy, iz);

            let added = hydrogenated.pickup_atoms(&atom.name);
            assert_eq!(
                added.len(),
                1,
                "Hydrogen {} should be present in hydrogenated component",
                atom.name
            );

            let dist = distance(&added[0].xyz, &ideal_pos);
            assert!(
                dist < 1e-6,
                "Hydrogen {} coordinate mismatch: actual {:?}, ideal {:?}, dist = {}",
                atom.name,
                added[0].xyz,
                ideal_pos,
                dist
            );
        }
    }
}

// 2. Real PDB fixture (1hls.pdb) regression test:
// For real structures without hydrogens, verify that added hydrogens have chemically
// reasonable covalent bond lengths to their parent heavy atoms (expected range ~0.95 to 1.15 Angstroms).
#[test]
fn test_hydrogenation_real_fixture_1hls_bond_lengths() {
    let db = CcdTemplateDb::global();
    let pdb_path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&pdb_path, None).expect("failed to load 1hls.pdb");
    let protein = pdb.get_atomgroup(None, None).expect("failed to get ag");

    // Test ALA B/14: rigid methyl sidechain, fits on N, CA, C, CB (O excluded automatically)
    let template = db.lookup("ALA").expect("Template for ALA must exist");
    let res = get_residue(&protein, "B", "14");
    let res_heavy = strip_hydrogens(res);

    let report = add_hydrogens_to_component(&res_heavy, template)
        .expect("Failed to hydrogenate residue B/14 (ALA)");

    assert!(
        report.get_number_of_atoms() > res_heavy.get_number_of_atoms(),
        "Hydrogens should have been added to residue B/14"
    );

    // Verify bond lengths for HA and all methyl hydrogens (HB1, HB2, HB3)
    let tested_hydrogens = ["HA", "HB1", "HB2", "HB3"];
    for h_name in tested_hydrogens {
        let parent_heavy = template
            .bonds
            .iter()
            .find_map(|(a1, a2, _)| {
                if a1 == h_name {
                    Some(a2.as_str())
                } else if a2 == h_name {
                    Some(a1.as_str())
                } else {
                    None
                }
            })
            .unwrap_or_else(|| panic!("Bond for {h_name} not found in template ALA"));

        let h_atoms = report.pickup_atoms(h_name);
        let heavy_atoms = report.pickup_atoms(parent_heavy);

        assert!(
            !h_atoms.is_empty(),
            "Hydrogen {h_name} should be present in ALA"
        );
        assert!(
            !heavy_atoms.is_empty(),
            "Parent heavy atom {parent_heavy} should be present in ALA"
        );

        let bond_len = distance(&h_atoms[0].xyz, &heavy_atoms[0].xyz);
        assert!(
            (0.95..=1.15).contains(&bond_len),
            "Bond length {parent_heavy}-{h_name} in ALA is {bond_len:.3} A, outside chemically expected range [0.95, 1.15] A"
        );
    }
}

// 3. Custom options: explicit heavy atom selection for flexible / aromatic sidechains:
// Demonstrates using `HydrogenationOptions` to specify fitting heavy atoms
// (e.g. backbone heavy atoms for VAL, aromatic ring heavy atoms for PHE).
#[test]
fn test_hydrogenation_options_custom_heavy_atoms() {
    let db = CcdTemplateDb::global();
    let pdb_path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&pdb_path, None).expect("failed to load 1hls.pdb");
    let protein = pdb.get_atomgroup(None, None).expect("failed to get ag");

    // Case A: Backbone fit for VAL A/3 -> CA-HA bond length
    {
        let template = db.lookup("VAL").expect("VAL template exists");
        let res = get_residue(&protein, "A", "3");
        let res_heavy = strip_hydrogens(res);

        let options = HydrogenationOptions {
            fit_heavy_atoms: Some(&["N", "CA", "C", "CB"]),
            auto_exclude_distorted_atoms: true,
        };

        let hydrogenated = add_hydrogens_to_component_with_options(&res_heavy, template, &options)
            .expect("VAL hydrogenation with custom backbone atoms should succeed");

        let ha = &hydrogenated.pickup_atoms("HA")[0].xyz;
        let ca = &hydrogenated.pickup_atoms("CA")[0].xyz;
        let dist_ha = distance(ha, ca);
        assert!(
            (0.95..=1.15).contains(&dist_ha),
            "Bond length CA-HA in VAL is {dist_ha:.3} A, expected [0.95, 1.15] A"
        );
    }

    // Case B: Aromatic ring fit for PHE B/24 -> ring hydrogens (HD1, HD2, HE1, HE2, HZ)
    {
        let template = db.lookup("PHE").expect("PHE template exists");
        let res = get_residue(&protein, "B", "24");
        let res_heavy = strip_hydrogens(res);

        let options = HydrogenationOptions {
            fit_heavy_atoms: Some(&["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]),
            auto_exclude_distorted_atoms: true,
        };

        let hydrogenated = add_hydrogens_to_component_with_options(&res_heavy, template, &options)
            .expect("PHE hydrogenation with ring atoms should succeed");

        let ring_hydrogens = [
            ("HD1", "CD1"),
            ("HD2", "CD2"),
            ("HE1", "CE1"),
            ("HE2", "CE2"),
            ("HZ", "CZ"),
        ];

        for (h_name, parent_heavy) in ring_hydrogens {
            let h_pos = &hydrogenated.pickup_atoms(h_name)[0].xyz;
            let heavy_pos = &hydrogenated.pickup_atoms(parent_heavy)[0].xyz;
            let dist_h = distance(h_pos, heavy_pos);
            assert!(
                (0.95..=1.15).contains(&dist_h),
                "Bond length {parent_heavy}-{h_name} in PHE is {dist_h:.3} A, expected [0.95, 1.15] A"
            );
        }
    }
}

// 4. Verification of Bug 2: fit_heavy_atoms explicitly naming 'O' still guards against distortion
// by default, unless auto_exclude_distorted_atoms is explicitly disabled.
#[test]
fn test_hydrogenation_guard_not_bypassed_by_explicit_fit_atoms() {
    let db = CcdTemplateDb::global();
    let pdb_path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&pdb_path, None).expect("failed to load 1hls.pdb");
    let protein = pdb.get_atomgroup(None, None).expect("failed to get ag");

    let template = db.lookup("ALA").expect("ALA template exists");
    let res = get_residue(&protein, "B", "14");
    let res_heavy = strip_hydrogens(res);

    // Case A: User explicitly provides ["N", "CA", "C", "O", "CB"] with default auto_exclude = true.
    // The distorted 'O' should still be excluded, yielding proper CA-HA bond length (~1.08 A).
    let options_safe = HydrogenationOptions {
        fit_heavy_atoms: Some(&["N", "CA", "C", "O", "CB"]),
        auto_exclude_distorted_atoms: true,
    };
    let h_safe = add_hydrogens_to_component_with_options(&res_heavy, template, &options_safe)
        .expect("Safe fit should succeed");
    let ca_ha_dist = distance(
        &h_safe.pickup_atoms("CA")[0].xyz,
        &h_safe.pickup_atoms("HA")[0].xyz,
    );
    assert!(
        (0.95..=1.15).contains(&ca_ha_dist),
        "With auto_exclude=true, CA-HA bond length should be valid ({ca_ha_dist:.3} A)"
    );

    // Case B: User forces auto_exclude = false with 'O' included.
    // The ~170 deg psi discrepancy distorts the fit, resulting in invalid CA-HA length (> 1.5 A).
    let options_forced = HydrogenationOptions {
        fit_heavy_atoms: Some(&["N", "CA", "C", "O", "CB"]),
        auto_exclude_distorted_atoms: false,
    };
    let h_forced = add_hydrogens_to_component_with_options(&res_heavy, template, &options_forced)
        .expect("Forced fit should succeed");
    let ca_ha_dist_forced = distance(
        &h_forced.pickup_atoms("CA")[0].xyz,
        &h_forced.pickup_atoms("HA")[0].xyz,
    );
    assert!(
        ca_ha_dist_forced > 1.5,
        "With auto_exclude=false and 'O' forced, CA-HA length should show distortion (> 1.5 A), got {ca_ha_dist_forced:.3} A"
    );
}

// 5. Verification of Bug 1: Collinear and degenerate heavy atoms must be rejected with an error.
#[test]
fn test_hydrogenation_collinear_heavy_atoms_rejected() {
    let db = CcdTemplateDb::global();
    let template = db.lookup("ALA").expect("ALA template exists");

    // Construct synthetic component where 3 heavy atoms are strictly collinear: (0,0,0), (1,0,0), (2,0,0)
    let mut collinear_group = AtomGroup::with_name("ALA");
    collinear_group.set_atom(
        "N",
        Atom::new_with_pos("N", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    collinear_group.set_atom(
        "CA",
        Atom::new_with_pos("C", Position::new(1.0, 0.0, 0.0)).unwrap(),
    );
    collinear_group.set_atom(
        "C",
        Atom::new_with_pos("C", Position::new(2.0, 0.0, 0.0)).unwrap(),
    );

    let err = add_hydrogens_to_component(&collinear_group, template);
    assert!(
        err.is_err(),
        "Collinear points must be rejected as rigid superposition cannot be uniquely determined"
    );
    let err_msg = err.unwrap_err().to_string();
    assert!(
        err_msg.contains("collinear or degenerate"),
        "Unexpected error message: {err_msg}"
    );

    // Direct Superposer check on collinear points
    let mut g1 = AtomGroup::new();
    let mut g2 = AtomGroup::new();
    for (i, name) in ["A1", "A2", "A3"].iter().enumerate() {
        let x = i as f64;
        let mut a1 = Atom::new_with_pos("C", Position::new(x, 0.0, 0.0)).unwrap();
        a1.name = name.to_string();
        let mut a2 = Atom::new_with_pos("C", Position::new(0.0, x, 0.0)).unwrap();
        a2.name = name.to_string();
        g1.set_atom(name, a1);
        g2.set_atom(name, a2);
    }
    let sp_err = Superposer::new(&g1, &g2);
    assert!(
        sp_err.is_err(),
        "Superposer::new must reject collinear points"
    );
}

// 6. Insufficient common heavy atoms error test:
// When fewer than MIN_SUPERPOSE_HEAVY_ATOMS (3) matching heavy atoms exist,
// the function must return an error and not perform silent fallback.
#[test]
fn test_hydrogenation_insufficient_heavy_atoms_error() {
    assert_eq!(MIN_SUPERPOSE_HEAVY_ATOMS, 3);
    let db = CcdTemplateDb::global();
    let template = db.lookup("ALA").expect("ALA template exists");

    // Case A: 0 heavy atoms
    let empty_group = AtomGroup::with_name("ALA");
    let err_0 = add_hydrogens_to_component(&empty_group, template);
    assert!(err_0.is_err());
    let err_msg_0 = err_0.unwrap_err().to_string();
    assert!(
        err_msg_0.contains("Insufficient common heavy atoms"),
        "Unexpected error message: {err_msg_0}"
    );

    // Case B: 1 heavy atom
    let mut group_1 = AtomGroup::with_name("ALA");
    group_1.set_atom(
        "CA",
        Atom::new_with_pos("C", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    let err_1 = add_hydrogens_to_component(&group_1, template);
    assert!(err_1.is_err());

    // Case C: 2 heavy atoms (below minimum 3)
    let mut group_2 = AtomGroup::with_name("ALA");
    group_2.set_atom(
        "CA",
        Atom::new_with_pos("C", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    group_2.set_atom(
        "N",
        Atom::new_with_pos("N", Position::new(1.0, 0.0, 0.0)).unwrap(),
    );
    let err_2 = add_hydrogens_to_component(&group_2, template);
    assert!(err_2.is_err());
    let err_msg_2 = err_2.unwrap_err().to_string();
    assert!(
        err_msg_2.contains("Insufficient common heavy atoms"),
        "Unexpected error message: {err_msg_2}"
    );

    // Case D: 3 non-collinear heavy atoms succeeds
    let mut group_3 = AtomGroup::with_name("ALA");
    let n = template.get_atom("N").unwrap().ideal_xyz.unwrap();
    let ca = template.get_atom("CA").unwrap().ideal_xyz.unwrap();
    let c = template.get_atom("C").unwrap().ideal_xyz.unwrap();
    group_3.set_atom(
        "N",
        Atom::new_with_pos("N", Position::new(n.0, n.1, n.2)).unwrap(),
    );
    group_3.set_atom(
        "CA",
        Atom::new_with_pos("C", Position::new(ca.0, ca.1, ca.2)).unwrap(),
    );
    group_3.set_atom(
        "C",
        Atom::new_with_pos("C", Position::new(c.0, c.1, c.2)).unwrap(),
    );
    let ok_3 = add_hydrogens_to_component(&group_3, template);
    assert!(
        ok_3.is_ok(),
        "3 non-collinear heavy atoms should satisfy MIN_SUPERPOSE_HEAVY_ATOMS: {:?}",
        ok_3.err()
    );
}

// 7. In-place modification and idempotency:
// Calling add_hydrogens_to_component_in_place a second time should add 0 hydrogens.
#[test]
fn test_hydrogenation_in_place_idempotency() {
    let db = CcdTemplateDb::global();
    let template = db.lookup("SER").expect("SER template exists");

    let mut ser = AtomGroup::with_name("SER");
    for atom in &template.atoms {
        if !atom.is_hydrogen() {
            let (x, y, z) = atom.ideal_xyz.unwrap();
            ser.set_atom(
                &atom.name,
                Atom::new_with_pos(&atom.element, Position::new(x, y, z)).unwrap(),
            );
        }
    }

    let report1 =
        add_hydrogens_to_component_in_place(&mut ser, template).expect("first pass should succeed");
    assert!(report1.added_hydrogens > 0);

    let count_after_first = ser.get_number_of_atoms();

    let report2 = add_hydrogens_to_component_in_place(&mut ser, template)
        .expect("second pass should succeed");
    assert_eq!(report2.added_hydrogens, 0);
    assert!(report2.added_atom_names.is_empty());
    assert_eq!(ser.get_number_of_atoms(), count_after_first);
}

// 8. Valence rule and parent heavy atom check (Review Round 4 Bug 1):
// - Free amino acid templates (e.g. ALA) contain an N-terminal free amine (N with H and H2)
//   and a C-terminal carboxylate (C with O, OXT, and HXT).
// - In an internal peptide residue (where OXT is absent):
//   1. The backbone amide nitrogen must only receive 'H', and secondary amine hydrogens ('H2', 'H3')
//      must NOT be added (avoiding over-protonation).
//   2. Hydrogens whose parent heavy atom is missing from the component (e.g. HXT bonded to OXT)
//      must NOT be added.
#[test]
fn test_hydrogenation_internal_peptide_valence_and_parent_heavy_check() {
    let db = CcdTemplateDb::global();
    let template = db.lookup("ALA").expect("ALA template exists");

    // Construct an internal ALA residue: has N, CA, C, O, CB (lacks OXT)
    let mut internal_ala = AtomGroup::with_name("ALA");
    for atom in &template.atoms {
        if !atom.is_hydrogen() && atom.name != "OXT" {
            let (x, y, z) = atom.ideal_xyz.unwrap();
            internal_ala.set_atom(
                &atom.name,
                Atom::new_with_pos(&atom.element, Position::new(x, y, z)).unwrap(),
            );
        }
    }

    assert!(internal_ala.has_atom("N"));
    assert!(!internal_ala.has_atom("OXT"));

    let hydrogenated = add_hydrogens_to_component(&internal_ala, template)
        .expect("Hydrogenation of internal ALA residue should succeed");

    // 1. Primary amide hydrogen 'H' must be added
    assert!(
        hydrogenated.has_atom("H"),
        "Backbone amide hydrogen 'H' should be added"
    );

    // 2. Secondary amine hydrogen 'H2' must NOT be added to an internal peptide residue
    assert!(
        !hydrogenated.has_atom("H2"),
        "Secondary amine hydrogen 'H2' must NOT be added to internal peptide residue (over-protonation)"
    );

    // 3. C-terminal carboxylate hydrogen 'HXT' must NOT be added when parent 'OXT' is missing
    assert!(
        !hydrogenated.has_atom("HXT"),
        "C-terminal hydrogen 'HXT' must NOT be added when parent 'OXT' is absent"
    );

    // 4. Sidechain and alpha hydrogens must be present
    assert!(hydrogenated.has_atom("HA"));
    assert!(hydrogenated.has_atom("HB1"));
    assert!(hydrogenated.has_atom("HB2"));
    assert!(hydrogenated.has_atom("HB3"));
}

// 9. Atomic rollback on error during in-place hydrogenation:
// If any hydrogen atom cannot be resolved (e.g. missing ideal_xyz),
// add_hydrogens_to_component_in_place must return Err AND leave the component
// completely unmodified (0 hydrogens partially added).
#[test]
fn test_hydrogenation_in_place_atomic_on_error() {
    let db = CcdTemplateDb::global();
    let template = db.lookup("ALA").expect("ALA template exists");

    // Create a modified template where the LAST hydrogen has None for ideal_xyz
    let mut bad_template = template.clone();
    let last_h_idx = bad_template
        .atoms
        .iter()
        .rposition(|a| a.is_hydrogen())
        .expect("ALA has hydrogens");
    bad_template.atoms[last_h_idx].ideal_xyz = None;

    // Prepare component with ALA heavy atoms
    let mut ala_heavy = AtomGroup::with_name("ALA");
    for atom in &template.atoms {
        if !atom.is_hydrogen() {
            let (x, y, z) = atom.ideal_xyz.unwrap();
            ala_heavy.set_atom(
                &atom.name,
                Atom::new_with_pos(&atom.element, Position::new(x, y, z)).unwrap(),
            );
        }
    }

    let initial_atom_count = ala_heavy.get_number_of_atoms();

    // Call in-place hydrogenation with the flawed template
    let res = add_hydrogens_to_component_in_place(&mut ala_heavy, &bad_template);
    assert!(res.is_err(), "Must fail when a hydrogen has no ideal_xyz");

    // Verify atomicity: NO hydrogens were added to ala_heavy
    assert_eq!(
        ala_heavy.get_number_of_atoms(),
        initial_atom_count,
        "Component must remain completely unmodified on failure"
    );
    for atom in &template.atoms {
        if atom.is_hydrogen() {
            assert!(
                !ala_heavy.has_atom(&atom.name),
                "Hydrogen '{}' should not exist in component after atomic failure",
                atom.name
            );
        }
    }
}

// 10. Regression test: H2 / H3 on backbone N must be unconditionally skipped in single-component
// hydrogenation (PR#38), whether the component is an isolated amino acid, C-terminal (with OXT),
// or internal. Only standard amide 'H' is added; full N-terminal amine capping/protonation (NH3+)
// is deferred to PR#39.
#[test]
fn test_hydrogenation_backbone_n_skips_h2_h3_unconditionally() {
    let db = CcdTemplateDb::global();
    let template = db.lookup("ALA").expect("ALA template exists");

    // Case A: Isolated / C-terminal ALA (has OXT)
    let mut c_term_ala = AtomGroup::with_name("ALA");
    for atom in &template.atoms {
        if !atom.is_hydrogen() {
            let (x, y, z) = atom.ideal_xyz.unwrap();
            c_term_ala.set_atom(
                &atom.name,
                Atom::new_with_pos(&atom.element, Position::new(x, y, z)).unwrap(),
            );
        }
    }
    assert!(c_term_ala.has_atom("OXT"));
    assert!(c_term_ala.has_atom("N"));

    let hydrogenated_c_term = add_hydrogens_to_component(&c_term_ala, template)
        .expect("Hydrogenation of C-terminal ALA should succeed");

    assert!(
        hydrogenated_c_term.has_atom("H"),
        "Primary backbone hydrogen 'H' must be added"
    );
    assert!(
        !hydrogenated_c_term.has_atom("H2"),
        "'H2' must NOT be added to backbone N even when OXT is present"
    );
    assert!(
        !hydrogenated_c_term.has_atom("H3"),
        "'H3' must NOT be added to backbone N"
    );
    assert!(
        hydrogenated_c_term.has_atom("HXT"),
        "'HXT' must be added when parent 'OXT' is present"
    );

    // Case B: Internal ALA (lacks OXT)
    let mut internal_ala = AtomGroup::with_name("ALA");
    for atom in &template.atoms {
        if !atom.is_hydrogen() && atom.name != "OXT" {
            let (x, y, z) = atom.ideal_xyz.unwrap();
            internal_ala.set_atom(
                &atom.name,
                Atom::new_with_pos(&atom.element, Position::new(x, y, z)).unwrap(),
            );
        }
    }
    assert!(!internal_ala.has_atom("OXT"));

    let hydrogenated_internal = add_hydrogens_to_component(&internal_ala, template)
        .expect("Hydrogenation of internal ALA should succeed");

    assert!(
        hydrogenated_internal.has_atom("H"),
        "Primary backbone hydrogen 'H' must be added"
    );
    assert!(
        !hydrogenated_internal.has_atom("H2"),
        "'H2' must NOT be added to internal peptide residue"
    );
    assert!(
        !hydrogenated_internal.has_atom("H3"),
        "'H3' must NOT be added to internal peptide residue"
    );
    assert!(
        !hydrogenated_internal.has_atom("HXT"),
        "'HXT' must NOT be added when parent 'OXT' is absent"
    );
}

// 11. Regression test: Hydrogens without bond information in template are safely skipped.
// If a hydrogen in template.atoms has no corresponding bond in template.bonds, it must NOT
// be added blindly, but skipped safely.
#[test]
fn test_hydrogenation_skips_hydrogen_with_missing_bond_info() {
    let db = CcdTemplateDb::global();
    let template = db.lookup("ALA").expect("ALA template exists");

    // Construct a template where the bond connecting 'HA' to 'CA' is removed from template.bonds
    let mut stripped_bonds_template = template.clone();
    stripped_bonds_template
        .bonds
        .retain(|(a1, a2, _)| !((a1 == "HA" && a2 == "CA") || (a1 == "CA" && a2 == "HA")));

    let mut ala_heavy = AtomGroup::with_name("ALA");
    for atom in &template.atoms {
        if !atom.is_hydrogen() && atom.name != "OXT" {
            let (x, y, z) = atom.ideal_xyz.unwrap();
            ala_heavy.set_atom(
                &atom.name,
                Atom::new_with_pos(&atom.element, Position::new(x, y, z)).unwrap(),
            );
        }
    }

    let hydrogenated = add_hydrogens_to_component(&ala_heavy, &stripped_bonds_template)
        .expect("Hydrogenation should succeed even if some hydrogens have no bond information");

    // 'HA' should have been safely skipped because its parent heavy atom could not be determined
    assert!(
        !hydrogenated.has_atom("HA"),
        "Hydrogen 'HA' without bond information in template must be safely skipped"
    );
    // Other hydrogens with valid bond info (e.g. 'H', 'HB1') should still be added
    assert!(hydrogenated.has_atom("H"));
    assert!(hydrogenated.has_atom("HB1"));
}
