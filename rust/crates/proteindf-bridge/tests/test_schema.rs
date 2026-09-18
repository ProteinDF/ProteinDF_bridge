// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::path::PathBuf;

use proteindf_bridge::atom::Atom;
use proteindf_bridge::atom_group::{AtomGroup, SchemaViolation};
use proteindf_bridge::format::{Format, Pdb};

fn fixture_path(filename: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join(filename)
}

#[test]
fn test_schema_level_checks_normal_hierarchy() {
    let mut root = AtomGroup::new();
    let mut model = AtomGroup::with_name("model_1");
    let mut chain = AtomGroup::with_name("A");
    let mut residue = AtomGroup::with_name("1");
    let mut atom = Atom::from_symbol("C").unwrap();
    atom.name = "CA".to_string();

    residue.set_atom("CA", atom);
    chain.set_group("1", residue);
    model.set_group("A", chain);
    root.set_group("model_1", model);

    // Root level: depth 0
    assert_eq!(root.path_depth(), 0);
    assert!(!root.is_model_level());
    assert!(!root.is_chain_level());
    assert!(!root.is_residue_level());

    // Model level: depth 1
    let m = root.get_group("model_1").unwrap();
    assert_eq!(m.path_depth(), 1);
    assert!(m.is_model_level());
    assert!(!m.is_chain_level());
    assert!(!m.is_residue_level());
    assert!(Format::is_protein(m));

    // Chain level: depth 2
    let c = m.get_group("A").unwrap();
    assert_eq!(c.path_depth(), 2);
    assert!(!c.is_model_level());
    assert!(c.is_chain_level());
    assert!(!c.is_residue_level());
    assert!(Format::is_chain(c));

    // Residue level: depth 3
    let r = c.get_group("1").unwrap();
    assert_eq!(r.path_depth(), 3);
    assert!(!r.is_model_level());
    assert!(!r.is_chain_level());
    assert!(r.is_residue_level());
    assert!(Format::is_residue(r));

    // Standard schema validation produces zero violations
    let violations = root.validate_schema();
    assert!(
        violations.is_empty(),
        "Expected no violations, got {:?}",
        violations
    );
}

#[test]
fn test_schema_real_pdb_fixture_1hls() {
    let path = fixture_path("1hls.pdb");
    let pdb = Pdb::from_file(&path, None).expect("Failed to load 1hls.pdb");
    let root = pdb
        .get_atomgroup(None, None)
        .expect("Failed to get atomgroup");

    // Real PDB loaded via Pdb parser must conform to the protein schema
    let violations = root.validate_schema();
    assert!(
        violations.is_empty(),
        "Real fixture 1hls.pdb should not have schema violations, but got: {:?}",
        violations
    );

    // Verify level checks on real fixture nodes
    for (_, model) in root.groups() {
        assert!(model.is_model_level());
        assert!(!model.is_chain_level());
        assert!(!model.is_residue_level());

        for (_, chain) in model.groups() {
            assert!(!chain.is_model_level());
            assert!(chain.is_chain_level());
            assert!(!chain.is_residue_level());

            for (_, residue) in chain.groups() {
                assert!(!residue.is_model_level());
                assert!(!residue.is_chain_level());
                assert!(residue.is_residue_level());
            }
        }
    }
}

#[test]
fn test_schema_violation_direct_atoms_in_chain() {
    let mut root = AtomGroup::new();
    let mut model = AtomGroup::with_name("model_1");
    let mut chain = AtomGroup::with_name("A");
    let mut residue = AtomGroup::with_name("1");

    let ca = Atom::from_symbol("C").unwrap();
    residue.set_atom("CA", ca);
    chain.set_group("1", residue);

    // Schema violation: direct water/HETATM atom attached directly to chain
    let mut water = Atom::from_symbol("O").unwrap();
    water.name = "O".to_string();
    chain.set_atom("HOH_100", water);

    model.set_group("A", chain);
    root.set_group("model_1", model);

    let chain_ref = root.get_group("model_1").unwrap().get_group("A").unwrap();
    // Positional check still says depth 2 (chain level)
    assert!(chain_ref.is_chain_level());
    // But structural check (Format::is_chain) fails because chain has direct atoms
    assert!(!Format::is_chain(chain_ref));

    let violations = root.validate_schema();
    assert_eq!(violations.len(), 1);
    match &violations[0] {
        SchemaViolation::DirectAtomsAtNonResidueLevel {
            path,
            depth,
            atom_keys,
        } => {
            assert_eq!(path, "/model_1/A/");
            assert_eq!(*depth, 2);
            assert_eq!(atom_keys, &vec!["HOH_100".to_string()]);
        }
        other => panic!("Unexpected violation type: {:?}", other),
    }

    // Check Display message formatting
    let msg = format!("{}", violations[0]);
    assert!(msg.contains("/model_1/A/"));
    assert!(msg.contains("depth 2"));
    assert!(msg.contains("HOH_100"));
}

#[test]
fn test_schema_violation_subgroup_in_residue_and_depth() {
    let mut root = AtomGroup::new();
    let mut model = AtomGroup::with_name("model_1");
    let mut chain = AtomGroup::with_name("A");
    let mut residue = AtomGroup::with_name("1");
    let mut child = AtomGroup::with_name("child_group");

    let atom = Atom::from_symbol("N").unwrap();
    child.set_atom("N1", atom);
    residue.set_group("child_group", child);
    chain.set_group("1", residue);
    model.set_group("A", chain);
    root.set_group("model_1", model);

    let violations = root.validate_schema();
    assert!(violations.iter().any(|v| matches!(
        v,
        SchemaViolation::SubgroupsInResidue { path, group_keys }
            if path == "/model_1/A/1/" && group_keys == &vec!["child_group".to_string()]
    )));
    assert!(violations.iter().any(|v| matches!(
        v,
        SchemaViolation::ExcessiveDepth { path, depth }
            if path == "/model_1/A/1/child_group/" && *depth == 4
    )));
}
