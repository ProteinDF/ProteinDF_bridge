// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::path::PathBuf;

use proteindf_bridge::atom::Atom;
use proteindf_bridge::atom_group::AtomGroup;
use proteindf_bridge::format::pdb::Pdb;
use proteindf_bridge::position::Position;
use proteindf_bridge::ssbond::SSBond;

fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

#[test]
fn test_ssbond_1hls_doctest() {
    // 1:1 reproduction of the doctest in proteindf_bridge/ssbond.py:
    // >>> tmp_pdb = Pdb('data/1hls.pdb')
    // >>> models = tmp_pdb.get_atomgroup()
    // >>> model = models.get_group('model_1')
    // >>> ssb = SSBond(model)
    // >>> print(ssb.get_bonds())
    // [('/model_1/A/6/', '/model_1/A/11/'), ('/model_1/A/7/', '/model_1/B/7/'), ('/model_1/A/20/', '/model_1/B/19/')]
    let path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&path, None).expect("failed to load 1hls.pdb");
    let models = pdb.get_atomgroup(None, None).expect("get_atomgroup failed");
    let model = models.get_group("model_1").expect("model_1 not found");

    let mut ssb = SSBond::new(model);
    let bonds = ssb.get_bonds();

    assert_eq!(bonds.len(), 3, "expected 3 disulfide bond pairs in 1hls");
    assert_eq!(
        bonds,
        &[
            ("/model_1/A/6/".to_string(), "/model_1/A/11/".to_string()),
            ("/model_1/A/7/".to_string(), "/model_1/B/7/".to_string()),
            ("/model_1/A/20/".to_string(), "/model_1/B/19/".to_string()),
        ]
    );

    // Verify distance values corresponding to doctest reference values:
    // /model_1/A/6/ <-> /model_1/A/11/: ~2.020336 Å
    // /model_1/A/7/ <-> /model_1/B/7/:  ~2.018259 Å
    // /model_1/A/20/ <-> /model_1/B/19/: ~2.017930 Å
    let sg_a6 = model["A"]["6"].get_atom("SG").unwrap();
    let sg_a11 = model["A"]["11"].get_atom("SG").unwrap();
    let dist1 = sg_a6.xyz.distance_from(&sg_a11.xyz);
    assert!((dist1 - 2.020336).abs() < 1e-4);

    let sg_a7 = model["A"]["7"].get_atom("SG").unwrap();
    let sg_b7 = model["B"]["7"].get_atom("SG").unwrap();
    let dist2 = sg_a7.xyz.distance_from(&sg_b7.xyz);
    assert!((dist2 - 2.018259).abs() < 1e-4);

    let sg_a20 = model["A"]["20"].get_atom("SG").unwrap();
    let sg_b19 = model["B"]["19"].get_atom("SG").unwrap();
    let dist3 = sg_a20.xyz.distance_from(&sg_b19.xyz);
    assert!((dist3 - 2.017930).abs() < 1e-4);
}

#[test]
fn test_ssbond_convenience_find_bonds() {
    let path = test_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&path, None).expect("failed to load 1hls.pdb");
    let models = pdb.get_atomgroup(None, None).expect("get_atomgroup failed");
    let model = models.get_group("model_1").expect("model_1 not found");

    let bonds = SSBond::find_bonds(model);
    assert_eq!(bonds.len(), 3);
}

#[test]
fn test_ssbond_threshold_and_cyx() {
    // Synthetic model with CYS and CYX
    let mut model = AtomGroup::new();
    model.set_path("/model_1/".to_string());

    let mut chain_a = AtomGroup::new();

    // Residue 1: CYS at (0.0, 0.0, 0.0)
    let mut res1 = AtomGroup::new();
    res1.name = "CYS".to_string();
    let mut sg1 = Atom::new_with_pos("S", Position::new(0.0, 0.0, 0.0)).unwrap();
    sg1.name = "SG".to_string();
    res1.set_atom("SG", sg1);
    chain_a.set_group("1", res1);

    // Residue 2: CYX at (2.30, 0.0, 0.0) -> distance 2.30 < 2.31 (should match)
    let mut res2 = AtomGroup::new();
    res2.name = "CYX".to_string();
    let mut sg2 = Atom::new_with_pos("S", Position::new(2.30, 0.0, 0.0)).unwrap();
    sg2.name = "SG".to_string();
    res2.set_atom("SG", sg2);
    chain_a.set_group("2", res2);

    // Residue 3: CYS at (5.0, 0.0, 0.0) -> distance to res2 is 2.70 > 2.31 (should NOT match)
    let mut res3 = AtomGroup::new();
    res3.name = "CYS".to_string();
    let mut sg3 = Atom::new_with_pos("S", Position::new(5.0, 0.0, 0.0)).unwrap();
    sg3.name = "SG".to_string();
    res3.set_atom("SG", sg3);
    chain_a.set_group("3", res3);

    // Residue 4: CYS without SG atom -> should be safely ignored
    let mut res4 = AtomGroup::new();
    res4.name = "CYS".to_string();
    let mut ca4 = Atom::new_with_pos("C", Position::new(0.0, 1.0, 0.0)).unwrap();
    ca4.name = "CA".to_string();
    res4.set_atom("CA", ca4);
    chain_a.set_group("4", res4);

    model.set_group("A", chain_a);

    let bonds = SSBond::find_bonds(&model);
    assert_eq!(bonds.len(), 1);
    assert_eq!(
        bonds[0],
        ("/model_1/A/1/".to_string(), "/model_1/A/2/".to_string())
    );
}
