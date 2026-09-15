// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::path::PathBuf;

use proteindf_bridge::atom::Atom;
use proteindf_bridge::atom_group::AtomGroup;
use proteindf_bridge::format::pdb::Pdb;
use proteindf_bridge::neutralize::Neutralize;
use proteindf_bridge::position::Position;

fn get_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

fn make_atom(name: &str, symbol: &str, pos: Position) -> Atom {
    let mut atom = Atom::new_with_pos(symbol, pos).unwrap();
    atom.name = name.to_string();
    atom
}

#[test]
fn test_exempt_list() {
    let mut protein = AtomGroup::with_name("protein");
    let mut model = AtomGroup::with_name("model_1");
    let mut chain_a = AtomGroup::with_name("A");

    let mut glu = AtomGroup::with_name("1");
    glu.name = "GLU".to_string();
    glu.set_atom("CD", make_atom("CD", "C", Position::new(0.0, 0.0, 0.0)));
    glu.set_atom("OE1", make_atom("OE1", "O", Position::new(0.0, 1.0, 0.0)));
    glu.set_atom("OE2", make_atom("OE2", "O", Position::new(1.0, 0.0, 0.0)));
    chain_a.set_group("1", glu);

    let mut lys = AtomGroup::with_name("2");
    lys.name = "LYS".to_string();
    lys.set_atom("NZ", make_atom("NZ", "N", Position::new(0.5, 2.0, 0.0)));
    lys.set_atom("HZ1", make_atom("HZ1", "H", Position::new(0.5, 2.5, 0.0)));
    lys.set_atom("HZ2", make_atom("HZ2", "H", Position::new(0.0, 2.0, 0.5)));
    lys.set_atom("HZ3", make_atom("HZ3", "H", Position::new(1.0, 2.0, 0.5)));
    chain_a.set_group("2", lys);

    model.set_group("A", chain_a);
    protein.set_group("model_1", model.clone());

    let exempt = Neutralize::exempt_list(&model);
    assert_eq!(exempt.len(), 2);
    assert_eq!(
        exempt,
        vec![
            ("A".to_string(), "1".to_string(), "GLU".to_string()),
            ("A".to_string(), "2".to_string(), "LYS".to_string()),
        ]
    );

    let neut = Neutralize::new(&protein).unwrap();
    assert_eq!(neut.neutralized().groups().count(), 1);
}

#[test]
fn test_non_destructive_input() {
    let mut protein = AtomGroup::with_name("protein");
    let mut model = AtomGroup::with_name("model_1");
    let mut chain_a = AtomGroup::with_name("A");

    let mut glu = AtomGroup::with_name("1");
    glu.name = "GLU".to_string();
    glu.set_atom("CD", make_atom("CD", "C", Position::new(0.0, 0.0, 0.0)));
    glu.set_atom("OE1", make_atom("OE1", "O", Position::new(0.0, 1.0, 0.0)));
    glu.set_atom("OE2", make_atom("OE2", "O", Position::new(1.0, 0.0, 0.0)));
    chain_a.set_group("1", glu);
    model.set_group("A", chain_a);
    protein.set_group("model_1", model);

    let initial_atom_count = protein.get_number_of_all_atoms();
    assert_eq!(initial_atom_count, 3);

    let neut = Neutralize::new(&protein).unwrap();
    // Verify original protein is not modified
    assert_eq!(protein.get_number_of_all_atoms(), initial_atom_count);
    // Neutralized object has additional ion atom (3 + 1 = 4)
    assert_eq!(neut.neutralized().get_number_of_all_atoms(), 4);
}

#[test]
fn test_neutralize_charged_residues() {
    let mut protein = AtomGroup::with_name("protein");
    let mut model = AtomGroup::with_name("model_1");
    let mut chain_a = AtomGroup::with_name("A");

    // ASP residue
    let mut asp = AtomGroup::with_name("1");
    asp.name = "ASP".to_string();
    asp.set_atom("CG", make_atom("CG", "C", Position::new(0.0, 0.0, 0.0)));
    asp.set_atom("OD1", make_atom("OD1", "O", Position::new(0.0, 1.0, 0.0)));
    asp.set_atom("OD2", make_atom("OD2", "O", Position::new(1.0, 0.0, 0.0)));
    chain_a.set_group("1", asp);

    // LYS residue
    let mut lys = AtomGroup::with_name("2");
    lys.name = "LYS".to_string();
    lys.set_atom("NZ", make_atom("NZ", "N", Position::new(10.0, 0.0, 0.0)));
    lys.set_atom("HZ1", make_atom("HZ1", "H", Position::new(10.0, 1.0, 0.0)));
    lys.set_atom("HZ2", make_atom("HZ2", "H", Position::new(10.0, 0.0, 1.0)));
    lys.set_atom("HZ3", make_atom("HZ3", "H", Position::new(11.0, 0.0, 0.0)));
    chain_a.set_group("2", lys);

    // ARG residue
    let mut arg = AtomGroup::with_name("3");
    arg.name = "ARG".to_string();
    arg.set_atom("CZ", make_atom("CZ", "C", Position::new(20.0, 0.0, 0.0)));
    arg.set_atom("NH1", make_atom("NH1", "N", Position::new(20.0, 1.0, 0.0)));
    arg.set_atom("NH2", make_atom("NH2", "N", Position::new(21.0, 0.0, 0.0)));
    chain_a.set_group("3", arg);

    model.set_group("A", chain_a);
    protein.set_group("model_1", model);

    let neut = Neutralize::new(&protein).unwrap();
    let neut_obj = neut.neutralized();

    // Check ASP got Na
    let m1 = neut_obj.get_group("model_1").unwrap();
    let ca = m1.get_group("A").unwrap();
    let neut_asp = ca.get_group("1").unwrap();
    assert!(neut_asp.has_atomkey("3_Na0"));
    let na = neut_asp.get_atom("3_Na0").unwrap();
    assert_eq!(na.name, "Na");

    // Check LYS got Cl
    let neut_lys = ca.get_group("2").unwrap();
    assert!(neut_lys.has_atomkey("4_Cl0"));
    let cl = neut_lys.get_atom("4_Cl0").unwrap();
    assert_eq!(cl.name, "Cl");

    // Check ARG got Cl
    let neut_arg = ca.get_group("3").unwrap();
    assert!(neut_arg.has_atomkey("3_Cl0"));
    let cl_arg = neut_arg.get_atom("3_Cl0").unwrap();
    assert_eq!(cl_arg.name, "Cl");
}

#[test]
fn test_neutralize_termini() {
    let mut protein = AtomGroup::with_name("protein");
    let mut model = AtomGroup::with_name("model_1");
    let mut chain_a = AtomGroup::with_name("A");

    // N-terminal residue with H3
    let mut nterm = AtomGroup::with_name("1");
    nterm.name = "ALA".to_string();
    nterm.set_atom("N", make_atom("N", "N", Position::new(0.0, 0.0, 0.0)));
    nterm.set_atom("H1", make_atom("H1", "H", Position::new(0.0, 1.0, 0.0)));
    nterm.set_atom("H2", make_atom("H2", "H", Position::new(1.0, 0.0, 0.0)));
    nterm.set_atom("H3", make_atom("H3", "H", Position::new(0.0, 0.0, 1.0)));
    chain_a.set_group("1", nterm);

    // C-terminal residue with OXT
    let mut cterm = AtomGroup::with_name("2");
    cterm.name = "ALA".to_string();
    cterm.set_atom("C", make_atom("C", "C", Position::new(10.0, 0.0, 0.0)));
    cterm.set_atom("O", make_atom("O", "O", Position::new(10.0, 1.0, 0.0)));
    cterm.set_atom("OXT", make_atom("OXT", "O", Position::new(11.0, 0.0, 0.0)));
    chain_a.set_group("2", cterm);

    model.set_group("A", chain_a);
    protein.set_group("model_1", model);

    let neut = Neutralize::new(&protein).unwrap();
    let neut_obj = neut.neutralized();

    let m1 = neut_obj.get_group("model_1").unwrap();
    let ca = m1.get_group("A").unwrap();
    let res1 = ca.get_group("1").unwrap();
    assert!(res1.has_atomkey("Cl0"));
    assert_eq!(res1.get_atom("Cl0").unwrap().name, "Cl");

    let res2 = ca.get_group("2").unwrap();
    assert!(res2.has_atomkey("Na0"));
    assert_eq!(res2.get_atom("Na0").unwrap().name, "Na");
}

#[test]
fn test_neutralize_fad() {
    let mut protein = AtomGroup::with_name("protein");
    let mut model = AtomGroup::with_name("model_1");
    let mut chain_a = AtomGroup::with_name("A");

    let mut fad = AtomGroup::with_name("1");
    fad.name = "FAD".to_string();
    fad.set_atom("P", make_atom("P", "P", Position::new(0.0, 0.0, 0.0)));
    fad.set_atom("O1P", make_atom("O1P", "O", Position::new(0.0, 1.0, 0.0)));
    fad.set_atom("O2P", make_atom("O2P", "O", Position::new(1.0, 0.0, 0.0)));
    fad.set_atom("PA", make_atom("PA", "P", Position::new(5.0, 0.0, 0.0)));
    fad.set_atom("O1A", make_atom("O1A", "O", Position::new(5.0, 1.0, 0.0)));
    fad.set_atom("O2A", make_atom("O2A", "O", Position::new(6.0, 0.0, 0.0)));
    chain_a.set_group("1", fad);

    model.set_group("A", chain_a);
    protein.set_group("model_1", model);

    let neut = Neutralize::new(&protein).unwrap();
    let neut_obj = neut.neutralized();

    let m1 = neut_obj.get_group("model_1").unwrap();
    let ca = m1.get_group("A").unwrap();
    let res = ca.get_group("1").unwrap();
    assert!(res.has_atomkey("3_Na10"));
    assert!(res.has_atomkey("3_Na20"));
    assert_eq!(res.get_atom("3_Na10").unwrap().name, "Na");
    assert_eq!(res.get_atom("3_Na20").unwrap().name, "Na");
}

#[test]
fn test_neutralize_real_fixture_1hls() {
    let pdb_path = get_data_dir().join("1hls.pdb");
    let pdb = Pdb::from_file(&pdb_path, None).expect("failed to load 1hls.pdb");
    let ag = pdb.get_atomgroup(None, None).expect("failed to get ag");
    assert_eq!(ag.get_number_of_all_atoms(), 782);

    let neut = Neutralize::new(&ag).unwrap();
    let neut_obj = neut.neutralized();
    assert_eq!(neut_obj.get_number_of_all_atoms(), 792);

    let m1 = neut_obj.get_group("model_1").expect("model_1 not found");

    // Expected ions and positions matching Python's output
    let expected_ions: Vec<(&str, &str, &str, &str, [f64; 3])> = vec![
        ("A", "1", "Cl0", "Cl", [-1.540482, 10.924491, 4.901284]),
        ("A", "4", "60_Na0", "Na", [3.254284, 13.207381, 8.759200]),
        (
            "A",
            "17",
            "253_Na0",
            "Na",
            [-1.673550, -12.096908, 7.090933],
        ),
        ("A", "21", "Na0", "Na", [-12.418245, -9.681051, 4.127976]),
        ("B", "1", "Cl0", "Cl", [5.662554, -8.224111, -7.883454]),
        (
            "B",
            "13",
            "514_Na0",
            "Na",
            [0.546786, -1.892307, -11.894061],
        ),
        (
            "B",
            "21",
            "627_Na0",
            "Na",
            [-14.276474, -8.451312, -11.152788],
        ),
        (
            "B",
            "22",
            "651_Cl0",
            "Cl",
            [-17.195667, -11.960994, 0.848280],
        ),
        (
            "B",
            "29",
            "769_Cl0",
            "Cl",
            [-7.185101, 19.671047, -1.362223],
        ),
        ("B", "30", "Na0", "Na", [-5.508727, 12.052101, -6.815907]),
    ];

    for (chain_key, res_key, ion_key, ion_name, expected_pos) in expected_ions {
        let chain = m1
            .get_group(chain_key)
            .unwrap_or_else(|| panic!("Chain {} not found", chain_key));
        let res = chain
            .get_group(res_key)
            .unwrap_or_else(|| panic!("Residue {}/{} not found", chain_key, res_key));
        let atom = res
            .get_atom(ion_key)
            .unwrap_or_else(|| panic!("Atom {} not found in {}/{}", ion_key, chain_key, res_key));
        assert_eq!(atom.name, ion_name);
        let pos = &atom.xyz;
        assert!(
            (pos.x - expected_pos[0]).abs() < 1e-4,
            "x mismatch for {}/{}/{}: expected {}, got {}",
            chain_key,
            res_key,
            ion_key,
            expected_pos[0],
            pos.x
        );
        assert!(
            (pos.y - expected_pos[1]).abs() < 1e-4,
            "y mismatch for {}/{}/{}: expected {}, got {}",
            chain_key,
            res_key,
            ion_key,
            expected_pos[1],
            pos.y
        );
        assert!(
            (pos.z - expected_pos[2]).abs() < 1e-4,
            "z mismatch for {}/{}/{}: expected {}, got {}",
            chain_key,
            res_key,
            ion_key,
            expected_pos[2],
            pos.z
        );
    }
}

#[test]
fn test_neutralize_real_fixture_2mgo() {
    let pdb_path = get_data_dir().join("2MGO.pdb");
    let pdb = Pdb::from_file(&pdb_path, None).expect("failed to load 2MGO.pdb");
    let ag = pdb.get_atomgroup(None, None).expect("failed to get ag");
    assert_eq!(ag.get_number_of_all_atoms(), 2680);

    let neut = Neutralize::new(&ag).unwrap();
    let neut_obj = neut.neutralized();
    // 20 models, each adding 2 ions (40 ions total)
    assert_eq!(neut_obj.get_number_of_all_atoms(), 2720);
}
