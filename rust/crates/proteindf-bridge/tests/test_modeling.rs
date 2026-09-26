// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::path::PathBuf;

use proteindf_bridge::atom::Atom;
use proteindf_bridge::atom_group::AtomGroup;
use proteindf_bridge::brd::load_atomgroup;
use proteindf_bridge::modeling::Modeling;
use proteindf_bridge::position::Position;

fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

fn make_atom(name: &str, symbol: &str, pos: Position) -> Atom {
    let mut atom = Atom::new_with_pos(symbol, pos).unwrap();
    atom.name = name.to_string();
    atom
}

#[test]
fn test_modeling_initialization() {
    let m = Modeling::new().expect("failed to create Modeling instance");
    for conf in Modeling::CONFORMERS {
        assert!(m.get_ace_ala_nme(conf).is_some());
        let ag = m.get_ace_ala_nme(conf).unwrap();
        assert_eq!(ag.get_number_of_groups(), 3);
    }

    let m_dir = Modeling::from_data_dir(test_data_dir()).expect("failed from_data_dir");
    for conf in Modeling::CONFORMERS {
        assert!(m_dir.get_ace_ala_nme(conf).is_some());
    }
}

#[test]
fn test_get_ace_and_nme_with_fixtures() {
    let m = Modeling::new().unwrap();

    for conf in Modeling::CONFORMERS {
        let path = test_data_dir().join(format!("ACE_ALA_NME_{conf}.brd"));
        let ag = load_atomgroup(&path).unwrap();
        let res2 = ag.get_group("2").unwrap();
        let res3 = ag.get_group("3").unwrap();

        // 1. With next_aa
        let ace = m.get_ACE(res2, Some(res3)).unwrap();
        assert_eq!(ace.get_number_of_atoms(), 6);
        assert_eq!(ace.path(), "/ACE/");

        let nme = m.get_NME(res2, Some(res3)).unwrap();
        assert_eq!(nme.get_number_of_atoms(), 6);
        assert_eq!(nme.path(), "/NME/");

        // 2. Without next_aa
        let ace_none = m.get_ACE(res2, None).unwrap();
        assert_eq!(ace_none.get_number_of_atoms(), 6);

        let nme_none = m.get_NME(res2, None).unwrap();
        assert_eq!(nme_none.get_number_of_atoms(), 6);
    }

    // Exact coordinate checks for trans1 against Python ProteinDF_bridge values
    let path = test_data_dir().join("ACE_ALA_NME_trans1.brd");
    let ag = load_atomgroup(&path).unwrap();
    let res2 = ag.get_group("2").unwrap();
    let res3 = ag.get_group("3").unwrap();

    let ace = m.get_ACE(res2, Some(res3)).unwrap();
    let ch3 = ace.get_atom("CH3").unwrap();
    assert!((ch3.xyz.x - (-3.3718)).abs() < 1e-3);
    assert!((ch3.xyz.y - (-1.1446)).abs() < 1e-3);
    assert!((ch3.xyz.z - 0.0463).abs() < 1e-3);

    let nme = m.get_NME(res2, Some(res3)).unwrap();
    let n = nme.get_atom("N").unwrap();
    assert!((n.xyz.x - 2.4756).abs() < 1e-3);
    assert!((n.xyz.y - 0.3160).abs() < 1e-3);
    assert!((n.xyz.z - (-0.4246)).abs() < 1e-3);
}

#[test]
fn test_get_ace_simple_and_nme_simple() {
    let m = Modeling::new().unwrap();

    let mut next_aa = AtomGroup::with_name("ALA");
    next_aa.set_atom(
        "CA",
        Atom::new_with_pos("C", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    next_aa.set_atom(
        "C",
        Atom::new_with_pos("C", Position::new(1.0, 0.0, 0.0)).unwrap(),
    );
    next_aa.set_atom(
        "O",
        Atom::new_with_pos("O", Position::new(1.0, 1.0, 0.0)).unwrap(),
    );
    next_aa.set_atom(
        "N",
        Atom::new_with_pos("N", Position::new(-1.0, 0.0, 0.0)).unwrap(),
    );
    next_aa.set_atom(
        "H",
        Atom::new_with_pos("H", Position::new(-1.0, 1.0, 0.0)).unwrap(),
    );

    let ace = m.get_ACE_simple(&next_aa).unwrap();
    assert_eq!(ace.get_number_of_atoms(), 6); // CA, C, O, H11, H12, H13
    assert!(ace.has_atom("CA"));
    assert!(ace.has_atom("C"));
    assert!(ace.has_atom("O"));
    assert!(ace.has_atom("H11"));
    assert!(ace.has_atom("H12"));
    assert!(ace.has_atom("H13"));

    let nme = m.get_NME_simple(&next_aa).unwrap();
    assert_eq!(nme.get_number_of_atoms(), 6); // CA, N, H, H11, H12, H13
    assert!(nme.has_atom("CA"));
    assert!(nme.has_atom("N"));
    assert!(nme.has_atom("H"));

    // Proline CD fallback in NME_simple
    let mut pro = AtomGroup::with_name("PRO");
    pro.set_atom(
        "CA",
        Atom::new_with_pos("C", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    pro.set_atom(
        "N",
        Atom::new_with_pos("N", Position::new(-1.0, 0.0, 0.0)).unwrap(),
    );
    pro.set_atom(
        "CD",
        Atom::new_with_pos("C", Position::new(-1.0, 1.0, 0.0)).unwrap(),
    );
    let nme_pro = m.get_NME_simple(&pro).unwrap();
    assert_eq!(nme_pro.get_number_of_atoms(), 6);
    assert_eq!(nme_pro.get_atom("H").unwrap().symbol().unwrap(), "H");
}

#[test]
fn test_neutralize_glu_asp_lys() {
    let m = Modeling::new().unwrap();

    // GLU
    let mut glu = AtomGroup::with_name("GLU");
    glu.set_atom("1_CD", make_atom("CD", "C", Position::new(0.0, 0.0, 0.0)));
    glu.set_atom("2_OE1", make_atom("OE1", "O", Position::new(1.0, 1.0, 0.0)));
    glu.set_atom(
        "3_OE2",
        make_atom("OE2", "O", Position::new(-1.0, 1.0, 0.0)),
    );

    let na = m.neutralize_GLU(&glu).unwrap();
    assert_eq!(na.get_number_of_atoms(), 1);
    let ion = na.get_atom("4_Na").unwrap();
    assert_eq!(ion.symbol().unwrap(), "Na");
    assert!((ion.xyz.x - 0.0).abs() < 1e-5);
    assert!((ion.xyz.y - 2.521).abs() < 1e-5);
    assert!((ion.xyz.z - 0.0).abs() < 1e-5);

    // ASP
    let mut asp = AtomGroup::with_name("ASP");
    asp.set_atom("1_CG", make_atom("CG", "C", Position::new(0.0, 0.0, 0.0)));
    asp.set_atom("2_OD1", make_atom("OD1", "O", Position::new(1.0, 1.0, 0.0)));
    asp.set_atom(
        "3_OD2",
        make_atom("OD2", "O", Position::new(-1.0, 1.0, 0.0)),
    );

    let na_asp = m.neutralize_ASP(&asp).unwrap();
    assert_eq!(na_asp.get_number_of_atoms(), 1);
    let ion_asp = na_asp.get_atom("4_Na").unwrap();
    assert!((ion_asp.xyz.y - 2.521).abs() < 1e-5);

    // LYS
    let mut lys = AtomGroup::with_name("LYS");
    lys.set_atom("1_NZ", make_atom("NZ", "N", Position::new(0.0, 0.0, 0.0)));
    lys.set_atom("2_HZ1", make_atom("HZ1", "H", Position::new(1.0, 0.0, 0.0)));
    lys.set_atom("3_HZ2", make_atom("HZ2", "H", Position::new(0.0, 1.0, 0.0)));
    lys.set_atom("4_HZ3", make_atom("HZ3", "H", Position::new(0.0, 0.0, 1.0)));

    let cl = m.neutralize_LYS(&lys).unwrap();
    assert_eq!(cl.get_number_of_atoms(), 1);
    let ion_cl = cl.get_atom("5_Cl").unwrap();
    assert_eq!(ion_cl.symbol().unwrap(), "Cl");
    let expected = 1.840015;
    assert!((ion_cl.xyz.x - expected).abs() < 1e-4);
    assert!((ion_cl.xyz.y - expected).abs() < 1e-4);
    assert!((ion_cl.xyz.z - expected).abs() < 1e-4);
}

#[test]
fn test_neutralize_arg_cases() {
    let m = Modeling::new().unwrap();

    let mut arg = AtomGroup::with_name("ARG");
    arg.set_atom("1_CZ", make_atom("CZ", "C", Position::new(0.0, 0.0, 0.0)));
    arg.set_atom("2_NH1", make_atom("NH1", "N", Position::new(1.0, 0.0, 0.0)));
    arg.set_atom("3_NH2", make_atom("NH2", "N", Position::new(0.0, 1.0, 0.0)));
    arg.set_atom(
        "4_HH11",
        make_atom("HH11", "H", Position::new(1.5, 0.5, 0.0)),
    );
    arg.set_atom(
        "5_HH12",
        make_atom("HH12", "H", Position::new(1.5, -0.5, 0.0)),
    );
    arg.set_atom(
        "6_HH21",
        make_atom("HH21", "H", Position::new(0.5, 1.5, 0.0)),
    );
    arg.set_atom(
        "7_HH22",
        make_atom("HH22", "H", Position::new(-0.5, 1.5, 0.0)),
    );

    // Case 0: center
    let cl0 = m.neutralize_ARG(&arg, 0).unwrap();
    let ion0 = cl0.get_atom("8_Cl").unwrap();
    let exp0 = 2.121320;
    assert!((ion0.xyz.x - exp0).abs() < 1e-4);
    assert!((ion0.xyz.y - exp0).abs() < 1e-4);
    assert!((ion0.xyz.z - 0.0).abs() < 1e-4);

    // Case 1: NH1 side
    let cl1 = m.neutralize_ARG(&arg, 1).unwrap();
    let ion1 = cl1.get_atom("8_Cl").unwrap();
    assert!((ion1.xyz.x - 3.0).abs() < 1e-4);
    assert!((ion1.xyz.y - 0.0).abs() < 1e-4);

    // Case 2: NH2 side
    let cl2 = m.neutralize_ARG(&arg, 2).unwrap();
    let ion2 = cl2.get_atom("8_Cl").unwrap();
    assert!((ion2.xyz.x - 0.0).abs() < 1e-4);
    assert!((ion2.xyz.y - 3.0).abs() < 1e-4);
}

#[test]
fn test_neutralize_nterm_and_cterm() {
    let m = Modeling::new().unwrap();

    // Normal N-term
    let mut ala = AtomGroup::with_name("ALA");
    ala.set_atom(
        "N",
        Atom::new_with_pos("N", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    ala.set_atom(
        "H1",
        Atom::new_with_pos("H", Position::new(1.0, 0.0, 0.0)).unwrap(),
    );
    ala.set_atom(
        "H2",
        Atom::new_with_pos("H", Position::new(0.0, 1.0, 0.0)).unwrap(),
    );
    ala.set_atom(
        "H3",
        Atom::new_with_pos("H", Position::new(0.0, 0.0, 1.0)).unwrap(),
    );
    let cl_nterm = m.neutralize_Nterm(&ala).unwrap();
    assert_eq!(cl_nterm.get_number_of_atoms(), 1);
    assert_eq!(cl_nterm.get_atom("Cl").unwrap().symbol().unwrap(), "Cl");

    // PRO N-term
    let mut pro = AtomGroup::with_name("PRO");
    pro.set_atom(
        "N",
        Atom::new_with_pos("N", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    pro.set_atom(
        "H2",
        Atom::new_with_pos("H", Position::new(1.0, 0.0, 0.0)).unwrap(),
    );
    pro.set_atom(
        "HXT",
        Atom::new_with_pos("H", Position::new(0.0, 1.0, 0.0)).unwrap(),
    );
    let cl_pro = m.neutralize_Nterm(&pro).unwrap();
    assert_eq!(cl_pro.get_number_of_atoms(), 1);

    // C-term
    let mut cterm = AtomGroup::with_name("ALA");
    cterm.set_atom(
        "C",
        Atom::new_with_pos("C", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    cterm.set_atom(
        "O",
        Atom::new_with_pos("O", Position::new(1.0, 0.0, 0.0)).unwrap(),
    );
    cterm.set_atom(
        "OXT",
        Atom::new_with_pos("O", Position::new(0.0, 1.0, 0.0)).unwrap(),
    );
    let na_cterm = m.neutralize_Cterm(&cterm).unwrap();
    assert_eq!(na_cterm.get_number_of_atoms(), 1);
    assert_eq!(na_cterm.get_atom("Na").unwrap().symbol().unwrap(), "Na");
}

#[test]
fn test_neutralize_fad_success_and_error() {
    let m = Modeling::new().unwrap();

    let mut fad = AtomGroup::with_name("FAD");
    fad.set_atom("1_P", make_atom("P", "P", Position::new(0.0, 0.0, 0.0)));
    fad.set_atom("2_O1P", make_atom("O1P", "O", Position::new(1.0, 0.0, 0.0)));
    fad.set_atom("3_O2P", make_atom("O2P", "O", Position::new(0.0, 1.0, 0.0)));
    fad.set_atom("4_PA", make_atom("PA", "P", Position::new(2.0, 0.0, 0.0)));
    fad.set_atom("5_O1A", make_atom("O1A", "O", Position::new(3.0, 0.0, 0.0)));
    fad.set_atom("6_O2A", make_atom("O2A", "O", Position::new(2.0, 1.0, 0.0)));

    let na_fad = m.neutralize_FAD(&fad).unwrap();
    assert_eq!(na_fad.get_number_of_atoms(), 2);
    let na1 = na_fad.get_atom("7_Na1").unwrap();
    let na2 = na_fad.get_atom("7_Na2").unwrap();
    assert!((na1.xyz.x - 1.9431).abs() < 1e-3);
    assert!((na1.xyz.y - 1.9431).abs() < 1e-3);
    assert!((na2.xyz.x - 3.9431).abs() < 1e-3);
    assert!((na2.xyz.y - 1.9431).abs() < 1e-3);

    // Missing O1P / OP1 raises error
    let mut broken_fad = AtomGroup::with_name("FAD");
    broken_fad.set_atom("1_P", make_atom("P", "P", Position::new(0.0, 0.0, 0.0)));
    broken_fad.set_atom("3_O2P", make_atom("O2P", "O", Position::new(0.0, 1.0, 0.0)));
    assert!(m.neutralize_FAD(&broken_fad).is_err());
}

#[test]
fn test_geometry_helpers() {
    let m = Modeling::new().unwrap();

    // arbitary_rotate_matrix: aligns vector a=(1,0,0) towards b=(0,1,0),
    // matching Python ProteinDF_bridge result: rotated a = (0.0, -1.0, 0.0)
    let a = Position::new(1.0, 0.0, 0.0);
    let b = Position::new(0.0, 1.0, 0.0);
    let rot = m.arbitary_rotate_matrix(a, b).unwrap();
    let mut test_p = a;
    test_p.rotate(&rot).unwrap();
    assert!((test_p.x - 0.0).abs() < 1e-5);
    assert!((test_p.y - (-1.0)).abs() < 1e-5);
    assert!((test_p.z - 0.0).abs() < 1e-5);

    // get_NH3
    let nh3 = m.get_NH3(std::f64::consts::FRAC_PI_2, 1.0).unwrap();
    assert_eq!(nh3.get_number_of_atoms(), 4);
    assert!(nh3.has_atom("N"));
    assert!(nh3.has_atom("H1"));
    assert!(nh3.has_atom("H2"));
    assert!(nh3.has_atom("H3"));

    // select_residues
    let mut chain = AtomGroup::with_name("A");
    for i in 1..=5 {
        let mut res = AtomGroup::with_name("ALA");
        res.set_atom(
            "CA",
            Atom::new_with_pos("C", Position::new(i as f64, 0.0, 0.0)).unwrap(),
        );
        chain.set_group(&i.to_string(), res);
    }
    let selected = m.select_residues(&chain, 2, 4);
    assert_eq!(selected.get_number_of_groups(), 0);
    assert_eq!(selected.get_number_of_atoms(), 1);
    assert_eq!(selected.get_atom("CA").unwrap().xyz.x, 4.0);

    // get_last_index
    let mut res = AtomGroup::new();
    res.set_atom(
        "12_CA",
        Atom::new_with_pos("C", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    res.set_atom(
        "42_CB",
        Atom::new_with_pos("C", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    res.set_atom(
        "3_N",
        Atom::new_with_pos("N", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    assert_eq!(m.get_last_index(&res), 42);
}

#[test]
fn test_get_ace_and_nme_resilient_to_conformer_failure() {
    let m = Modeling::new().unwrap();

    // Construct a minimal residue where exactly 3 atoms match (N, CA, C).
    // The atoms are positioned non-collinearly so valid superposition can occur,
    // but without CB, HA, O, ensuring only 3 points are matched.
    // In some reference conformers or geometries, if 3 points are collinear/degenerate,
    // Superposer::new will return Err. The loop in get_ACE / get_NME should continue
    // and successfully find a matching conformer rather than terminating early with `?`.
    let mut min_res = AtomGroup::with_name("GLY");
    min_res.set_atom("N", make_atom("N", "N", Position::new(-1.0, 1.0, 0.0)));
    min_res.set_atom("CA", make_atom("CA", "C", Position::new(0.0, 0.0, 0.0)));
    min_res.set_atom("C", make_atom("C", "C", Position::new(1.0, 1.0, 0.0)));

    let ace = m.get_ACE(&min_res, None);
    assert!(
        ace.is_ok(),
        "get_ACE should succeed when valid conformers exist: {:?}",
        ace.err()
    );
    let ace_grp = ace.unwrap();
    assert_eq!(ace_grp.get_number_of_atoms(), 6);

    let nme = m.get_NME(&min_res, None);
    assert!(
        nme.is_ok(),
        "get_NME should succeed when valid conformers exist: {:?}",
        nme.err()
    );
    let nme_grp = nme.unwrap();
    assert_eq!(nme_grp.get_number_of_atoms(), 6);
}

#[test]
fn test_get_ace_and_nme_all_conformers_fail_error_context() {
    let m = Modeling::new().unwrap();

    // Construct a residue where atoms are collinear (N, CA, C on a single line)
    // and fewer than 3 non-collinear atoms exist. Superposer::new will fail with
    // a collinearity error on all conformers.
    let mut collinear_res = AtomGroup::with_name("GLY");
    collinear_res.set_atom("N", make_atom("N", "N", Position::new(0.0, 0.0, 0.0)));
    collinear_res.set_atom("CA", make_atom("CA", "C", Position::new(1.0, 0.0, 0.0)));
    collinear_res.set_atom("C", make_atom("C", "C", Position::new(2.0, 0.0, 0.0)));

    let ace_res = m.get_ACE(&collinear_res, None);
    assert!(
        ace_res.is_err(),
        "Must fail when all conformers are collinear"
    );
    let ace_err = ace_res.err().unwrap().to_string();
    assert!(
        ace_err.contains("no matching ACE conformer found"),
        "Error should state no matching conformer: {ace_err}"
    );
    assert!(
        ace_err.contains("collinear") || ace_err.contains("trans1"),
        "Error must contain underlying failure details: {ace_err}"
    );

    let nme_res = m.get_NME(&collinear_res, None);
    assert!(
        nme_res.is_err(),
        "Must fail when all conformers are collinear"
    );
    let nme_err = nme_res.err().unwrap().to_string();
    assert!(
        nme_err.contains("no matching NME conformer found"),
        "Error should state no matching conformer: {nme_err}"
    );
    assert!(
        nme_err.contains("collinear") || nme_err.contains("trans1"),
        "Error must contain underlying failure details: {nme_err}"
    );
}
