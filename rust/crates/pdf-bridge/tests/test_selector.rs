// Copyright (C) 2014 The ProteinDF development team.
// see also AUTHORS and README if provided.
//
// This file is a part of the ProteinDF software package.
//
// The ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

use pdf_bridge::atom::Atom;
use pdf_bridge::atom_group::AtomGroup;
use pdf_bridge::position::Position;
use pdf_bridge::selector::*;

fn build_fixture_group1() -> AtomGroup {
    let mut group1 = AtomGroup::new();

    let mut atom11 = Atom::new_with_pos("C", Position::new(0.0, 0.0, 0.0)).unwrap();
    atom11.name = "C1".to_string();
    let mut atom12 = Atom::new_with_pos("H", Position::new(1.0, 0.0, 0.0)).unwrap();
    atom12.name = "H1".to_string();
    let mut atom13 = Atom::new_with_pos("N", Position::new(2.0, 0.0, 0.0)).unwrap();
    atom13.name = "N1".to_string();
    let mut subgrp1 = AtomGroup::new();
    subgrp1.set_atom("C1", atom11);
    subgrp1.set_atom("H1", atom12);
    subgrp1.set_atom("N1", atom13);
    group1.set_group("sub1", subgrp1);

    let mut atom21 = Atom::new_with_pos("C", Position::new(0.0, 1.0, 0.0)).unwrap();
    atom21.name = "C1".to_string();
    let mut atom22 = Atom::new_with_pos("H", Position::new(0.0, 2.0, 0.0)).unwrap();
    atom22.name = "H1".to_string();
    let mut atom23 = Atom::new_with_pos("N", Position::new(0.0, 3.0, 0.0)).unwrap();
    atom23.name = "N1".to_string();
    let mut subgrp2 = AtomGroup::new();
    subgrp2.set_atom("C1", atom21);
    subgrp2.set_atom("H1", atom22);
    subgrp2.set_atom("N1", atom23);
    group1.set_group("sub2", subgrp2);

    group1
}

#[test]
fn test_range1() {
    let group1 = build_fixture_group1();
    let sel_range = Select_Range::from_str("0.0 0.0 0.0", 0.1).expect("valid range");
    let sel = group1.select(&sel_range);

    assert_eq!(sel.get_number_of_all_atoms(), 1);
    let c1 = sel["sub1"].get_atom("C1").expect("C1 in sub1");
    assert_eq!(c1.symbol().unwrap(), "C");
    assert_eq!(c1.path, "/sub1/C1");
}

#[test]
fn test_range2() {
    let group1 = build_fixture_group1();
    let sel_range = Select_Range::from_str("0.0 0.0 0.0", 1.1).expect("valid range");
    let sel = group1.select(&sel_range);

    assert_eq!(sel.get_number_of_all_atoms(), 3);
    let sub1_c1 = sel["sub1"].get_atom("C1").expect("C1 in sub1");
    assert_eq!(sub1_c1.symbol().unwrap(), "C");
    assert_eq!(sub1_c1.path, "/sub1/C1");

    let sub1_h1 = sel["sub1"].get_atom("H1").expect("H1 in sub1");
    assert_eq!(sub1_h1.path, "/sub1/H1");

    let sub2_c1 = sel["sub2"].get_atom("C1").expect("C1 in sub2");
    assert_eq!(sub2_c1.path, "/sub2/C1");
}

#[test]
fn test_select_symbol() {
    let group1 = build_fixture_group1();

    let sel_c = SelectSymbol::new("C");
    let res_c = group1.select(&sel_c);
    assert_eq!(res_c.get_number_of_all_atoms(), 2);

    // Case-insensitive
    let sel_c_lower = Select_Symbol::new("c");
    let res_c_lower = group1.select(&sel_c_lower);
    assert_eq!(res_c_lower.get_number_of_all_atoms(), 2);

    let sel_n = SelectSymbol::new("N");
    let res_n = group1.select(&sel_n);
    assert_eq!(res_n.get_number_of_all_atoms(), 2);
}

#[test]
fn test_select_name() {
    let group1 = build_fixture_group1();

    // Select by atom name
    let sel_c1 = SelectName::new("C1");
    let res = group1.select(&sel_c1);
    assert_eq!(res.get_number_of_all_atoms(), 2);

    // Trimming
    let sel_trimmed = Select_Name::new("  C1  ");
    let res_trimmed = group1.select(&sel_trimmed);
    assert_eq!(res_trimmed.get_number_of_all_atoms(), 2);

    // Select whole group by group name
    let mut root = AtomGroup::new();
    let mut g1 = build_fixture_group1();
    g1.name = "my_protein".to_string();
    root.set_group("prot", g1);

    let sel_grp = SelectName::new("my_protein");
    let res_grp = root.select(&sel_grp);
    assert_eq!(res_grp.get_number_of_all_atoms(), 6);
}

#[test]
fn test_select_path_simple() {
    let group1 = build_fixture_group1();

    let sel = SelectPathSimple::new("/sub1/H1");
    let res = group1.select(&sel);
    assert_eq!(res.get_number_of_all_atoms(), 1);
    assert_eq!(res["sub1"].get_atom("H1").unwrap().path, "/sub1/H1");

    let sel_none = Select_Path_simple::new("/nonexistent/path");
    let res_none = group1.select(&sel_none);
    assert_eq!(res_none.get_number_of_all_atoms(), 0);
}

#[test]
fn test_select_path_wildcard() {
    let group1 = build_fixture_group1();

    // Matches everything under /sub1/
    let sel_sub1 = SelectPathWildcard::new("/sub1/*").expect("valid pattern");
    let res_sub1 = group1.select(&sel_sub1);
    assert_eq!(res_sub1.get_number_of_all_atoms(), 3);

    // Matches C1 in any subgroup
    let sel_c1 = Select_Path_wildcard::new("/*/C1").expect("valid pattern");
    let res_c1 = group1.select(&sel_c1);
    assert_eq!(res_c1.get_number_of_all_atoms(), 2);

    // Deprecated Select_Path wrapper
    let sel_dep = Select_Path::new("/sub2/*", true).expect("valid pattern");
    let res_dep = group1.select(&sel_dep);
    assert_eq!(res_dep.get_number_of_all_atoms(), 3);
}

#[test]
fn test_select_path_regex() {
    let group1 = build_fixture_group1();

    // Regex matching any /sub[12]/N1
    let sel = SelectPathRegex::new(r"^/sub[12]/N1$").expect("valid regex");
    let res = group1.select(&sel);
    assert_eq!(res.get_number_of_all_atoms(), 2);

    // Partial match regex
    let sel_partial = Select_PathRegex::new("sub1").expect("valid regex");
    let res_partial = group1.select(&sel_partial);
    assert_eq!(res_partial.get_number_of_all_atoms(), 3);
}

#[test]
fn test_select_atom() {
    let group1 = build_fixture_group1();

    // Probe atom near sub1/C1 (0.0, 0.0, 0.0) with small distance
    let probe_c = Atom::new_with_pos("C", Position::new(0.05, 0.0, 0.0)).unwrap();
    let sel_close = SelectAtom::new(probe_c.clone(), 0.1);
    let res = group1.select(&sel_close);
    assert_eq!(res.get_number_of_all_atoms(), 1);
    assert_eq!(res["sub1"].get_atom("C1").unwrap().path, "/sub1/C1");

    // Same position but different element symbol (N instead of C) -> should NOT match
    let probe_n = Atom::new_with_pos("N", Position::new(0.05, 0.0, 0.0)).unwrap();
    let sel_diff_elem = Select_Atom::new(probe_n, 0.1);
    let res_diff = group1.select(&sel_diff_elem);
    assert_eq!(res_diff.get_number_of_all_atoms(), 0);
}

#[test]
fn test_select_atomgroup() {
    let group1 = build_fixture_group1();

    // Reference group containing 2 atoms matching sub1/C1 and sub2/N1
    let mut ref_ag = AtomGroup::new();
    ref_ag.set_atom(
        "ref1",
        Atom::new_with_pos("C", Position::new(0.0, 0.0, 0.0)).unwrap(),
    );
    ref_ag.set_atom(
        "ref2",
        Atom::new_with_pos("N", Position::new(0.0, 3.0, 0.0)).unwrap(),
    );

    let sel = SelectAtomGroup::with_default_range(&ref_ag);
    let res = group1.select(&sel);
    assert_eq!(res.get_number_of_all_atoms(), 2);
    assert!(res["sub1"].has_atom("C1"));
    assert!(res["sub2"].has_atom("N1"));
}

#[test]
fn test_error_handling() {
    let err_regex = SelectPathRegex::new("(?P<invalid");
    assert!(err_regex.is_err());

    let err_pos = SelectRange::from_str("invalid position string", 1.0);
    assert!(err_pos.is_err());
}
