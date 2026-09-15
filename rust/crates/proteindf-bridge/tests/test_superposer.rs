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

use proteindf_bridge::atom::Atom;
use proteindf_bridge::atom_group::AtomGroup;
use proteindf_bridge::format::Pdb;
use proteindf_bridge::matrix::Matrix;
use proteindf_bridge::position::Position;
use proteindf_bridge::superposer::Superposer;
use std::path::PathBuf;

fn make_atom(name: &str, xyz: Position) -> Atom {
    let mut a = Atom::new_with_pos("C", xyz).unwrap();
    a.name = name.to_string();
    a
}

fn create_test_atom_groups() -> (AtomGroup, AtomGroup) {
    let mut ag1 = AtomGroup::new();
    ag1.name = "mol1".to_string();
    ag1.set_atom("A1", make_atom("A1", Position::new(1.0, 1.0, 1.0)));
    ag1.set_atom("A2", make_atom("A2", Position::new(1.0, -1.0, -1.0)));
    ag1.set_atom("A3", make_atom("A3", Position::new(-1.0, 1.0, -1.0)));
    ag1.set_atom("A4", make_atom("A4", Position::new(-1.0, -1.0, 1.0)));

    // Translated: (x+1.0, y+2.0, z+3.0)
    let mut ag2 = AtomGroup::new();
    ag2.name = "mol2".to_string();
    ag2.set_atom("A1", make_atom("A1", Position::new(2.0, 3.0, 4.0)));
    ag2.set_atom("A2", make_atom("A2", Position::new(2.0, 1.0, 2.0)));
    ag2.set_atom("A3", make_atom("A3", Position::new(0.0, 3.0, 2.0)));
    ag2.set_atom("A4", make_atom("A4", Position::new(0.0, 1.0, 4.0)));

    (ag1, ag2)
}

#[test]
fn test_superposer_rmsd_and_centers() {
    let (ag1, ag2) = create_test_atom_groups();
    let sp = Superposer::new(&ag1, &ag2).unwrap();

    assert_eq!(sp.num_of_positions(), 4);
    assert!((sp.rmsd() - 0.0).abs() < 1e-5);

    let c1 = sp.center1();
    assert!((c1.x - 0.0).abs() < 1e-6);
    assert!((c1.y - 0.0).abs() < 1e-6);
    assert!((c1.z - 0.0).abs() < 1e-6);

    let c2 = sp.center2();
    assert!((c2.x - 1.0).abs() < 1e-6);
    assert!((c2.y - 2.0).abs() < 1e-6);
    assert!((c2.z - 3.0).abs() < 1e-6);

    // Rotation matrix should be Identity (within numerical precision)
    let rot = sp.rotation_mat();
    for r in 0..3 {
        for c in 0..3 {
            let expected = if r == c { 1.0 } else { 0.0 };
            assert!(
                (rot.get(r, c).unwrap() - expected).abs() < 1e-6,
                "rot({}, {}) = {}, expected {}",
                r,
                c,
                rot.get(r, c).unwrap(),
                expected
            );
        }
    }
}

#[test]
fn test_superposer_superimpose() {
    let (ag1, ag2) = create_test_atom_groups();
    let sp = Superposer::new(&ag1, &ag2).unwrap();

    let superimposed = sp.superimpose(&ag1).unwrap();
    assert_eq!(superimposed.get_number_of_atoms(), 4);

    for key in &["A1", "A2", "A3", "A4"] {
        let p1 = superimposed.get_atom(key).unwrap().xyz;
        let p2 = ag2.get_atom(key).unwrap().xyz;
        let dist = p1.distance_from(&p2);
        assert!(dist < 1e-4, "dist for {} was {}", key, dist);
    }
}

#[test]
fn test_superposer_rotation() {
    let (ag1, _) = create_test_atom_groups();

    // Rotate ag1 90 degrees around Z: (x, y, z) -> (-y, x, z) + (5.0, 2.0, -3.0)
    let mut ag2 = AtomGroup::new();
    ag2.name = "mol2".to_string();
    ag2.set_atom(
        "A1",
        make_atom("A1", Position::new(-1.0 + 5.0, 1.0 + 2.0, 1.0 - 3.0)),
    );
    ag2.set_atom(
        "A2",
        make_atom("A2", Position::new(1.0 + 5.0, 1.0 + 2.0, -1.0 - 3.0)),
    );
    ag2.set_atom(
        "A3",
        make_atom("A3", Position::new(-1.0 + 5.0, -1.0 + 2.0, -1.0 - 3.0)),
    );
    ag2.set_atom(
        "A4",
        make_atom("A4", Position::new(1.0 + 5.0, -1.0 + 2.0, 1.0 - 3.0)),
    );

    let sp = Superposer::new(&ag1, &ag2).unwrap();
    assert!(sp.rmsd() < 1e-10);

    // Rotation matrix should be:
    // [ 0 -1  0]
    // [ 1  0  0]
    // [ 0  0  1]
    let rot = sp.rotation_mat();
    let expected = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
    for (r, row) in expected.iter().enumerate() {
        for (c, &exp_val) in row.iter().enumerate() {
            assert!(
                (rot.get(r, c).unwrap() - exp_val).abs() < 1e-6,
                "rot({}, {}) = {}, expected {}",
                r,
                c,
                rot.get(r, c).unwrap(),
                exp_val
            );
        }
    }

    let superimposed = sp.superimpose(&ag1).unwrap();
    for key in &["A1", "A2", "A3", "A4"] {
        let p1 = superimposed.get_atom(key).unwrap().xyz;
        let p2 = ag2.get_atom(key).unwrap().xyz;
        assert!(p1.distance_from(&p2) < 1e-4);
    }
}

#[test]
fn test_superposer_with_real_pdb_1hls() {
    let pdb_path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/1hls.pdb");
    let pdb = Pdb::from_file(&pdb_path, None).unwrap();
    let ag1 = pdb.get_atomgroup(None, None).unwrap();

    // 1. Translation test: shift by (10.0, -5.0, 3.0)
    let mut ag2 = ag1.clone();
    ag2.shift_by(Position::new(10.0, -5.0, 3.0));

    let sp = Superposer::new(&ag1, &ag2).unwrap();
    assert_eq!(sp.num_of_positions(), 782);
    assert!(sp.rmsd() < 1e-10);

    let rot = sp.rotation_mat();
    for r in 0..3 {
        for c in 0..3 {
            let expected = if r == c { 1.0 } else { 0.0 };
            assert!((rot.get(r, c).unwrap() - expected).abs() < 1e-6);
        }
    }

    // 2. Rotation test: rotate by 90 degrees around Z and shift
    let mut rot_z = Matrix::new(3, 3);
    rot_z.set(0, 1, -1.0);
    rot_z.set(1, 0, 1.0);
    rot_z.set(2, 2, 1.0);

    let mut ag3 = ag1.clone();
    ag3.rotate(&rot_z).unwrap();
    ag3.shift_by(Position::new(5.0, -3.0, 2.0));

    let sp_rot = Superposer::new(&ag1, &ag3).unwrap();
    assert_eq!(sp_rot.num_of_positions(), 782);
    assert!(sp_rot.rmsd() < 1e-10);

    let rot3 = sp_rot.rotation_mat();
    let expected = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
    for (r, row) in expected.iter().enumerate() {
        for (c, &exp_val) in row.iter().enumerate() {
            assert!((rot3.get(r, c).unwrap() - exp_val).abs() < 1e-6);
        }
    }
}

#[test]
fn test_superposer_no_common_atoms_error() {
    let mut ag1 = AtomGroup::new();
    ag1.set_atom("A1", make_atom("A1", Position::new(1.0, 1.0, 1.0)));

    let mut ag2 = AtomGroup::new();
    ag2.set_atom("B1", make_atom("B1", Position::new(2.0, 2.0, 2.0)));

    let result = Superposer::new(&ag1, &ag2);
    assert!(result.is_err());
}
