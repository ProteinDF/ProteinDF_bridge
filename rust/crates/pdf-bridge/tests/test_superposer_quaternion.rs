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
use pdf_bridge::format::Pdb;
use pdf_bridge::position::Position;
use pdf_bridge::superposer::Superposer;
use pdf_bridge::{SuperposerQuaternion, Superposer_quaternion};
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
fn test_quaternion_rmsd_and_calc() {
    let (ag1, ag2) = create_test_atom_groups();
    let sq = SuperposerQuaternion::new(&ag1, &ag2).unwrap();

    let rmsd = sq.rmsd();
    assert!((rmsd - 0.0).abs() < 1e-5);
    assert!((sq.calc() - rmsd).abs() < 1e-12);
}

#[test]
fn test_quaternion_rotation_mat() {
    let (ag1, ag2) = create_test_atom_groups();
    // Test alias Superposer_quaternion
    let sq = Superposer_quaternion::new(&ag1, &ag2).unwrap();

    let rot_mat = sq.rotation_mat();
    assert_eq!(rot_mat.rows(), 3);
    assert_eq!(rot_mat.cols(), 3);

    for r in 0..3 {
        for c in 0..3 {
            let expected = if r == c { 1.0 } else { 0.0 };
            assert!(
                (rot_mat.get(r, c).unwrap() - expected).abs() < 1e-6,
                "rot({}, {}) = {}, expected {}",
                r,
                c,
                rot_mat.get(r, c).unwrap(),
                expected
            );
        }
    }
}

#[test]
fn test_quaternion_superimpose() {
    let (ag1, ag2) = create_test_atom_groups();
    let sq = SuperposerQuaternion::new(&ag1, &ag2).unwrap();

    let superimposed = sq.superimpose(&ag1).unwrap();
    assert_eq!(superimposed.get_number_of_atoms(), 4);

    for key in &["A1", "A2", "A3", "A4"] {
        let p1 = superimposed.get_atom(key).unwrap().xyz;
        let p2 = ag2.get_atom(key).unwrap().xyz;
        let dist = p1.distance_from(&p2);
        assert!(dist < 1e-4, "dist for {} was {}", key, dist);
    }
}

#[test]
fn test_quaternion_with_rotation_and_agreement_with_kabsch() {
    // Non-symmetric positions with 90-degree rotation around X + translation
    let mut ag1 = AtomGroup::new();
    ag1.set_atom("A1", make_atom("A1", Position::new(1.2, 2.3, 3.4)));
    ag1.set_atom("A2", make_atom("A2", Position::new(4.5, 1.1, 0.2)));
    ag1.set_atom("A3", make_atom("A3", Position::new(0.1, 5.6, 2.7)));
    ag1.set_atom("A4", make_atom("A4", Position::new(3.3, 0.4, 6.1)));

    // Rotate 90 deg around X: (x, y, z) -> (x, -z, y) + (2.0, 3.0, 4.0)
    let mut ag2 = AtomGroup::new();
    ag2.set_atom(
        "A1",
        make_atom("A1", Position::new(1.2 + 2.0, -3.4 + 3.0, 2.3 + 4.0)),
    );
    ag2.set_atom(
        "A2",
        make_atom("A2", Position::new(4.5 + 2.0, -0.2 + 3.0, 1.1 + 4.0)),
    );
    ag2.set_atom(
        "A3",
        make_atom("A3", Position::new(0.1 + 2.0, -2.7 + 3.0, 5.6 + 4.0)),
    );
    ag2.set_atom(
        "A4",
        make_atom("A4", Position::new(3.3 + 2.0, -6.1 + 3.0, 0.4 + 4.0)),
    );

    let sq = SuperposerQuaternion::new(&ag1, &ag2).unwrap();
    let sp = Superposer::new(&ag1, &ag2).unwrap();

    // Both methods should achieve RMSD near 0.0
    assert!(sq.rmsd() < 1e-10, "sq.rmsd was {}", sq.rmsd());
    assert!(sp.rmsd() < 1e-10, "sp.rmsd was {}", sp.rmsd());
    assert!(
        (sq.rmsd() - sp.rmsd()).abs() < 1e-10,
        "RMSD difference between quaternion and Kabsch: {}",
        (sq.rmsd() - sp.rmsd()).abs()
    );

    // Rotation matrices should agree
    let rot_q = sq.rotation_mat();
    let rot_k = sp.rotation_mat();
    for r in 0..3 {
        for c in 0..3 {
            assert!(
                (rot_q.get(r, c).unwrap() - rot_k.get(r, c).unwrap()).abs() < 1e-6,
                "Rotation matrix mismatch at ({}, {}): quaternion={}, kabsch={}",
                r,
                c,
                rot_q.get(r, c).unwrap(),
                rot_k.get(r, c).unwrap()
            );
        }
    }
}

#[test]
fn test_quaternion_with_real_pdb_1hls() {
    let pdb_path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/1hls.pdb");
    let pdb = Pdb::from_file(&pdb_path, None).unwrap();
    let ag1 = pdb.get_atomgroup(None, None).unwrap();

    let mut ag2 = ag1.clone();
    ag2.shift_by(Position::new(10.0, -5.0, 3.0));

    let sq = SuperposerQuaternion::new(&ag1, &ag2).unwrap();
    assert_eq!(sq.num_of_positions(), 782);
    assert!(sq.rmsd() < 1e-10);

    let rot = sq.rotation_mat();
    for r in 0..3 {
        for c in 0..3 {
            let expected = if r == c { 1.0 } else { 0.0 };
            assert!((rot.get(r, c).unwrap() - expected).abs() < 1e-6);
        }
    }
}

#[test]
fn test_quaternion_no_common_atoms_error() {
    let mut ag1 = AtomGroup::new();
    ag1.set_atom("A1", make_atom("A1", Position::new(1.0, 1.0, 1.0)));

    let mut ag2 = AtomGroup::new();
    ag2.set_atom("B1", make_atom("B1", Position::new(2.0, 2.0, 2.0)));

    let result = SuperposerQuaternion::new(&ag1, &ag2);
    assert!(result.is_err());
}
