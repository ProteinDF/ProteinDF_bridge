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

use crate::atom_group::AtomGroup;
use crate::error::{BridgeError, Result};
use crate::matrix::{Matrix, SymmetricMatrix};
use crate::position::Position;
use crate::vector::Vector;

/// Structure for superposing two AtomGroups using Kearsley's quaternion method (1989).
///
/// Ports `proteindf_bridge.superposer_quaternion.Superposer_quaternion`.
///
/// # Note on Divergence from Python
/// In Python's `matrix.py`, `SymmetricMatrix` did not override `add()`, causing
/// `B.add(row, col, value)` to add only to the upper-triangular storage (`row < col`),
/// while `numpy.linalg.eigh(..., 'L')` only read the lower-triangular entries (which remained 0.0).
/// In this Rust implementation, `SymmetricMatrix::add` properly adds to both symmetric entries,
/// allowing the quaternion algorithm to correctly compute the optimal rotation matrix and minimal
/// RMSD for arbitrary 3D rotations, matching the mathematical definition (Kearsley 1989) and
/// agreeing with the Kabsch algorithm (`Superposer`).
#[derive(Debug, Clone)]
pub struct SuperposerQuaternion {
    atomgroup1: AtomGroup,
    atomgroup2: AtomGroup,
    positions1: Vec<Position>,
    positions2: Vec<Position>,
    center1: Position,
    center2: Position,
    r_a: Vec<Position>,
    r_b: Vec<Position>,
    va: Vec<Position>,
    vb: Vec<Position>,
    mat_b: SymmetricMatrix,
    eigval: Vector,
    eigvec: Matrix,
    mat_r: Matrix,
    rmsd: f64,
}

impl SuperposerQuaternion {
    /// Creates a new `SuperposerQuaternion` from two `AtomGroup` instances and performs superimposition.
    pub fn new(atomgroup1: &AtomGroup, atomgroup2: &AtomGroup) -> Result<Self> {
        let (positions1, positions2) = Self::match_positions(atomgroup1, atomgroup2);
        let num_of_positions = positions1.len();
        if num_of_positions != positions2.len() {
            return Err(BridgeError::value_error(
                "positions",
                "Matched positions between atom groups must have equal length",
            ));
        }
        if num_of_positions == 0 {
            return Err(BridgeError::value_error(
                "positions",
                "No common atom positions found to superpose",
            ));
        }

        let center1 = Self::calc_center(&positions1);
        let center2 = Self::calc_center(&positions2);

        let r_a = Self::shift_positions(&positions1, center1);
        let r_b = Self::shift_positions(&positions2, center2);

        let va = Self::make_va(&r_a, &r_b);
        let vb = Self::make_vb(&r_a, &r_b);

        let mat_b = Self::make_b(&va, &vb);

        let (eigval, eigvec) = mat_b.eig()?;

        // Row 0 corresponds to the smallest eigenvalue, representing the optimal quaternion
        let q = eigvec.get_row_vector(0)?;
        let mat_r = Self::make_r(&q)?;

        let rmsd = Self::calc_rmsd(&r_a, &r_b, &mat_r)?;

        Ok(Self {
            atomgroup1: atomgroup1.clone(),
            atomgroup2: atomgroup2.clone(),
            positions1,
            positions2,
            center1,
            center2,
            r_a,
            r_b,
            va,
            vb,
            mat_b,
            eigval,
            eigvec,
            mat_r,
            rmsd,
        })
    }

    /// Returns a reference to the first atom group.
    pub fn atomgroup1(&self) -> &AtomGroup {
        &self.atomgroup1
    }

    /// Returns a reference to the second atom group.
    pub fn atomgroup2(&self) -> &AtomGroup {
        &self.atomgroup2
    }

    /// Number of matched common points between the two atom groups.
    pub fn num_of_positions(&self) -> usize {
        self.positions1.len()
    }

    /// First atom group's matched positions.
    pub fn positions1(&self) -> &[Position] {
        &self.positions1
    }

    /// Second atom group's matched positions.
    pub fn positions2(&self) -> &[Position] {
        &self.positions2
    }

    /// Centroid of positions in the first atom group.
    pub fn center1(&self) -> Position {
        self.center1
    }

    /// Centroid of positions in the second atom group.
    pub fn center2(&self) -> Position {
        self.center2
    }

    /// Centroid-shifted positions of the first atom group (r_A).
    pub fn r_a(&self) -> &[Position] {
        &self.r_a
    }

    /// Centroid-shifted positions of the second atom group (r_B).
    pub fn r_b(&self) -> &[Position] {
        &self.r_b
    }

    /// Sum vector va = r_B + r_A.
    pub fn va(&self) -> &[Position] {
        &self.va
    }

    /// Difference vector vb = r_B - r_A.
    pub fn vb(&self) -> &[Position] {
        &self.vb
    }

    /// The 4x4 symmetric matrix B.
    pub fn mat_b(&self) -> &SymmetricMatrix {
        &self.mat_b
    }

    /// Eigenvalues of matrix B (in ascending order).
    pub fn eigval(&self) -> &Vector {
        &self.eigval
    }

    /// Eigenvectors of matrix B (row i corresponds to eigenvalue i).
    pub fn eigvec(&self) -> &Matrix {
        &self.eigvec
    }

    /// 3x3 optimal rotation matrix R.
    pub fn mat_r(&self) -> &Matrix {
        &self.mat_r
    }

    /// Alias for `mat_r()` matching Python's `rotation_mat` property.
    pub fn rotation_mat(&self) -> &Matrix {
        &self.mat_r
    }

    /// Root-mean-square deviation (RMSD).
    pub fn rmsd(&self) -> f64 {
        self.rmsd
    }

    /// Compute and return the RMSD (matching Python's `calc()` method).
    pub fn calc(&self) -> f64 {
        self.rmsd
    }

    /// Superimposes the given `atomgroup` by shifting by `-center1`, rotating by `mat_r`,
    /// and shifting by `+center2`.
    pub fn superimpose(&self, atomgroup: &AtomGroup) -> Result<AtomGroup> {
        let mut answer = atomgroup.clone();
        answer.shift_by(-self.center1);
        answer.rotate(&self.mat_r)?;
        answer.shift_by(self.center2);
        Ok(answer)
    }

    /// Recursively matches common atom positions between two atom groups.
    fn match_positions(
        atomgroup1: &AtomGroup,
        atomgroup2: &AtomGroup,
    ) -> (Vec<Position>, Vec<Position>) {
        let mut positions1 = Vec::new();
        let mut positions2 = Vec::new();

        for (key, ag1) in atomgroup1.groups() {
            if let Some(ag2) = atomgroup2.get_group(key) {
                let (p1, p2) = Self::match_positions(ag1, ag2);
                positions1.extend(p1);
                positions2.extend(p2);
            }
        }

        for (key, atom1) in atomgroup1.atoms() {
            if let Some(atom2) = atomgroup2.get_atom(key) {
                positions1.push(atom1.xyz);
                positions2.push(atom2.xyz);
            }
        }

        (positions1, positions2)
    }

    /// Computes the centroid of positions.
    fn calc_center(positions: &[Position]) -> Position {
        if positions.is_empty() {
            return Position::default();
        }
        let mut sum = Position::default();
        for p in positions {
            sum += *p;
        }
        sum / (positions.len() as f64)
    }

    /// Shifts positions by subtracting `center`.
    fn shift_positions(positions: &[Position], center: Position) -> Vec<Position> {
        positions.iter().map(|p| *p - center).collect()
    }

    /// Makes va = r_B + r_A.
    fn make_va(r_a: &[Position], r_b: &[Position]) -> Vec<Position> {
        r_a.iter().zip(r_b.iter()).map(|(a, b)| *b + *a).collect()
    }

    /// Makes vb = r_B - r_A.
    fn make_vb(r_a: &[Position], r_b: &[Position]) -> Vec<Position> {
        r_a.iter().zip(r_b.iter()).map(|(a, b)| *b - *a).collect()
    }

    /// Constructs the 4x4 symmetric matrix B according to Kearsley (1989).
    fn make_b(va: &[Position], vb: &[Position]) -> SymmetricMatrix {
        let n = va.len();
        let mut b = SymmetricMatrix::new(4);
        for i in 0..n {
            let a = va[i];
            let b_pos = vb[i];
            let ax = a.x;
            let ay = a.y;
            let az = a.z;
            let bx = b_pos.x;
            let by = b_pos.y;
            let bz = b_pos.z;

            b.add(0, 0, bx * bx + by * by + bz * bz);
            b.add(0, 1, az * by - ay * bz);
            b.add(0, 2, -az * bx + ax * bz);
            b.add(0, 3, ay * bx - ax * by);
            b.add(1, 1, bx * bx + ay * ay + az * az);
            b.add(1, 2, bx * by - ax * ay);
            b.add(1, 3, bx * bz - ax * az);
            b.add(2, 2, ax * ax + by * by + az * az);
            b.add(2, 3, by * bz - ay * az);
            b.add(3, 3, ax * ax + ay * ay + bz * bz);
        }

        let scale = 1.0 / ((n * n) as f64);
        for r in 0..4 {
            for c in 0..=r {
                if let Ok(v) = b.get(r, c) {
                    b.set(r, c, v * scale);
                }
            }
        }
        b
    }

    /// Constructs 3x3 rotation matrix R from quaternion vector q.
    fn make_r(q: &Vector) -> Result<Matrix> {
        if q.len() != 4 {
            return Err(BridgeError::value_error(
                "q",
                "Quaternion vector must have length 4",
            ));
        }
        let q0 = q.get(0)?;
        let q1 = q.get(1)?;
        let q2 = q.get(2)?;
        let q3 = q.get(3)?;

        let mut r = Matrix::new(3, 3);
        r.set(0, 0, 2.0 * q0 * q0 + 2.0 * q1 * q1 - 1.0);
        r.set(0, 1, 2.0 * q1 * q2 - 2.0 * q0 * q3);
        r.set(0, 2, 2.0 * q1 * q3 + 2.0 * q0 * q2);

        r.set(1, 0, 2.0 * q1 * q2 + 2.0 * q0 * q3);
        r.set(1, 1, 2.0 * q0 * q0 + 2.0 * q2 * q2 - 1.0);
        r.set(1, 2, 2.0 * q2 * q3 - 2.0 * q0 * q1);

        r.set(2, 0, 2.0 * q1 * q3 - 2.0 * q0 * q2);
        r.set(2, 1, 2.0 * q2 * q3 + 2.0 * q0 * q1);
        r.set(2, 2, 2.0 * q0 * q0 + 2.0 * q3 * q3 - 1.0);

        Ok(r)
    }

    /// Calculates RMSD between r_A and r_B given rotation matrix R.
    fn calc_rmsd(r_a: &[Position], r_b: &[Position], mat_r: &Matrix) -> Result<f64> {
        let n = r_a.len();
        let mut msd = 0.0;
        for i in 0..n {
            let mut rotated_a = r_a[i];
            rotated_a.rotate(mat_r)?;
            let diff = r_b[i] - rotated_a;
            msd += diff.square_distance_from(&Position::default());
        }
        msd /= n as f64;
        Ok(msd.sqrt())
    }
}
