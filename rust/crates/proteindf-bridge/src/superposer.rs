// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom_group::AtomGroup;
use crate::error::{BridgeError, Result};
use crate::matrix::Matrix;
use crate::position::Position;
use crate::vector::Vector;

/// Structure for superposing two AtomGroups using the Kabsch algorithm.
///
/// Ports `proteindf_bridge.superposer.Superposer`.
#[derive(Debug, Clone)]
pub struct Superposer {
    atom_group1: AtomGroup,
    atom_group2: AtomGroup,
    num_of_positions: usize,
    positions1: Vec<Position>,
    positions2: Vec<Position>,
    center1: Position,
    center2: Position,
    shift_positions1: Vec<Position>,
    shift_positions2: Vec<Position>,
    rotation_mat: Matrix,
    update_positions1: Vec<Position>,
    update_positions2: Vec<Position>,
    rmsd: f64,
}

impl Superposer {
    /// Creates a new `Superposer` from two `AtomGroup` instances and computes RMSD and rotation matrix.
    pub fn new(atom_group1: &AtomGroup, atom_group2: &AtomGroup) -> Result<Self> {
        let (positions1, positions2) = Self::match_positions(atom_group1, atom_group2);
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

        let shift_positions1 = Self::shift_positions(&positions1, center1);
        let shift_positions2 = Self::shift_positions(&positions2, center2);

        let rotation_mat =
            Self::calc_rotation_matrix(num_of_positions, &shift_positions1, &shift_positions2)?;

        // Rotate update_positions1 by rotation_mat
        let mut update_positions1 = shift_positions1.clone();
        for pos in &mut update_positions1 {
            pos.rotate(&rotation_mat)?;
        }

        let update_positions2 = shift_positions2.clone();

        let rmsd = Self::calc_rmsd(&update_positions1, &update_positions2);

        Ok(Self {
            atom_group1: atom_group1.clone(),
            atom_group2: atom_group2.clone(),
            num_of_positions,
            positions1,
            positions2,
            center1,
            center2,
            shift_positions1,
            shift_positions2,
            rotation_mat,
            update_positions1,
            update_positions2,
            rmsd,
        })
    }

    /// Returns a reference to the first atom group.
    pub fn atom_group1(&self) -> &AtomGroup {
        &self.atom_group1
    }

    /// Returns a reference to the second atom group.
    pub fn atom_group2(&self) -> &AtomGroup {
        &self.atom_group2
    }

    /// Number of matched common points between the two atom groups.
    pub fn num_of_positions(&self) -> usize {
        self.num_of_positions
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

    /// Centroid-shifted positions of the first atom group.
    pub fn shift_positions1(&self) -> &[Position] {
        &self.shift_positions1
    }

    /// Centroid-shifted positions of the second atom group.
    pub fn shift_positions2(&self) -> &[Position] {
        &self.shift_positions2
    }

    /// 3x3 rotation matrix calculated to minimize RMSD.
    pub fn rotation_mat(&self) -> &Matrix {
        &self.rotation_mat
    }

    /// Rotated and shifted positions of the first atom group.
    pub fn update_positions1(&self) -> &[Position] {
        &self.update_positions1
    }

    /// Shifted positions of the second atom group (reference).
    pub fn update_positions2(&self) -> &[Position] {
        &self.update_positions2
    }

    /// Root-mean-square deviation (RMSD) between the two groups.
    pub fn rmsd(&self) -> f64 {
        self.rmsd
    }

    /// Superimposes the given `atomgroup` by shifting by `-center1`, rotating by `rotation_mat`,
    /// and shifting by `+center2`.
    pub fn superimpose(&self, atomgroup: &AtomGroup) -> Result<AtomGroup> {
        let mut answer = atomgroup.clone();
        answer.shift_by(-self.center1);
        answer.rotate(&self.rotation_mat)?;
        answer.shift_by(self.center2);
        Ok(answer)
    }

    /// Recursively matches common atom positions between two atom groups.
    fn match_positions(
        atom_group1: &AtomGroup,
        atom_group2: &AtomGroup,
    ) -> (Vec<Position>, Vec<Position>) {
        let mut positions1 = Vec::new();
        let mut positions2 = Vec::new();

        for (key, ag1) in atom_group1.groups() {
            if let Some(ag2) = atom_group2.get_group(key) {
                let (p1, p2) = Self::match_positions(ag1, ag2);
                positions1.extend(p1);
                positions2.extend(p2);
            }
        }

        for (key, atom1) in atom_group1.atoms() {
            if let Some(atom2) = atom_group2.get_atom(key) {
                positions1.push(atom1.xyz);
                positions2.push(atom2.xyz);
            }
        }

        (positions1, positions2)
    }

    /// Computes the geometric center (centroid) of positions.
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

    /// Shifts all positions by subtracting `center`.
    fn shift_positions(positions: &[Position], center: Position) -> Vec<Position> {
        positions.iter().map(|p| *p - center).collect()
    }

    /// Computes the optimal 3x3 rotation matrix using the Kabsch algorithm.
    fn calc_rotation_matrix(
        num_of_points: usize,
        positions1: &[Position],
        positions2: &[Position],
    ) -> Result<Matrix> {
        let mut r = Matrix::new(3, 3);

        // r_ij = Sum_over_k { p2(k, i) * p1(k, j) }
        for k in 0..num_of_points {
            let x1 = positions1[k].x;
            let y1 = positions1[k].y;
            let z1 = positions1[k].z;
            let x2 = positions2[k].x;
            let y2 = positions2[k].y;
            let z2 = positions2[k].z;

            r.add(0, 0, x2 * x1);
            r.add(0, 1, x2 * y1);
            r.add(0, 2, x2 * z1);
            r.add(1, 0, y2 * x1);
            r.add(1, 1, y2 * y1);
            r.add(1, 2, y2 * z1);
            r.add(2, 0, z2 * x1);
            r.add(2, 1, z2 * y1);
            r.add(2, 2, z2 * z1);
        }

        let tr = r.transpose();

        let trr = &tr * &r;
        let trr_sym = trr.get_symmetric_matrix()?;

        let (_eigval, eigvec) = trr_sym.eig()?;

        // NOTE: In the original Python implementation:
        // `a = self._make_right_handed(eigvec)` is called with `eigvec`.
        let a = Self::make_right_handed(&eigvec)?;

        let mut b = Matrix::new(3, 3);
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    let v = r.get(j, k)? * a.get(i, k)?;
                    b.add(i, j, v);
                }
            }

            // Normalize row b[i]
            let mut w = 0.0;
            for j in 0..3 {
                let v = b.get(i, j)?;
                w += v * v;
            }

            if w > 0.0 {
                let t = (1.0 / w).sqrt();
                for j in 0..3 {
                    let v = b.get(i, j)?;
                    b.set(i, j, v * t);
                }
            }
        }

        // rotation matrix r_ij = b_ki * a_kj
        Self::set_rotation(&a, &b)
    }

    /// Makes the coordinate system right-handed.
    fn make_right_handed(mat: &Matrix) -> Result<Matrix> {
        let mut v1 = Vector::new(3);
        let mut v2 = Vector::new(3);
        for i in 0..3 {
            v1.set(i, mat.get(0, i)?);
            v2.set(i, mat.get(1, i)?);
        }

        let v3 = Self::calc_vector_product(&v1, &v2)?;

        let mut answer = Matrix::new(3, 3);
        for i in 0..3 {
            answer.set(0, i, v1.get(i)?);
            answer.set(1, i, v2.get(i)?);
            answer.set(2, i, v3.get(i)?);
        }

        Ok(answer)
    }

    /// Computes the cross product of two 3D vectors.
    fn calc_vector_product(v1: &Vector, v2: &Vector) -> Result<Vector> {
        let mut v3 = Vector::new(3);
        v3.set(0, v1.get(1)? * v2.get(2)? - v1.get(2)? * v2.get(1)?);
        v3.set(1, v1.get(2)? * v2.get(0)? - v1.get(0)? * v2.get(2)?);
        v3.set(2, v1.get(0)? * v2.get(1)? - v1.get(1)? * v2.get(0)?);
        Ok(v3)
    }

    /// Computes rotation matrix R_ij = sum_k { b_ki * a_kj }.
    fn set_rotation(a: &Matrix, b: &Matrix) -> Result<Matrix> {
        let mut r = Matrix::new(3, 3);
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    let v = b.get(k, i)? * a.get(k, j)?;
                    r.add(i, j, v);
                }
            }
        }
        Ok(r)
    }

    /// Computes RMSD between two sets of positions.
    fn calc_rmsd(positions1: &[Position], positions2: &[Position]) -> f64 {
        let n = positions1.len().min(positions2.len());
        if n == 0 {
            return 0.0;
        }
        let mut msd = 0.0;
        for i in 0..n {
            msd += positions1[i].square_distance_from(&positions2[i]);
        }
        (msd / (n as f64)).sqrt()
    }
}
