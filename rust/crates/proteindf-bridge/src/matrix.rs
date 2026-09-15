// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::fmt;
use std::ops::{AddAssign, Mul, MulAssign, Sub, SubAssign};

use crate::error::{BridgeError, Result};
use crate::vector::Vector;

/// General 2D numerical matrix corresponding to `proteindf_bridge.matrix.Matrix`.
#[derive(Debug, Clone)]
pub struct Matrix {
    rows: usize,
    cols: usize,
    data: Vec<f64>,
}

impl Matrix {
    /// Creates a new matrix of size `rows` x `cols`, initialized to 0.0.
    pub fn new(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            data: vec![0.0; rows * cols],
        }
    }

    /// Creates a 1x1 matrix initialized to 0.0 (default constructor in Python).
    pub fn default_1x1() -> Self {
        Self::new(1, 1)
    }

    /// Creates a matrix from a nested slice `&[&[f64]]`.
    pub fn from_nested_slice(slice: &[&[f64]]) -> Self {
        let rows = slice.len();
        assert!(rows > 0, "Matrix must have at least one row");
        let cols = slice[0].len();
        assert!(cols > 0, "Matrix must have at least one column");

        let mut data = Vec::with_capacity(rows * cols);
        for row in slice {
            assert_eq!(row.len(), cols, "All rows must have the same length");
            data.extend_from_slice(row);
        }

        Self { rows, cols, data }
    }

    /// Creates a matrix from flat data with dimensions.
    pub fn from_vec(rows: usize, cols: usize, data: Vec<f64>) -> Self {
        assert_eq!(
            rows * cols,
            data.len(),
            "Data length does not match matrix dimensions"
        );
        Self { rows, cols, data }
    }

    /// Returns the number of rows.
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Returns the number of columns.
    pub fn cols(&self) -> usize {
        self.cols
    }

    #[inline]
    fn index_of(&self, row: usize, col: usize) -> usize {
        row * self.cols + col
    }

    /// Gets the element at (row, col).
    pub fn get(&self, row: usize, col: usize) -> Result<f64> {
        if row >= self.rows {
            return Err(BridgeError::value_error(
                format!("row: {}", row),
                format!("out of range in row: 0 <= {} < {}", row, self.rows),
            ));
        }
        if col >= self.cols {
            return Err(BridgeError::value_error(
                format!("col: {}", col),
                format!("out of range in col: 0 <= {} < {}", col, self.cols),
            ));
        }
        Ok(self.data[self.index_of(row, col)])
    }

    /// Sets the element at (row, col).
    pub fn set(&mut self, row: usize, col: usize, value: f64) {
        assert!(
            row < self.rows && col < self.cols,
            "Index out of bounds: ({}, {}) for shape ({}, {})",
            row,
            col,
            self.rows,
            self.cols
        );
        let idx = self.index_of(row, col);
        self.data[idx] = value;
    }

    /// Adds `value` to the element at (row, col).
    pub fn add(&mut self, row: usize, col: usize, value: f64) {
        assert!(
            row < self.rows && col < self.cols,
            "Index out of bounds: ({}, {}) for shape ({}, {})",
            row,
            col,
            self.rows,
            self.cols
        );
        let idx = self.index_of(row, col);
        self.data[idx] += value;
    }

    /// Returns the transpose of this matrix.
    pub fn transpose(&self) -> Self {
        let mut result = Matrix::new(self.cols, self.rows);
        for r in 0..self.rows {
            for c in 0..self.cols {
                result.set(c, r, self.data[self.index_of(r, c)]);
            }
        }
        result
    }

    /// Selects a sub-block of this matrix.
    pub fn select(
        &self,
        start_row: usize,
        start_col: usize,
        end_row: usize,
        end_col: usize,
    ) -> Result<Self> {
        if start_row >= end_row
            || end_row > self.rows
            || start_col >= end_col
            || end_col > self.cols
        {
            return Err(BridgeError::value_error(
                format!(
                    "({}, {}) to ({}, {})",
                    start_row, start_col, end_row, end_col
                ),
                "invalid sub-matrix bounds",
            ));
        }
        let sub_rows = end_row - start_row;
        let sub_cols = end_col - start_col;
        let mut sub = Matrix::new(sub_rows, sub_cols);
        for r in 0..sub_rows {
            for c in 0..sub_cols {
                sub.set(r, c, self.get(start_row + r, start_col + c)?);
            }
        }
        Ok(sub)
    }

    /// Extracts row `row` as a `Vector`.
    pub fn get_row_vector(&self, row: usize) -> Result<Vector> {
        if row >= self.rows {
            return Err(BridgeError::value_error(
                row.to_string(),
                "row index out of range",
            ));
        }
        let mut v = Vector::new(self.cols);
        for c in 0..self.cols {
            v.set(c, self.data[self.index_of(row, c)]);
        }
        Ok(v)
    }

    /// Extracts column `col` as a `Vector`.
    pub fn get_col_vector(&self, col: usize) -> Result<Vector> {
        if col >= self.cols {
            return Err(BridgeError::value_error(
                col.to_string(),
                "col index out of range",
            ));
        }
        let mut v = Vector::new(self.rows);
        for r in 0..self.rows {
            v.set(r, self.data[self.index_of(r, col)]);
        }
        Ok(v)
    }

    /// Returns the maximum value in the matrix.
    pub fn max(&self) -> f64 {
        self.data.iter().copied().fold(f64::NEG_INFINITY, f64::max)
    }

    /// Returns the minimum value in the matrix.
    pub fn min(&self) -> f64 {
        self.data.iter().copied().fold(f64::INFINITY, f64::min)
    }

    /// Resizes the matrix in-place. Existing values are preserved where indices overlap.
    pub fn resize(&mut self, new_rows: usize, new_cols: usize) {
        let mut new_data = vec![0.0; new_rows * new_cols];
        let min_r = std::cmp::min(self.rows, new_rows);
        let min_c = std::cmp::min(self.cols, new_cols);
        for r in 0..min_r {
            for c in 0..min_c {
                new_data[r * new_cols + c] = self.data[self.index_of(r, c)];
            }
        }
        self.rows = new_rows;
        self.cols = new_cols;
        self.data = new_data;
    }

    /// Computes the inverse of a square matrix using Gauss-Jordan elimination with partial pivoting.
    pub fn inverse(&self) -> Result<Matrix> {
        if self.rows != self.cols {
            return Err(BridgeError::value_error(
                format!("({}, {})", self.rows, self.cols),
                "Cannot invert non-square matrix",
            ));
        }
        let n = self.rows;
        let mut aug = vec![0.0; n * 2 * n];
        let aug_cols = 2 * n;

        for r in 0..n {
            for c in 0..n {
                aug[r * aug_cols + c] = self.data[self.index_of(r, c)];
            }
            aug[r * aug_cols + (n + r)] = 1.0;
        }

        for i in 0..n {
            // Find pivot
            let mut max_row = i;
            let mut max_val = aug[i * aug_cols + i].abs();
            for r in (i + 1)..n {
                let val = aug[r * aug_cols + i].abs();
                if val > max_val {
                    max_val = val;
                    max_row = r;
                }
            }

            if max_val < 1e-15 {
                return Err(BridgeError::value_error(
                    "Matrix",
                    "Matrix is singular and cannot be inverted",
                ));
            }

            // Swap rows
            if max_row != i {
                for c in 0..aug_cols {
                    aug.swap(i * aug_cols + c, max_row * aug_cols + c);
                }
            }

            // Scale pivot row
            let pivot = aug[i * aug_cols + i];
            for c in 0..aug_cols {
                aug[i * aug_cols + c] /= pivot;
            }

            // Eliminate other rows
            for r in 0..n {
                if r != i {
                    let factor = aug[r * aug_cols + i];
                    for c in 0..aug_cols {
                        aug[r * aug_cols + c] -= factor * aug[i * aug_cols + c];
                    }
                }
            }
        }

        let mut inv = Matrix::new(n, n);
        for r in 0..n {
            for c in 0..n {
                inv.set(r, c, aug[r * aug_cols + (n + c)]);
            }
        }
        Ok(inv)
    }

    /// Converts this square matrix to a `SymmetricMatrix`.
    pub fn get_symmetric_matrix(&self) -> Result<SymmetricMatrix> {
        if self.rows != self.cols {
            return Err(BridgeError::value_error(
                format!("({}, {})", self.rows, self.cols),
                "Cannot convert non-square Matrix to SymmetricMatrix",
            ));
        }
        let mut sm = SymmetricMatrix::new(self.rows);
        for r in 0..self.rows {
            for c in 0..=r {
                sm.set(r, c, self.get(r, c)?);
            }
        }
        Ok(sm)
    }

    /// Returns a flat vector of all elements (row-major).
    pub fn to_vec(&self) -> Vec<f64> {
        self.data.clone()
    }
}

// Equality with tolerance matching Python's `math.fabs(a - b) < 1.0e-5`.
impl PartialEq for Matrix {
    fn eq(&self, other: &Self) -> bool {
        if self.rows != other.rows || self.cols != other.cols {
            return false;
        }
        self.data
            .iter()
            .zip(other.data.iter())
            .all(|(&a, &b)| (a - b).abs() <= 1.0e-5)
    }
}

impl std::ops::Add for &Matrix {
    type Output = Matrix;
    fn add(self, rhs: Self) -> Self::Output {
        assert_eq!(
            (self.rows, self.cols),
            (rhs.rows, rhs.cols),
            "Matrix dimensions must match for addition"
        );
        Matrix {
            rows: self.rows,
            cols: self.cols,
            data: self
                .data
                .iter()
                .zip(rhs.data.iter())
                .map(|(&a, &b)| a + b)
                .collect(),
        }
    }
}

impl std::ops::Add for Matrix {
    type Output = Matrix;
    fn add(self, rhs: Self) -> Self::Output {
        &self + &rhs
    }
}

impl AddAssign<&Matrix> for Matrix {
    fn add_assign(&mut self, rhs: &Matrix) {
        assert_eq!(
            (self.rows, self.cols),
            (rhs.rows, rhs.cols),
            "Matrix dimensions must match for add_assign"
        );
        for (a, &b) in self.data.iter_mut().zip(rhs.data.iter()) {
            *a += b;
        }
    }
}

impl Sub for &Matrix {
    type Output = Matrix;
    fn sub(self, rhs: Self) -> Self::Output {
        assert_eq!(
            (self.rows, self.cols),
            (rhs.rows, rhs.cols),
            "Matrix dimensions must match for subtraction"
        );
        Matrix {
            rows: self.rows,
            cols: self.cols,
            data: self
                .data
                .iter()
                .zip(rhs.data.iter())
                .map(|(&a, &b)| a - b)
                .collect(),
        }
    }
}

impl Sub for Matrix {
    type Output = Matrix;
    fn sub(self, rhs: Self) -> Self::Output {
        &self - &rhs
    }
}

impl SubAssign<&Matrix> for Matrix {
    fn sub_assign(&mut self, rhs: &Matrix) {
        assert_eq!(
            (self.rows, self.cols),
            (rhs.rows, rhs.cols),
            "Matrix dimensions must match for sub_assign"
        );
        for (a, &b) in self.data.iter_mut().zip(rhs.data.iter()) {
            *a -= b;
        }
    }
}

// Matrix * scalar
impl Mul<f64> for &Matrix {
    type Output = Matrix;
    fn mul(self, rhs: f64) -> Self::Output {
        Matrix {
            rows: self.rows,
            cols: self.cols,
            data: self.data.iter().map(|&x| x * rhs).collect(),
        }
    }
}

impl Mul<f64> for Matrix {
    type Output = Matrix;
    fn mul(self, rhs: f64) -> Self::Output {
        &self * rhs
    }
}

impl MulAssign<f64> for Matrix {
    fn mul_assign(&mut self, rhs: f64) {
        for x in self.data.iter_mut() {
            *x *= rhs;
        }
    }
}

// Matrix * Matrix
impl Mul<&Matrix> for &Matrix {
    type Output = Matrix;
    fn mul(self, rhs: &Matrix) -> Self::Output {
        assert_eq!(
            self.cols, rhs.rows,
            "Matrix multiplication inner dimensions must match: {} != {}",
            self.cols, rhs.rows
        );
        let mut result = Matrix::new(self.rows, rhs.cols);
        for i in 0..self.rows {
            for k in 0..self.cols {
                let a = self.data[self.index_of(i, k)];
                for j in 0..rhs.cols {
                    let out_idx = result.index_of(i, j);
                    result.data[out_idx] += a * rhs.data[rhs.index_of(k, j)];
                }
            }
        }
        result
    }
}

impl Mul<Matrix> for Matrix {
    type Output = Matrix;
    fn mul(self, rhs: Matrix) -> Self::Output {
        &self * &rhs
    }
}

// Matrix * Vector -> Vector
impl Mul<&Vector> for &Matrix {
    type Output = Vector;
    fn mul(self, rhs: &Vector) -> Self::Output {
        assert_eq!(
            self.cols,
            rhs.len(),
            "Matrix cols must match Vector length for multiplication: {} != {}",
            self.cols,
            rhs.len()
        );
        let mut result = Vector::new(self.rows);
        for i in 0..self.rows {
            let mut sum = 0.0;
            for j in 0..self.cols {
                sum += self.data[self.index_of(i, j)] * rhs[j];
            }
            result.set(i, sum);
        }
        result
    }
}

impl Mul<Vector> for Matrix {
    type Output = Vector;
    fn mul(self, rhs: Vector) -> Self::Output {
        &self * &rhs
    }
}

impl fmt::Display for Matrix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for order in (0..self.cols).step_by(10) {
            write!(f, "       ")?;
            for j in order..std::cmp::min(order + 10, self.cols) {
                write!(f, "   {:5} th", j + 1)?;
            }
            writeln!(f, "\n  ---- -----------")?;
            for i in 0..self.rows {
                write!(f, " {:5} ", i + 1)?;
                for j in order..std::cmp::min(order + 10, self.cols) {
                    write!(f, " {:10.6}", self.get(i, j).unwrap_or(0.0))?;
                }
                writeln!(f)?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

/// Symmetric matrix corresponding to `proteindf_bridge.matrix.SymmetricMatrix`.
#[derive(Debug, Clone)]
pub struct SymmetricMatrix {
    dim: usize,
    matrix: Matrix,
}

impl SymmetricMatrix {
    /// Creates a new symmetric matrix of dimension `dim` x `dim`, initialized to 0.0.
    pub fn new(dim: usize) -> Self {
        Self {
            dim,
            matrix: Matrix::new(dim, dim),
        }
    }

    /// Returns the dimension of the matrix.
    pub fn dim(&self) -> usize {
        self.dim
    }

    /// Gets the element at (row, col), automatically handling symmetry (r < c swapped to c, r).
    pub fn get(&self, mut row: usize, mut col: usize) -> Result<f64> {
        if row < col {
            std::mem::swap(&mut row, &mut col);
        }
        self.matrix.get(row, col)
    }

    /// Sets the element at (row, col) and its symmetric counterpart.
    pub fn set(&mut self, mut row: usize, mut col: usize, value: f64) {
        if row < col {
            std::mem::swap(&mut row, &mut col);
        }
        self.matrix.set(row, col, value);
        self.matrix.set(col, row, value);
    }

    /// Adds `value` to the element at (row, col) and its symmetric counterpart.
    pub fn add(&mut self, mut row: usize, mut col: usize, value: f64) {
        if row < col {
            std::mem::swap(&mut row, &mut col);
        }
        self.matrix.add(row, col, value);
        if row != col {
            self.matrix.add(col, row, value);
        }
    }

    /// Converts this to a general `Matrix`.
    pub fn get_general_matrix(&self) -> Matrix {
        self.matrix.clone()
    }

    /// Resizes the symmetric matrix in-place.
    pub fn resize(&mut self, new_dim: usize) {
        let mut new_sm = SymmetricMatrix::new(new_dim);
        let min_d = std::cmp::min(self.dim, new_dim);
        for r in 0..min_d {
            for c in 0..=r {
                if let Ok(v) = self.get(r, c) {
                    new_sm.set(r, c, v);
                }
            }
        }
        *self = new_sm;
    }

    /// Returns the raw packed lower triangular data ('SP' format) of length `dim * (dim + 1) / 2`.
    pub fn get_raw_data(&self) -> Vec<f64> {
        let mut data = Vec::with_capacity(self.dim * (self.dim + 1) / 2);
        for r in 0..self.dim {
            for c in 0..=r {
                data.push(self.get(r, c).unwrap_or(0.0));
            }
        }
        data
    }

    /// Computes eigenvalues and eigenvectors using the Jacobi eigenvalue algorithm.
    ///
    /// Returns `(eigenvalues, eigenvectors)` where eigenvalues are in ascending order,
    /// and row `i` of the returned `Matrix` is the eigenvector corresponding to `eigenvalues[i]`
    /// (matching Python's `numpy.linalg.eigh(A)` followed by transpose).
    pub fn eig(&self) -> Result<(Vector, Matrix)> {
        let n = self.dim;
        if n == 0 {
            return Ok((Vector::new(0), Matrix::new(0, 0)));
        }
        if n == 1 {
            let val = self.get(0, 0)?;
            let mut vec_mat = Matrix::new(1, 1);
            vec_mat.set(0, 0, 1.0);
            return Ok((Vector::from_slice(&[val]), vec_mat));
        }

        // Work on a copy of the matrix
        let mut a = vec![0.0; n * n];
        for r in 0..n {
            for c in 0..n {
                a[r * n + c] = self.get(r, c)?;
            }
        }

        // V initialized to identity matrix
        let mut v = vec![0.0; n * n];
        for i in 0..n {
            v[i * n + i] = 1.0;
        }

        let max_iterations = 100 * n * n;
        for _ in 0..max_iterations {
            // Find largest off-diagonal element
            let mut max_val = 0.0;
            let mut p = 0;
            let mut q = 1;
            for r in 0..n {
                for c in (r + 1)..n {
                    let val = a[r * n + c].abs();
                    if val > max_val {
                        max_val = val;
                        p = r;
                        q = c;
                    }
                }
            }

            if max_val < 1e-15 {
                break;
            }

            // Compute Jacobi rotation angle
            let app = a[p * n + p];
            let aqq = a[q * n + q];
            let apq = a[p * n + q];

            let theta = 0.5 * (aqq - app) / apq;
            let t = if theta >= 0.0 {
                1.0 / (theta + (1.0 + theta * theta).sqrt())
            } else {
                -1.0 / (-theta + (1.0 + theta * theta).sqrt())
            };
            let c = 1.0 / (1.0 + t * t).sqrt();
            let s = t * c;

            // Apply rotation to A
            let tau = s / (1.0 + c);
            a[p * n + p] -= t * apq;
            a[q * n + q] += t * apq;
            a[p * n + q] = 0.0;
            a[q * n + p] = 0.0;

            for r in 0..n {
                if r != p && r != q {
                    let arp = a[r * n + p];
                    let arq = a[r * n + q];
                    a[r * n + p] = arp - s * (arq + tau * arp);
                    a[p * n + r] = a[r * n + p];
                    a[r * n + q] = arq + s * (arp - tau * arq);
                    a[q * n + r] = a[r * n + q];
                }
            }

            // Accumulate transformations into V
            for r in 0..n {
                let vrp = v[r * n + p];
                let vrq = v[r * n + q];
                v[r * n + p] = vrp - s * (vrq + tau * vrp);
                v[r * n + q] = vrq + s * (vrp - tau * vrq);
            }
        }

        // Extract eigenvalues
        let mut eigen_pairs: Vec<(f64, Vec<f64>)> = (0..n)
            .map(|col| {
                let val = a[col * n + col];
                let vec: Vec<f64> = (0..n).map(|row| v[row * n + col]).collect();
                (val, vec)
            })
            .collect();

        // Sort eigenvalues in ascending order
        eigen_pairs.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));

        let mut eigenvalues = Vector::new(n);
        let mut eigenvectors = Matrix::new(n, n);
        for (i, (val, vec)) in eigen_pairs.into_iter().enumerate() {
            eigenvalues.set(i, val);
            // Python transpose puts eigenvector i as row i
            for (j, &coord) in vec.iter().enumerate() {
                eigenvectors.set(i, j, coord);
            }
        }

        Ok((eigenvalues, eigenvectors))
    }
}

impl PartialEq for SymmetricMatrix {
    fn eq(&self, other: &Self) -> bool {
        if self.dim != other.dim {
            return false;
        }
        for r in 0..self.dim {
            for c in 0..=r {
                let a = self.get(r, c).unwrap_or(0.0);
                let b = other.get(r, c).unwrap_or(0.0);
                if (a - b).abs() > 1.0e-5 {
                    return false;
                }
            }
        }
        true
    }
}

/// Creates an identity matrix of size `dim` x `dim`.
pub fn identity_matrix(dim: usize) -> SymmetricMatrix {
    let mut mat = SymmetricMatrix::new(dim);
    for i in 0..dim {
        mat.set(i, i, 1.0);
    }
    mat
}

#[cfg(test)]
mod tests {
    use super::*;

    // Ported from tests/test_matrix.py
    #[test]
    fn test_init_and_shape() {
        let m = Matrix::new(2, 3);
        assert_eq!(m.rows(), 2);
        assert_eq!(m.cols(), 3);
    }

    #[test]
    fn test_init_with_nested_list() {
        let data: &[&[f64]] = &[&[1.0, 2.0], &[3.0, 4.0]];
        let m = Matrix::from_nested_slice(data);
        assert_eq!(m.rows(), 2);
        assert_eq!(m.cols(), 2);
        assert!((m.get(0, 1).unwrap() - 2.0).abs() < 1e-10);
        assert!((m.get(1, 0).unwrap() - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_set_and_get() {
        let mut m = Matrix::new(2, 2);
        m.set(0, 1, 5.5);
        assert!((m.get(0, 1).unwrap() - 5.5).abs() < 1e-10);
        m.add(0, 1, 2.0);
        assert!((m.get(0, 1).unwrap() - 7.5).abs() < 1e-10);
    }

    #[test]
    fn test_multiplication() {
        let c = Matrix::from_nested_slice(&[&[7.0, 4.0, -1.0], &[3.0, 0.0, 5.0]]);
        let d =
            Matrix::from_nested_slice(&[&[8.0, 4.0, 2.0], &[1.0, 3.0, -6.0], &[-7.0, 0.0, 5.0]]);
        let cd = &c * &d;
        let expected = Matrix::from_nested_slice(&[&[67.0, 40.0, -15.0], &[-11.0, 12.0, 31.0]]);
        assert_eq!(cd, expected);
    }

    #[test]
    fn test_vector_multiplication() {
        let m = Matrix::from_nested_slice(&[&[1.0, 2.0], &[3.0, 4.0]]);
        let v = Vector::from_slice(&[1.0, 1.0]);
        let res = &m * &v;
        assert!((res[0] - 3.0).abs() < 1e-10);
        assert!((res[1] - 7.0).abs() < 1e-10);
    }

    #[test]
    fn test_inverse() {
        let m = Matrix::from_nested_slice(&[&[4.0, 7.0], &[2.0, 6.0]]);
        let inv = m.inverse().unwrap();
        let eye = &m * &inv;
        let expected_eye = Matrix::from_nested_slice(&[&[1.0, 0.0], &[0.0, 1.0]]);
        assert_eq!(eye, expected_eye);
    }

    // SymmetricMatrix tests
    #[test]
    fn test_symmetric_indexing() {
        let mut sm = SymmetricMatrix::new(3);
        assert_eq!(sm.dim(), 3);
        sm.set(0, 1, 2.5);
        assert!((sm.get(1, 0).unwrap() - 2.5).abs() < 1e-10);
        assert!((sm.get(0, 1).unwrap() - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_get_raw_data() {
        let mut sm = SymmetricMatrix::new(3);
        sm.set(0, 0, 1.0);
        sm.set(1, 0, 2.0);
        sm.set(1, 1, 3.0);
        let raw = sm.get_raw_data();
        assert_eq!(raw.len(), 6);
        assert_eq!(raw, vec![1.0, 2.0, 3.0, 0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_eigenvalues() {
        let mut sm = SymmetricMatrix::new(2);
        sm.set(0, 0, 2.0);
        sm.set(1, 1, 2.0);
        sm.set(0, 1, 1.0);
        let (eigvals, _eigvecs) = sm.eig().unwrap();
        assert_eq!(eigvals.len(), 2);
        assert!((eigvals[0] - 1.0).abs() < 1e-10);
        assert!((eigvals[1] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_identity_matrix() {
        let id = identity_matrix(3);
        assert_eq!(id.get(0, 0).unwrap(), 1.0);
        assert_eq!(id.get(1, 1).unwrap(), 1.0);
        assert_eq!(id.get(0, 1).unwrap(), 0.0);
    }
}
