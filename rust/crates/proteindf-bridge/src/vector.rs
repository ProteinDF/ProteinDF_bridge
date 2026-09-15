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

use std::fmt;
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::error::{BridgeError, Result};

/// A 1D numerical vector corresponding to `proteindf_bridge.vector.Vector`.
#[derive(Debug, Clone, PartialEq)]
pub struct Vector {
    data: Vec<f64>,
}

impl Vector {
    /// Creates a new vector of the given size, initialized to 0.0.
    pub fn new(size: usize) -> Self {
        Self {
            data: vec![0.0; size],
        }
    }

    /// Creates a vector from a slice of `f64`.
    pub fn from_slice(slice: &[f64]) -> Self {
        Self {
            data: slice.to_vec(),
        }
    }

    /// Creates a vector from a `Vec<f64>`.
    pub fn from_vec(vec: Vec<f64>) -> Self {
        Self { data: vec }
    }

    /// Returns the number of elements in the vector.
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Returns true if the vector contains no elements.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Returns the number of elements (alias for `len`).
    pub fn size(&self) -> usize {
        self.len()
    }

    /// Resizes the vector in-place. If expanded, new elements are initialized to 0.0.
    pub fn resize(&mut self, new_size: usize) {
        self.data.resize(new_size, 0.0);
    }

    /// Gets the value at the given index.
    pub fn get(&self, index: usize) -> Result<f64> {
        self.data
            .get(index)
            .copied()
            .ok_or_else(|| BridgeError::value_error(index.to_string(), "index out of range"))
    }

    /// Sets the value at the given index.
    pub fn set(&mut self, index: usize, value: f64) {
        assert!(index < self.len(), "index out of range: {}", index);
        self.data[index] = value;
    }

    /// Returns the maximum value in the vector.
    pub fn max(&self) -> f64 {
        assert!(!self.data.is_empty(), "cannot compute max of empty Vector");
        self.data.iter().copied().fold(f64::NEG_INFINITY, f64::max)
    }

    /// Returns the minimum value in the vector.
    pub fn min(&self) -> f64 {
        assert!(!self.data.is_empty(), "cannot compute min of empty Vector");
        self.data.iter().copied().fold(f64::INFINITY, f64::min)
    }

    /// Returns a new vector containing the absolute value of each element.
    pub fn abs(&self) -> Self {
        Self {
            data: self.data.iter().map(|&x| x.abs()).collect(),
        }
    }

    /// Computes the dot product with another vector.
    pub fn dot(&self, other: &Vector) -> f64 {
        assert_eq!(
            self.len(),
            other.len(),
            "Vector lengths must match for dot product: {} != {}",
            self.len(),
            other.len()
        );
        self.data
            .iter()
            .zip(other.data.iter())
            .map(|(&a, &b)| a * b)
            .sum()
    }

    /// Returns the underlying elements as a `Vec<f64>`.
    pub fn to_vec(&self) -> Vec<f64> {
        self.data.clone()
    }

    /// Returns a slice of the elements.
    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }

    /// Returns the indices that would sort this vector in ascending order.
    pub fn argsort(&self) -> Vec<usize> {
        let mut indices: Vec<usize> = (0..self.len()).collect();
        indices.sort_by(|&a, &b| {
            self.data[a]
                .partial_cmp(&self.data[b])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        indices
    }

    /// Returns a new vector with elements in reverse order.
    pub fn flip(&self) -> Self {
        let mut rev = self.data.clone();
        rev.reverse();
        Self { data: rev }
    }
}

impl From<Vec<f64>> for Vector {
    fn from(vec: Vec<f64>) -> Self {
        Self::from_vec(vec)
    }
}

impl From<&[f64]> for Vector {
    fn from(slice: &[f64]) -> Self {
        Self::from_slice(slice)
    }
}

impl Index<usize> for Vector {
    type Output = f64;
    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl IndexMut<usize> for Vector {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index]
    }
}

impl Add for Vector {
    type Output = Vector;
    fn add(self, rhs: Self) -> Self::Output {
        &self + &rhs
    }
}

impl Add for &Vector {
    type Output = Vector;
    fn add(self, rhs: Self) -> Self::Output {
        assert_eq!(self.len(), rhs.len(), "Vector lengths must match for add");
        Vector {
            data: self
                .data
                .iter()
                .zip(rhs.data.iter())
                .map(|(&a, &b)| a + b)
                .collect(),
        }
    }
}

impl AddAssign<&Vector> for Vector {
    fn add_assign(&mut self, rhs: &Vector) {
        assert_eq!(
            self.len(),
            rhs.len(),
            "Vector lengths must match for add_assign"
        );
        for (a, &b) in self.data.iter_mut().zip(rhs.data.iter()) {
            *a += b;
        }
    }
}

impl AddAssign for Vector {
    fn add_assign(&mut self, rhs: Vector) {
        *self += &rhs;
    }
}

impl Sub for Vector {
    type Output = Vector;
    fn sub(self, rhs: Self) -> Self::Output {
        &self - &rhs
    }
}

impl Sub for &Vector {
    type Output = Vector;
    fn sub(self, rhs: Self) -> Self::Output {
        assert_eq!(self.len(), rhs.len(), "Vector lengths must match for sub");
        Vector {
            data: self
                .data
                .iter()
                .zip(rhs.data.iter())
                .map(|(&a, &b)| a - b)
                .collect(),
        }
    }
}

impl SubAssign<&Vector> for Vector {
    fn sub_assign(&mut self, rhs: &Vector) {
        assert_eq!(
            self.len(),
            rhs.len(),
            "Vector lengths must match for sub_assign"
        );
        for (a, &b) in self.data.iter_mut().zip(rhs.data.iter()) {
            *a -= b;
        }
    }
}

impl SubAssign for Vector {
    fn sub_assign(&mut self, rhs: Vector) {
        *self -= &rhs;
    }
}

// Vector * scalar
impl Mul<f64> for Vector {
    type Output = Vector;
    fn mul(self, rhs: f64) -> Self::Output {
        &self * rhs
    }
}

impl Mul<f64> for &Vector {
    type Output = Vector;
    fn mul(self, rhs: f64) -> Self::Output {
        Vector {
            data: self.data.iter().map(|&x| x * rhs).collect(),
        }
    }
}

impl Mul<Vector> for f64 {
    type Output = Vector;
    fn mul(self, rhs: Vector) -> Self::Output {
        &rhs * self
    }
}

impl Mul<&Vector> for f64 {
    type Output = Vector;
    fn mul(self, rhs: &Vector) -> Self::Output {
        rhs * self
    }
}

impl MulAssign<f64> for Vector {
    fn mul_assign(&mut self, rhs: f64) {
        for x in self.data.iter_mut() {
            *x *= rhs;
        }
    }
}

// Vector * Vector -> f64 (dot product, matching Python `Vector.__mul__`)
impl Mul<&Vector> for &Vector {
    type Output = f64;
    fn mul(self, rhs: &Vector) -> Self::Output {
        self.dot(rhs)
    }
}

impl Mul<Vector> for Vector {
    type Output = f64;
    fn mul(self, rhs: Vector) -> Self::Output {
        self.dot(&rhs)
    }
}

impl Neg for Vector {
    type Output = Vector;
    fn neg(self) -> Self::Output {
        -&self
    }
}

impl Neg for &Vector {
    type Output = Vector;
    fn neg(self) -> Self::Output {
        Vector {
            data: self.data.iter().map(|&x| -x).collect(),
        }
    }
}

impl fmt::Display for Vector {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let len = self.len();
        for order in (0..len).step_by(10) {
            writeln!(f)?;
            for j in order..std::cmp::min(order + 10, len) {
                write!(f, "   {:5} th", j + 1)?;
            }
            writeln!(f)?;
            for _ in order..std::cmp::min(order + 10, len) {
                write!(f, "-----------")?;
            }
            writeln!(f, "----")?;
            for j in order..std::cmp::min(order + 10, len) {
                write!(f, " {:10.6}", self.data[j])?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Ported from tests/test_vector.py
    #[test]
    fn test_init_with_size() {
        let v = Vector::new(5);
        assert_eq!(v.len(), 5);
        for i in 0..5 {
            assert_eq!(v[i], 0.0);
        }
    }

    #[test]
    fn test_init_with_list() {
        let v = Vector::from_slice(&[1.0, 2.0, 3.0]);
        assert_eq!(v.len(), 3);
        assert_eq!(v[0], 1.0);
        assert_eq!(v[1], 2.0);
        assert_eq!(v[2], 3.0);
    }

    #[test]
    fn test_set_and_get() {
        let mut v = Vector::new(3);
        v.set(1, 4.5);
        assert_eq!(v.get(1).unwrap(), 4.5);
        v[2] = 9.0;
        assert_eq!(v[2], 9.0);
    }

    #[test]
    fn test_resize() {
        let mut v = Vector::from_slice(&[1.0, 2.0, 3.0]);
        v.resize(5);
        assert_eq!(v.len(), 5);
        assert_eq!(v[0], 1.0);
        assert_eq!(v[3], 0.0);
    }

    #[test]
    fn test_dot_and_abs() {
        let v1 = Vector::from_slice(&[3.0, 4.0]);
        let v1_abs = v1.abs();
        assert_eq!(v1_abs.to_vec(), vec![3.0, 4.0]);

        let v2 = Vector::from_slice(&[1.0, 2.0]);
        let v3 = Vector::from_slice(&[3.0, 4.0]);
        let dot = &v2 * &v3; // Vector multiplication computes dot product
        assert!((dot - 11.0).abs() < 1e-10);
    }

    #[test]
    fn test_arithmetic_operations() {
        let v1 = Vector::from_slice(&[1.0, 2.0, 3.0]);
        let v2 = Vector::from_slice(&[4.0, 5.0, 6.0]);

        let v_add = &v1 + &v2;
        assert_eq!(v_add.to_vec(), vec![5.0, 7.0, 9.0]);

        let v_sub = &v2 - &v1;
        assert_eq!(v_sub.to_vec(), vec![3.0, 3.0, 3.0]);

        let v_mul = &v1 * 2.0;
        assert_eq!(v_mul.to_vec(), vec![2.0, 4.0, 6.0]);

        let v_neg = -&v1;
        assert_eq!(v_neg.to_vec(), vec![-1.0, -2.0, -3.0]);
    }
}
