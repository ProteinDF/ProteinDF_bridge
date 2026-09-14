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
use std::ops::{
    Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign,
};
use std::str::FromStr;

use crate::error::{BridgeError, Result};
use crate::matrix::Matrix;
use crate::vector::Vector;

pub const DEFAULT_EPSILON: f64 = 1.0e-5;

/// 3D coordinate vector corresponding to `proteindf_bridge.position.Position`.
#[derive(Debug, Clone, Copy)]
pub struct Position {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub epsilon: f64,
}

impl Default for Position {
    fn default() -> Self {
        Self::new(0.0, 0.0, 0.0)
    }
}

impl Position {
    /// Creates a new `Position` with the given coordinates and default epsilon (1.0e-5).
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self {
            x,
            y,
            z,
            epsilon: DEFAULT_EPSILON,
        }
    }

    /// Creates a `Position` from a slice of at least 3 elements.
    pub fn from_slice(slice: &[f64]) -> Result<Self> {
        if slice.len() != 3 {
            return Err(BridgeError::input_error(
                format!("{:?}", slice),
                "position requires exactly 3 elements",
            ));
        }
        Ok(Self::new(slice[0], slice[1], slice[2]))
    }

    /// Creates a `Position` from a `Vector`.
    pub fn from_vector(v: &Vector) -> Result<Self> {
        Self::from_slice(v.as_slice())
    }

    /// Returns the coordinates as an array `[x, y, z]`.
    pub fn xyz(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }

    /// Returns the raw data as `[x, y, z]`.
    pub fn get_raw_data(&self) -> [f64; 3] {
        self.xyz()
    }

    /// Updates coordinates to match another position.
    pub fn move_to(&mut self, other: Position) {
        self.x = other.x;
        self.y = other.y;
        self.z = other.z;
    }

    /// Computes the squared Euclidean distance to another position.
    pub fn square_distance_from(&self, other: &Position) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        dx * dx + dy * dy + dz * dz
    }

    /// Computes the Euclidean distance to another position.
    pub fn distance_from(&self, other: &Position) -> f64 {
        self.square_distance_from(other).sqrt()
    }

    /// Normalizes this position vector in-place.
    /// Returns `BridgeError::ValueError` if the vector length is near zero (< 1e-15).
    pub fn norm(&mut self) -> Result<&mut Self> {
        let len = self.length();
        if len < 1.0e-15 {
            return Err(BridgeError::value_error(
                "Position.norm()",
                "cannot normalize zero vector",
            ));
        }
        self.x /= len;
        self.y /= len;
        self.z /= len;
        Ok(self)
    }

    /// Computes the Euclidean norm (length) of the vector.
    pub fn length(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Rotates this position by a 3x3 `Matrix`.
    pub fn rotate(&mut self, mat: &Matrix) -> Result<()> {
        if mat.rows() != 3 || mat.cols() != 3 {
            return Err(BridgeError::value_error(
                format!("({}, {})", mat.rows(), mat.cols()),
                "Rotation matrix must be 3x3",
            ));
        }
        let v = Vector::from_slice(&self.xyz());
        let rotated = mat * &v;
        self.x = rotated[0];
        self.y = rotated[1];
        self.z = rotated[2];
        Ok(())
    }

    /// Computes the dot product with another position.
    pub fn dot(&self, other: &Position) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Computes the cross product with another position (`self x other`).
    pub fn cross(&self, other: &Position) -> Position {
        Position::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }
}

impl FromStr for Position {
    type Err = BridgeError;

    /// Parses a string of numbers (space- or comma-separated, e.g. "1.0 2.0 -3.0" or "1.0, 2.0, -3.0").
    fn from_str(s: &str) -> Result<Self> {
        let cleaned = s.replace(',', " ");
        let tokens: Vec<&str> = cleaned.split_whitespace().collect();
        if tokens.len() < 3 {
            return Err(BridgeError::input_error(
                s,
                "expected at least 3 coordinates",
            ));
        }
        let x = tokens[0].parse::<f64>().map_err(|_| {
            BridgeError::input_error(tokens[0], "failed to parse coordinate as f64")
        })?;
        let y = tokens[1].parse::<f64>().map_err(|_| {
            BridgeError::input_error(tokens[1], "failed to parse coordinate as f64")
        })?;
        let z = tokens[2].parse::<f64>().map_err(|_| {
            BridgeError::input_error(tokens[2], "failed to parse coordinate as f64")
        })?;

        Ok(Position::new(x, y, z))
    }
}

impl Index<usize> for Position {
    type Output = f64;
    fn index(&self, index: usize) -> &Self::Output {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Position index out of bounds: {}", index),
        }
    }
}

impl IndexMut<usize> for Position {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => panic!("Position index out of bounds: {}", index),
        }
    }
}

// Equality with epsilon matching Python's `self.distance_from(rhs) < self.epsilon`.
impl PartialEq for Position {
    fn eq(&self, other: &Self) -> bool {
        self.distance_from(other) < self.epsilon
    }
}

impl Add for Position {
    type Output = Position;
    fn add(self, rhs: Self) -> Self::Output {
        Position::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl Add for &Position {
    type Output = Position;
    fn add(self, rhs: Self) -> Self::Output {
        Position::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl AddAssign for Position {
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl Sub for Position {
    type Output = Position;
    fn sub(self, rhs: Self) -> Self::Output {
        Position::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl Sub for &Position {
    type Output = Position;
    fn sub(self, rhs: Self) -> Self::Output {
        Position::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl SubAssign for Position {
    fn sub_assign(&mut self, rhs: Self) {
        self.x -= rhs.x;
        self.y -= rhs.y;
        self.z -= rhs.z;
    }
}

// Position * scalar
impl Mul<f64> for Position {
    type Output = Position;
    fn mul(self, rhs: f64) -> Self::Output {
        Position::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

impl Mul<f64> for &Position {
    type Output = Position;
    fn mul(self, rhs: f64) -> Self::Output {
        Position::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

impl Mul<Position> for f64 {
    type Output = Position;
    fn mul(self, rhs: Position) -> Self::Output {
        rhs * self
    }
}

impl Mul<&Position> for f64 {
    type Output = Position;
    fn mul(self, rhs: &Position) -> Self::Output {
        rhs * self
    }
}

impl MulAssign<f64> for Position {
    fn mul_assign(&mut self, rhs: f64) {
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
    }
}

// Position * Position -> f64 (dot product, matching Python's `Position.__mul__`)
impl Mul<Position> for Position {
    type Output = f64;
    fn mul(self, rhs: Position) -> Self::Output {
        self.dot(&rhs)
    }
}

impl Mul<&Position> for &Position {
    type Output = f64;
    fn mul(self, rhs: &Position) -> Self::Output {
        self.dot(rhs)
    }
}

impl Div<f64> for Position {
    type Output = Position;
    fn div(self, rhs: f64) -> Self::Output {
        Position::new(self.x / rhs, self.y / rhs, self.z / rhs)
    }
}

impl Div<f64> for &Position {
    type Output = Position;
    fn div(self, rhs: f64) -> Self::Output {
        Position::new(self.x / rhs, self.y / rhs, self.z / rhs)
    }
}

impl DivAssign<f64> for Position {
    fn div_assign(&mut self, rhs: f64) {
        self.x /= rhs;
        self.y /= rhs;
        self.z /= rhs;
    }
}

impl Neg for Position {
    type Output = Position;
    fn neg(self) -> Self::Output {
        Position::new(-self.x, -self.y, -self.z)
    }
}

impl Neg for &Position {
    type Output = Position;
    fn neg(self) -> Self::Output {
        Position::new(-self.x, -self.y, -self.z)
    }
}

impl fmt::Display for Position {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({:10.6}, {:10.6}, {:10.6})", self.x, self.y, self.z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Ported from tests/test_position.py
    #[test]
    fn test_init1() {
        let pos1 = Position::default();
        assert!((pos1.x - 0.0).abs() < 1e-10);
        assert!((pos1.y - 0.0).abs() < 1e-10);
        assert!((pos1.z - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_init2() {
        let pos1 = Position::from_slice(&[1.0, 2.0, -3.0]).unwrap();
        assert!((pos1.x - 1.0).abs() < 1e-10);
        assert!((pos1.y - 2.0).abs() < 1e-10);
        assert!((pos1.z - (-3.0)).abs() < 1e-10);
    }

    #[test]
    fn test_init3() {
        let pos1: Position = "1.0  2.0 -3.0".parse().unwrap();
        assert!((pos1.x - 1.0).abs() < 1e-10);
        assert!((pos1.y - 2.0).abs() < 1e-10);
        assert!((pos1.z - (-3.0)).abs() < 1e-10);

        let pos2: Position = "1.0,  2.0,-3.0".parse().unwrap();
        assert!((pos2.x - 1.0).abs() < 1e-10);
        assert!((pos2.y - 2.0).abs() < 1e-10);
        assert!((pos2.z - (-3.0)).abs() < 1e-10);
    }

    #[test]
    fn test_norm_zero_vector() {
        let mut pos = Position::new(0.0, 0.0, 0.0);
        assert!(matches!(pos.norm(), Err(BridgeError::ValueError { .. })));
    }

    // Doctests from proteindf_bridge/position.py
    #[test]
    fn test_doctest_operations() {
        let mut p = Position::new(0.0, 1.0, 2.0);
        assert!((p.length() - 2.23606).abs() < 1e-4);

        p.norm().unwrap();
        let sqrt5 = 5.0_f64.sqrt();
        assert_eq!(p, Position::new(0.0, 1.0 / sqrt5, 2.0 / sqrt5));

        p.move_to(Position::new(3.0, 4.0, 5.0));
        assert_eq!(p, Position::new(3.0, 4.0, 5.0));

        let tmp = p * 2.0;
        assert_eq!(tmp, Position::new(6.0, 8.0, 10.0));

        let tmp_neg = -1.0 * p;
        assert_eq!(tmp_neg, Position::new(-3.0, -4.0, -5.0));

        let mut a = Position::new(1.0, 2.0, 3.0);
        let b = Position::new(2.0, 3.0, 4.0);
        assert!(((a * b) - 20.0).abs() < 1e-10);

        let sum = a + b;
        assert_eq!(sum, Position::new(3.0, 5.0, 7.0));

        a += b;
        assert_eq!(a, Position::new(3.0, 5.0, 7.0));

        let diff = a - b;
        assert_eq!(diff, Position::new(1.0, 2.0, 3.0));

        a -= b;
        assert_eq!(a, Position::new(1.0, 2.0, 3.0));

        assert!((a.dot(&b) - 20.0).abs() < 1e-10);

        let cross = a.cross(&b);
        assert_eq!(cross, Position::new(-1.0, 2.0, -1.0));
    }

    // From proteindf_bridge/position_test.py
    #[test]
    fn test_position_test_py() {
        let p = Position::from_slice(&[5.0, 3.0, -1.2]).unwrap();
        assert!((p.x - 5.0).abs() < 1e-5);
    }
}
