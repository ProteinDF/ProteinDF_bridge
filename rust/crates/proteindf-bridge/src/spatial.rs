// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::collections::HashMap;

use crate::position::Position;

/// A uniform grid spatial partition for 3D coordinates.
///
/// Divides 3D space into cubic cells of side length `cell_size`.
/// Enables efficient $O(N)$ spatial querying such as neighbor pair enumeration
/// without building hierarchical tree structures.
#[derive(Debug, Clone)]
pub struct CellList {
    cell_size: f64,
    grid: HashMap<(i32, i32, i32), Vec<usize>>,
    positions: Vec<Position>,
}

impl CellList {
    /// Creates a new `CellList` from a slice of 3D positions and a cubic cell size.
    ///
    /// # Panics
    /// Panics if `cell_size <= 0.0` or not finite.
    pub fn new(positions: &[Position], cell_size: f64) -> Self {
        assert!(
            cell_size > 0.0 && cell_size.is_finite(),
            "cell_size must be positive and finite"
        );

        let mut grid: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();
        for (idx, pos) in positions.iter().enumerate() {
            let key = Self::cell_key(pos, cell_size);
            grid.entry(key).or_default().push(idx);
        }

        Self {
            cell_size,
            grid,
            positions: positions.to_vec(),
        }
    }

    /// Returns the cell size.
    pub fn cell_size(&self) -> f64 {
        self.cell_size
    }

    /// Returns the number of positions indexed in this cell list.
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// Returns `true` if the cell list contains no positions.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    /// Computes the integer grid coordinate for a given position.
    #[inline]
    fn cell_key(pos: &Position, cell_size: f64) -> (i32, i32, i32) {
        let ix = (pos.x / cell_size).floor() as i32;
        let iy = (pos.y / cell_size).floor() as i32;
        let iz = (pos.z / cell_size).floor() as i32;
        (ix, iy, iz)
    }

    /// Iterates over all candidate neighbor pairs `(i, j)` where `i < j` within the
    /// same or adjacent 26 cells, invoking `callback(i, j, distance)`.
    ///
    /// The callback receives:
    /// - `i`: Index of the first position
    /// - `j`: Index of the second position, always satisfying `i < j` (no self pairs, no duplicate reverse pairs)
    /// - `distance`: Euclidean distance between `positions[i]` and `positions[j]`
    ///
    /// Only pairs with `distance <= max_radius` are reported.
    ///
    /// # Note
    /// For completeness, `max_radius` should satisfy `max_radius <= self.cell_size`.
    /// If `max_radius > self.cell_size`, pairs beyond the 27 neighboring cells might not be checked.
    pub fn for_each_neighbor_pair<F>(&self, max_radius: f64, mut callback: F)
    where
        F: FnMut(usize, usize, f64),
    {
        if self.positions.is_empty() {
            return;
        }

        let max_r2 = max_radius * max_radius;

        // 13 relative offsets for half-neighborhood.
        // Using dictionary-order positive directions (z > 0, or z == 0 && y > 0, or z == 0 && y == 0 && x > 0)
        // ensures that each adjacent cell pair is checked exactly once with i < j ordering.
        const HALF_NEIGHBORHOOD: [(i32, i32, i32); 13] = [
            (1, 0, 0),
            (-1, 1, 0),
            (0, 1, 0),
            (1, 1, 0),
            (-1, -1, 1),
            (0, -1, 1),
            (1, -1, 1),
            (-1, 0, 1),
            (0, 0, 1),
            (1, 0, 1),
            (-1, 1, 1),
            (0, 1, 1),
            (1, 1, 1),
        ];

        for (&(cx, cy, cz), cell_atoms) in &self.grid {
            // 1. Within the same cell: pairs (i, j) with i < j
            let n = cell_atoms.len();
            for a in 0..n {
                let idx_i = cell_atoms[a];
                let pos_i = &self.positions[idx_i];
                for &idx_j in &cell_atoms[a + 1..n] {
                    let pos_j = &self.positions[idx_j];
                    let dx = pos_i.x - pos_j.x;
                    let dy = pos_i.y - pos_j.y;
                    let dz = pos_i.z - pos_j.z;
                    let d2 = dx * dx + dy * dy + dz * dz;
                    if d2 <= max_r2 {
                        let (i, j) = if idx_i < idx_j {
                            (idx_i, idx_j)
                        } else {
                            (idx_j, idx_i)
                        };
                        callback(i, j, d2.sqrt());
                    }
                }
            }

            // 2. Between this cell and forward neighbor cells
            for &(dx, dy, dz) in &HALF_NEIGHBORHOOD {
                let neighbor_key = (cx + dx, cy + dy, cz + dz);
                if let Some(neighbor_atoms) = self.grid.get(&neighbor_key) {
                    for &idx_i in cell_atoms {
                        let pos_i = &self.positions[idx_i];
                        for &idx_j in neighbor_atoms {
                            let pos_j = &self.positions[idx_j];
                            let dx = pos_i.x - pos_j.x;
                            let dy = pos_i.y - pos_j.y;
                            let dz = pos_i.z - pos_j.z;
                            let d2 = dx * dx + dy * dy + dz * dz;
                            if d2 <= max_r2 {
                                let (i, j) = if idx_i < idx_j {
                                    (idx_i, idx_j)
                                } else {
                                    (idx_j, idx_i)
                                };
                                callback(i, j, d2.sqrt());
                            }
                        }
                    }
                }
            }
        }
    }

    /// Queries all pairs `(i, j)` where `i < j` with distance `<= radius`.
    ///
    /// Returns a `Vec` of `(i, j, distance)` tuples.
    pub fn query_pairs_within(&self, radius: f64) -> Vec<(usize, usize, f64)> {
        let mut pairs = Vec::new();
        self.for_each_neighbor_pair(radius, |i, j, d| {
            pairs.push((i, j, d));
        });
        pairs
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cell_list_basic() {
        let positions = vec![
            Position::new(0.0, 0.0, 0.0),
            Position::new(1.0, 0.0, 0.0),
            Position::new(5.0, 0.0, 0.0),
            Position::new(10.0, 0.0, 0.0),
        ];

        let cl = CellList::new(&positions, 3.0);
        assert_eq!(cl.len(), 4);
        assert!(!cl.is_empty());

        let pairs = cl.query_pairs_within(2.0);
        // Only (0, 1) has distance <= 2.0 (distance 1.0)
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].0, 0);
        assert_eq!(pairs[0].1, 1);
        assert!((pairs[0].2 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_cell_list_empty() {
        let cl = CellList::new(&[], 5.0);
        assert_eq!(cl.len(), 0);
        assert!(cl.is_empty());
        assert!(cl.query_pairs_within(5.0).is_empty());
    }

    #[test]
    fn test_cell_list_against_brute_force() {
        // Generate a grid of points
        let mut positions = Vec::new();
        for x in 0..5 {
            for y in 0..5 {
                for z in 0..5 {
                    positions.push(Position::new(
                        x as f64 * 1.5,
                        y as f64 * 1.5,
                        z as f64 * 1.5,
                    ));
                }
            }
        }

        let radius = 2.0;
        let cl = CellList::new(&positions, radius);
        let mut cell_pairs = cl.query_pairs_within(radius);
        cell_pairs.sort_by_key(|&(i, j, _)| (i, j));

        // Brute force pairs
        let mut brute_pairs = Vec::new();
        for i in 0..positions.len() {
            for j in i + 1..positions.len() {
                let d = positions[i].distance_from(&positions[j]);
                if d <= radius {
                    brute_pairs.push((i, j, d));
                }
            }
        }
        brute_pairs.sort_by_key(|&(i, j, _)| (i, j));

        assert_eq!(cell_pairs.len(), brute_pairs.len());
        for (cp, bp) in cell_pairs.iter().zip(brute_pairs.iter()) {
            assert_eq!(cp.0, bp.0);
            assert_eq!(cp.1, bp.1);
            assert!((cp.2 - bp.2).abs() < 1e-10);
        }
    }
}
