// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom::Atom;
use crate::atom_group::AtomGroup;
use crate::brd::sort_nicely;
use crate::error::{BridgeError, Result};
use crate::matrix::SymmetricMatrix;
use crate::position::Position;

/// Default maximum distance between carbon atom and ring centroid (4.5 Å).
pub const DEFAULT_MAX_DISTANCE: f64 = 4.5;

/// Default maximum angle between centroid-carbon vector and ring normal (40.0°).
pub const DEFAULT_MAX_ANGLE_DEG: f64 = 40.0;

/// Definition of an aromatic ring by its residue name and constituent atom names.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AromaticRingDef {
    pub residue_name: &'static str,
    pub atom_names: &'static [&'static str],
}

/// Static table of aromatic ring definitions.
pub const AROMATIC_RINGS: &[AromaticRingDef] = &[
    AromaticRingDef {
        residue_name: "PHE",
        atom_names: &["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    },
    AromaticRingDef {
        residue_name: "TYR",
        atom_names: &["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    },
    AromaticRingDef {
        residue_name: "HIS",
        atom_names: &["CG", "ND1", "CD2", "CE1", "NE2"],
    },
    AromaticRingDef {
        residue_name: "TRP",
        atom_names: &["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    },
];

/// Represents an aromatic ring with its centroid and normal vector.
#[derive(Debug, Clone, PartialEq)]
pub struct AromaticRing {
    pub residue_path: String,
    pub residue_name: String,
    pub center: Position,
    pub normal: Position,
}

impl AromaticRing {
    pub fn new(
        residue_path: impl Into<String>,
        residue_name: impl Into<String>,
        center: Position,
        normal: Position,
    ) -> Self {
        Self {
            residue_path: residue_path.into(),
            residue_name: residue_name.into(),
            center,
            normal,
        }
    }
}

/// Represents a detected CH-pi interaction between a carbon atom and an aromatic ring.
#[derive(Debug, Clone, PartialEq)]
pub struct ChPiInteraction {
    pub carbon_path: String,
    pub carbon_atom: String,
    pub ring_path: String,
    pub ring_residue: String,
    pub distance: f64,
    pub angle: f64,
}

impl ChPiInteraction {
    pub fn new(
        carbon_path: impl Into<String>,
        carbon_atom: impl Into<String>,
        ring_path: impl Into<String>,
        ring_residue: impl Into<String>,
        distance: f64,
        angle: f64,
    ) -> Self {
        Self {
            carbon_path: carbon_path.into(),
            carbon_atom: carbon_atom.into(),
            ring_path: ring_path.into(),
            ring_residue: ring_residue.into(),
            distance,
            angle,
        }
    }
}

/// Calculates ring centroid and unit normal vector by least-squares plane fitting.
///
/// Uses Jacobi eigenvalue decomposition of the centered covariance matrix (equivalent to SVD).
pub fn calc_ring_center_and_normal(coords: &[Position]) -> Result<(Position, Position)> {
    let n = coords.len();
    if n < 3 {
        return Err(BridgeError::InputError {
            expr: "calc_ring_center_and_normal".to_string(),
            msg: "At least 3 atoms required to define a ring plane".to_string(),
        });
    }

    // 1. Centroid
    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    let mut sum_z = 0.0;
    for p in coords {
        sum_x += p.x;
        sum_y += p.y;
        sum_z += p.z;
    }
    let center = Position::new(sum_x / n as f64, sum_y / n as f64, sum_z / n as f64);

    // 2. Covariance matrix M = sum d_k * d_k^T
    let mut m00 = 0.0;
    let mut m01 = 0.0;
    let mut m02 = 0.0;
    let mut m11 = 0.0;
    let mut m12 = 0.0;
    let mut m22 = 0.0;

    for p in coords {
        let dx = p.x - center.x;
        let dy = p.y - center.y;
        let dz = p.z - center.z;
        m00 += dx * dx;
        m01 += dx * dy;
        m02 += dx * dz;
        m11 += dy * dy;
        m12 += dy * dz;
        m22 += dz * dz;
    }

    let mut cov = SymmetricMatrix::new(3);
    cov.set(0, 0, m00);
    cov.set(0, 1, m01);
    cov.set(0, 2, m02);
    cov.set(1, 1, m11);
    cov.set(1, 2, m12);
    cov.set(2, 2, m22);

    // Eigenvalues are sorted in ascending order; row 0 of eigenvectors matrix is the normal
    let (_eigvals, eigvecs) = cov.eig()?;
    let nx = eigvecs.get(0, 0)?;
    let ny = eigvecs.get(0, 1)?;
    let nz = eigvecs.get(0, 2)?;
    let norm = (nx * nx + ny * ny + nz * nz).sqrt();
    if norm < 1e-12 {
        return Err(BridgeError::InputError {
            expr: "calc_ring_center_and_normal".to_string(),
            msg: "Degenerate ring coordinates, cannot compute normal".to_string(),
        });
    }

    let normal = Position::new(nx / norm, ny / norm, nz / norm);
    Ok((center, normal))
}

/// Computes centroid and unit normal vector for an aromatic ring in a residue.
/// Returns `None` if any of the required ring atoms is missing.
pub fn calc_ring_geometry(res: &AtomGroup, ring_atoms: &[&str]) -> Option<(Position, Position)> {
    let mut coords = Vec::with_capacity(ring_atoms.len());
    for &aname in ring_atoms {
        let atom = res.get_atom(aname)?;
        coords.push(atom.xyz);
    }
    calc_ring_center_and_normal(&coords).ok()
}

fn collect_leaf_residues<'a>(group: &'a AtomGroup, out: &mut Vec<&'a AtomGroup>) {
    if group.atoms().next().is_some() && group.groups().next().is_none() {
        out.push(group);
    } else {
        let mut keys = group.get_group_list();
        sort_nicely(&mut keys);
        for k in keys {
            if let Some(child) = group.get_group(&k) {
                collect_leaf_residues(child, out);
            }
        }
    }
}

fn is_carbon_atom(atom: &Atom) -> bool {
    atom.atomic_number() == 6
}

/// Detects CH-pi interactions using default thresholds (distance <= 4.5 Å, angle <= 40.0°).
pub fn calc_ch_pi_interactions(root: &AtomGroup) -> Vec<ChPiInteraction> {
    calc_ch_pi_interactions_with_thresholds(root, DEFAULT_MAX_DISTANCE, DEFAULT_MAX_ANGLE_DEG)
}

/// Detects CH-pi interactions with custom distance and angle thresholds.
///
/// Rules:
/// - Aromatic rings: PHE, TYR (6-membered ring), HIS (5-membered ring), TRP (6-membered ring).
/// - Carbon atoms: any atom in other residues with atomic number 6 / symbol 'C'.
/// - Intra-residue pairs (carbon in the same residue as the ring) are excluded.
/// - Distance: distance between carbon and ring centroid <= `max_distance`.
/// - Angle: angle between centroid-to-carbon vector and ring unit normal <= `max_angle_deg`.
///   (Uses absolute value of dot product to be independent of normal sign).
pub fn calc_ch_pi_interactions_with_thresholds(
    root: &AtomGroup,
    max_distance: f64,
    max_angle_deg: f64,
) -> Vec<ChPiInteraction> {
    let mut residues = Vec::new();
    collect_leaf_residues(root, &mut residues);

    // 1. Identify all aromatic rings
    let mut rings = Vec::new();
    for res in &residues {
        let rname = res.name.trim();
        for def in AROMATIC_RINGS {
            if def.residue_name == rname {
                if let Some((center, normal)) = calc_ring_geometry(res, def.atom_names) {
                    rings.push(AromaticRing::new(
                        res.path().to_string(),
                        rname,
                        center,
                        normal,
                    ));
                }
            }
        }
    }

    // 2. Identify all carbon atoms
    struct CarbonCandidate<'a> {
        res: &'a AtomGroup,
        atom_name: String,
        pos: Position,
    }

    let mut carbons = Vec::new();
    for res in &residues {
        for (_, atom) in res.atoms() {
            if is_carbon_atom(atom) {
                carbons.push(CarbonCandidate {
                    res,
                    atom_name: atom.name.clone(),
                    pos: atom.xyz,
                });
            }
        }
    }

    let mut interactions = Vec::new();

    // 3. Evaluate pairs
    for ring in &rings {
        for carbon in &carbons {
            // Exclude same residue
            if ring.residue_path == carbon.res.path() {
                continue;
            }

            let v = carbon.pos - ring.center;
            let distance = v.length();
            if distance > max_distance || distance < 1e-6 {
                continue;
            }

            let v_unit = v / distance;
            let normal_unit = ring.normal / ring.normal.length();
            let cos_val = v_unit.dot(&normal_unit).abs().clamp(-1.0, 1.0);
            let angle = cos_val.acos().to_degrees();

            if angle <= max_angle_deg {
                interactions.push(ChPiInteraction::new(
                    carbon.res.path().to_string(),
                    carbon.atom_name.clone(),
                    ring.residue_path.clone(),
                    ring.residue_name.clone(),
                    distance,
                    angle,
                ));
            }
        }
    }

    interactions
}
