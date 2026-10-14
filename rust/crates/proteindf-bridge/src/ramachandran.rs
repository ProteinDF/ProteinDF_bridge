// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom_group::AtomGroup;
use crate::position::{dihedral_angle, Position};

/// Represents the Ramachandran backbone dihedral angles (phi and psi) for a residue.
#[derive(Debug, Clone, PartialEq)]
pub struct RamachandranAngle {
    pub residue_key: String,
    pub residue_name: String,
    pub phi: Option<f64>,
    pub psi: Option<f64>,
}

impl RamachandranAngle {
    /// Creates a new `RamachandranAngle` record.
    pub fn new(
        residue_key: impl Into<String>,
        residue_name: impl Into<String>,
        phi: Option<f64>,
        psi: Option<f64>,
    ) -> Self {
        Self {
            residue_key: residue_key.into(),
            residue_name: residue_name.into(),
            phi,
            psi,
        }
    }
}

/// Maximum distance in Angstroms between C(i-1) and N(i) to be considered a continuous peptide bond.
pub const PEPTIDE_BOND_MAX_DISTANCE: f64 = 2.0;

/// Calculates backbone dihedral angles (phi and psi) in degrees for all residues in a protein chain.
///
/// - Residues missing any of the essential backbone atoms (`N`, `CA`, `C`) are safely skipped.
/// - Adjacent residues are considered connected only if sequential in order and the C(i-1)–N(i)
///   distance does not exceed `PEPTIDE_BOND_MAX_DISTANCE` (2.0 Å).
/// - For the first residue in a contiguous segment, `phi` is `None` (no preceding `C` atom).
/// - For the last residue in a contiguous segment, `psi` is `None` (no following `N` atom).
pub fn calc_phi_psi(chain: &AtomGroup) -> Vec<RamachandranAngle> {
    let mut res_keys = chain.get_group_list();
    crate::brd::sort_nicely(&mut res_keys);
    let mut results = Vec::new();

    struct BackboneInfo {
        orig_idx: usize,
        key: String,
        name: String,
        n_pos: Position,
        ca_pos: Position,
        c_pos: Position,
    }

    let mut valid_residues: Vec<BackboneInfo> = Vec::new();
    for (orig_idx, key) in res_keys.iter().enumerate() {
        if let Some(res) = chain.get_group(key) {
            let n = res.get_atom("N");
            let ca = res.get_atom("CA");
            let c = res.get_atom("C");

            if let (Some(n), Some(ca), Some(c)) = (n, ca, c) {
                valid_residues.push(BackboneInfo {
                    orig_idx,
                    key: key.clone(),
                    name: res.name.clone(),
                    n_pos: n.xyz,
                    ca_pos: ca.xyz,
                    c_pos: c.xyz,
                });
            }
        }
    }

    for i in 0..valid_residues.len() {
        let curr = &valid_residues[i];

        // phi is defined if the immediately preceding residue in the chain sequence exists
        // and the C(i-1)–N(i) distance does not exceed PEPTIDE_BOND_MAX_DISTANCE.
        let phi = if i > 0
            && valid_residues[i - 1].orig_idx == curr.orig_idx - 1
            && valid_residues[i - 1].c_pos.distance_from(&curr.n_pos) <= PEPTIDE_BOND_MAX_DISTANCE
        {
            Some(dihedral_angle(
                &valid_residues[i - 1].c_pos,
                &curr.n_pos,
                &curr.ca_pos,
                &curr.c_pos,
            ))
        } else {
            None
        };

        // psi is defined if the immediately following residue in the chain sequence exists
        // and the C(i)–N(i+1) distance does not exceed PEPTIDE_BOND_MAX_DISTANCE.
        let psi = if i + 1 < valid_residues.len()
            && valid_residues[i + 1].orig_idx == curr.orig_idx + 1
            && curr.c_pos.distance_from(&valid_residues[i + 1].n_pos) <= PEPTIDE_BOND_MAX_DISTANCE
        {
            Some(dihedral_angle(
                &curr.n_pos,
                &curr.ca_pos,
                &curr.c_pos,
                &valid_residues[i + 1].n_pos,
            ))
        } else {
            None
        };

        results.push(RamachandranAngle::new(&curr.key, &curr.name, phi, psi));
    }

    results
}
