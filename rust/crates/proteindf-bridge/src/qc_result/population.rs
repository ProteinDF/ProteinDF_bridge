// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

//! Calculation of orbital-resolved Mulliken atomic populations.
//!
//! Mulliken population analysis partitions total electron density across atoms.
//! For molecular orbital (MO) $i$, the gross population on atom $A$ is given by:
//!
//! $$q_A^{(i)} = \sum_{\mu \in A} C_{\mu i} (\mathbf{S} \mathbf{C})_{\mu i}$$
//!
//! where:
//! - $C_{\mu i}$ is the molecular orbital coefficient for atomic orbital $\mu$ in MO $i$.
//! - $\mathbf{S}$ is the atomic orbital overlap matrix.
//! - The summation runs over all basis functions $\mu$ belonging to atom $A$.
//!
//! ### Properties
//! - Sum over all atoms in any given orbital $i$ is exactly 1.0 (by MO orthonormality $\mathbf{C}_i^T \mathbf{S} \mathbf{C}_i = 1$).
//! - Individual atomic contributions can be positive or negative.
//! - For dummy atoms (such as 'X' without basis functions), the contribution is 0.0.
//! - The total gross Mulliken population on atom $A$ across all occupied orbitals is:
//!   $$Q_A = \sum_i \mathrm{occ}_i \cdot q_A^{(i)}$$
//!   and the net Mulliken atomic charge is $Z_A - Q_A$, where $Z_A$ is the charge column
//!   of `molecule/xyzc` (`QcAtom::charge`): the nuclear charge for real atoms, and the
//!   assigned point charge for dummy atoms (e.g. 1.0 for `X` in ProteinDF_test `CH3COOX_HF`).

use crate::error::{BridgeError, Result};
use crate::matrix::SymmetricMatrix;
use crate::qc_result::types::{MoCoefficients, OrbitalInfo};

/// Orbital-resolved Mulliken atomic contribution for a specific molecular orbital.
#[derive(Debug, Clone, PartialEq)]
pub struct OrbitalContribution {
    /// Zero-based molecular orbital index $i$.
    pub orbital_index: usize,
    /// List of atomic contributions `(atom_index, contribution)` in the order of
    /// atoms in the molecule, length $N_{\text{atom}}$.
    pub atom_contributions: Vec<(usize, f64)>,
}

/// Calculates orbital-resolved Mulliken atomic contributions for all molecular orbitals.
///
/// # Arguments
/// * `mo_coeffs` - MO coefficient matrix $\mathbf{C}$ ($N_{\text{AO}} \times N_{\text{MO}}$).
/// * `s_matrix` - Overlap matrix $\mathbf{S}$ ($N_{\text{AO}} \times N_{\text{AO}}$).
/// * `orb_info` - Mapping between atomic orbitals and atoms.
/// * `num_atoms` - Total number of atoms in the molecule ($N_{\text{atom}}$).
///
/// # Returns
/// A vector of [`OrbitalContribution`] of length $N_{\text{MO}}$, one for each molecular orbital.
///
/// # Errors
/// Returns [`BridgeError::ValueError`] if matrix or mapping dimensions do not match:
/// - `s_matrix.dim() != mo_coeffs.n_ao()`
/// - `orb_info.num_of_aos != mo_coeffs.n_ao()`
/// - Any AO's `atom_index >= num_atoms`
pub fn calc_orbital_mulliken_contributions(
    mo_coeffs: &MoCoefficients,
    s_matrix: &SymmetricMatrix,
    orb_info: &OrbitalInfo,
    num_atoms: usize,
) -> Result<Vec<OrbitalContribution>> {
    let n_ao = mo_coeffs.n_ao();
    let n_mo = mo_coeffs.n_mo();

    if s_matrix.dim() != n_ao {
        return Err(BridgeError::value_error(
            format!("s_matrix dim {}", s_matrix.dim()),
            format!("Overlap matrix dimension must match MO AO count {n_ao}"),
        ));
    }

    if orb_info.num_of_aos != n_ao {
        return Err(BridgeError::value_error(
            format!("orb_info AOs {}", orb_info.num_of_aos),
            format!("OrbitalInfo AO count must match MO AO count {n_ao}"),
        ));
    }

    // Verify atom indices in orb_info
    for ao in &orb_info.orbitals {
        if ao.atom_index >= num_atoms {
            return Err(BridgeError::value_error(
                format!("atom_index {}", ao.atom_index),
                format!("AO atom_index must be less than num_atoms {num_atoms}"),
            ));
        }
    }

    // Precompute (SC)_{mu, i} = sum_{nu} S_{mu, nu} * C_{nu, i}
    // Shape: n_ao x n_mo
    let mut sc = vec![0.0; n_ao * n_mo];
    for mu in 0..n_ao {
        for i in 0..n_mo {
            let mut sum = 0.0;
            for nu in 0..n_ao {
                // S is symmetric; get(mu, nu) handles symmetry internally.
                let s_val = s_matrix.get(mu, nu)?;
                let c_val = mo_coeffs.matrix.get(nu, i)?;
                sum += s_val * c_val;
            }
            sc[mu * n_mo + i] = sum;
        }
    }

    let mut result = Vec::with_capacity(n_mo);

    for i in 0..n_mo {
        let mut atom_contrib_map = vec![0.0; num_atoms];

        // For each AO mu, add C_{mu, i} * (SC)_{mu, i} to its parent atom
        for ao in &orb_info.orbitals {
            let mu = ao.ao_index;
            let c_mu_i = mo_coeffs.matrix.get(mu, i)?;
            let sc_mu_i = sc[mu * n_mo + i];
            let term = c_mu_i * sc_mu_i;
            atom_contrib_map[ao.atom_index] += term;
        }

        let atom_contributions: Vec<(usize, f64)> =
            atom_contrib_map.into_iter().enumerate().collect();

        result.push(OrbitalContribution {
            orbital_index: i,
            atom_contributions,
        });
    }

    Ok(result)
}
