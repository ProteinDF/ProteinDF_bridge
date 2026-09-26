// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

//! Protein backbone amide and N-terminal hydrogen addition.
//!
//! Provides functions to construct missing backbone hydrogens for protein residues
//! based on polymer connectivity:
//! - **Internal peptide backbone amide hydrogen (`H`)**: Constructed from the preceding
//!   residue's carbonyl carbon (`C`), the current residue's amide nitrogen (`N`), and
//!   alpha carbon (`CA`) assuming standard planar trans-peptide geometry ($N-H = 1.01\ \text{Å}$).
//!   Delegates to [`crate::hydrogen_bond::calc_pseudo_hydrogen`].
//!   Proline residues (`PRO`) lack backbone amide hydrogens and are skipped.
//! - **N-terminal ammonium hydrogens (`H1`, `H2`, `H3`)**: For true N-terminal residues
//!   where no preceding residue exists, standard tetrahedral ammonium geometry is generated
//!   using [`crate::Modeling::get_NH3`] aligned against the $N-CA$ bond.
//!
//! See `RUST_PORT_SPEC.md` §3.16 for full design details and literature sources.

use crate::atom::Atom;
use crate::atom_group::AtomGroup;
use crate::error::{BridgeError, Result};
use crate::hydrogen_bond::calc_pseudo_hydrogen;
use crate::hydrogenation::HydrogenationReport;
use crate::modeling::Modeling;
use crate::position::Position;

/// Standard bond length for backbone amide N-H in Angstroms.
///
/// Source: Engh, R. A. & Huber, R. (1991). Accurate bond and angle parameters for X-ray
/// protein structure refinement. Acta Cryst. A47, 392-400. Also matches the DSSP /
/// Kabsch & Sander (1983) electrostatic model used across structural bioinformatics.
pub const STANDARD_AMIDE_NH_BOND_LENGTH: f64 = 1.01;

/// Half cone angle from CA->N axis for tetrahedral N-terminal ammonium hydrogens (rad).
///
/// In [`Modeling::get_NH3`], the `angle` parameter specifies the cone angle from the +Z axis,
/// which is rotated to align with the $CA \to N$ vector. The angle between the generated N-H bonds
/// and the $N \to CA$ bond (back towards CA) is therefore $\pi - \text{angle}$.
/// To obtain the ideal tetrahedral angle of $\arccos(-1/3) \approx 109.47^\circ$, the cone angle
/// must satisfy $\cos(\text{angle}) = -\cos(109.47^\circ) = 1/3$, yielding $\arccos(1/3) \approx 70.53^\circ$.
pub const NH3_TETRAHEDRAL_HALF_ANGLE: f64 = 1.230_959_417_340_774_7; // (1.0_f64 / 3.0).acos()

/// Calculates the position of the peptide backbone amide hydrogen (`H`) for a residue
/// given the preceding residue's carbonyl carbon `C`, the current residue's `N`, and `CA`.
///
/// Assumes planar trans-peptide bond geometry where the $N-H$ bond bisects the angle
/// formed by the $C_{prev}-N$ vector and the $CA-N$ vector in the $C_{prev}-N-CA$ plane.
/// Delegates directly to [`calc_pseudo_hydrogen`].
///
/// Returns `None` if the input coordinates are degenerate (e.g. $N == C_{prev}$ or $N == CA$)
/// or collinear.
pub fn build_backbone_amide_hydrogen(
    c_prev: &Position,
    n_curr: &Position,
    ca_curr: &Position,
) -> Option<Position> {
    calc_pseudo_hydrogen(c_prev, n_curr, ca_curr)
}

/// Builds N-terminal hydrogens (`H1`, `H2`, `H3`) for an N-terminal residue
/// using [`Modeling::get_NH3`] aligned along the $CA \to N$ vector.
///
/// For proline (`PRO`), whose ring nitrogen is a secondary amine at the N-terminus,
/// two hydrogens (`H1`, `H2`) are generated.
pub fn build_nterm_hydrogens(residue: &AtomGroup) -> Result<Vec<(String, Atom)>> {
    let n = residue
        .get_atom("N")
        .ok_or_else(|| BridgeError::input_error("residue", "missing 'N' atom for N-terminus"))?;
    let ca = residue
        .get_atom("CA")
        .ok_or_else(|| BridgeError::input_error("residue", "missing 'CA' atom for N-terminus"))?;

    let modeling = Modeling::new()?;

    // Standard tetrahedral NH3 geometry from Modeling::get_NH3:
    // angle = arccos(1/3) (~70.53 deg) so that the resulting H-N-CA angle is tetrahedral (109.47 deg).
    // length = 1.0 A
    let nh3 = modeling.get_NH3(NH3_TETRAHEDRAL_HALF_ANGLE, 1.0)?;

    // Align NH3 local frame: In get_NH3, N is at (0,0,0) and the 3 hydrogens point generally towards +Z.
    // In the residue, the ammonium group points away from CA (vector CA -> N).
    let target_dir = n.xyz - ca.xyz;
    let ref_dir = Position::new(0.0, 0.0, 1.0);

    let rot = modeling.arbitary_rotate_matrix(target_dir, ref_dir)?;

    let mut result = Vec::new();
    let is_pro = residue.name == "PRO";

    let h_names: &[&str] = if is_pro {
        &["H1", "H2"]
    } else {
        &["H1", "H2", "H3"]
    };

    for &h_name in h_names {
        if let Some(src_h) = nh3.get_atom(h_name) {
            let mut pos = src_h.xyz;
            pos.rotate(&rot)?;
            pos += n.xyz;

            let mut atom = Atom::new_with_pos("H", pos)?;
            atom.name = h_name.to_string();
            result.push((h_name.to_string(), atom));
        }
    }

    Ok(result)
}

/// Adds missing backbone hydrogens to a single protein residue.
///
/// # Arguments
/// * `residue` - The residue [`AtomGroup`] to modify.
/// * `prev_residue` - Optional reference to the preceding residue in the chain.
///   - If `Some(prev)`: Treats `residue` as an internal/C-terminal residue connected to `prev`.
///     Adds peptide backbone amide hydrogen (`H`) using [`build_backbone_amide_hydrogen`].
///     If the residue is Proline (`PRO`), no backbone hydrogen is added.
///   - If `None`: Treats `residue` as an N-terminal residue. Adds N-terminal ammonium
///     hydrogens (`H1`, `H2`, `H3`, or `H1`, `H2` for PRO) using [`build_nterm_hydrogens`].
///
/// # Errors
/// Returns an error if required heavy atoms (`N`, `CA`, or `prev.C`) are missing,
/// have degenerate coordinates, or if the $C_{prev}-N$ distance is outside plausible peptide bond bounds (0.8 - 2.5 Å).
pub fn add_backbone_hydrogens_to_residue_in_place(
    residue: &mut AtomGroup,
    prev_residue: Option<&AtomGroup>,
) -> Result<HydrogenationReport> {
    let mut staged = Vec::new();

    if let Some(prev) = prev_residue {
        // Internal peptide residue
        // Proline has a tertiary amine ring nitrogen in peptide chains and carries NO amide hydrogen.
        if residue.name == "PRO" {
            return Ok(HydrogenationReport {
                added_hydrogens: 0,
                added_atom_names: Vec::new(),
            });
        }

        // Only add 'H' if not already present
        if !residue.has_atom("H") {
            let n = residue.get_atom("N").ok_or_else(|| {
                BridgeError::input_error("residue", "missing 'N' atom for backbone amide hydrogen")
            })?;
            let ca = residue.get_atom("CA").ok_or_else(|| {
                BridgeError::input_error("residue", "missing 'CA' atom for backbone amide hydrogen")
            })?;
            let c_prev = prev.get_atom("C").ok_or_else(|| {
                BridgeError::input_error(
                    "prev_residue",
                    "missing 'C' atom on preceding residue for backbone amide hydrogen",
                )
            })?;

            // Guard against disconnected residues or chain breaks
            let c_n_dist = (n.xyz - c_prev.xyz).length();
            if !(0.8..=2.5).contains(&c_n_dist) {
                return Err(BridgeError::value_error(
                    "coordinates",
                    format!(
                        "peptide bond C_prev-N distance {c_n_dist:.3} A is out of plausible range (0.8 - 2.5 A)"
                    ),
                ));
            }

            let pos =
                build_backbone_amide_hydrogen(&c_prev.xyz, &n.xyz, &ca.xyz).ok_or_else(|| {
                    BridgeError::value_error(
                        "coordinates",
                        "degenerate backbone geometry (C_prev, N, CA are collinear or coincident)",
                    )
                })?;

            let mut atom = Atom::new_with_pos("H", pos)?;
            atom.name = "H".to_string();
            staged.push(("H".to_string(), atom));
        }
    } else {
        // True N-terminal residue
        // For N-terminal PRO (secondary amine), purge any spurious existing H3
        if residue.name == "PRO" && residue.has_atom("H3") {
            residue.remove_atom("H3");
        }

        let nterm_hydrogens = build_nterm_hydrogens(residue)?;
        for (h_name, atom) in nterm_hydrogens {
            if !residue.has_atom(&h_name) {
                staged.push((h_name, atom));
            }
        }
    }

    let added_atom_names: Vec<String> = staged.iter().map(|(name, _)| name.clone()).collect();
    let added_hydrogens = staged.len();

    // Apply atomically
    for (name, atom) in staged {
        residue.set_atom(&name, atom);
    }

    Ok(HydrogenationReport {
        added_hydrogens,
        added_atom_names,
    })
}

/// Out-of-place variant of [`add_backbone_hydrogens_to_residue_in_place`].
pub fn add_backbone_hydrogens_to_residue(
    residue: &AtomGroup,
    prev_residue: Option<&AtomGroup>,
) -> Result<AtomGroup> {
    let mut copy = residue.clone();
    add_backbone_hydrogens_to_residue_in_place(&mut copy, prev_residue)?;
    Ok(copy)
}
