// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom_group::AtomGroup;
use crate::brd::sort_nicely;
use crate::position::Position;

/// Represents a detected backbone hydrogen bond from a donor residue (N-H)
/// to an acceptor residue (C=O).
#[derive(Debug, Clone, PartialEq)]
pub struct HydrogenBond {
    pub donor_residue_key: String,
    pub acceptor_residue_key: String,
    pub energy: f64,
}

impl HydrogenBond {
    /// Creates a new `HydrogenBond` record.
    pub fn new(
        donor_residue_key: impl Into<String>,
        acceptor_residue_key: impl Into<String>,
        energy: f64,
    ) -> Self {
        Self {
            donor_residue_key: donor_residue_key.into(),
            acceptor_residue_key: acceptor_residue_key.into(),
            energy,
        }
    }
}

/// Calculates the pseudo amide hydrogen coordinate for residue `i` given
/// the C atom of residue `i-1`, the N atom of residue `i`, and the CA atom of residue `i`.
///
/// Formula:
/// ```text
/// vec_cn  = normalize(N(i) - C(i-1))
/// vec_can = normalize(N(i) - CA(i))
/// vec_nh  = normalize(vec_cn + vec_can)
/// H(i)    = N(i) + 1.01 * vec_nh
/// ```
pub fn calc_pseudo_hydrogen(
    c_prev: &Position,
    n_curr: &Position,
    ca_curr: &Position,
) -> Option<Position> {
    let diff_cn = *n_curr - *c_prev;
    let diff_can = *n_curr - *ca_curr;

    let len_cn = diff_cn.length();
    let len_can = diff_can.length();
    if len_cn < 1.0e-12 || len_can < 1.0e-12 {
        return None;
    }

    let vec_cn = diff_cn / len_cn;
    let vec_can = diff_can / len_can;

    let sum_vec = vec_cn + vec_can;
    let len_sum = sum_vec.length();
    if len_sum < 1.0e-12 {
        return None;
    }

    let vec_nh = sum_vec / len_sum;
    Some(*n_curr + vec_nh * 1.01)
}

/// Calculates the Kabsch-Sander electrostatic interaction energy (in kcal/mol)
/// between an amide donor group (N, H) and a carbonyl acceptor group (C, O).
///
/// Formula:
/// ```text
/// E = q1 * q2 * (1/r(O_a, N_d) + 1/r(C_a, H_d) - 1/r(O_a, H_d) - 1/r(C_a, N_d)) * 332.0
/// ```
/// where `q1 * q2 = 0.084` (q1 = 0.42, q2 = 0.20, factor 332.0).
pub fn calc_kabsch_sander_energy(
    n_d: &Position,
    h_d: &Position,
    c_a: &Position,
    o_a: &Position,
) -> f64 {
    let r_on = o_a.distance_from(n_d).max(1.0e-6);
    let r_ch = c_a.distance_from(h_d).max(1.0e-6);
    let r_oh = o_a.distance_from(h_d).max(1.0e-6);
    let r_cn = c_a.distance_from(n_d).max(1.0e-6);

    0.084 * (1.0 / r_on + 1.0 / r_ch - 1.0 / r_oh - 1.0 / r_cn) * 332.0
}

/// Detects backbone hydrogen bonds within a protein chain using the Kabsch-Sander electrostatic model.
///
/// - Residues missing any of `N`, `CA`, `C`, `O` are safely skipped.
/// - Residue keys are sorted using `sort_nicely` to guarantee sequence order independent of insertion order.
/// - A hydrogen bond is detected if:
///   1. `E(donor, acceptor) < -0.5` kcal/mol
///   2. `|donor_index - acceptor_index| > 2` (excluding trivial local interactions `|d - a| <= 2`)
///   3. The donor residue has an immediately preceding residue in the chain to compute pseudo-H.
pub fn calc_backbone_hbonds(chain: &AtomGroup) -> Vec<HydrogenBond> {
    let mut res_keys = chain.get_group_list();
    sort_nicely(&mut res_keys);

    struct ResidueBackbone {
        orig_idx: usize,
        key: String,
        n: Position,
        ca: Position,
        c: Position,
        o: Position,
    }

    let mut residues: Vec<ResidueBackbone> = Vec::new();
    for (orig_idx, key) in res_keys.iter().enumerate() {
        if let Some(res) = chain.get_group(key) {
            let n = res.get_atom("N");
            let ca = res.get_atom("CA");
            let c = res.get_atom("C");
            let o = res.get_atom("O");

            if let (Some(n), Some(ca), Some(c), Some(o)) = (n, ca, c, o) {
                residues.push(ResidueBackbone {
                    orig_idx,
                    key: key.clone(),
                    n: n.xyz,
                    ca: ca.xyz,
                    c: c.xyz,
                    o: o.xyz,
                });
            }
        }
    }

    // Compute pseudo-H for each residue if its immediate predecessor exists in the chain
    let mut pseudo_h_list: Vec<Option<Position>> = Vec::with_capacity(residues.len());
    for i in 0..residues.len() {
        let h = if i > 0 && residues[i - 1].orig_idx == residues[i].orig_idx - 1 {
            calc_pseudo_hydrogen(&residues[i - 1].c, &residues[i].n, &residues[i].ca)
        } else {
            None
        };
        pseudo_h_list.push(h);
    }

    let mut hbonds = Vec::new();

    // Iterate over all donor residues
    for (d_idx, d_res) in residues.iter().enumerate() {
        let h_d = match &pseudo_h_list[d_idx] {
            Some(h) => h,
            None => continue, // Cannot act as donor if pseudo-H cannot be computed
        };

        // Iterate over all acceptor residues
        for (a_idx, a_res) in residues.iter().enumerate() {
            // Condition 2: |d - a| > 2 (exclude trivial local interactions |d - a| <= 2)
            if (d_idx as isize - a_idx as isize).abs() <= 2 {
                continue;
            }

            let energy = calc_kabsch_sander_energy(&d_res.n, h_d, &a_res.c, &a_res.o);

            // Condition 1: E < -0.5 kcal/mol
            if energy < -0.5 {
                hbonds.push(HydrogenBond::new(&d_res.key, &a_res.key, energy));
            }
        }
    }

    hbonds
}
