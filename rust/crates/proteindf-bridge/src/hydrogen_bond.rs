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
        for a_res in &residues {
            // Condition 2: |d - a| > 2 (exclude trivial local interactions |d - a| <= 2)
            if (d_res.orig_idx as isize - a_res.orig_idx as isize).abs() <= 2 {
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

/// Definition of a sidechain donor/acceptor atom type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SidechainAtomType {
    pub residue_name: &'static str,
    pub atom_name: &'static str,
    pub is_donor: bool,
    pub is_acceptor: bool,
}

/// Static table of standard amino acid sidechain donor and acceptor atoms.
pub const SIDECHAIN_ATOM_TYPES: &[SidechainAtomType] = &[
    SidechainAtomType {
        residue_name: "SER",
        atom_name: "OG",
        is_donor: true,
        is_acceptor: true,
    },
    SidechainAtomType {
        residue_name: "THR",
        atom_name: "OG1",
        is_donor: true,
        is_acceptor: true,
    },
    SidechainAtomType {
        residue_name: "TYR",
        atom_name: "OH",
        is_donor: true,
        is_acceptor: true,
    },
    SidechainAtomType {
        residue_name: "ASN",
        atom_name: "ND2",
        is_donor: true,
        is_acceptor: false,
    },
    SidechainAtomType {
        residue_name: "ASN",
        atom_name: "OD1",
        is_donor: false,
        is_acceptor: true,
    },
    SidechainAtomType {
        residue_name: "GLN",
        atom_name: "NE2",
        is_donor: true,
        is_acceptor: false,
    },
    SidechainAtomType {
        residue_name: "GLN",
        atom_name: "OE1",
        is_donor: false,
        is_acceptor: true,
    },
    SidechainAtomType {
        residue_name: "HIS",
        atom_name: "ND1",
        is_donor: true,
        is_acceptor: true,
    },
    SidechainAtomType {
        residue_name: "HIS",
        atom_name: "NE2",
        is_donor: true,
        is_acceptor: true,
    },
];

/// Represents a detected sidechain hydrogen bond between a donor and acceptor atom.
#[derive(Debug, Clone, PartialEq)]
pub struct SidechainHydrogenBond {
    pub donor_path: String,
    pub donor_atom: String,
    pub acceptor_path: String,
    pub acceptor_atom: String,
    pub distance: f64,
    pub angle: Option<f64>,
}

impl SidechainHydrogenBond {
    pub fn new(
        donor_path: impl Into<String>,
        donor_atom: impl Into<String>,
        acceptor_path: impl Into<String>,
        acceptor_atom: impl Into<String>,
        distance: f64,
        angle: Option<f64>,
    ) -> Self {
        Self {
            donor_path: donor_path.into(),
            donor_atom: donor_atom.into(),
            acceptor_path: acceptor_path.into(),
            acceptor_atom: acceptor_atom.into(),
            distance,
            angle,
        }
    }
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

fn is_hydrogen_atom(symbol: &str, name: &str) -> bool {
    symbol.eq_ignore_ascii_case("H")
        || name.starts_with('H')
        || name.starts_with("1H")
        || name.starts_with("2H")
        || name.starts_with("3H")
}

/// Detects sidechain hydrogen bonds in the structure (within a chain or across chains)
/// using default heavy-atom-only mode (distances < 3.5 Å, angles not evaluated).
pub fn calc_sidechain_hbonds(root: &AtomGroup) -> Vec<SidechainHydrogenBond> {
    calc_sidechain_hbonds_with_options(root, false)
}

/// Detects sidechain hydrogen bonds in the structure with optional explicit hydrogen angle evaluation.
///
/// Rules:
/// - Candidate atoms include sidechain donors/acceptors from `SIDECHAIN_ATOM_TYPES`
///   plus backbone N (donor) and backbone O (acceptor).
/// - Backbone-backbone interactions are excluded (handled by `calc_backbone_hbonds`).
/// - Pairs within the same residue (`donor_path == acceptor_path`) are excluded.
/// - Distance requirement: `distance < 3.5 Å`.
/// - Explicit hydrogen mode (`explicit_h_mode = true`): if the donor atom has a bonded
///   hydrogen (distance < 1.3 Å) in the residue, the angle `D-H...A` must be > 120°.
///   If no bonded hydrogen exists or `explicit_h_mode = false`, `angle` is `None`.
pub fn calc_sidechain_hbonds_with_options(
    root: &AtomGroup,
    explicit_h_mode: bool,
) -> Vec<SidechainHydrogenBond> {
    let mut residues = Vec::new();
    collect_leaf_residues(root, &mut residues);

    struct DonorCandidate<'a> {
        res: &'a AtomGroup,
        atom_name: &'static str,
        pos: Position,
        is_sidechain: bool,
        bonded_hydrogens: Vec<Position>,
    }

    struct AcceptorCandidate<'a> {
        res: &'a AtomGroup,
        atom_name: &'static str,
        pos: Position,
        is_sidechain: bool,
    }

    let mut donors: Vec<DonorCandidate> = Vec::new();
    let mut acceptors: Vec<AcceptorCandidate> = Vec::new();

    for res in &residues {
        let rname = res.name.trim();

        // Collect all hydrogen positions in this residue
        let mut h_positions = Vec::new();
        for (_, atom) in res.atoms() {
            if is_hydrogen_atom(atom.symbol().unwrap_or(""), &atom.name) {
                h_positions.push(atom.xyz);
            }
        }

        // 1. Check backbone N (donor)
        if let Some(n_atom) = res.get_atom("N") {
            let n_pos = n_atom.xyz;
            let bonded_h: Vec<Position> = h_positions
                .iter()
                .copied()
                .filter(|h_pos| h_pos.distance_from(&n_pos) < 1.3)
                .collect();
            donors.push(DonorCandidate {
                res,
                atom_name: "N",
                pos: n_pos,
                is_sidechain: false,
                bonded_hydrogens: bonded_h,
            });
        }

        // 2. Check backbone O (acceptor)
        if let Some(o_atom) = res.get_atom("O") {
            acceptors.push(AcceptorCandidate {
                res,
                atom_name: "O",
                pos: o_atom.xyz,
                is_sidechain: false,
            });
        }

        // 3. Check sidechain donors and acceptors
        for entry in SIDECHAIN_ATOM_TYPES {
            if entry.residue_name == rname {
                if let Some(atom) = res.get_atom(entry.atom_name) {
                    let pos = atom.xyz;
                    if entry.is_donor {
                        let bonded_h: Vec<Position> = h_positions
                            .iter()
                            .copied()
                            .filter(|h_pos| h_pos.distance_from(&pos) < 1.3)
                            .collect();
                        donors.push(DonorCandidate {
                            res,
                            atom_name: entry.atom_name,
                            pos,
                            is_sidechain: true,
                            bonded_hydrogens: bonded_h,
                        });
                    }
                    if entry.is_acceptor {
                        acceptors.push(AcceptorCandidate {
                            res,
                            atom_name: entry.atom_name,
                            pos,
                            is_sidechain: true,
                        });
                    }
                }
            }
        }
    }

    let mut hbonds = Vec::new();

    for d in &donors {
        for a in &acceptors {
            // Exclude same residue
            if d.res.path() == a.res.path() {
                continue;
            }

            // Exclude backbone-backbone pairs
            if !d.is_sidechain && !a.is_sidechain {
                continue;
            }

            let distance = d.pos.distance_from(&a.pos);
            if distance >= 3.5 {
                continue;
            }

            // Explicit hydrogen check
            let mut angle = None;
            if explicit_h_mode && !d.bonded_hydrogens.is_empty() {
                let mut satisfied = false;
                let mut best_angle = 0.0;
                for h_pos in &d.bonded_hydrogens {
                    let v_hd = d.pos - *h_pos;
                    let v_ha = a.pos - *h_pos;
                    let len_hd = v_hd.length();
                    let len_ha = v_ha.length();
                    if len_hd > 1e-6 && len_ha > 1e-6 {
                        let cos_theta = (v_hd.dot(&v_ha) / (len_hd * len_ha)).clamp(-1.0, 1.0);
                        let deg = cos_theta.acos().to_degrees();
                        if deg > 120.0 {
                            satisfied = true;
                            if deg > best_angle {
                                best_angle = deg;
                            }
                        }
                    }
                }

                if !satisfied {
                    continue;
                }
                angle = Some(best_angle);
            }

            hbonds.push(SidechainHydrogenBond::new(
                d.res.path().to_string(),
                d.atom_name,
                a.res.path().to_string(),
                a.atom_name,
                distance,
                angle,
            ));
        }
    }

    hbonds
}
