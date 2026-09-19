// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::collections::HashMap;
use std::fmt;

use crate::atom_group::AtomGroup;
use crate::brd::sort_nicely;
use crate::hydrogen_bond::calc_backbone_hbonds;

/// 3-state secondary structure code (DSSP-style simplified classification).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SsCode {
    /// Helix (3-10, alpha, or pi helix)
    Helix,
    /// Extended strand (participating in beta-ladder / bridge)
    Strand,
    /// Loop / coil
    Loop,
}

impl SsCode {
    /// Returns the single-letter code character ('H', 'E', or '-').
    pub fn as_char(&self) -> char {
        match self {
            Self::Helix => 'H',
            Self::Strand => 'E',
            Self::Loop => '-',
        }
    }
}

impl fmt::Display for SsCode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_char())
    }
}

/// Secondary structure assignment for a residue.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SecondaryStructure {
    pub residue_key: String,
    pub residue_name: String,
    pub code: SsCode,
}

impl SecondaryStructure {
    pub fn new(
        residue_key: impl Into<String>,
        residue_name: impl Into<String>,
        code: SsCode,
    ) -> Self {
        Self {
            residue_key: residue_key.into(),
            residue_name: residue_name.into(),
            code,
        }
    }
}

/// Calculates 3-state secondary structure assignments for all residues in a chain.
///
/// Uses the Kabsch-Sander hydrogen bond network detected by `calc_backbone_hbonds`
/// and identifies:
/// - n-turns (n = 3, 4, 5) where HBond(i+n, i) is present
/// - n-helices formed by two consecutive n-turns (4-helix prioritized over 3- and 5-helix)
/// - parallel and antiparallel beta-bridges forming ladders
pub fn calc_secondary_structure(chain: &AtomGroup) -> Vec<SecondaryStructure> {
    let mut res_keys = chain.get_group_list();
    sort_nicely(&mut res_keys);
    let l = res_keys.len();

    if l == 0 {
        return Vec::new();
    }

    // Map residue key to index
    let key_to_idx: HashMap<&str, usize> = res_keys
        .iter()
        .enumerate()
        .map(|(idx, key)| (key.as_str(), idx))
        .collect();

    // Get backbone hydrogen bonds
    let hbonds = calc_backbone_hbonds(chain);

    // Build hbmap[acceptor_idx][donor_idx]:
    // true if residue d_idx acts as donor to residue a_idx
    let mut hbmap = vec![vec![false; l]; l];
    for hb in &hbonds {
        if let (Some(&d_idx), Some(&a_idx)) = (
            key_to_idx.get(hb.donor_residue_key.as_str()),
            key_to_idx.get(hb.acceptor_residue_key.as_str()),
        ) {
            hbmap[a_idx][d_idx] = true;
        }
    }

    // Identify n-turns: turn_n[i] is true if HBond(donor=i+n, acceptor=i)
    let get_turns = |n: usize| -> Vec<bool> {
        let mut turns = vec![false; l];
        for i in 0..l {
            if i + n < l && hbmap[i][i + n] {
                turns[i] = true;
            }
        }
        turns
    };

    let turn3 = get_turns(3);
    let turn4 = get_turns(4);
    let turn5 = get_turns(5);

    // Helical core segments (two consecutive n-turns):
    // h_n[i] is true if turn_n[i-1] && turn_n[i]
    let get_consecutive_turns = |turns: &[bool], n: usize| -> Vec<bool> {
        let mut h = vec![false; l];
        for i in 1..l {
            if i + n < l && turns[i - 1] && turns[i] {
                h[i] = true;
            }
        }
        h
    };

    let h3 = get_consecutive_turns(&turn3, 3);
    let h4 = get_consecutive_turns(&turn4, 4);
    let h5 = get_consecutive_turns(&turn5, 5);

    // Helix 4 has priority
    let mut helix4 = vec![false; l];
    for i in 0..l {
        if h4[i] {
            for offset in 0..4 {
                if i + offset < l {
                    helix4[i + offset] = true;
                }
            }
        }
    }

    // Mask h3 and h5 where helix4 is present or immediately following
    let mut h3_masked = vec![false; l];
    for i in 0..l {
        if h3[i] && !helix4[i] && !(i + 1 < l && helix4[i + 1]) {
            h3_masked[i] = true;
        }
    }

    let mut h5_masked = vec![false; l];
    for i in 0..l {
        if h5[i] && !helix4[i] && !(i + 1 < l && helix4[i + 1]) {
            h5_masked[i] = true;
        }
    }

    let mut helix3 = vec![false; l];
    for i in 0..l {
        if h3_masked[i] {
            for offset in 0..3 {
                if i + offset < l {
                    helix3[i + offset] = true;
                }
            }
        }
    }

    let mut helix5 = vec![false; l];
    for i in 0..l {
        if h5_masked[i] {
            for offset in 0..5 {
                if i + offset < l {
                    helix5[i + offset] = true;
                }
            }
        }
    }

    // Combined helix
    let mut is_helix = vec![false; l];
    for i in 0..l {
        if helix4[i] || helix3[i] || helix5[i] {
            is_helix[i] = true;
        }
    }

    // Identify bridges (ladders)
    let mut is_strand = vec![false; l];
    if l >= 3 {
        for i in 1..(l - 1) {
            for j in 1..(l - 1) {
                if (i as isize - j as isize).abs() <= 2 {
                    continue;
                }

                // Parallel bridge:
                // [HBond(j, i-1) && HBond(i+1, j)] || [HBond(i, j-1) && HBond(j+1, i)]
                let p1 = hbmap[i - 1][j] && hbmap[j][i + 1];
                let p2 = hbmap[j - 1][i] && hbmap[i][j + 1];

                // Antiparallel bridge:
                // [HBond(j, i) && HBond(i, j)] || [HBond(j+1, i-1) && HBond(i+1, j-1)]
                let a1 = hbmap[i][j] && hbmap[j][i];
                let a2 = hbmap[i - 1][j + 1] && hbmap[j - 1][i + 1];

                if p1 || p2 || a1 || a2 {
                    is_strand[i] = true;
                    is_strand[j] = true;
                }
            }
        }
    }

    // Final assignment
    let mut results = Vec::with_capacity(l);
    for (i, key) in res_keys.iter().enumerate() {
        let code = if is_helix[i] {
            SsCode::Helix
        } else if is_strand[i] {
            SsCode::Strand
        } else {
            SsCode::Loop
        };

        let residue_name = chain
            .get_group(key)
            .map(|r| r.name.clone())
            .unwrap_or_default();

        results.push(SecondaryStructure::new(key.clone(), residue_name, code));
    }

    results
}

/// Applies 3-state secondary structure assignments to each residue in a chain.
///
/// Calls [`calc_secondary_structure`] on the chain and writes the resulting
/// [`SsCode`] directly into each residue group's `secondary_structure` field.
pub fn apply_secondary_structure(chain: &mut AtomGroup) {
    let assignments = calc_secondary_structure(chain);
    for assignment in assignments {
        if let Some(res) = chain.get_group_mut(&assignment.residue_key) {
            res.secondary_structure = Some(assignment.code);
        }
    }
}
