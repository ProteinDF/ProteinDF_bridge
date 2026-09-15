// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom_group::AtomGroup;

/// Amino acid identifier utility corresponding to `proteindf_bridge.aminoacid.AminoAcid`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct AminoAcid;

impl AminoAcid {
    /// List of recognizable amino acid 3-letter codes.
    pub const AA_LIST: &'static [&'static str] = &[
        "ALA", "ASX", "ASN", "ASP", "CYS", "CYX", "GLU", "PHE", "GLY", "HIS", "HIE", "HIP", "ILE",
        "LYS", "LEU", "MET", "PRO", "GLN", "ARG", "SER", "THR", "SEC", "VAL", "TRP", "XAA", "TYR",
        "GLX",
    ];

    /// Returns `true` if the given `AtomGroup`'s name represents an amino acid.
    pub fn is_aminoacid(atomgroup: &AtomGroup) -> bool {
        Self::is_aminoacid_name(&atomgroup.name)
    }

    /// Returns `true` if the given residue name represents an amino acid.
    pub fn is_aminoacid_name(name: &str) -> bool {
        let trimmed = name.trim();
        Self::AA_LIST.contains(&trimmed)
    }

    /// Returns the static slice of supported amino acid codes.
    pub fn aa_list() -> &'static [&'static str] {
        Self::AA_LIST
    }
}
