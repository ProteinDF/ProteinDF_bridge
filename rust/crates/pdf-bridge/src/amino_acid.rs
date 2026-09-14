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
