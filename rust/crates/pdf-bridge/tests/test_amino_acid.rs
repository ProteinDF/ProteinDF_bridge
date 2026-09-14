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

use pdf_bridge::amino_acid::AminoAcid;
use pdf_bridge::atom_group::AtomGroup;

#[test]
fn test_is_aminoacid() {
    let mut ag_ala = AtomGroup::new();
    ag_ala.name = "ALA".to_string();
    assert!(AminoAcid::is_aminoacid(&ag_ala));

    let mut ag_other = AtomGroup::new();
    ag_other.name = "HOH".to_string();
    assert!(!AminoAcid::is_aminoacid(&ag_other));
}

#[test]
fn test_all_amino_acids() {
    for &code in AminoAcid::aa_list() {
        let mut ag = AtomGroup::new();
        ag.name = code.to_string();
        assert!(
            AminoAcid::is_aminoacid(&ag),
            "expected {code} to be recognized as an amino acid"
        );
        assert!(AminoAcid::is_aminoacid_name(code));
    }
    assert_eq!(AminoAcid::aa_list().len(), 27);
}

#[test]
fn test_non_amino_acids() {
    let non_aa = ["HOH", "WAT", "LIG", "ATP", "NA", "CL", "", "123"];
    for &code in &non_aa {
        let mut ag = AtomGroup::new();
        ag.name = code.to_string();
        assert!(
            !AminoAcid::is_aminoacid(&ag),
            "expected {code} to not be recognized as an amino acid"
        );
        assert!(!AminoAcid::is_aminoacid_name(code));
    }
}

#[test]
fn test_trimmed_name() {
    assert!(AminoAcid::is_aminoacid_name("  ALA  "));
    assert!(AminoAcid::is_aminoacid_name("GLY\n"));
    assert!(!AminoAcid::is_aminoacid_name("ala")); // Case-sensitive uppercase
}
