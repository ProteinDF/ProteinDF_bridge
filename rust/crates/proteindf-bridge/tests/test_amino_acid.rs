// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use proteindf_bridge::amino_acid::AminoAcid;
use proteindf_bridge::atom_group::AtomGroup;

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
