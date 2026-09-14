// Copyright (C) 2019 The ProteinDF development team.
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

use std::path::PathBuf;

use pdf_bridge::format::SimpleMmcif;

#[test]
fn test_load_mmcif_basic() {
    let mmcif_content = r#"data_test_entry
_entry.id test_entry
_cell.length_a 10.0
_cell.length_b 20.0
_cell.length_c 30.0
"#;

    let mut cif = SimpleMmcif::new();
    cif.load_from_str(mmcif_content).unwrap();

    assert!(cif.data().contains_key("data_test_entry"));
    let block = cif.get_data_block("data_test_entry").unwrap();
    assert_eq!(
        block.key_values.get("_entry.id").map(|s| s.as_str()),
        Some("test_entry")
    );
    assert_eq!(
        block.key_values.get("_cell.length_a").map(|s| s.as_str()),
        Some("10.0")
    );
    assert_eq!(
        block.key_values.get("_cell.length_b").map(|s| s.as_str()),
        Some("20.0")
    );
    assert_eq!(
        block.key_values.get("_cell.length_c").map(|s| s.as_str()),
        Some("30.0")
    );
}

#[test]
fn test_ala_ccd_atomgroup() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/ALA.cif");
    let cif = SimpleMmcif::from_file(&path).expect("failed to load ALA.cif");

    let mol_names = cif.get_molecule_names();
    assert_eq!(mol_names, vec!["data_ALA".to_string()]);

    let ag = cif
        .get_atomgroup("data_ALA")
        .expect("failed to get ALA atomgroup");

    // Check molecule/residue name
    assert_eq!(ag.name, "ALA");

    // Check atom count: 13 atoms
    assert_eq!(ag.get_number_of_atoms(), 13);

    // Check expected coordinates for key atoms (from pdbx_model_Cartn_*_ideal)
    let n = ag.get_atom("N").expect("N atom missing");
    assert_eq!(n.name, "N");
    assert_eq!(n.symbol().unwrap(), "N");
    assert!((n.xyz.x - (-0.966)).abs() < 1e-4);
    assert!((n.xyz.y - 0.493).abs() < 1e-4);
    assert!((n.xyz.z - 1.500).abs() < 1e-4);

    let ca = ag.get_atom("CA").expect("CA atom missing");
    assert_eq!(ca.name, "CA");
    assert_eq!(ca.symbol().unwrap(), "C");
    assert!((ca.xyz.x - 0.257).abs() < 1e-4);
    assert!((ca.xyz.y - 0.418).abs() < 1e-4);
    assert!((ca.xyz.z - 0.692).abs() < 1e-4);

    let c = ag.get_atom("C").expect("C atom missing");
    assert_eq!(c.name, "C");
    assert_eq!(c.symbol().unwrap(), "C");
    assert!((c.xyz.x - (-0.094)).abs() < 1e-4);
    assert!((c.xyz.y - 0.017).abs() < 1e-4);
    assert!((c.xyz.z - (-0.716)).abs() < 1e-4);

    let o = ag.get_atom("O").expect("O atom missing");
    assert_eq!(o.name, "O");
    assert_eq!(o.symbol().unwrap(), "O");
    assert!((o.xyz.x - (-1.056)).abs() < 1e-4);
    assert!((o.xyz.y - (-0.682)).abs() < 1e-4);
    assert!((o.xyz.z - (-0.923)).abs() < 1e-4);

    let cb = ag.get_atom("CB").expect("CB atom missing");
    assert_eq!(cb.symbol().unwrap(), "C");
    assert!((cb.xyz.x - 1.204).abs() < 1e-4);
    assert!((cb.xyz.y - (-0.620)).abs() < 1e-4);
    assert!((cb.xyz.z - 1.296).abs() < 1e-4);

    let oxt = ag.get_atom("OXT").expect("OXT atom missing");
    assert_eq!(oxt.symbol().unwrap(), "O");
    assert!((oxt.xyz.x - 0.661).abs() < 1e-4);
    assert!((oxt.xyz.y - 0.439).abs() < 1e-4);
    assert!((oxt.xyz.z - (-1.742)).abs() < 1e-4);

    // Check bonds: 12 bonds total
    let bonds = ag.bonds();
    assert_eq!(bonds.len(), 12);

    // Verify bond orders: C=O is bond order 2, all other 11 bonds are order 1
    let mut double_bonds = 0;
    let mut single_bonds = 0;
    for bond in bonds {
        let p1 = bond.atom1_path.trim_start_matches('/');
        let p2 = bond.atom2_path.trim_start_matches('/');
        if bond.order == 2 {
            double_bonds += 1;
            // Should be between C and O
            assert!(
                (p1 == "C" && p2 == "O") || (p1 == "O" && p2 == "C"),
                "Unexpected double bond between {} and {}",
                p1,
                p2
            );
        } else if bond.order == 1 {
            single_bonds += 1;
        }
    }
    assert_eq!(double_bonds, 1);
    assert_eq!(single_bonds, 11);
}

#[test]
fn test_deuterium_and_coordinate_fallback() {
    let mmcif_content = r#"data_D_test
_chem_comp.id D_test
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.model_Cartn_x
_chem_comp_atom.model_Cartn_y
_chem_comp_atom.model_Cartn_z
D_test D1 D 1.5 2.5 3.5
D_test C1 C 0.0 0.0 0.0
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
D_test D1 C1 SING
"#;

    let mut cif = SimpleMmcif::new();
    cif.load_from_str(mmcif_content).unwrap();
    let ag = cif.get_atomgroup("data_D_test").unwrap();

    assert_eq!(ag.name, "D_test");
    assert_eq!(ag.get_number_of_atoms(), 2);

    let d1 = ag.get_atom("D1").unwrap();
    // Deuterium "D" should be mapped to Hydrogen "H" (symbol H, atomic number 1)
    assert_eq!(d1.symbol().unwrap(), "H");
    assert_eq!(d1.atomic_number(), 1);
    // Should fallback to model_Cartn_* coordinates (1.5, 2.5, 3.5)
    assert!((d1.xyz.x - 1.5).abs() < 1e-6);
    assert!((d1.xyz.y - 2.5).abs() < 1e-6);
    assert!((d1.xyz.z - 3.5).abs() < 1e-6);

    let bonds = ag.bonds();
    assert_eq!(bonds.len(), 1);
    assert_eq!(bonds[0].order, 1);
}

#[test]
fn test_error_handling() {
    let cif = SimpleMmcif::new();
    let result = cif.get_atomgroup("non_existent");
    assert!(result.is_err());
}
