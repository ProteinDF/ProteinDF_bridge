// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::path::PathBuf;

use proteindf_bridge::format::SimpleMol2;

fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

#[test]
fn test_load_sample_mol2_fixture() {
    let fixture_path = test_data_dir().join("sample.mol2");
    let mol2 = SimpleMol2::from_file(&fixture_path).expect("failed to load sample.mol2");
    let ag = mol2.get_atomgroup();

    assert_eq!(ag.name, "ethanol");
    assert_eq!(ag.get_number_of_all_atoms(), 9);

    // Verify atoms and properties
    let c1 = ag.get_atom("C1").expect("C1 not found");
    assert_eq!(c1.symbol().unwrap(), "C");
    assert!((c1.xyz.x - 0.0).abs() < 1e-4);
    assert!((c1.charge - (-0.10)).abs() < 1e-4);

    let c2 = ag.get_atom("C2").expect("C2 not found");
    assert_eq!(c2.symbol().unwrap(), "C");
    assert!((c2.xyz.x - 1.52).abs() < 1e-4);
    assert!((c2.charge - 0.10).abs() < 1e-4);

    let o1 = ag.get_atom("O1").expect("O1 not found");
    assert_eq!(o1.symbol().unwrap(), "O");
    assert!((o1.xyz.x - 2.00).abs() < 1e-4);
    assert!((o1.xyz.y - 1.30).abs() < 1e-4);
    assert!((o1.charge - (-0.60)).abs() < 1e-4);

    let h6 = ag.get_atom("H6").expect("H6 not found");
    assert_eq!(h6.symbol().unwrap(), "H");
    assert!((h6.xyz.x - 2.95).abs() < 1e-4);
    assert!((h6.charge - 0.40).abs() < 1e-4);

    // Verify bonds
    let mut ag_mut = ag.clone();
    let bonds = ag_mut.get_bond_list();
    assert_eq!(bonds.len(), 8);

    // All bonds in ethanol are single bonds (order 1)
    for b in &bonds {
        assert_eq!(b.order, 1);
    }

    // Verify connectivity: resolve bonds and check endpoints
    let has_c1_c2 = bonds.iter().any(|b| {
        let (a1, a2) = ag.resolve_bond(b).expect("bond endpoints must resolve");
        (a1.name == "C1" && a2.name == "C2") || (a1.name == "C2" && a2.name == "C1")
    });
    assert!(has_c1_c2, "C1-C2 bond must exist");

    let has_c2_o1 = bonds.iter().any(|b| {
        let (a1, a2) = ag.resolve_bond(b).expect("bond endpoints must resolve");
        (a1.name == "C2" && a2.name == "O1") || (a1.name == "O1" && a2.name == "C2")
    });
    assert!(has_c2_o1, "C2-O1 bond must exist");

    let has_o1_h6 = bonds.iter().any(|b| {
        let (a1, a2) = ag.resolve_bond(b).expect("bond endpoints must resolve");
        (a1.name == "O1" && a2.name == "H6") || (a1.name == "H6" && a2.name == "O1")
    });
    assert!(has_o1_h6, "O1-H6 bond must exist");
}
