// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use proteindf_bridge::error::BridgeError;
use proteindf_bridge::format::cube::{Cube, CubeUnit, BOHR_TO_ANGSTROM};
use std::fs;

#[test]
fn test_cube_basic_bohr() {
    let cube_str = "\
Test Cube File Title
Total SCF Density
    2    0.000000    0.500000    1.000000
    2    0.500000    0.000000    0.000000
    3    0.000000    0.400000    0.000000
    4    0.000000    0.000000    0.300000
    6    6.000000    0.000000    0.000000    0.000000
    8    8.000000    1.200000    0.000000    0.000000
  0.01  0.02  0.03  0.04  0.05  0.06
  0.07  0.08  0.09  0.10  0.11  0.12
  0.13  0.14  0.15  0.16  0.17  0.18
  0.19  0.20  0.21  0.22  0.23  0.24
";

    let cube = Cube::parse_str(cube_str).unwrap();

    // Headers
    assert_eq!(cube.comment1(), "Test Cube File Title");
    assert_eq!(cube.comment2(), "Total SCF Density");
    assert_eq!(cube.origin(), [0.0, 0.5, 1.0]);
    assert_eq!(cube.unit(), CubeUnit::Bohr);
    assert_eq!(cube.n_grid(), [2, 3, 4]);
    assert_eq!(cube.spacing()[0], [0.5, 0.0, 0.0]);
    assert_eq!(cube.spacing()[1], [0.0, 0.4, 0.0]);
    assert_eq!(cube.spacing()[2], [0.0, 0.0, 0.3]);

    // Atoms
    assert_eq!(cube.num_atoms(), 2);
    let atoms = cube.atoms();
    assert_eq!(atoms[0].atomic_number, 6);
    assert_eq!(atoms[0].charge, 6.0);
    assert_eq!(atoms[0].x, 0.0);
    assert_eq!(atoms[0].y, 0.0);
    assert_eq!(atoms[0].z, 0.0);

    assert_eq!(atoms[1].atomic_number, 8);
    assert_eq!(atoms[1].charge, 8.0);
    assert_eq!(atoms[1].x, 1.2);
    assert_eq!(atoms[1].y, 0.0);
    assert_eq!(atoms[1].z, 0.0);

    // Molecular orbital check
    assert!(!cube.is_molecular_orbital());
    assert_eq!(cube.orbitals(), None);
    assert_eq!(cube.num_datasets(), 1);

    // Grid positions
    let pos_000 = cube.grid_point_position(0, 0, 0).unwrap();
    assert_eq!(pos_000, [0.0, 0.5, 1.0]);

    let pos_123 = cube.grid_point_position(1, 2, 3).unwrap();
    assert!((pos_123[0] - (0.0 + 1.0 * 0.5)).abs() < 1e-9);
    assert!((pos_123[1] - (0.5 + 2.0 * 0.4)).abs() < 1e-9);
    assert!((pos_123[2] - (1.0 + 3.0 * 0.3)).abs() < 1e-9);

    assert!(cube.grid_point_position(2, 0, 0).is_err());
    assert!(cube.grid_point_position(0, 3, 0).is_err());
    assert!(cube.grid_point_position(0, 0, 4).is_err());

    // Volumetric values
    // Flat index: (ix * 3 + iy) * 4 + iz
    // (0, 0, 0) -> idx 0 -> 0.01
    assert!((cube.value(0, 0, 0).unwrap() - 0.01).abs() < 1e-9);
    // (0, 0, 1) -> idx 1 -> 0.02
    assert!((cube.value(0, 0, 1).unwrap() - 0.02).abs() < 1e-9);
    // (0, 1, 0) -> idx 4 -> 0.05
    assert!((cube.value(0, 1, 0).unwrap() - 0.05).abs() < 1e-9);
    // (1, 2, 3) -> idx 23 -> 0.24
    assert!((cube.value(1, 2, 3).unwrap() - 0.24).abs() < 1e-9);

    let raw_data = cube.data(0).unwrap();
    assert_eq!(raw_data.len(), 24);
    assert!((raw_data[0] - 0.01).abs() < 1e-9);
    assert!((raw_data[23] - 0.24).abs() < 1e-9);

    // AtomGroup extraction (Bohr -> Angstrom conversion)
    let ag = cube.get_atomgroup().unwrap();
    assert_eq!(ag.get_number_of_atoms(), 2);
    let a0 = ag.get_atom("0").unwrap();
    assert_eq!(a0.symbol().unwrap(), "C");
    assert_eq!(a0.xyz.x, 0.0);
    assert_eq!(a0.xyz.y, 0.0);
    assert_eq!(a0.xyz.z, 0.0);

    let a1 = ag.get_atom("1").unwrap();
    assert_eq!(a1.symbol().unwrap(), "O");
    assert!((a1.xyz.x - (1.2 * BOHR_TO_ANGSTROM)).abs() < 1e-9);
    assert_eq!(a1.xyz.y, 0.0);
    assert_eq!(a1.xyz.z, 0.0);
    assert_eq!(a1.charge, 8.0);
}

#[test]
fn test_cube_molecular_orbitals() {
    // 2x2x2 grid = 8 grid points. 2 MOs -> 8 * 2 = 16 values.
    // NATOMS = -2 (negative indicates MO output)
    let cube_str = "\
MO Cube File
Alpha MOs
   -2    0.000000    0.000000    0.000000
    2    1.000000    0.000000    0.000000
    2    0.000000    1.000000    0.000000
    2    0.000000    0.000000    1.000000
    1    1.000000    0.000000    0.000000   -0.700000
    1    1.000000    0.000000    0.000000    0.700000
    2   10   11
  1.0  2.0    1.1  2.1
  1.2  2.2    1.3  2.3
  1.4  2.4    1.5  2.5
  1.6  2.6    1.7  2.7
";

    let cube = Cube::parse_str(cube_str).unwrap();

    assert!(cube.is_molecular_orbital());
    assert_eq!(cube.num_atoms(), 2);
    assert_eq!(cube.orbitals(), Some([10, 11].as_slice()));
    assert_eq!(cube.num_datasets(), 2);

    let mo1 = cube.data(0).unwrap();
    let mo2 = cube.data(1).unwrap();
    assert_eq!(mo1.len(), 8);
    assert_eq!(mo2.len(), 8);

    assert!((mo1[0] - 1.0).abs() < 1e-9);
    assert!((mo2[0] - 2.0).abs() < 1e-9);

    assert!((mo1[1] - 1.1).abs() < 1e-9);
    assert!((mo2[1] - 2.1).abs() < 1e-9);

    assert!((mo1[7] - 1.7).abs() < 1e-9);
    assert!((mo2[7] - 2.7).abs() < 1e-9);

    // Test orbital_value
    assert!((cube.orbital_value(0, 0, 0, 0).unwrap() - 1.0).abs() < 1e-9);
    assert!((cube.orbital_value(1, 0, 0, 0).unwrap() - 2.0).abs() < 1e-9);
    assert!((cube.orbital_value(0, 1, 1, 1).unwrap() - 1.7).abs() < 1e-9);
    assert!((cube.orbital_value(1, 1, 1, 1).unwrap() - 2.7).abs() < 1e-9);
}

#[test]
fn test_cube_angstrom_units() {
    // Negative grid count indicates units are in Angstroms
    let cube_str = "\
Angstrom Cube File
Generated by tool
    1    1.000000    2.000000    3.000000
   -2    0.500000    0.000000    0.000000
   -2    0.000000    0.500000    0.000000
   -2    0.000000    0.000000    0.500000
    6    6.000000    1.500000    2.500000    3.500000
  1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0
";

    let cube = Cube::parse_str(cube_str).unwrap();

    assert_eq!(cube.unit(), CubeUnit::Angstrom);
    assert_eq!(cube.n_grid(), [2, 2, 2]);

    let ag = cube.get_atomgroup().unwrap();
    let atom = ag.get_atom("0").unwrap();
    // In Angstrom mode, coordinates must be untouched (no conversion)
    assert!((atom.xyz.x - 1.5).abs() < 1e-9);
    assert!((atom.xyz.y - 2.5).abs() < 1e-9);
    assert!((atom.xyz.z - 3.5).abs() < 1e-9);
}

#[test]
fn test_cube_format_errors() {
    // 1. Missing header lines
    let err1 = Cube::parse_str("only one line").unwrap_err();
    match &err1 {
        BridgeError::InputError { expr, .. } => {
            assert!(expr.contains("line 2"), "Expected line 2 in {}", expr);
        }
        _ => panic!("Expected BridgeError::InputError, got {:?}", err1),
    }

    // 2. Line 3 insufficient fields
    let err2 = Cube::parse_str("title\nsubtitle\n2 0.0 0.0\n").unwrap_err();
    match &err2 {
        BridgeError::InputError { expr, .. } => {
            assert!(expr.contains("line 3"), "Expected line 3 in {}", expr);
        }
        _ => panic!("Expected BridgeError::InputError, got {:?}", err2),
    }

    // 3. Line 3 invalid natoms
    let err3 = Cube::parse_str("title\nsubtitle\nnot_an_int 0.0 0.0 0.0\n").unwrap_err();
    match &err3 {
        BridgeError::InputError { expr, .. } => {
            assert!(expr.contains("line 3"), "Expected line 3 in {}", expr);
        }
        _ => panic!("Expected BridgeError::InputError, got {:?}", err3),
    }

    // 4. Line 4 grid axis count 0
    let err4 = Cube::parse_str("title\nsubtitle\n1 0.0 0.0 0.0\n0 0.1 0.0 0.0\n").unwrap_err();
    match &err4 {
        BridgeError::InputError { expr, .. } => {
            assert!(expr.contains("line 4"), "Expected line 4 in {}", expr);
        }
        _ => panic!("Expected BridgeError::InputError, got {:?}", err4),
    }

    // 5. Unexpected EOF at atom line
    let err5 = Cube::parse_str(
        "\
title
subtitle
1 0.0 0.0 0.0
1 0.1 0.0 0.0
1 0.0 0.1 0.0
1 0.0 0.0 0.1
",
    )
    .unwrap_err();
    match &err5 {
        BridgeError::InputError { expr, .. } => {
            assert!(expr.contains("line 7"), "Expected line 7 in {}", expr);
        }
        _ => panic!("Expected BridgeError::InputError, got {:?}", err5),
    }

    // 6. Invalid float in atom line
    let err6 = Cube::parse_str(
        "\
title
subtitle
1 0.0 0.0 0.0
1 0.1 0.0 0.0
1 0.0 0.1 0.0
1 0.0 0.0 0.1
6 6.0 0.0 not_a_float 0.0
1.0
",
    )
    .unwrap_err();
    match &err6 {
        BridgeError::InputError { expr, .. } => {
            assert!(expr.contains("line 7"), "Expected line 7 in {}", expr);
        }
        _ => panic!("Expected BridgeError::InputError, got {:?}", err6),
    }

    // 7. Insufficient volumetric data
    let err7 = Cube::parse_str(
        "\
title
subtitle
1 0.0 0.0 0.0
2 0.1 0.0 0.0
2 0.0 0.1 0.0
2 0.0 0.0 0.1
6 6.0 0.0 0.0 0.0
1.0 2.0 3.0
",
    )
    .unwrap_err();
    match &err7 {
        BridgeError::InputError { msg, .. } => {
            assert!(
                msg.contains("insufficient volumetric values"),
                "Expected error in msg: {}",
                msg
            );
        }
        _ => panic!("Expected BridgeError::InputError, got {:?}", err7),
    }

    // 8. Too many volumetric values
    let err8 = Cube::parse_str(
        "\
title
subtitle
1 0.0 0.0 0.0
1 0.1 0.0 0.0
1 0.0 0.1 0.0
1 0.0 0.0 0.1
6 6.0 0.0 0.0 0.0
1.0 2.0
",
    )
    .unwrap_err();
    match &err8 {
        BridgeError::InputError { expr, msg } => {
            assert!(expr.contains("line 8"), "Expected line 8 in {}", expr);
            assert!(
                msg.contains("too many volumetric values"),
                "Expected message: {}",
                msg
            );
        }
        _ => panic!("Expected BridgeError::InputError, got {:?}", err8),
    }

    // 9. Non-numeric value in volumetric data
    let err9 = Cube::parse_str(
        "\
title
subtitle
1 0.0 0.0 0.0
1 0.1 0.0 0.0
1 0.0 0.1 0.0
1 0.0 0.0 0.1
6 6.0 0.0 0.0 0.0
invalid_float
",
    )
    .unwrap_err();
    match &err9 {
        BridgeError::InputError { expr, msg } => {
            assert!(expr.contains("line 8"), "Expected line 8 in {}", expr);
            assert!(
                msg.contains("invalid float value"),
                "Expected message: {}",
                msg
            );
        }
        _ => panic!("Expected BridgeError::InputError, got {:?}", err9),
    }
}

#[test]
fn test_cube_file_io() {
    let cube_str = "\
File IO Test
Subtitle
    1    0.0 0.0 0.0
    1    1.0 0.0 0.0
    1    0.0 1.0 0.0
    1    0.0 0.0 1.0
    1    1.0 0.0 0.0 0.0
  42.0
";
    let temp_dir = std::env::temp_dir();
    let temp_file = temp_dir.join(format!("test_cube_{}.cube", std::process::id()));

    fs::write(&temp_file, cube_str).unwrap();

    let cube = Cube::from_file(&temp_file).unwrap();
    assert_eq!(cube.comment1(), "File IO Test");
    assert_eq!(cube.value(0, 0, 0).unwrap(), 42.0);

    let _ = fs::remove_file(&temp_file);
}
