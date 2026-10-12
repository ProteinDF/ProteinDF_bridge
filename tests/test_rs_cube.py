#!/usr/bin/env python
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: The ProteinDF development team
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Unit tests for Gaussian CUBE format parser via Python bindings.
"""

import os
import tempfile
import unittest

import proteindf_bridge.rs as rs_br

BOHR_TO_ANGSTROM = 0.529177210903


class TestRsCube(unittest.TestCase):
    def test_cube_basic_bohr(self):
        cube_str = """\
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
"""
        cube = rs_br.Cube()
        cube.parse_str(cube_str)

        # Header properties
        self.assertEqual(cube.comment1, "Test Cube File Title")
        self.assertEqual(cube.comment2, "Total SCF Density")
        self.assertEqual(cube.origin, (0.0, 0.5, 1.0))
        self.assertEqual(cube.unit, "bohr")
        self.assertEqual(cube.n_grid, (2, 3, 4))

        sp = cube.spacing
        self.assertEqual(sp[0], (0.5, 0.0, 0.0))
        self.assertEqual(sp[1], (0.0, 0.4, 0.0))
        self.assertEqual(sp[2], (0.0, 0.0, 0.3))

        # Atoms
        self.assertEqual(cube.num_of_atoms, 2)
        self.assertEqual(cube.num_atoms, 2)

        # Molecular orbital attributes
        self.assertFalse(cube.is_molecular_orbital)
        self.assertIsNone(cube.orbitals)
        self.assertEqual(cube.num_datasets, 1)

        # Grid positions
        pos_000 = cube.grid_point_position(0, 0, 0)
        self.assertAlmostEqual(pos_000[0], 0.0)
        self.assertAlmostEqual(pos_000[1], 0.5)
        self.assertAlmostEqual(pos_000[2], 1.0)

        pos_123 = cube.grid_point_position(1, 2, 3)
        self.assertAlmostEqual(pos_123[0], 0.0 + 1.0 * 0.5)
        self.assertAlmostEqual(pos_123[1], 0.5 + 2.0 * 0.4)
        self.assertAlmostEqual(pos_123[2], 1.0 + 3.0 * 0.3)

        # Out of bounds grid position
        with self.assertRaises(Exception):
            cube.grid_point_position(2, 0, 0)

        # Values
        self.assertAlmostEqual(cube.value(0, 0, 0), 0.01)
        self.assertAlmostEqual(cube.value(0, 0, 1), 0.02)
        self.assertAlmostEqual(cube.value(0, 1, 0), 0.05)
        self.assertAlmostEqual(cube.value(1, 2, 3), 0.24)

        data = cube.data()
        self.assertEqual(len(data), 24)
        self.assertAlmostEqual(data[0], 0.01)
        self.assertAlmostEqual(data[23], 0.24)

        # AtomGroup extraction (Bohr to Angstrom conversion)
        ag = cube.get_atomgroup()
        self.assertEqual(ag.get_number_of_atoms(), 2)

        a0 = ag.get_atom("0")
        self.assertEqual(a0.symbol, "C")
        self.assertAlmostEqual(a0.xyz.x, 0.0)
        self.assertAlmostEqual(a0.xyz.y, 0.0)
        self.assertAlmostEqual(a0.xyz.z, 0.0)
        self.assertAlmostEqual(a0.charge, 6.0)

        a1 = ag.get_atom("1")
        self.assertEqual(a1.symbol, "O")
        self.assertAlmostEqual(a1.xyz.x, 1.2 * BOHR_TO_ANGSTROM, places=6)
        self.assertAlmostEqual(a1.xyz.y, 0.0)
        self.assertAlmostEqual(a1.xyz.z, 0.0)
        self.assertAlmostEqual(a1.charge, 8.0)

        # Test alias get_atom_group
        ag_alias = cube.get_atom_group()
        self.assertEqual(ag_alias.get_number_of_atoms(), 2)

    def test_cube_molecular_orbitals(self):
        cube_str = """\
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
"""
        cube = rs_br.Cube()
        cube.parse_str(cube_str)

        self.assertTrue(cube.is_molecular_orbital)
        self.assertEqual(cube.num_of_atoms, 2)
        self.assertEqual(cube.orbitals, [10, 11])
        self.assertEqual(cube.num_datasets, 2)

        mo1 = cube.data(0)
        mo2 = cube.data(1)
        self.assertEqual(len(mo1), 8)
        self.assertEqual(len(mo2), 8)

        self.assertAlmostEqual(mo1[0], 1.0)
        self.assertAlmostEqual(mo2[0], 2.0)
        self.assertAlmostEqual(mo1[1], 1.1)
        self.assertAlmostEqual(mo2[1], 2.1)
        self.assertAlmostEqual(mo1[7], 1.7)
        self.assertAlmostEqual(mo2[7], 2.7)

        self.assertAlmostEqual(cube.value(0, 0, 0, 0), 1.0)
        self.assertAlmostEqual(cube.value(0, 0, 0, 1), 2.0)
        self.assertAlmostEqual(cube.value(1, 1, 1, 0), 1.7)
        self.assertAlmostEqual(cube.value(1, 1, 1, 1), 2.7)

    def test_cube_angstrom_units(self):
        cube_str = """\
Angstrom Cube File
Generated by tool
    1    1.000000    2.000000    3.000000
   -2    0.500000    0.000000    0.000000
   -2    0.000000    0.500000    0.000000
   -2    0.000000    0.000000    0.500000
    6    6.000000    1.500000    2.500000    3.500000
  1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0
"""
        cube = rs_br.Cube()
        cube.parse_str(cube_str)

        self.assertEqual(cube.unit, "angstrom")
        self.assertEqual(cube.n_grid, (2, 2, 2))

        ag = cube.get_atomgroup()
        atom = ag.get_atom("0")
        # In Angstrom mode, coordinates must be untouched
        self.assertAlmostEqual(atom.xyz.x, 1.5)
        self.assertAlmostEqual(atom.xyz.y, 2.5)
        self.assertAlmostEqual(atom.xyz.z, 3.5)

    def test_cube_format_errors(self):
        # 1. Missing header lines
        cube = rs_br.Cube()
        with self.assertRaises((rs_br.BrInputError, rs_br.BrError)) as ctx:
            cube.parse_str("only one line")
        self.assertIn("line 2", str(ctx.exception))

        # 2. Line 3 insufficient fields
        with self.assertRaises((rs_br.BrInputError, rs_br.BrError)) as ctx:
            cube.parse_str("title\nsubtitle\n2 0.0 0.0\n")
        self.assertIn("line 3", str(ctx.exception))

        # 3. Line 3 invalid natoms
        with self.assertRaises((rs_br.BrInputError, rs_br.BrError)) as ctx:
            cube.parse_str("title\nsubtitle\nnot_an_int 0.0 0.0 0.0\n")
        self.assertIn("line 3", str(ctx.exception))

        # 4. Line 4 grid axis count 0
        with self.assertRaises((rs_br.BrInputError, rs_br.BrError)) as ctx:
            cube.parse_str("title\nsubtitle\n1 0.0 0.0 0.0\n0 0.1 0.0 0.0\n")
        self.assertIn("line 4", str(ctx.exception))

        # 5. Unexpected EOF at atom line
        with self.assertRaises((rs_br.BrInputError, rs_br.BrError)) as ctx:
            cube.parse_str("""\
title
subtitle
1 0.0 0.0 0.0
1 0.1 0.0 0.0
1 0.0 0.1 0.0
1 0.0 0.0 0.1
""")
        self.assertIn("line 7", str(ctx.exception))

        # 6. Invalid float in atom line
        with self.assertRaises((rs_br.BrInputError, rs_br.BrError)) as ctx:
            cube.parse_str("""\
title
subtitle
1 0.0 0.0 0.0
1 0.1 0.0 0.0
1 0.0 0.1 0.0
1 0.0 0.0 0.1
6 6.0 0.0 not_a_float 0.0
1.0
""")
        self.assertIn("line 7", str(ctx.exception))

        # 7. Insufficient volumetric data
        with self.assertRaises((rs_br.BrInputError, rs_br.BrError)) as ctx:
            cube.parse_str("""\
title
subtitle
1 0.0 0.0 0.0
2 0.1 0.0 0.0
2 0.0 0.1 0.0
2 0.0 0.0 0.1
6 6.0 0.0 0.0 0.0
1.0 2.0 3.0
""")
        self.assertIn("insufficient volumetric values", str(ctx.exception))

        # 8. Too many volumetric values
        with self.assertRaises((rs_br.BrInputError, rs_br.BrError)) as ctx:
            cube.parse_str("""\
title
subtitle
1 0.0 0.0 0.0
1 0.1 0.0 0.0
1 0.0 0.1 0.0
1 0.0 0.0 0.1
6 6.0 0.0 0.0 0.0
1.0 2.0
""")
        self.assertIn("line 8", str(ctx.exception))
        self.assertIn("too many volumetric values", str(ctx.exception))

        # 9. Non-numeric value in volumetric data
        with self.assertRaises((rs_br.BrInputError, rs_br.BrError)) as ctx:
            cube.parse_str("""\
title
subtitle
1 0.0 0.0 0.0
1 0.1 0.0 0.0
1 0.0 0.1 0.0
1 0.0 0.0 0.1
6 6.0 0.0 0.0 0.0
invalid_float
""")
        self.assertIn("line 8", str(ctx.exception))
        self.assertIn("invalid float value", str(ctx.exception))

    def test_cube_file_io(self):
        cube_str = """\
File IO Test
Subtitle
    1    0.0 0.0 0.0
    1    1.0 0.0 0.0
    1    0.0 1.0 0.0
    1    0.0 0.0 1.0
    1    1.0 0.0 0.0 0.0
  42.0
"""
        with tempfile.NamedTemporaryFile(suffix=".cube", delete=False, mode="w") as f:
            f.write(cube_str)
            temp_path = f.name

        try:
            cube = rs_br.Cube(temp_path)
            self.assertEqual(cube.comment1, "File IO Test")
            self.assertEqual(cube.comment2, "Subtitle")
            self.assertAlmostEqual(cube.value(0, 0, 0), 42.0)

            ag = cube.get_atomgroup()
            self.assertEqual(ag.get_number_of_atoms(), 1)
            self.assertEqual(ag.get_atom("0").symbol, "H")
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)


if __name__ == "__main__":
    unittest.main()
