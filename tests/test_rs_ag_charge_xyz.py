#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test suite validating rs.AtomGroup charge, nuclei_charge, real_nuclei_charge,
get_atom_kinds, and get_xyz methods.
"""

import glob
import os
import unittest

try:
    import pytest
    approx = pytest.approx
except ImportError:
    import math

    class approx:
        def __init__(self, expected, rel=1e-6, abs=1e-12):
            self.expected = expected
            self.rel = rel
            self.abs = abs

        def __eq__(self, actual):
            return math.isclose(actual, self.expected, rel_tol=self.rel, abs_tol=self.abs)

        def __repr__(self):
            return f"approx({self.expected})"

import proteindf_bridge.rs as rs


class TestRsAgChargeXyz(unittest.TestCase):
    """Test suite for rs.AtomGroup charge, atomic species, and XYZ generation."""

    def setUp(self):
        self.data_dir = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "proteindf_bridge", "data")
        )

    def test_brd_files_charge_and_kinds(self):
        """Criterion 1: Verify charge, nuclei_charge, real_nuclei_charge, and get_atom_kinds()

        against expected values for all .brd files in proteindf_bridge/data/.
        """
        brd_files = sorted(glob.glob(os.path.join(self.data_dir, "*.brd")))
        self.assertGreater(len(brd_files), 0, "No .brd files found in data directory")

        expected_brd = {
            "ACE_ALA_NME.brd": (0.0, 78.0, 78.0, {"C", "H", "N", "O"}),
            "ACE_ALA_NME_cis1.brd": (0.0, 78.0, 78.0, {"C", "H", "N", "O"}),
            "ACE_ALA_NME_cis2.brd": (0.0, 78.0, 78.0, {"C", "H", "N", "O"}),
            "ACE_ALA_NME_trans1.brd": (0.0, 78.0, 78.0, {"C", "H", "N", "O"}),
            "ACE_ALA_NME_trans2.brd": (0.0, 78.0, 78.0, {"C", "H", "N", "O"}),
            "NML.brd": (0.0, 40.0, 40.0, {"C", "H", "N", "O"}),
            "NML_trans.brd": (0.0, 40.0, 40.0, {"C", "H", "N", "O"}),
        }

        for brd_path in brd_files:
            fname = os.path.basename(brd_path)
            with self.subTest(file=fname):
                rs_ag = rs.load_atomgroup(brd_path)
                exp_charge, exp_nuclei, exp_real_nuclei, exp_kinds = expected_brd[fname]

                self.assertEqual(rs_ag.charge, approx(exp_charge))
                self.assertEqual(rs_ag.nuclei_charge, approx(exp_nuclei))
                self.assertEqual(
                    rs_ag.real_nuclei_charge, approx(exp_real_nuclei)
                )
                self.assertEqual(rs_ag.get_atom_kinds(), exp_kinds)

    def test_2mgo_pdb_charge_and_kinds(self):
        """Criterion 1: Verify charge, nuclei_charge, real_nuclei_charge, and get_atom_kinds()

        against expected values for 2MGO.pdb.
        """
        pdb_path = os.path.join(self.data_dir, "2MGO.pdb")
        self.assertTrue(os.path.exists(pdb_path), "2MGO.pdb fixture not found")

        rs_ag = rs.Pdb(pdb_path).get_atomgroup()

        self.assertEqual(rs_ag.charge, approx(0.0))
        self.assertEqual(rs_ag.nuclei_charge, approx(10720.0))
        self.assertEqual(rs_ag.real_nuclei_charge, approx(10720.0))
        self.assertEqual(rs_ag.get_atom_kinds(), {"C", "H", "N", "O", "S"})

    def test_dummy_atoms_tree(self):
        """Criterion 2: Verify dummy atoms (atomic_number 0, with charge) in small tree."""
        # Rust tree
        rs_root = rs.AtomGroup("root")
        rs_sub = rs.AtomGroup("sub")
        rs_a1 = rs.Atom(symbol="C")
        rs_a1.charge = -0.5
        rs_d1 = rs.Atom(symbol="X")
        rs_d1.charge = 1.5
        rs_sub.set_atom("C1", rs_a1)
        rs_sub.set_atom("X1", rs_d1)
        rs_root.set_group("sub", rs_sub)

        rs_a2 = rs.Atom(symbol="H")
        rs_a2.charge = 0.2
        rs_d2 = rs.Atom(symbol="X")
        rs_d2.charge = -0.8
        rs_root.set_atom("H1", rs_a2)
        rs_root.set_atom("X2", rs_d2)

        self.assertIn("X", rs_root.get_atom_kinds())

        # Verify expected exact values
        # charge: -0.5 + 1.5 + 0.2 + (-0.8) = 0.4
        self.assertEqual(rs_root.charge, approx(0.4))
        # nuclei_charge: 6.0 + 1.5 (dummy charge) + 1.0 + (-0.8 dummy charge) = 7.7
        self.assertEqual(rs_root.nuclei_charge, approx(7.7))
        # real_nuclei_charge: 6.0 + 1.0 = 7.0
        self.assertEqual(rs_root.real_nuclei_charge, approx(7.0))
        # atom kinds: C, H, X
        self.assertEqual(rs_root.get_atom_kinds(), {"C", "H", "X"})

    def test_get_xyz_flat(self):
        """Criterion 3: get_xyz() for flat molecule returns formatted XYZ string."""
        rs_ag = rs.AtomGroup("flat")
        rs_a1 = rs.Atom(symbol="C", xyz=[1.0, 2.0, 3.0])
        rs_a2 = rs.Atom(symbol="H", xyz=[0.5, -1.2, 3.456])
        rs_ag.set_atom("C1", rs_a1)
        rs_ag.set_atom("H1", rs_a2)

        expected_xyz = (
            "2\n"
            "# \n"
            f"{'C':>2} {1.0:8.3f} {2.0:8.3f} {3.0:8.3f}\n"
            f"{'H':>2} {0.5:8.3f} {-1.2:8.3f} {3.456:8.3f}\n"
        )
        self.assertEqual(rs_ag.get_xyz(), expected_xyz)

    def test_get_xyz_nested(self):
        """Criterion 3: get_xyz() for nested groups outputs header exactly once

        with child group atoms preceding parent atoms.
        """
        rs_root = rs.AtomGroup("root")
        rs_sub1 = rs.AtomGroup("sub1")
        rs_sub2 = rs.AtomGroup("sub2")
        rs_a1 = rs.Atom(symbol="C", xyz=[1.0, 2.0, 3.0])
        rs_a2 = rs.Atom(symbol="O", xyz=[4.0, 5.0, 6.0])
        rs_a3 = rs.Atom(symbol="N", xyz=[7.0, 8.0, 9.0])

        rs_sub1.set_atom("C1", rs_a1)
        rs_sub2.set_atom("O1", rs_a2)
        rs_root.set_group("g1", rs_sub1)
        rs_root.set_group("g2", rs_sub2)
        rs_root.set_atom("N1", rs_a3)

        xyz_str = rs_root.get_xyz()
        lines = xyz_str.splitlines()

        # Total 3 atoms: 1 atom count line + 1 comment line + 3 atom lines = 5 lines
        self.assertEqual(len(lines), 5)
        self.assertEqual(lines[0], "3")
        self.assertEqual(lines[1], "# ")

        # Children atoms come first, then own atoms
        self.assertTrue(lines[2].startswith(" C "))
        self.assertTrue(lines[3].startswith(" O "))
        self.assertTrue(lines[4].startswith(" N "))

        # Verify formatting "%2s %8.3f %8.3f %8.3f"
        self.assertEqual(lines[2], f"{'C':>2} {1.0:8.3f} {2.0:8.3f} {3.0:8.3f}")
        self.assertEqual(lines[3], f"{'O':>2} {4.0:8.3f} {5.0:8.3f} {6.0:8.3f}")
        self.assertEqual(lines[4], f"{'N':>2} {7.0:8.3f} {8.0:8.3f} {9.0:8.3f}")


if __name__ == "__main__":
    unittest.main()
