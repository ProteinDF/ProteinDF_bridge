#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Validation test suite for Phase 7 (PR#21):
Verifying Ramachandran calculation and dihedral angle PyO3 bindings (proteindf_bridge_rs).
"""

import math
import os
import unittest

from proteindf_bridge.biopdb import Pdb as PyPdb
import proteindf_bridge_rs as rs_br

DATA_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "proteindf_bridge", "data")
)


class TestRsPhase7(unittest.TestCase):

    # --------------------------------------------------------------------------
    # 1. Dihedral Angle Geometric Sanity
    # --------------------------------------------------------------------------
    def test_dihedral_angle_function_and_method(self):
        p1 = rs_br.Position(0.0, 1.0, 0.0)
        p2 = rs_br.Position(0.0, 0.0, 0.0)
        p3 = rs_br.Position(1.0, 0.0, 0.0)

        # Cis (0 degrees)
        p4_cis = rs_br.Position(1.0, 1.0, 0.0)
        angle_cis_func = rs_br.dihedral_angle(p1, p2, p3, p4_cis)
        angle_cis_method = p1.dihedral(p2, p3, p4_cis)
        self.assertAlmostEqual(angle_cis_func, 0.0, places=7)
        self.assertAlmostEqual(angle_cis_method, 0.0, places=7)

        # Trans (180 degrees)
        p4_trans = rs_br.Position(1.0, -1.0, 0.0)
        angle_trans_func = rs_br.dihedral_angle(p1, p2, p3, p4_trans)
        angle_trans_method = p1.dihedral(p2, p3, p4_trans)
        self.assertAlmostEqual(abs(angle_trans_func), 180.0, places=7)
        self.assertAlmostEqual(abs(angle_trans_method), 180.0, places=7)

        # +90 degrees (rotated towards -Z in right-handed system looking along +X)
        p4_pos90 = rs_br.Position(1.0, 0.0, -1.0)
        angle_pos90 = rs_br.dihedral_angle(p1, p2, p3, p4_pos90)
        self.assertAlmostEqual(angle_pos90, 90.0, places=7)

        # -90 degrees (rotated towards +Z in right-handed system looking along +X)
        p4_neg90 = rs_br.Position(1.0, 0.0, 1.0)
        angle_neg90 = rs_br.dihedral_angle(p1, p2, p3, p4_neg90)
        self.assertAlmostEqual(angle_neg90, -90.0, places=7)

        # Flexible argument extraction: accept lists and strings
        angle_list = rs_br.dihedral_angle(
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
        )
        self.assertAlmostEqual(angle_list, 0.0, places=7)

        # Degenerate collinear case: returns 0.0 without crashing
        angle_collinear = rs_br.dihedral_angle(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        )
        self.assertAlmostEqual(angle_collinear, 0.0, places=7)

    # --------------------------------------------------------------------------
    # 2. Real PDB Fixture Validation (1hls.pdb model_1 / chain A)
    # --------------------------------------------------------------------------
    def test_ramachandran_1hls_ground_truth(self):
        pdb_path = os.path.join(DATA_DIR, "1hls.pdb")
        pdb = PyPdb(pdb_path)
        py_ag = pdb.get_atomgroup()

        # Load via Rust Pdb as well to get PyAtomGroup
        rs_pdb = rs_br.Pdb(pdb_path)
        rs_ag = rs_pdb.get_atomgroup()

        chain_a = rs_ag.get_group("model_1").get_group("A")
        angles = rs_br.calc_phi_psi(chain_a)

        self.assertEqual(len(angles), 21)

        # First residue (GLY 1): phi is None, psi is not None
        first = angles[0]
        self.assertEqual(first.residue_key, "1")
        self.assertEqual(first.residue_name, "GLY")
        self.assertIsNone(first.phi)
        self.assertIsNotNone(first.psi)

        # Last residue (ASN 21): phi is not None, psi is None
        last = angles[20]
        self.assertEqual(last.residue_key, "21")
        self.assertEqual(last.residue_name, "ASN")
        self.assertIsNotNone(last.phi)
        self.assertIsNone(last.psi)

        # Verification against independently computed ground truth in docs/rust-port-handoff.md:
        # Residue 4 (GLU): phi=70.5993, psi=3.3945
        res4 = next(a for a in angles if a.residue_key == "4")
        self.assertEqual(res4.residue_name, "GLU")
        self.assertAlmostEqual(res4.phi, 70.5993, places=3)
        self.assertAlmostEqual(res4.psi, 3.3945, places=3)

        # Residue 5 (GLN): phi=121.6311, psi=26.5618
        res5 = next(a for a in angles if a.residue_key == "5")
        self.assertEqual(res5.residue_name, "GLN")
        self.assertAlmostEqual(res5.phi, 121.6311, places=3)
        self.assertAlmostEqual(res5.psi, 26.5618, places=3)

        # Residue 10 (ILE): phi=84.5507, psi=-99.6314
        res10 = next(a for a in angles if a.residue_key == "10")
        self.assertEqual(res10.residue_name, "ILE")
        self.assertAlmostEqual(res10.phi, 84.5507, places=3)
        self.assertAlmostEqual(res10.psi, -99.6314, places=3)

        # __repr__ test
        self.assertIn("RamachandranAngle", repr(res4))
        self.assertIn("residue_key='4'", repr(res4))

    # --------------------------------------------------------------------------
    # 3. Deliberately Unsorted Chain Keys
    # --------------------------------------------------------------------------
    def test_ramachandran_unsorted_chain_keys(self):
        def make_residue(key, name, offset):
            res = rs_br.AtomGroup(key)
            res.name = name
            res.set_atom("N", rs_br.Atom(symbol="N", name="N", position=[offset, 0.0, 0.0]))
            res.set_atom("CA", rs_br.Atom(symbol="C", name="CA", position=[offset + 1.0, 0.0, 0.0]))
            res.set_atom("C", rs_br.Atom(symbol="C", name="C", position=[offset + 1.5, 1.0, 0.0]))
            return res

        # Sorted chain
        chain_sorted = rs_br.AtomGroup("A")
        chain_sorted.set_group("1", make_residue("1", "ALA", 0.0))
        chain_sorted.set_group("2", make_residue("2", "GLY", 3.0))
        chain_sorted.set_group("3", make_residue("3", "VAL", 6.0))
        chain_sorted.set_group("4", make_residue("4", "LEU", 9.0))
        chain_sorted.set_group("10", make_residue("10", "ILE", 12.0))

        # Scrambled chain
        chain_scrambled = rs_br.AtomGroup("A")
        chain_scrambled.set_group("3", make_residue("3", "VAL", 6.0))
        chain_scrambled.set_group("10", make_residue("10", "ILE", 12.0))
        chain_scrambled.set_group("1", make_residue("1", "ALA", 0.0))
        chain_scrambled.set_group("4", make_residue("4", "LEU", 9.0))
        chain_scrambled.set_group("2", make_residue("2", "GLY", 3.0))

        angles_sorted = rs_br.calc_phi_psi(chain_sorted)
        angles_scrambled = rs_br.calc_phi_psi(chain_scrambled)

        self.assertEqual(len(angles_sorted), 5)
        self.assertEqual(len(angles_scrambled), 5)

        for s, sc in zip(angles_sorted, angles_scrambled):
            self.assertEqual(s.residue_key, sc.residue_key)
            self.assertEqual(s.residue_name, sc.residue_name)
            if s.phi is not None and sc.phi is not None:
                self.assertAlmostEqual(s.phi, sc.phi, places=7)
            else:
                self.assertEqual(s.phi, sc.phi)
            if s.psi is not None and sc.psi is not None:
                self.assertAlmostEqual(s.psi, sc.psi, places=7)
            else:
                self.assertEqual(s.psi, sc.psi)

    # --------------------------------------------------------------------------
    # 4. Safe Skipping of Incomplete Backbone Residues
    # --------------------------------------------------------------------------
    def test_ramachandran_missing_atoms_safe_skip(self):
        chain = rs_br.AtomGroup("A")

        # Res 1: complete
        res1 = rs_br.AtomGroup("1")
        res1.name = "ALA"
        res1.set_atom("N", rs_br.Atom(symbol="N", name="N", position=[0.0, 0.0, 0.0]))
        res1.set_atom("CA", rs_br.Atom(symbol="C", name="CA", position=[1.0, 0.0, 0.0]))
        res1.set_atom("C", rs_br.Atom(symbol="C", name="C", position=[1.5, 1.0, 0.0]))
        chain.set_group("1", res1)

        # Res 2: missing C
        res2 = rs_br.AtomGroup("2")
        res2.name = "GLY"
        res2.set_atom("N", rs_br.Atom(symbol="N", name="N", position=[2.5, 1.0, 0.0]))
        res2.set_atom("CA", rs_br.Atom(symbol="C", name="CA", position=[3.0, 2.0, 0.0]))
        chain.set_group("2", res2)

        # Res 3: complete
        res3 = rs_br.AtomGroup("3")
        res3.name = "VAL"
        res3.set_atom("N", rs_br.Atom(symbol="N", name="N", position=[4.5, 1.0, 0.0]))
        res3.set_atom("CA", rs_br.Atom(symbol="C", name="CA", position=[5.5, 1.0, 0.0]))
        res3.set_atom("C", rs_br.Atom(symbol="C", name="C", position=[6.0, 2.0, 0.0]))
        chain.set_group("3", res3)

        angles = rs_br.calc_phi_psi(chain)
        # Res 2 was skipped
        self.assertEqual(len(angles), 2)
        self.assertEqual(angles[0].residue_key, "1")
        self.assertEqual(angles[1].residue_key, "3")

        # Res 1: phi is None; psi is None because res 2 lacks C
        self.assertIsNone(angles[0].phi)
        self.assertIsNone(angles[0].psi)

        # Res 3: phi is None because res 2 is broken; psi is None because it is the last residue
        self.assertIsNone(angles[1].phi)
        self.assertIsNone(angles[1].psi)


if __name__ == "__main__":
    unittest.main()
