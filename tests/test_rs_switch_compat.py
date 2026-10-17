#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tests for rs switch compatibility:
Ensuring rs.Position, rs.Atom, and rs.AtomGroup match pure Python behaviors.
"""

import math
import unittest
import numpy as np

from proteindf_bridge.position import Position as PyPos
from proteindf_bridge.atom import Atom as PyAtom
from proteindf_bridge.atomgroup import AtomGroup as PyAtomGroup
from proteindf_bridge.matrix import Matrix as PyMat
import proteindf_bridge.rs as rs_br


class TestRsSwitchCompat(unittest.TestCase):

    # --------------------------------------------------------------------------
    # 1. rotate with pure Python Matrix and rs.Matrix
    # --------------------------------------------------------------------------
    def _create_rot_matrices(self):
        # 90-degree rotation around z-axis
        # [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
        py_mat = PyMat([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        rs_mat = rs_br.Matrix([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        return py_mat, rs_mat

    def test_position_rotate(self):
        py_mat, rs_mat = self._create_rot_matrices()

        pos_rs_with_py_mat = rs_br.Position(1.0, 2.0, 3.0)
        pos_rs_with_py_mat.rotate(py_mat)

        pos_rs_with_rs_mat = rs_br.Position(1.0, 2.0, 3.0)
        pos_rs_with_rs_mat.rotate(rs_mat)

        self.assertAlmostEqual(pos_rs_with_py_mat.x, pos_rs_with_rs_mat.x, places=6)
        self.assertAlmostEqual(pos_rs_with_py_mat.y, pos_rs_with_rs_mat.y, places=6)
        self.assertAlmostEqual(pos_rs_with_py_mat.z, pos_rs_with_rs_mat.z, places=6)

        pos_py = PyPos([1.0, 2.0, 3.0])
        pos_py.rotate(py_mat)
        self.assertAlmostEqual(pos_rs_with_py_mat.x, pos_py.x, places=6)
        self.assertAlmostEqual(pos_rs_with_py_mat.y, pos_py.y, places=6)
        self.assertAlmostEqual(pos_rs_with_py_mat.z, pos_py.z, places=6)

        # Non-3x3 matrix should raise ValueError
        bad_mat_py = PyMat([[1.0, 0.0], [0.0, 1.0]])
        bad_mat_rs = rs_br.Matrix([[1.0, 0.0], [0.0, 1.0]])
        p = rs_br.Position(1.0, 2.0, 3.0)
        with self.assertRaises(ValueError):
            p.rotate(bad_mat_py)
        with self.assertRaises(ValueError):
            p.rotate(bad_mat_rs)

    def test_atom_rotate(self):
        py_mat, rs_mat = self._create_rot_matrices()

        atom_rs_with_py_mat = rs_br.Atom(xyz=[1.0, 2.0, 3.0])
        atom_rs_with_py_mat.rotate(py_mat)

        atom_rs_with_rs_mat = rs_br.Atom(xyz=[1.0, 2.0, 3.0])
        atom_rs_with_rs_mat.rotate(rs_mat)

        self.assertAlmostEqual(atom_rs_with_py_mat.xyz.x, atom_rs_with_rs_mat.xyz.x, places=6)
        self.assertAlmostEqual(atom_rs_with_py_mat.xyz.y, atom_rs_with_rs_mat.xyz.y, places=6)
        self.assertAlmostEqual(atom_rs_with_py_mat.xyz.z, atom_rs_with_rs_mat.xyz.z, places=6)

        atom_py = PyAtom(xyz=[1.0, 2.0, 3.0])
        atom_py.rotate(py_mat)
        self.assertAlmostEqual(atom_rs_with_py_mat.xyz.x, atom_py.xyz.x, places=6)
        self.assertAlmostEqual(atom_rs_with_py_mat.xyz.y, atom_py.xyz.y, places=6)
        self.assertAlmostEqual(atom_rs_with_py_mat.xyz.z, atom_py.xyz.z, places=6)

        bad_mat_py = PyMat([[1.0, 0.0], [0.0, 1.0]])
        a = rs_br.Atom(xyz=[1.0, 2.0, 3.0])
        with self.assertRaises(ValueError):
            a.rotate(bad_mat_py)

    def test_atomgroup_rotate(self):
        py_mat, rs_mat = self._create_rot_matrices()

        # Build rs atom group with nested groups and atoms
        ag_py_mat = rs_br.AtomGroup()
        ag_py_mat.set_atom("a1", rs_br.Atom(xyz=[1.0, 2.0, 3.0]))
        sub_py_mat = rs_br.AtomGroup()
        sub_py_mat.set_atom("a2", rs_br.Atom(xyz=[4.0, 5.0, 6.0]))
        ag_py_mat.set_group("sub", sub_py_mat)

        ag_rs_mat = rs_br.AtomGroup()
        ag_rs_mat.set_atom("a1", rs_br.Atom(xyz=[1.0, 2.0, 3.0]))
        sub_rs_mat = rs_br.AtomGroup()
        sub_rs_mat.set_atom("a2", rs_br.Atom(xyz=[4.0, 5.0, 6.0]))
        ag_rs_mat.set_group("sub", sub_rs_mat)

        ag_py_mat.rotate(py_mat)
        ag_rs_mat.rotate(rs_mat)

        self.assertAlmostEqual(ag_py_mat["a1"].xyz.x, ag_rs_mat["a1"].xyz.x, places=6)
        self.assertAlmostEqual(ag_py_mat["a1"].xyz.y, ag_rs_mat["a1"].xyz.y, places=6)
        self.assertAlmostEqual(ag_py_mat["a1"].xyz.z, ag_rs_mat["a1"].xyz.z, places=6)

        self.assertAlmostEqual(ag_py_mat["sub"]["a2"].xyz.x, ag_rs_mat["sub"]["a2"].xyz.x, places=6)
        self.assertAlmostEqual(ag_py_mat["sub"]["a2"].xyz.y, ag_rs_mat["sub"]["a2"].xyz.y, places=6)
        self.assertAlmostEqual(ag_py_mat["sub"]["a2"].xyz.z, ag_rs_mat["sub"]["a2"].xyz.z, places=6)

        # Compare with pure Python AtomGroup
        ag_pure_py = PyAtomGroup()
        ag_pure_py.set_atom("a1", PyAtom(xyz=[1.0, 2.0, 3.0]))
        sub_pure_py = PyAtomGroup()
        sub_pure_py.set_atom("a2", PyAtom(xyz=[4.0, 5.0, 6.0]))
        ag_pure_py.set_group("sub", sub_pure_py)
        ag_pure_py.rotate(py_mat)

        self.assertAlmostEqual(ag_py_mat["a1"].xyz.x, ag_pure_py["a1"].xyz.x, places=6)
        self.assertAlmostEqual(ag_py_mat["a1"].xyz.y, ag_pure_py["a1"].xyz.y, places=6)
        self.assertAlmostEqual(ag_py_mat["a1"].xyz.z, ag_pure_py["a1"].xyz.z, places=6)
        self.assertAlmostEqual(ag_py_mat["sub"]["a2"].xyz.x, ag_pure_py["sub"]["a2"].xyz.x, places=6)
        self.assertAlmostEqual(ag_py_mat["sub"]["a2"].xyz.y, ag_pure_py["sub"]["a2"].xyz.y, places=6)
        self.assertAlmostEqual(ag_py_mat["sub"]["a2"].xyz.z, ag_pure_py["sub"]["a2"].xyz.z, places=6)

        bad_mat_py = PyMat([[1.0, 0.0], [0.0, 1.0]])
        with self.assertRaises(ValueError):
            ag_py_mat.rotate(bad_mat_py)

    # --------------------------------------------------------------------------
    # 2. Position len, iter, multiplication with numbers.Real, __setstate__
    # --------------------------------------------------------------------------
    def test_position_len_iter(self):
        p_rs = rs_br.Position(1.5, 2.5, 3.5)
        p_py = PyPos([1.5, 2.5, 3.5])

        self.assertEqual(len(p_rs), len(p_py))
        self.assertEqual(len(p_rs), 3)

        self.assertEqual(list(p_rs), list(p_py))
        self.assertEqual(list(p_rs), [1.5, 2.5, 3.5])

    def test_position_multiplication(self):
        p_rs = rs_br.Position(1.0, 2.0, 3.0)
        p_py = PyPos([1.0, 2.0, 3.0])

        # p * 2
        res_rs = p_rs * 2
        res_py = p_py * 2
        self.assertEqual(list(res_rs), list(res_py))

        # 2 * p
        res_rs = 2 * p_rs
        res_py = 2 * p_py
        self.assertEqual(list(res_rs), list(res_py))

        # np.float64(2) * p
        np_val = np.float64(2.0)
        res_rs = np_val * p_rs
        res_py = np_val * p_py
        self.assertEqual(list(res_rs), list(res_py))

        # p * np.float64(2)
        res_rs = p_rs * np_val
        res_py = p_py * np_val
        self.assertEqual(list(res_rs), list(res_py))

        # int multiplication
        res_rs = p_rs * 3
        res_py = p_py * 3
        self.assertEqual(list(res_rs), list(res_py))
        res_rs = 3 * p_rs
        res_py = 3 * p_py
        self.assertEqual(list(res_rs), list(res_py))

        # rmul delegating to other.__mul__
        class CustomMultiplier:
            def __mul__(self, other):
                return "custom_result"

        custom = CustomMultiplier()
        self.assertEqual(custom * p_rs, "custom_result")

        # Unsupported type raises TypeError
        with self.assertRaises(TypeError):
            _ = "string" * p_rs
        with self.assertRaises(TypeError):
            _ = p_rs * "string"


    def test_position_setstate(self):
        p = rs_br.Position()

        # Valid list and tuple
        p.__setstate__([1.0, 2.0, 3.0])
        self.assertEqual(list(p), [1.0, 2.0, 3.0])

        p.__setstate__((4.0, 5.0, 6.0))
        self.assertEqual(list(p), [4.0, 5.0, 6.0])

        # Non-list and non-tuple should raise TypeError
        with self.assertRaises(TypeError):
            p.__setstate__("1.0 2.0 3.0")

        with self.assertRaises(TypeError):
            p.__setstate__({"x": 1.0, "y": 2.0, "z": 3.0})

        with self.assertRaises(TypeError):
            p.__setstate__(123)

    # --------------------------------------------------------------------------
    # 3. Atom.charge float conversion and ValueError
    # --------------------------------------------------------------------------
    def test_atom_charge(self):
        atom = rs_br.Atom()
        atom.charge = "0.5"
        self.assertEqual(atom.charge, 0.5)

        atom.charge = "-1.25"
        self.assertEqual(atom.charge, -1.25)

        atom.charge = 2
        self.assertEqual(atom.charge, 2.0)

        atom.charge = 0.75
        self.assertEqual(atom.charge, 0.75)

        # Invalid string must raise ValueError
        with self.assertRaises(ValueError):
            atom.charge = "invalid_number"

    # --------------------------------------------------------------------------
    # 4. AtomGroup.set_atom with slashes
    # --------------------------------------------------------------------------
    def test_atomgroup_set_atom_with_slashes(self):
        ag = rs_br.AtomGroup()
        a = rs_br.Atom(name="C3", symbol="C")
        ag.set_atom("/group_A/group_B/C3", a)

        # Pure Python equivalent
        ag_py = PyAtomGroup()
        a_py = PyAtom(name="C3", symbol="C")
        ag_py.set_atom("/group_A/group_B/C3", a_py)

        self.assertEqual(ag.get_number_of_atoms(), ag_py.get_number_of_atoms())
        self.assertEqual(ag.get_number_of_atoms(), 0)

        self.assertEqual(ag.get_number_of_all_atoms(), ag_py.get_number_of_all_atoms())
        self.assertEqual(ag.get_number_of_all_atoms(), 1)

        self.assertTrue(ag["group_A"]["group_B"].has_atomkey("C3"))
        self.assertEqual(
            ag["group_A"]["group_B"].has_atomkey("C3"),
            ag_py["group_A"]["group_B"].has_atomkey("C3"),
        )

        # Path check
        self.assertEqual(ag["group_A"]["group_B"]["C3"].path, "/group_A/group_B/C3")


if __name__ == "__main__":
    unittest.main()
