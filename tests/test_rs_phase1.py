#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Comparison test suite for Phase 5 (PR#14):
Validating proteindf_bridge_rs (PyO3 Rust bindings) against pure Python proteindf_bridge.
"""

import math
import unittest
import numpy as np

# Pure Python implementation
import proteindf_bridge as py_br
from proteindf_bridge.periodictable import PeriodicTable as PyPT
from proteindf_bridge.vector import Vector as PyVec
from proteindf_bridge.matrix import Matrix as PyMat, SymmetricMatrix as PySymMat
from proteindf_bridge.position import Position as PyPos
from proteindf_bridge.atom import Atom as PyAtom
from proteindf_bridge.bond import Bond as PyBond
from proteindf_bridge.atomgroup import AtomGroup as PyAtomGroup
from proteindf_bridge.error import BrError as PyBrError, BrInputError as PyBrInputError, BrValueError as PyBrValueError

# Rust PyO3 bindings
import proteindf_bridge_rs as rs_br


class TestRsPhase1(unittest.TestCase):

    # --------------------------------------------------------------------------
    # 1. Exceptions
    # --------------------------------------------------------------------------
    def test_exceptions(self):
        self.assertTrue(issubclass(rs_br.BrError, Exception))
        self.assertTrue(issubclass(rs_br.BrInputError, rs_br.BrError))
        self.assertTrue(issubclass(rs_br.BrValueError, rs_br.BrError))

    # --------------------------------------------------------------------------
    # 2. PeriodicTable
    # --------------------------------------------------------------------------
    def test_periodic_table(self):
        rs_pt = rs_br.PeriodicTable
        # Symbols & Atomic numbers
        for num in [1, 6, 7, 8, 12, 20, 26, 29]:
            self.assertEqual(rs_pt.get_symbol(num), PyPT.get_symbol(num))
            sym = PyPT.get_symbol(num)
            self.assertEqual(rs_pt.get_atomic_number(sym), PyPT.get_atomic_number(sym))

        # Weights & VDW radii
        for sym in ["H", "C", "N", "O", "S", "Fe"]:
            self.assertAlmostEqual(rs_pt.atomic_weight(sym), PyPT.atomic_weight(sym), places=5)
            self.assertAlmostEqual(rs_pt.get_weight(sym), PyPT.atomic_weight(sym), places=5)
            self.assertAlmostEqual(rs_pt.vdw(sym), PyPT.vdw(sym), places=5)
            self.assertAlmostEqual(rs_pt.get_vdw(sym), PyPT.vdw(sym), places=5)

        # Total number of atoms
        self.assertEqual(rs_pt.get_num_of_atoms(), PyPT.get_num_of_atoms())

        # Error cases raise BrValueError
        with self.assertRaises(rs_br.BrValueError):
            rs_pt.get_symbol(999)
        with self.assertRaises(rs_br.BrValueError):
            rs_pt.get_atomic_number("UnknownElement")

    # --------------------------------------------------------------------------
    # 3. Vector
    # --------------------------------------------------------------------------
    def test_vector(self):
        py_v = PyVec([1.0, -2.5, 3.0, 4.2])
        rs_v = rs_br.Vector([1.0, -2.5, 3.0, 4.2])

        self.assertEqual(len(rs_v), len(py_v))
        self.assertEqual(rs_v.size(), py_v.size())
        for i in range(len(py_v)):
            self.assertAlmostEqual(rs_v[i], py_v[i], places=6)
            self.assertAlmostEqual(rs_v.get(i), py_v.get(i), places=6)

        # max, min, abs
        self.assertAlmostEqual(rs_v.max, py_v.max, places=6)
        self.assertAlmostEqual(rs_v.min, py_v.min, places=6)
        rs_abs = rs_v.abs()
        py_abs = py_v.abs()
        self.assertEqual(rs_abs.to_list(), py_abs.to_list())

        # set / indexing
        rs_v[1] = 5.5
        py_v[1] = 5.5
        self.assertAlmostEqual(rs_v[1], py_v[1], places=6)

        # addition / subtraction
        py_v2 = PyVec([0.5, 1.0, -1.0, 2.0])
        rs_v2 = rs_br.Vector([0.5, 1.0, -1.0, 2.0])
        rs_add = rs_v + rs_v2
        py_add = py_v + py_v2
        self.assertEqual(rs_add.to_list(), py_add.to_list())

        rs_sub = rs_v - rs_v2
        py_sub = py_v - py_v2
        self.assertEqual(rs_sub.to_list(), py_sub.to_list())

        # scalar mul & dot
        rs_smul = rs_v * 2.0
        py_smul = py_v * 2.0
        self.assertEqual(rs_smul.to_list(), py_smul.to_list())

        rs_dot = rs_v * rs_v2
        py_dot = py_v * py_v2
        self.assertAlmostEqual(rs_dot, py_dot, places=6)
        self.assertAlmostEqual(rs_v.dot(rs_v2), py_v.data.dot(py_v2.data), places=6)

        # flip & argsort
        rs_flip = rs_v.flip()
        py_flip = py_v.flip()
        self.assertEqual(rs_flip.to_list(), py_flip.to_list())

        rs_sort = rs_v.argsort()
        py_sort = py_v.argsort()
        self.assertEqual(rs_sort.to_list(), py_sort.to_list())

        # resize
        rs_v.resize(6)
        py_v.resize(6)
        self.assertEqual(rs_v.to_list(), py_v.to_list())

    # --------------------------------------------------------------------------
    # 4. Matrix & SymmetricMatrix
    # --------------------------------------------------------------------------
    def test_matrix(self):
        py_m = PyMat([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        rs_m = rs_br.Matrix([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])

        self.assertEqual(rs_m.rows, py_m.rows)
        self.assertEqual(rs_m.cols, py_m.cols)
        for r in range(py_m.rows):
            for c in range(py_m.cols):
                self.assertAlmostEqual(rs_m.get(r, c), py_m.get(r, c), places=6)

        # transpose (in-place in Python proteindf_bridge)
        rs_m_t = rs_br.Matrix([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        py_m_t = PyMat([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        rs_t = rs_m_t.transpose()
        py_t = py_m_t.transpose()
        self.assertEqual(rs_t.rows, py_t.rows)
        self.assertEqual(rs_t.cols, py_t.cols)
        self.assertEqual(rs_t.to_list(), py_t._data.tolist())
        self.assertEqual(rs_m_t.rows, 3)
        self.assertEqual(py_m_t.rows, 3)

        # Matrix * Matrix
        py_m2 = PyMat([[1.0, 0.0], [2.0, -1.0], [0.0, 3.0]])
        rs_m2 = rs_br.Matrix([[1.0, 0.0], [2.0, -1.0], [0.0, 3.0]])
        rs_prod = rs_m * rs_m2
        py_prod = py_m * py_m2
        for r in range(rs_prod.rows):
            for c in range(rs_prod.cols):
                self.assertAlmostEqual(rs_prod.get(r, c), py_prod.get(r, c), places=6)

        # Matrix * Vector
        py_v = PyVec([1.0, 2.0, 3.0])
        rs_v = rs_br.Vector([1.0, 2.0, 3.0])
        rs_mv = rs_m * rs_v
        py_mv = py_m * py_v
        self.assertEqual(rs_mv.to_list(), py_mv.to_list())

        # select sub-matrix
        rs_sel = rs_m.select(0, 1, 2, 3)
        py_sel = py_m.select(0, 1, 2, 3)
        self.assertEqual(rs_sel.to_list(), py_sel._data.tolist())

        # inverse
        py_sq = PyMat([[4.0, 7.0], [2.0, 6.0]])
        rs_sq = rs_br.Matrix([[4.0, 7.0], [2.0, 6.0]])
        rs_inv = rs_sq.inverse()
        py_inv = py_sq.inverse()
        for r in range(2):
            for c in range(2):
                self.assertAlmostEqual(rs_inv.get(r, c), py_inv.get(r, c), places=5)

    def test_symmetric_matrix(self):
        py_sm = PySymMat(3)
        rs_sm = rs_br.SymmetricMatrix(3)

        self.assertEqual(rs_sm.dim, py_sm.dim)
        py_sm.set(0, 1, 3.0)
        rs_sm.set(0, 1, 3.0)
        self.assertAlmostEqual(rs_sm.get(0, 1), 3.0)
        self.assertAlmostEqual(rs_sm.get(1, 0), 3.0)
        self.assertAlmostEqual(py_sm.get(0, 1), 3.0)
        self.assertAlmostEqual(py_sm.get(1, 0), 3.0)

        # add
        py_sm.add(0, 2, 4.0)
        rs_sm.add(0, 2, 4.0)
        self.assertAlmostEqual(rs_sm.get(0, 2), 4.0)
        self.assertAlmostEqual(rs_sm.get(2, 0), 4.0)

        # eig
        # Create a symmetric matrix with known eigenvalues
        # A = [[2, 1], [1, 2]], eigenvalues: 1.0, 3.0
        rs_sym2 = rs_br.SymmetricMatrix([[2.0, 1.0], [1.0, 2.0]])
        py_sym2 = PySymMat(2)
        py_sym2.set(0, 0, 2.0)
        py_sym2.set(1, 1, 2.0)
        py_sym2.set(0, 1, 1.0)

        rs_w, rs_v = rs_sym2.eig()
        py_w, py_v = py_sym2.eig()

        # Check eigenvalues match
        for i in range(2):
            self.assertAlmostEqual(rs_w[i], py_w[i], places=5)

        # Check A * v_i = lambda_i * v_i
        for i in range(2):
            val = rs_w[i]
            vec = rs_v.get_row_vector(i)
            # A * vec
            A_gen = rs_sym2.get_general_matrix()
            Av = A_gen * vec
            scaled = vec * val
            for j in range(2):
                self.assertAlmostEqual(Av[j], scaled[j], places=5)

    # --------------------------------------------------------------------------
    # 5. Position
    # --------------------------------------------------------------------------
    def test_position(self):
        py_p1 = PyPos([1.0, 2.0, 3.0])
        rs_p1 = rs_br.Position([1.0, 2.0, 3.0])

        self.assertAlmostEqual(rs_p1.x, py_p1.x, places=6)
        self.assertAlmostEqual(rs_p1.y, py_p1.y, places=6)
        self.assertAlmostEqual(rs_p1.z, py_p1.z, places=6)

        # parse from string
        rs_from_str = rs_br.Position("1.5, -2.5, 3.5")
        py_from_str = PyPos("1.5, -2.5, 3.5")
        self.assertAlmostEqual(rs_from_str.x, py_from_str.x, places=6)
        self.assertAlmostEqual(rs_from_str.y, py_from_str.y, places=6)
        self.assertAlmostEqual(rs_from_str.z, py_from_str.z, places=6)

        # distance
        py_p2 = PyPos([4.0, 6.0, 3.0])
        rs_p2 = rs_br.Position([4.0, 6.0, 3.0])
        self.assertAlmostEqual(rs_p1.distance_from(rs_p2), py_p1.distance_from(py_p2), places=6)
        self.assertAlmostEqual(rs_p1.square_distance_from(rs_p2), py_p1.square_distance_from(py_p2), places=6)

        # norm / length
        self.assertAlmostEqual(abs(rs_p1), abs(py_p1), places=6)
        rs_p1_norm = rs_br.Position(rs_p1)
        rs_p1_norm.norm()
        self.assertAlmostEqual(abs(rs_p1_norm), 1.0, places=6)

        # zero vector norm error
        zero_p = rs_br.Position([0.0, 0.0, 0.0])
        with self.assertRaises(rs_br.BrValueError):
            zero_p.norm()

        # dot & cross
        self.assertAlmostEqual(rs_p1.dot(rs_p2), py_p1.dot(py_p2), places=6)
        rs_cross = rs_p1.cross(rs_p2)
        py_cross = py_p1.cross(py_p2)
        self.assertAlmostEqual(rs_cross.x, py_cross.x, places=6)
        self.assertAlmostEqual(rs_cross.y, py_cross.y, places=6)
        self.assertAlmostEqual(rs_cross.z, py_cross.z, places=6)

        # rotate
        rot = rs_br.Matrix([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        py_rot = PyMat([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        rs_p_rot = rs_br.Position([1.0, 0.0, 0.0])
        py_p_rot = PyPos([1.0, 0.0, 0.0])
        rs_p_rot.rotate(rot)
        py_p_rot.rotate(py_rot)
        self.assertAlmostEqual(rs_p_rot.x, py_p_rot.x, places=6)
        self.assertAlmostEqual(rs_p_rot.y, py_p_rot.y, places=6)
        self.assertAlmostEqual(rs_p_rot.z, py_p_rot.z, places=6)

    # --------------------------------------------------------------------------
    # 6. Atom
    # --------------------------------------------------------------------------
    def test_atom(self):
        py_atom = PyAtom(symbol="C", name="CA", xyz=[1.0, 2.0, 3.0], charge=-0.2)
        rs_atom = rs_br.Atom(symbol="C", name="CA", xyz=[1.0, 2.0, 3.0], charge=-0.2)

        self.assertEqual(rs_atom.atomic_number, py_atom.atomic_number)
        self.assertEqual(rs_atom.symbol, py_atom.symbol)
        self.assertEqual(rs_atom.name, py_atom.name)
        self.assertAlmostEqual(rs_atom.charge, py_atom.charge, places=6)
        self.assertEqual(rs_atom.is_real, py_atom.is_real)
        self.assertAlmostEqual(rs_atom.vdw, py_atom.vdw, places=5)
        self.assertAlmostEqual(rs_atom.weight(), py_atom.weight(), places=5)

        # move_to and shift_by chaining
        rs_atom.shift_by([1.0, 1.0, 1.0]).move_to([5.0, 5.0, 5.0])
        self.assertAlmostEqual(rs_atom.xyz.x, 5.0, places=6)

        # raw data roundtrip
        raw = rs_atom.get_raw_data()
        rs_atom2 = rs_br.Atom()
        rs_atom2.set_by_raw_data(raw)
        self.assertEqual(rs_atom.atomic_number, rs_atom2.atomic_number)
        self.assertEqual(rs_atom.name, rs_atom2.name)
        self.assertAlmostEqual(rs_atom.xyz.x, rs_atom2.xyz.x, places=6)

    # --------------------------------------------------------------------------
    # 7. AtomGroup & Bond
    # --------------------------------------------------------------------------
    def test_atom_group_hierarchy_and_bonds(self):
        # Build hierarchy: root -> chainA -> res1 -> atoms
        py_root = PyAtomGroup(name="mol")
        rs_root = rs_br.AtomGroup(name="mol")

        # Residue 1
        py_res1 = PyAtomGroup(name="ALA")
        rs_res1 = rs_br.AtomGroup(name="ALA")

        py_n = PyAtom(symbol="N", name="N", xyz=[0.0, 0.0, 0.0])
        rs_n = rs_br.Atom(symbol="N", name="N", xyz=[0.0, 0.0, 0.0])
        py_ca = PyAtom(symbol="C", name="CA", xyz=[1.45, 0.0, 0.0])
        rs_ca = rs_br.Atom(symbol="C", name="CA", xyz=[1.45, 0.0, 0.0])
        py_c = PyAtom(symbol="C", name="C", xyz=[2.0, 1.4, 0.0])
        rs_c = rs_br.Atom(symbol="C", name="C", xyz=[2.0, 1.4, 0.0])

        for k, a in [("N", py_n), ("CA", py_ca), ("C", py_c)]:
            py_res1.set_atom(k, a)
        for k, a in [("N", rs_n), ("CA", rs_ca), ("C", rs_c)]:
            rs_res1.set_atom(k, a)

        py_root.set_group("ALA_1", py_res1)
        rs_root.set_group("ALA_1", rs_res1)

        # Counts
        self.assertEqual(rs_root.get_number_of_all_atoms(), py_root.get_number_of_all_atoms())
        self.assertEqual(rs_root.get_number_of_groups(), py_root.get_number_of_groups())

        # Path updating
        rs_ca_ret = rs_root.get_group("ALA_1").get_atom("CA")
        self.assertEqual(rs_ca_ret.path, "/ALA_1/CA")

        # Center & Box
        rs_center = rs_root.center()
        py_center = py_root.center()
        self.assertAlmostEqual(rs_center.x, py_center.x, places=5)
        self.assertAlmostEqual(rs_center.y, py_center.y, places=5)
        self.assertAlmostEqual(rs_center.z, py_center.z, places=5)

        rs_bmin, rs_bmax = rs_root.box()
        py_bmin, py_bmax = py_root.box()
        self.assertAlmostEqual(rs_bmin.x, py_bmin.x, places=5)
        self.assertAlmostEqual(rs_bmax.x, py_bmax.x, places=5)

        # Formula & counts
        self.assertEqual(rs_root.get_formula(), py_root.get_formula())
        self.assertEqual(rs_root.get_atom_kinds_count(), py_root.get_atom_kinds_count())

        # Bond setup
        py_bond = PyBond()
        py_bond.setup(py_root)

        rs_bond = rs_br.Bond()
        rs_bond.setup(rs_root)

        self.assertEqual(rs_root.get_number_of_bonds(), py_root.get_number_of_bonds())
        rs_bonds = rs_root.get_bond_list()
        py_bonds = py_root.get_bond_list()
        self.assertEqual(len(rs_bonds), len(py_bonds))

    # --------------------------------------------------------------------------
    # 8. AtomGroup Set Operations (&, |, ^)
    # --------------------------------------------------------------------------
    def test_atom_group_set_ops(self):
        rs_g1 = rs_br.AtomGroup(name="g1")
        rs_g1.set_atom("C1", rs_br.Atom("C", name="C1"))
        rs_g1.set_atom("O1", rs_br.Atom("O", name="O1"))

        rs_g2 = rs_br.AtomGroup(name="g2")
        rs_g2.set_atom("C1", rs_br.Atom("C", name="C1"))
        rs_g2.set_atom("N1", rs_br.Atom("N", name="N1"))

        # Intersection &
        inter = rs_g1 & rs_g2
        self.assertEqual(inter.get_number_of_atoms(), 1)
        self.assertTrue(inter.has_atom("C1"))

        # Union |
        union = rs_g1 | rs_g2
        self.assertEqual(union.get_number_of_atoms(), 3)
        self.assertTrue(union.has_atom("C1"))
        self.assertTrue(union.has_atom("O1"))
        self.assertTrue(union.has_atom("N1"))

        # Symmetric difference ^
        diff = rs_g1 ^ rs_g2
        self.assertEqual(diff.get_number_of_atoms(), 2)
        self.assertTrue(diff.has_atom("O1"))
        self.assertTrue(diff.has_atom("N1"))
        self.assertFalse(diff.has_atom("C1"))

    # --------------------------------------------------------------------------
    # 9. AtomGroup Accessors and Manipulation
    # --------------------------------------------------------------------------
    def test_atom_group_methods(self):
        root = rs_br.AtomGroup(name="protein")
        chain = rs_br.AtomGroup(name="A")
        atom = rs_br.Atom(symbol="C", name="CA", xyz=[1.0, 2.0, 3.0])

        chain.set_atom("CA", atom)
        root.set_group("chainA", chain)

        # has / get / del
        self.assertTrue(root.has_group("chainA"))
        self.assertTrue(chain.has_atom("CA"))
        self.assertFalse(root.has_atom("CA"))

        # atoms() and groups()
        atom_entries = chain.atoms()
        self.assertEqual(len(atom_entries), 1)
        self.assertEqual(atom_entries[0][0], "CA")
        self.assertEqual(atom_entries[0][1].name, "CA")

        group_entries = root.groups()
        self.assertEqual(len(group_entries), 1)
        self.assertEqual(group_entries[0][0], "chainA")
        self.assertEqual(group_entries[0][1].name, "A")

        # index access root["chainA"]
        self.assertEqual(root["chainA"].name, "A")
        with self.assertRaises(KeyError):
            _ = root["nonexistent"]

        # del_atom / del_group
        removed_atom = chain.del_atom("CA")
        self.assertIsNotNone(removed_atom)
        self.assertFalse(chain.has_atom("CA"))

        removed_group = root.del_group("chainA")
        self.assertIsNotNone(removed_group)
        self.assertFalse(root.has_group("chainA"))


if __name__ == "__main__":
    unittest.main()
