#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Validation test suite for Phase 1 (PR#14):
Validating proteindf_bridge.rs (PyO3 Rust bindings) core types, periodic table, linear algebra, and data structures.
"""

import unittest

# Pure Python implementation
from proteindf_bridge.periodictable import PeriodicTable as PyPT
from proteindf_bridge.vector import Vector as PyVec
from proteindf_bridge.matrix import Matrix as PyMat, SymmetricMatrix as PySymMat

# Rust PyO3 bindings
import proteindf_bridge.rs as rs_br


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
        rs_p1 = rs_br.Position([1.0, 2.0, 3.0])

        self.assertAlmostEqual(rs_p1.x, 1.0, places=6)
        self.assertAlmostEqual(rs_p1.y, 2.0, places=6)
        self.assertAlmostEqual(rs_p1.z, 3.0, places=6)

        # parse from string
        rs_from_str = rs_br.Position("1.5, -2.5, 3.5")
        self.assertAlmostEqual(rs_from_str.x, 1.5, places=6)
        self.assertAlmostEqual(rs_from_str.y, -2.5, places=6)
        self.assertAlmostEqual(rs_from_str.z, 3.5, places=6)

        # distance
        rs_p2 = rs_br.Position([4.0, 6.0, 3.0])
        self.assertAlmostEqual(rs_p1.distance_from(rs_p2), 5.0, places=6)
        self.assertAlmostEqual(rs_p1.square_distance_from(rs_p2), 25.0, places=6)

        # norm / length
        self.assertAlmostEqual(abs(rs_p1), 3.7416573867739413, places=6)
        rs_p1_norm = rs_br.Position(rs_p1)
        rs_p1_norm.norm()
        self.assertAlmostEqual(abs(rs_p1_norm), 1.0, places=6)

        # zero vector norm error
        zero_p = rs_br.Position([0.0, 0.0, 0.0])
        with self.assertRaises(rs_br.BrValueError):
            zero_p.norm()

        # dot & cross
        self.assertAlmostEqual(rs_p1.dot(rs_p2), 25.0, places=6)
        rs_cross = rs_p1.cross(rs_p2)
        self.assertAlmostEqual(rs_cross.x, -12.0, places=6)
        self.assertAlmostEqual(rs_cross.y, 9.0, places=6)
        self.assertAlmostEqual(rs_cross.z, -2.0, places=6)

        # rotate
        rot = rs_br.Matrix([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        rs_p_rot = rs_br.Position([1.0, 0.0, 0.0])
        rs_p_rot.rotate(rot)
        self.assertAlmostEqual(rs_p_rot.x, 0.0, places=6)
        self.assertAlmostEqual(rs_p_rot.y, 1.0, places=6)
        self.assertAlmostEqual(rs_p_rot.z, 0.0, places=6)

    # --------------------------------------------------------------------------
    # 6. Atom
    # --------------------------------------------------------------------------
    def test_atom(self):
        rs_atom = rs_br.Atom(symbol="C", name="CA", xyz=[1.0, 2.0, 3.0], charge=-0.2)

        self.assertEqual(rs_atom.atomic_number, 6)
        self.assertEqual(rs_atom.symbol, "C")
        self.assertEqual(rs_atom.name, "CA")
        self.assertAlmostEqual(rs_atom.charge, -0.2, places=6)
        self.assertEqual(rs_atom.is_real, True)
        self.assertAlmostEqual(rs_atom.vdw, 1.7, places=5)
        self.assertAlmostEqual(rs_atom.weight(), 12.01, places=5)

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
        rs_root = rs_br.AtomGroup(name="mol")

        # Residue 1
        rs_res1 = rs_br.AtomGroup(name="ALA")

        rs_n = rs_br.Atom(symbol="N", name="N", xyz=[0.0, 0.0, 0.0])
        rs_ca = rs_br.Atom(symbol="C", name="CA", xyz=[1.45, 0.0, 0.0])
        rs_c = rs_br.Atom(symbol="C", name="C", xyz=[2.0, 1.4, 0.0])

        for k, a in [("N", rs_n), ("CA", rs_ca), ("C", rs_c)]:
            rs_res1.set_atom(k, a)

        rs_root.set_group("ALA_1", rs_res1)

        # Counts
        self.assertEqual(rs_root.get_number_of_all_atoms(), 3)
        self.assertEqual(rs_root.get_number_of_groups(), 1)

        # Path updating
        rs_ca_ret = rs_root.get_group("ALA_1").get_atom("CA")
        self.assertEqual(rs_ca_ret.path, "/ALA_1/CA")

        # Center & Box
        rs_center = rs_root.center()
        self.assertAlmostEqual(rs_center.x, 1.15, places=5)
        self.assertAlmostEqual(rs_center.y, 1.4 / 3.0, places=5)
        self.assertAlmostEqual(rs_center.z, 0.0, places=5)

        rs_bmin, rs_bmax = rs_root.box()
        self.assertAlmostEqual(rs_bmin.x, 0.0, places=5)
        self.assertAlmostEqual(rs_bmin.y, 0.0, places=5)
        self.assertAlmostEqual(rs_bmin.z, 0.0, places=5)
        self.assertAlmostEqual(rs_bmax.x, 2.0, places=5)
        self.assertAlmostEqual(rs_bmax.y, 1.4, places=5)
        self.assertAlmostEqual(rs_bmax.z, 0.0, places=5)

        # Formula & counts
        self.assertEqual(rs_root.get_formula(), "C2N")
        self.assertEqual(rs_root.formula(), "C2N")
        self.assertEqual(rs_root.get_atom_kinds_count(), {"C": 2, "N": 1})

        # Bond setup
        rs_bond = rs_br.Bond()
        rs_bond.setup(rs_root)

        self.assertEqual(rs_root.get_number_of_bonds(), 0)
        rs_bonds = rs_root.get_bond_list()
        self.assertEqual(len(rs_bonds), 2)
        bond_pairs = {frozenset([b[0], b[1]]) for b in rs_bonds}
        expected_pairs = {
            frozenset(["/ALA_1/N", "/ALA_1/CA"]),
            frozenset(["/ALA_1/CA", "/ALA_1/C"]),
        }
        self.assertEqual(bond_pairs, expected_pairs)
        self.assertEqual(
            rs_bonds,
            [("/ALA_1/N", "/ALA_1/CA", 1), ("/ALA_1/CA", "/ALA_1/C", 1)],
        )

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

        # Polymorphic index access: group by key, atom by key, group by name, atom by name
        # 1. Group by key
        self.assertIsInstance(root["chainA"], rs_br.AtomGroup)
        self.assertEqual(root["chainA"].name, "A")

        # 2. Atom by key
        self.assertIsInstance(chain["CA"], rs_br.Atom)
        self.assertEqual(chain["CA"].name, "CA")

        # 3. Group by name (e.g. key is "grp_key", name is "A")
        sub_group = rs_br.AtomGroup(name="RES_NAME")
        root.set_group("res_key_1", sub_group)
        self.assertEqual(root["res_key_1"].name, "RES_NAME")
        self.assertEqual(root["RES_NAME"].name, "RES_NAME")

        # 4. Atom by name (e.g. key is "atm_1", name is "CB")
        cb_atom = rs_br.Atom(symbol="C", name="CB", xyz=[2.0, 3.0, 4.0])
        sub_group.set_atom("atm_1", cb_atom)
        self.assertEqual(sub_group["atm_1"].name, "CB")
        self.assertEqual(sub_group["CB"].name, "CB")

        # 5. Non-existent key raises KeyError
        with self.assertRaises(KeyError):
            _ = root["nonexistent"]
        with self.assertRaises(KeyError):
            _ = sub_group["nonexistent"]

        # 6. __setitem__ for both Atom and AtomGroup
        test_grp = rs_br.AtomGroup(name="test")
        test_grp["new_atom"] = rs_br.Atom(symbol="O", name="O1")
        self.assertTrue(test_grp.has_atom("new_atom"))
        self.assertIsInstance(test_grp["new_atom"], rs_br.Atom)

        test_grp["new_sub"] = rs_br.AtomGroup(name="sub")
        self.assertTrue(test_grp.has_group("new_sub"))
        self.assertIsInstance(test_grp["new_sub"], rs_br.AtomGroup)

        # del_atom / del_group
        removed_atom = chain.del_atom("CA")
        self.assertIsNotNone(removed_atom)
        self.assertFalse(chain.has_atom("CA"))

        removed_group = root.del_group("chainA")
        self.assertIsNotNone(removed_group)
        self.assertFalse(root.has_group("chainA"))

    def test_atom_group_formula_hill(self):
        """Test AtomGroup.get_formula() and AtomGroup.formula() Hill system formatting."""
        # 1. Empty group
        empty = rs_br.AtomGroup(name="empty")
        self.assertEqual(empty.get_formula(), "")
        self.assertEqual(empty.formula(), "")

        # 2. Carbon with H (methane: CH4)
        ch4 = rs_br.AtomGroup(name="methane")
        ch4.set_atom("C1", rs_br.Atom(symbol="C"))
        for i in range(1, 5):
            ch4.set_atom(f"H{i}", rs_br.Atom(symbol="H"))
        self.assertEqual(ch4.get_formula(), "CH4")
        self.assertEqual(ch4.formula(), "CH4")

        # 3. Carbon with H and other elements (ethanol: C2H6O)
        c2h6o = rs_br.AtomGroup(name="ethanol")
        c2h6o.set_atom("C1", rs_br.Atom(symbol="C"))
        c2h6o.set_atom("C2", rs_br.Atom(symbol="C"))
        for i in range(1, 7):
            c2h6o.set_atom(f"H{i}", rs_br.Atom(symbol="H"))
        c2h6o.set_atom("O1", rs_br.Atom(symbol="O"))
        self.assertEqual(c2h6o.get_formula(), "C2H6O")
        self.assertEqual(c2h6o.formula(), "C2H6O")

        # 4. Carbon without H (carbon dioxide: CO2)
        co2 = rs_br.AtomGroup(name="co2")
        co2.set_atom("C1", rs_br.Atom(symbol="C"))
        co2.set_atom("O1", rs_br.Atom(symbol="O"))
        co2.set_atom("O2", rs_br.Atom(symbol="O"))
        self.assertEqual(co2.get_formula(), "CO2")
        self.assertEqual(co2.formula(), "CO2")

        # 5. Non-carbon (water: H2O)
        h2o = rs_br.AtomGroup(name="water")
        h2o.set_atom("H1", rs_br.Atom(symbol="H"))
        h2o.set_atom("H2", rs_br.Atom(symbol="H"))
        h2o.set_atom("O1", rs_br.Atom(symbol="O"))
        self.assertEqual(h2o.get_formula(), "H2O")
        self.assertEqual(h2o.formula(), "H2O")

        # 6. Non-carbon alphabetical (sodium chloride: ClNa)
        nacl = rs_br.AtomGroup(name="nacl")
        nacl.set_atom("Na1", rs_br.Atom(symbol="Na"))
        nacl.set_atom("Cl1", rs_br.Atom(symbol="Cl"))
        self.assertEqual(nacl.get_formula(), "ClNa")
        self.assertEqual(nacl.formula(), "ClNa")

        # 7. With unknown dummy atom X
        grp_x = rs_br.AtomGroup(name="with_x")
        grp_x.set_atom("C1", rs_br.Atom(symbol="C"))
        grp_x.set_atom("X1", rs_br.Atom(symbol="X"))
        grp_x.set_atom("X2", rs_br.Atom(symbol="X"))
        self.assertEqual(grp_x.get_formula(), "CX2")
        self.assertEqual(grp_x.formula(), "CX2")


if __name__ == "__main__":
    unittest.main()
