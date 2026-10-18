#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test suite validating rs.AtomGroup aggregation methods:
- get_number_of_all_atoms
- sum_of_atomic_number
- get_atom_kinds_count
- get_formula / formula
- weight
- center
- box / get_box
"""

import os
import unittest
import proteindf_bridge.rs as rs


class TestRsAgAggregates(unittest.TestCase):
    """Test suite validating rs.AtomGroup aggregation methods."""

    def setUp(self):
        self.data_dir = os.path.abspath(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                "rust",
                "crates",
                "proteindf-bridge",
                "tests",
                "data",
            )
        )
        self.pdb_2fb4 = os.path.join(self.data_dir, "2FB4.pdb")
        self.pdb_2mgo = os.path.join(self.data_dir, "2MGO.pdb")

    def test_2fb4_aggregates(self):
        """Validate aggregation methods on 2FB4.pdb (3417 atoms)."""
        self.assertTrue(os.path.exists(self.pdb_2fb4), f"2FB4.pdb not found at {self.pdb_2fb4}")
        ag = rs.Pdb(self.pdb_2fb4).get_atomgroup()

        # get_number_of_all_atoms
        self.assertEqual(ag.get_number_of_all_atoms(), 3417)

        # sum_of_atomic_number
        self.assertAlmostEqual(ag.sum_of_atomic_number(), 22759.0, places=5)

        # get_atom_kinds_count
        expected_kinds = {"C": 2072, "N": 561, "O": 768, "S": 16}
        self.assertEqual(ag.get_atom_kinds_count(), expected_kinds)

        # get_formula and formula
        expected_formula = "C2072N561O768S16"
        self.assertEqual(ag.get_formula(), expected_formula)
        self.assertEqual(ag.formula(), expected_formula)

        # weight
        self.assertAlmostEqual(ag.weight, 45545.45000000009, places=4)

        # center
        c = ag.center()
        self.assertAlmostEqual(c.x, 38.896935323383104, places=5)
        self.assertAlmostEqual(c.y, 72.71806292069068, places=5)
        self.assertAlmostEqual(c.z, 29.23772022241732, places=5)

        # box and get_box
        bmin, bmax = ag.box()
        self.assertAlmostEqual(bmin.x, 9.259, places=5)
        self.assertAlmostEqual(bmin.y, 34.066, places=5)
        self.assertAlmostEqual(bmin.z, 0.839, places=5)
        self.assertAlmostEqual(bmax.x, 69.869, places=5)
        self.assertAlmostEqual(bmax.y, 113.35, places=5)
        self.assertAlmostEqual(bmax.z, 56.707, places=5)

        gbmin, gbmax = ag.get_box()
        self.assertAlmostEqual(gbmin.x, 9.259, places=5)
        self.assertAlmostEqual(gbmin.y, 34.066, places=5)
        self.assertAlmostEqual(gbmin.z, 0.839, places=5)
        self.assertAlmostEqual(gbmax.x, 69.869, places=5)
        self.assertAlmostEqual(gbmax.y, 113.35, places=5)
        self.assertAlmostEqual(gbmax.z, 56.707, places=5)

    def test_2mgo_aggregates(self):
        """Validate aggregation methods on 2MGO.pdb (2680 atoms)."""
        self.assertTrue(os.path.exists(self.pdb_2mgo), f"2MGO.pdb not found at {self.pdb_2mgo}")
        ag = rs.Pdb(self.pdb_2mgo).get_atomgroup()

        # get_number_of_all_atoms
        self.assertEqual(ag.get_number_of_all_atoms(), 2680)

        # sum_of_atomic_number
        self.assertAlmostEqual(ag.sum_of_atomic_number(), 10720.0, places=5)

        # get_atom_kinds_count
        expected_kinds = {"C": 860, "H": 1300, "N": 220, "O": 260, "S": 40}
        self.assertEqual(ag.get_atom_kinds_count(), expected_kinds)

        # get_formula and formula
        expected_formula = "C860H1300N220O260S40"
        self.assertEqual(ag.get_formula(), expected_formula)
        self.assertEqual(ag.formula(), expected_formula)

        # weight
        self.assertAlmostEqual(ag.weight, 20164.00000000008, places=4)

        # center
        c = ag.center()
        self.assertAlmostEqual(c.x, -0.25817649253731345, places=5)
        self.assertAlmostEqual(c.y, -0.13786604477611944, places=5)
        self.assertAlmostEqual(c.z, 0.3300723880597014, places=5)

        # box and get_box
        bmin, bmax = ag.box()
        self.assertAlmostEqual(bmin.x, -9.846, places=5)
        self.assertAlmostEqual(bmin.y, -9.485, places=5)
        self.assertAlmostEqual(bmin.z, -9.729, places=5)
        self.assertAlmostEqual(bmax.x, 9.269, places=5)
        self.assertAlmostEqual(bmax.y, 11.487, places=5)
        self.assertAlmostEqual(bmax.z, 11.773, places=5)

        gbmin, gbmax = ag.get_box()
        self.assertAlmostEqual(gbmin.x, -9.846, places=5)
        self.assertAlmostEqual(gbmin.y, -9.485, places=5)
        self.assertAlmostEqual(gbmin.z, -9.729, places=5)
        self.assertAlmostEqual(gbmax.x, 9.269, places=5)
        self.assertAlmostEqual(gbmax.y, 11.487, places=5)
        self.assertAlmostEqual(gbmax.z, 11.773, places=5)

    def test_empty_atomgroup_aggregates(self):
        """Validate aggregation methods on an empty AtomGroup."""
        ag = rs.AtomGroup("empty")

        self.assertEqual(ag.get_number_of_all_atoms(), 0)
        self.assertAlmostEqual(ag.sum_of_atomic_number(), 0.0, places=5)
        self.assertEqual(ag.get_atom_kinds_count(), {})
        self.assertEqual(ag.get_formula(), "")
        self.assertEqual(ag.formula(), "")
        self.assertAlmostEqual(ag.weight, 0.0, places=5)

        c = ag.center()
        self.assertAlmostEqual(c.x, 0.0, places=5)
        self.assertAlmostEqual(c.y, 0.0, places=5)
        self.assertAlmostEqual(c.z, 0.0, places=5)

        bmin, bmax = ag.box()
        self.assertAlmostEqual(bmin.x, 0.0, places=5)
        self.assertAlmostEqual(bmin.y, 0.0, places=5)
        self.assertAlmostEqual(bmin.z, 0.0, places=5)
        self.assertAlmostEqual(bmax.x, 0.0, places=5)
        self.assertAlmostEqual(bmax.y, 0.0, places=5)
        self.assertAlmostEqual(bmax.z, 0.0, places=5)

        gbmin, gbmax = ag.get_box()
        self.assertAlmostEqual(gbmin.x, 0.0, places=5)
        self.assertAlmostEqual(gbmin.y, 0.0, places=5)
        self.assertAlmostEqual(gbmin.z, 0.0, places=5)
        self.assertAlmostEqual(gbmax.x, 0.0, places=5)
        self.assertAlmostEqual(gbmax.y, 0.0, places=5)
        self.assertAlmostEqual(gbmax.z, 0.0, places=5)

    def test_single_atom_aggregates(self):
        """Validate aggregation methods on an AtomGroup containing a single atom."""
        ag = rs.AtomGroup("single")
        a = rs.Atom(name="CA", symbol="C", xyz=rs.Position(1.0, 2.0, 3.0), atomic_number=6, weight=12.011)
        ag.set_atom("A1", a)

        self.assertEqual(ag.get_number_of_all_atoms(), 1)
        self.assertAlmostEqual(ag.sum_of_atomic_number(), 6.0, places=5)
        self.assertEqual(ag.get_atom_kinds_count(), {"C": 1})
        self.assertEqual(ag.get_formula(), "C")
        self.assertEqual(ag.formula(), "C")
        self.assertAlmostEqual(ag.weight, 12.01, places=4)

        c = ag.center()
        self.assertAlmostEqual(c.x, 1.0, places=5)
        self.assertAlmostEqual(c.y, 2.0, places=5)
        self.assertAlmostEqual(c.z, 3.0, places=5)

        bmin, bmax = ag.box()
        self.assertAlmostEqual(bmin.x, 1.0, places=5)
        self.assertAlmostEqual(bmin.y, 2.0, places=5)
        self.assertAlmostEqual(bmin.z, 3.0, places=5)
        self.assertAlmostEqual(bmax.x, 1.0, places=5)
        self.assertAlmostEqual(bmax.y, 2.0, places=5)
        self.assertAlmostEqual(bmax.z, 3.0, places=5)

        gbmin, gbmax = ag.get_box()
        self.assertAlmostEqual(gbmin.x, 1.0, places=5)
        self.assertAlmostEqual(gbmin.y, 2.0, places=5)
        self.assertAlmostEqual(gbmin.z, 3.0, places=5)
        self.assertAlmostEqual(gbmax.x, 1.0, places=5)
        self.assertAlmostEqual(gbmax.y, 2.0, places=5)
        self.assertAlmostEqual(gbmax.z, 3.0, places=5)

    def test_hierarchy_and_weightless_atom(self):
        """Validate aggregation methods with nested groups, sorting, and dummy atoms."""
        root = rs.AtomGroup("root")
        root.sort_groups = "nice"
        root.sort_atoms = "nice"

        sub1 = rs.AtomGroup("sub1")
        sub2 = rs.AtomGroup("sub2")

        # Unknown/dummy atom with atomic_number=0 (weight should be 0.0)
        dummy = rs.Atom(name="XX", symbol="X", xyz=rs.Position(0.0, 0.0, 0.0), atomic_number=0)
        sub1.set_atom("atom_2", dummy)

        # Real atom
        ca = rs.Atom(name="CA", symbol="C", xyz=rs.Position(10.0, 20.0, 30.0), atomic_number=6)
        sub1.set_atom("atom_1", ca)

        # Atom in sub2
        n = rs.Atom(name="N", symbol="N", xyz=rs.Position(20.0, 40.0, 60.0), atomic_number=7)
        sub2.set_atom("atom_1", n)

        root.set_group("grp2", sub2)
        root.set_group("grp1", sub1)

        # 3 atoms total: C, N, X
        self.assertEqual(root.get_number_of_all_atoms(), 3)
        self.assertAlmostEqual(root.sum_of_atomic_number(), 13.0, places=5)  # 6 + 7 + 0

        # Kinds: X has symbol X if PeriodicTable has it or symbol() ok
        kinds = root.get_atom_kinds_count()
        self.assertEqual(kinds.get("C"), 1)
        self.assertEqual(kinds.get("N"), 1)
        self.assertEqual(kinds.get("X"), 1)

        # Formula: C is present, so C then others (N), then X at end => "CNX"
        self.assertEqual(root.get_formula(), "CNX")
        self.assertEqual(root.formula(), "CNX")

        # Weight: C(12.01) + N(14.01) + X(0.0)
        expected_weight = 12.01 + 14.01
        self.assertAlmostEqual(root.weight, expected_weight, places=3)

        # Center: weighted sum of coordinates
        # sub1: (10, 20, 30) + (0, 0, 0) => avg = (5, 10, 15), count=2
        # sub2: (20, 40, 60) => avg = (20, 40, 60), count=1
        # root center: (5*2 + 20*1)/3 = 30/3 = 10, y: (10*2 + 40*1)/3 = 20, z: (15*2 + 60*1)/3 = 30
        c = root.center()
        self.assertAlmostEqual(c.x, 10.0, places=5)
        self.assertAlmostEqual(c.y, 20.0, places=5)
        self.assertAlmostEqual(c.z, 30.0, places=5)

        # Box: min=(0, 0, 0), max=(20, 40, 60)
        bmin, bmax = root.box()
        self.assertAlmostEqual(bmin.x, 0.0, places=5)
        self.assertAlmostEqual(bmin.y, 0.0, places=5)
        self.assertAlmostEqual(bmin.z, 0.0, places=5)
        self.assertAlmostEqual(bmax.x, 20.0, places=5)
        self.assertAlmostEqual(bmax.y, 40.0, places=5)
        self.assertAlmostEqual(bmax.z, 60.0, places=5)


if __name__ == "__main__":
    unittest.main()
