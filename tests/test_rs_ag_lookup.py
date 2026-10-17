#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test suite validating rs.AtomGroup lookup and selection methods (TASK_rs-ag-lookup).
"""

import glob
import os
import unittest

import proteindf_bridge.rs as rs


class TestRsAgLookup(unittest.TestCase):
    """Test suite for rs.AtomGroup lookup, selection, and hierarchy navigation."""

    def setUp(self):
        self.data_dir = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "proteindf_bridge", "data")
        )

    def test_brd_files_lookup(self):
        """Criterion 1: Verify get_atom_keys, has_atomname, has_groupname, and pickup_atoms

        against expected structure and queries for all .brd files in proteindf_bridge/data/.
        """
        brd_files = sorted(glob.glob(os.path.join(self.data_dir, "*.brd")))
        self.assertGreater(len(brd_files), 0, "No .brd files found in data directory")

        for brd_path in brd_files:
            fname = os.path.basename(brd_path)
            with self.subTest(file=fname):
                rs_ag = rs.load_atomgroup(brd_path)
                if fname.startswith("ACE_ALA_NME"):
                    self.assertEqual(rs_ag.get_group_list(), ["1", "2", "3"])
                    self.assertEqual(rs_ag.get_atom_keys(), [])
                    self.assertEqual(rs_ag.get_number_of_all_atoms(), 22)
                    self.assertTrue(rs_ag.has_group("1"))
                    self.assertTrue(rs_ag.has_groupname("ACE"))
                    self.assertTrue(rs_ag.has_groupname("ALA"))
                    self.assertTrue(rs_ag.has_groupname("NME"))
                    self.assertFalse(rs_ag.has_groupname("__NON_EXISTENT__"))
                    self.assertEqual(len(rs_ag.pickup_atoms("CA")), 1)
                    self.assertEqual(len(rs_ag.pickup_atoms("C")), 2)
                    self.assertEqual(len(rs_ag.pickup_atoms("N")), 2)
                    self.assertEqual(len(rs_ag.pickup_atoms("non_existent")), 0)
                elif fname == "NML.brd":
                    self.assertEqual(rs_ag.get_group_list(), [])
                    self.assertEqual(len(rs_ag.get_atom_keys()), 12)
                    self.assertEqual(rs_ag.get_number_of_all_atoms(), 12)
                    self.assertTrue(rs_ag.has_atomkey("0"))
                    self.assertTrue(rs_ag.has_atomname("C1"))
                    self.assertTrue(rs_ag.has_atomname("  C1  "))
                    self.assertFalse(rs_ag.has_atomname("__NON_EXISTENT__"))
                    self.assertEqual(len(rs_ag.pickup_atoms("C1")), 1)
                    self.assertEqual(len(rs_ag.pickup_atoms("non_existent")), 0)
                elif fname == "NML_trans.brd":
                    self.assertEqual(rs_ag.get_group_list(), [])
                    self.assertEqual(len(rs_ag.get_atom_keys()), 12)
                    self.assertEqual(rs_ag.get_number_of_all_atoms(), 12)
                    self.assertTrue(rs_ag.has_atomkey("0"))
                    self.assertTrue(rs_ag.has_atomname(""))
                    self.assertFalse(rs_ag.has_atomname("__NON_EXISTENT__"))
                    self.assertEqual(len(rs_ag.pickup_atoms("non_existent")), 0)

    def test_2mgo_pdb_lookup(self):
        """Criterion 1: Verify get_atom_keys, has_atomname, has_groupname, and pickup_atoms

        against expected counts and sample atoms for 2MGO.pdb.
        """
        pdb_path = os.path.join(self.data_dir, "2MGO.pdb")
        self.assertTrue(os.path.exists(pdb_path), "2MGO.pdb fixture not found")

        rs_ag = rs.Pdb(pdb_path).get_atomgroup()

        self.assertEqual(rs_ag.get_atom_keys(), [])
        self.assertEqual(rs_ag.get_number_of_all_atoms(), 2680)
        self.assertEqual(len(rs_ag.get_group_list()), 20)
        self.assertTrue(rs_ag.has_groupname("model_1"))
        self.assertTrue(rs_ag.has_groupname("model_20"))
        self.assertFalse(rs_ag.has_groupname("__NON_EXISTENT_GROUP__"))

        m1 = rs_ag.get_group("model_1")
        self.assertTrue(m1.has_groupname("A"))
        chain_a = m1.get_group("A")
        self.assertTrue(chain_a.has_group("1"))
        self.assertTrue(chain_a.has_groupname("CYS"))
        res1 = chain_a.get_group("1")
        self.assertTrue(res1.has_atomname("CA"))
        self.assertTrue(res1.has_atomname("  CA  "))
        self.assertFalse(res1.has_atomname("__NON_EXISTENT_ATOM__"))

        expected_queries = {
            "CA": (180, "/model_1/A/1/2_CA", "CA", [3.791, 1.923, 0.075]),
            "CB": (160, "/model_1/A/1/5_CB", "CB", [3.435, 1.946, -1.418]),
            "N": (180, "/model_1/A/1/1_N", "N", [4.874, 2.855, 0.366]),
            "C": (180, "/model_1/A/1/3_C", "C", [4.157, 0.518, 0.5]),
            "O": (180, "/model_1/A/1/4_O", "O", [5.232, 0.286, 1.045]),
            "SG": (40, "/model_1/A/1/6_SG", "SG", [2.665, 3.502, -1.99]),
            "1": (0, None, None, None),
            "10": (0, None, None, None),
            "100": (0, None, None, None),
            "non_existent": (0, None, None, None),
        }

        for query, (exp_count, exp_path, exp_name, exp_pos) in expected_queries.items():
            picked = rs_ag.pickup_atoms(query)
            self.assertEqual(
                len(picked),
                exp_count,
                f"Mismatch in count of picked atoms for query '{query}'",
            )
            if exp_count > 0:
                first = picked[0]
                self.assertEqual(first.path, exp_path)
                self.assertEqual(first.name, exp_name)
                self.assertAlmostEqual(first.position.x, exp_pos[0], places=3)
                self.assertAlmostEqual(first.position.y, exp_pos[1], places=3)
                self.assertAlmostEqual(first.position.z, exp_pos[2], places=3)

    def test_identity_and_mutation(self):
        """Criterion 2: Return values of pickup_atoms and get_family are the objects in the tree.

        Check identity with `is` and mutation of atom coordinates.
        """
        ag = rs.AtomGroup("root")
        chain = rs.AtomGroup("A")
        res = rs.AtomGroup("1")
        atm = rs.Atom(name="CA", xyz=[1.0, 2.0, 3.0])

        res.set_atom("1", atm)
        chain.set_group("1", res)
        ag.set_group("A", chain)

        chain_group = ag.get_group("A")
        res_group = chain_group.get_group("1")

        # 1. get_family returns the exact group object in the tree
        target_group = ag.get_family("/A/1/")
        self.assertIsNotNone(target_group)
        self.assertIs(target_group, res_group)

        # Identity from a sub-level
        self.assertIs(chain_group.get_family("/A/1/"), res_group)
        self.assertIs(res_group.get_family("/A/"), chain_group)
        self.assertIs(res_group.get_family("/"), ag)

        # 2. pickup_atoms returns the exact atom object in the tree
        picked = ag.pickup_atoms("CA")
        self.assertEqual(len(picked), 1)
        picked_atom = picked[0]
        tree_atom = res_group.get_atom("1")
        self.assertIs(picked_atom, tree_atom)

        # 3. Mutating the returned atom reflects in the tree
        picked_atom.position.x += 10.0
        picked_atom.position.y += 20.0
        picked_atom.position.z += 30.0

        self.assertAlmostEqual(tree_atom.position.x, 11.0)
        self.assertAlmostEqual(tree_atom.position.y, 22.0)
        self.assertAlmostEqual(tree_atom.position.z, 33.0)

    def test_get_family_navigation_and_edge_cases(self):
        """Criterion 3: Verify get_family navigation:

        - Descending to descendants
        - Moving to sibling groups via parent
        - Returning None for non-existent paths
        - Paths missing trailing '/'
        - Overlapping prefix keys like '/1/' and '/10/'
        """
        root = rs.AtomGroup("root")
        g1_init = rs.AtomGroup("grp1")
        g10_init = rs.AtomGroup("grp10")
        subA1 = rs.AtomGroup("subA")
        subB1 = rs.AtomGroup("subB")
        subA10 = rs.AtomGroup("subA")

        g1_init.set_group("subA", subA1)
        g1_init.set_group("subB", subB1)
        g10_init.set_group("subA", subA10)

        root.set_group("1", g1_init)
        root.set_group("10", g10_init)

        # Retrieve the groups from within the tree
        g1 = root.get_group("1")
        g10 = root.get_group("10")
        g1_subA = g1.get_group("subA")
        g1_subB = g1.get_group("subB")
        g10_subA = g10.get_group("subA")

        # Verify paths
        self.assertEqual(root.path, "/")
        self.assertEqual(g1.path, "/1/")
        self.assertEqual(g10.path, "/10/")
        self.assertEqual(g1_subA.path, "/1/subA/")
        self.assertEqual(g1_subB.path, "/1/subB/")
        self.assertEqual(g10_subA.path, "/10/subA/")

        # (a) Descending to descendants
        self.assertIs(root.get_family("/1/"), g1)
        self.assertIs(root.get_family("/10/"), g10)
        self.assertIs(root.get_family("/1/subA/"), g1_subA)
        self.assertIs(root.get_family("/1/subB/"), g1_subB)
        self.assertIs(root.get_family("/10/subA/"), g10_subA)

        # (b) Moving to sibling groups via parent
        self.assertIs(g1_subA.get_family("/1/subB/"), g1_subB)
        self.assertIs(g1_subB.get_family("/1/subA/"), g1_subA)
        self.assertIs(g1_subA.get_family("/10/subA/"), g10_subA)
        self.assertIs(g10_subA.get_family("/1/subA/"), g1_subA)
        self.assertIs(g1.get_family("/10/"), g10)
        self.assertIs(g10.get_family("/1/"), g1)

        # (c) Returning None for non-existent path
        self.assertIsNone(root.get_family("/non_existent/"))
        self.assertIsNone(root.get_family("/1/non_existent/"))
        self.assertIsNone(g1_subA.get_family("/1/non_existent/"))
        self.assertIsNone(g1_subA.get_family("/unknown_root/"))

        # (d) Missing trailing '/'
        self.assertIs(root.get_family("/1"), g1)
        self.assertIs(root.get_family("/10"), g10)
        self.assertIs(root.get_family("/1/subA"), g1_subA)
        self.assertIs(g1_subA.get_family("/1/subB"), g1_subB)
        self.assertIs(g1_subA.get_family("/10/subA"), g10_subA)

        # (e) Overlapping prefix keys: '/1/' vs '/10/'
        self.assertIs(root.get_family("/1/"), g1)
        self.assertIsNot(root.get_family("/1/"), g10)
        self.assertIs(root.get_family("/10/"), g10)
        self.assertIsNot(root.get_family("/10/"), g1)

        self.assertIs(root.get_family("/1"), g1)
        self.assertIsNot(root.get_family("/1"), g10)
        self.assertIs(root.get_family("/10"), g10)
        self.assertIsNot(root.get_family("/10"), g1)

        # Check self path match
        self.assertIs(root.get_family("/"), root)
        self.assertIs(g1.get_family("/1/"), g1)
        self.assertIs(g1.get_family("/1"), g1)

    def test_divide_path(self):
        """Target 4: divide_path(path) static method splits on '/' and filters empty elements."""
        cases = [
            ("/a/b/c/", ["a", "b", "c"]),
            ("///a///b//c///", ["a", "b", "c"]),
            ("", []),
            ("/", []),
            ("///", []),
            ("a", ["a"]),
            ("a/b", ["a", "b"]),
            ("/model/chain/res/atom", ["model", "chain", "res", "atom"]),
        ]
        for path_str, expected in cases:
            # Static call on class
            self.assertEqual(rs.AtomGroup.divide_path(path_str), expected)

        # Call on instance
        ag = rs.AtomGroup("test")
        self.assertEqual(ag.divide_path("/x/y/z/"), ["x", "y", "z"])


if __name__ == "__main__":
    unittest.main()
