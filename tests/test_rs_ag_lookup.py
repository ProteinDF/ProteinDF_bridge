#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test suite validating rs.AtomGroup lookup and selection methods (TASK_rs-ag-lookup).
"""

import glob
import os
import unittest

from proteindf_bridge.atom import Atom as PyAtom
from proteindf_bridge.atomgroup import AtomGroup as PyAtomGroup
from proteindf_bridge.biopdb import Pdb as PyPdb
from proteindf_bridge.functions import load_atomgroup as py_load_atomgroup
from proteindf_bridge.position import Position as PyPosition
import proteindf_bridge.rs as rs


class TestRsAgLookup(unittest.TestCase):
    """Test suite for rs.AtomGroup lookup, selection, and hierarchy navigation."""

    def setUp(self):
        self.data_dir = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "proteindf_bridge", "data")
        )

    def _compare_group_lookups(self, py_ag, rs_ag):
        """Helper to recursively compare get_atom_keys, has_atomname, and has_groupname."""
        # 1. get_atom_keys
        py_keys = py_ag.get_atom_keys()
        rs_keys = rs_ag.get_atom_keys()
        self.assertEqual(rs_keys, py_keys)

        # 2. has_atomname
        # Test existing atom names
        for _, py_atm in py_ag.atoms():
            atm_name = py_atm.name
            self.assertTrue(rs_ag.has_atomname(atm_name))
            # Test with leading and trailing whitespaces
            self.assertTrue(rs_ag.has_atomname(f"  {atm_name}  "))
            self.assertEqual(rs_ag.has_atomname(atm_name), py_ag.has_atomname(atm_name))
            self.assertEqual(
                rs_ag.has_atomname(f"  {atm_name}  "),
                py_ag.has_atomname(f"  {atm_name}  "),
            )

        # Test non-existent atom name
        self.assertFalse(rs_ag.has_atomname("__NON_EXISTENT_ATOM__"))
        self.assertEqual(
            rs_ag.has_atomname("__NON_EXISTENT_ATOM__"),
            py_ag.has_atomname("__NON_EXISTENT_ATOM__"),
        )

        # 3. has_groupname
        for _, py_sub in py_ag.groups():
            sub_name = py_sub.name
            self.assertTrue(rs_ag.has_groupname(sub_name))
            self.assertEqual(rs_ag.has_groupname(sub_name), py_ag.has_groupname(sub_name))

        self.assertFalse(rs_ag.has_groupname("__NON_EXISTENT_GROUP__"))
        self.assertEqual(
            rs_ag.has_groupname("__NON_EXISTENT_GROUP__"),
            py_ag.has_groupname("__NON_EXISTENT_GROUP__"),
        )

        # Recurse for child groups
        for (k_py, sub_py), (k_rs, sub_rs) in zip(py_ag.groups(), rs_ag.groups()):
            self.assertEqual(k_py, k_rs)
            self._compare_group_lookups(sub_py, sub_rs)

    def _compare_pickup_atoms(self, py_ag, rs_ag, query_keys_or_names):
        """Helper to compare pickup_atoms results by atom path and name."""
        for query in query_keys_or_names:
            py_picked = py_ag.pickup_atoms(query)
            rs_picked = rs_ag.pickup_atoms(query)

            self.assertEqual(
                len(rs_picked),
                len(py_picked),
                f"Mismatch in count of picked atoms for query '{query}'",
            )
            for py_atm, rs_atm in zip(py_picked, rs_picked):
                self.assertEqual(rs_atm.path, py_atm.path)
                self.assertEqual(rs_atm.name, py_atm.name)

    def test_brd_files_lookup(self):
        """Criterion 1: Verify get_atom_keys, has_atomname, has_groupname, and pickup_atoms

        against pure-Python for all .brd files in proteindf_bridge/data/.
        """
        brd_files = sorted(glob.glob(os.path.join(self.data_dir, "*.brd")))
        self.assertGreater(len(brd_files), 0, "No .brd files found in data directory")

        for brd_path in brd_files:
            with self.subTest(file=os.path.basename(brd_path)):
                py_ag = py_load_atomgroup(brd_path)
                rs_ag = rs.load_atomgroup(brd_path)

                self._compare_group_lookups(py_ag, rs_ag)

                # Pickup atoms test
                sample_queries = ["CA", "C", "N", "O", "H", "1", "2", "non_existent"]
                all_py_atoms = py_ag.get_atom_list()
                if all_py_atoms:
                    sample_queries.append(all_py_atoms[0].name)
                self._compare_pickup_atoms(py_ag, rs_ag, sample_queries)

    def test_2mgo_pdb_lookup(self):
        """Criterion 1: Verify get_atom_keys, has_atomname, has_groupname, and pickup_atoms

        against pure-Python for 2MGO.pdb.
        """
        pdb_path = os.path.join(self.data_dir, "2MGO.pdb")
        self.assertTrue(os.path.exists(pdb_path), "2MGO.pdb fixture not found")

        py_ag = PyPdb(pdb_path).get_atomgroup()
        rs_ag = rs.Pdb(pdb_path).get_atomgroup()

        self._compare_group_lookups(py_ag, rs_ag)

        sample_queries = [
            "CA",
            "CB",
            "N",
            "C",
            "O",
            "SG",
            "1",
            "10",
            "100",
            "non_existent",
        ]
        self._compare_pickup_atoms(py_ag, rs_ag, sample_queries)

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
            self.assertEqual(PyAtomGroup.divide_path(path_str), expected)
            self.assertEqual(
                rs.AtomGroup.divide_path(path_str),
                PyAtomGroup.divide_path(path_str),
            )

        # Call on instance
        ag = rs.AtomGroup("test")
        py_ag = PyAtomGroup("test")
        self.assertEqual(ag.divide_path("/x/y/z/"), ["x", "y", "z"])
        self.assertEqual(ag.divide_path("/x/y/z/"), py_ag.divide_path("/x/y/z/"))


if __name__ == "__main__":
    unittest.main()
