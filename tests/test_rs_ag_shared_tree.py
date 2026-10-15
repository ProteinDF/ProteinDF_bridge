#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test suite validating rs.AtomGroup shared tree behavior (SPEC 11 Stage 2).
"""

import copy
import gc
import os
import unittest
import weakref

import proteindf_bridge.rs as rs


class TestRsAgSharedTree(unittest.TestCase):
    """Test suite for rs.AtomGroup reference tree behavior."""

    def _build_test_tree(self):
        root = rs.AtomGroup(name="root")
        model = rs.AtomGroup(name="model_1")
        chain = rs.AtomGroup(name="A")
        res = rs.AtomGroup(name="1")

        atom_ca = rs.Atom(symbol="C", name="CA", xyz=[1.0, 2.0, 3.0])
        atom_n = rs.Atom(symbol="N", name="N", xyz=[0.0, 0.0, 0.0])

        res.set_atom("CA", atom_ca)
        res.set_atom("N", atom_n)
        chain.set_group("1", res)
        model.set_group("A", chain)
        root.set_group("model_1", model)
        return root

    def test_shared_tree_mutations_persist_in_parent(self):
        """Criterion 1: ag.get_group('A').get_group('1').get_atom('CA').position.x = 9 persists."""
        ag = rs.AtomGroup(name="root")
        chain = rs.AtomGroup(name="A")
        res = rs.AtomGroup(name="1")
        ca = rs.Atom(symbol="C", name="CA", xyz=[1.0, 2.0, 3.0])
        res.set_atom("CA", ca)
        chain.set_group("1", res)
        ag.set_group("A", chain)

        # 1. Modify atom position through chained get_group
        ag.get_group("A").get_group("1").get_atom("CA").position.x = 9.0
        self.assertEqual(ag.get_group("A").get_group("1").get_atom("CA").position.x, 9.0)

        # 2. Modify group name in for-loop
        for k, g in ag.groups():
            g.name = "X"
        self.assertEqual(ag.get_group("A").name, "X")

        # 3. Identity checks
        self.assertIs(ag.get_group("A"), ag.get_group("A"))
        self.assertIs(ag.get_group("A").parent, ag)
        self.assertIs(ag.get_group("A").get_group("1").parent, ag.get_group("A"))

        # 4. set_group copies the inserted group, mutation on original does not affect tree
        external_res = rs.AtomGroup(name="2")
        ag.get_group("A").set_group("2", external_res)
        external_res.name = "MUTATED"
        self.assertEqual(ag.get_group("A").get_group("2").name, "2")

        # 5. Detached child has parent is None and can be used
        detached = ag.get_group("A").del_group("2")
        self.assertIsNotNone(detached)
        self.assertIsNone(detached.parent)
        self.assertFalse(ag.get_group("A").has_groupkey("2"))
        # Detached group remains usable
        self.assertEqual(detached.name, "2")
        self.assertEqual(detached.path, "/A/2/")

    def test_parent_path_semantics(self):
        """Path and depth semantics: parented group raises ValueError on set_path."""
        ag = rs.AtomGroup(name="root")
        chain = rs.AtomGroup(name="A")
        ag.set_group("A", chain)

        child = ag.get_group("A")
        self.assertEqual(child.path, "/A/")
        self.assertEqual(child.path_depth(), 1)

        # Setting path on a group with a parent raises ValueError
        with self.assertRaises(ValueError):
            child.path = "/custom/"

        # Setting path on root (no parent) is allowed
        ag.path = "/my_root"
        self.assertEqual(ag.path, "/my_root/")
        self.assertEqual(child.path, "/my_root/A/")

    def test_update_from_core_preserves_object_identity_and_reflects_results(self):
        """Criterion 2: Objects retrieved before update_from_core are identical and reflect results."""
        root = rs.AtomGroup(name="root")
        res = rs.AtomGroup(name="ALA_1")
        n = rs.Atom(symbol="N", name="N", xyz=[0.0, 0.0, 0.0])
        ca = rs.Atom(symbol="C", name="CA", xyz=[1.45, 0.0, 0.0])
        res.set_atom("N", n)
        res.set_atom("CA", ca)
        root.set_group("ALA_1", res)

        # Retrieve child before Bond.setup
        res_before = root.get_group("ALA_1")
        ca_before = res_before.get_atom("CA")

        # Execute Bond setup (which runs to_core -> setup_heuristic -> update_from_core)
        bond = rs.Bond()
        bond.setup(root)

        # Child and atom must be the exact same objects
        res_after = root.get_group("ALA_1")
        ca_after = res_after.get_atom("CA")
        self.assertIs(res_before, res_after)
        self.assertIs(ca_before, ca_after)

        # Results (bonds) are reflected
        bonds = root.get_bond_list()
        self.assertGreater(len(bonds), 0)

    def test_child_retains_parent_and_path_without_root_variable(self):
        """Review fix 1: Child retains parent and path even when root is not held in variable."""
        data_dir = os.path.join(
            os.path.dirname(__file__), "..", "rust", "crates", "proteindf-bridge", "tests", "data"
        )
        pdb_path = os.path.join(data_dir, "1hls.pdb")
        pdb = rs.Pdb()
        pdb.load(pdb_path)

        # Retrieve child group directly without storing root in variable
        chain_a = pdb.get_atomgroup().get_group("model_1").get_group("A")
        del pdb
        gc.collect()

        # Path of chain, parent, residue, and atom
        self.assertEqual(chain_a.path, "/model_1/A/")
        self.assertIsNotNone(chain_a.parent)
        self.assertEqual(chain_a.parent.path, "/model_1/")
        first_res = chain_a.get_group("1")
        self.assertEqual(first_res.path, "/model_1/A/1/")
        first_atom = first_res.get_atom("1_N")
        self.assertEqual(first_atom.path, "/model_1/A/1/1_N")

        # Detached child has parent is None, but retains previous path
        detached_res = chain_a.del_group("1")
        self.assertIsNotNone(detached_res)
        self.assertIsNone(detached_res.parent)
        self.assertEqual(detached_res.path, "/model_1/A/1/")
        self.assertEqual(detached_res.get_atom("1_N").path, "/model_1/A/1/1_N")

    def test_cyclic_reference_collected_by_gc(self):
        """Review fix 1 & Criterion 3: Tree with cyclic references is collected by gc."""
        def create_and_observe():
            root = rs.AtomGroup(name="root")
            child = rs.AtomGroup(name="child")
            root.set_group("child", child)
            return weakref.ref(root), weakref.ref(root.get_group("child"))

        w_root, w_child = create_and_observe()
        self.assertIsNotNone(w_root())
        self.assertIsNotNone(w_child())

        # Collect cyclic references
        gc.collect()

        self.assertIsNone(w_root())
        self.assertIsNone(w_child())


if __name__ == "__main__":
    unittest.main()
