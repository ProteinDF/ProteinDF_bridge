#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test suite validating rs.Atom shared position / force references,
and copy / deepcopy / pickle roundtrip for rs.Position and rs.Atom.
"""

import copy
import pickle
import unittest
import proteindf_bridge.rs as rs


class TestRsAtomSharedPosition(unittest.TestCase):
    """Test suite for rs.Atom shared position references and copy/pickle support."""

    def test_atom_position_and_force_mutations_persist(self):
        """Verify that mutations to a.position and a.force persist in the atom."""
        a = rs.Atom(symbol="C", name="CA", xyz=[1.0, 2.0, 3.0], force=[0.0, 0.0, 0.0])

        # 1. a.position is a.position, a.xyz is a.xyz, and a.position is a.xyz
        self.assertIs(a.position, a.position)
        self.assertIs(a.xyz, a.xyz)
        self.assertIs(a.position, a.xyz)
        self.assertIs(a.force, a.force)

        # 2. a.position.x = 9 persists in atom
        a.position.x = 9.0
        self.assertEqual(a.position.x, 9.0)
        self.assertEqual(a.xyz.x, 9.0)

        # 3. a.force.y = 1 persists in atom
        a.force.y = 1.0
        self.assertEqual(a.force.y, 1.0)

        # 4. p = Position(1,2,3); a.position = p; p.x = 5 leaves a.position.x == 1
        p = rs.Position(1.0, 2.0, 3.0)
        a.position = p
        self.assertIsNot(a.position, p)
        p.x = 5.0
        self.assertEqual(a.position.x, 1.0)

        # Same for force
        f = rs.Position(4.0, 5.0, 6.0)
        a.force = f
        self.assertIsNot(a.force, f)
        f.y = 9.0
        self.assertEqual(a.force.y, 5.0)

    def test_position_copy_deepcopy_pickle(self):
        """Verify copy, deepcopy, and pickle roundtrip for rs.Position."""
        p = rs.Position(1.5, 2.5, 3.5)

        for copy_func, desc in [
            (copy.copy, "copy"),
            (copy.deepcopy, "deepcopy"),
            (lambda obj: pickle.loads(pickle.dumps(obj)), "pickle"),
        ]:
            with self.subTest(mode=desc):
                p2 = copy_func(p)
                self.assertIsNot(p, p2)
                self.assertEqual((p2.x, p2.y, p2.z), (1.5, 2.5, 3.5))
                p2.x = 99.0
                self.assertEqual(p.x, 1.5)

    def test_atom_copy_deepcopy_pickle(self):
        """Verify copy, deepcopy, and pickle roundtrip for rs.Atom."""
        a = rs.Atom(
            symbol="N",
            name="N1",
            label="LABEL1",
            charge=-0.5,
            path="/A/1/N",
            position=[1.0, 2.0, 3.0],
            force=[0.1, 0.2, 0.3],
        )

        for copy_func, desc in [
            (copy.copy, "copy"),
            (copy.deepcopy, "deepcopy"),
            (lambda obj: pickle.loads(pickle.dumps(obj)), "pickle"),
        ]:
            with self.subTest(mode=desc):
                a2 = copy_func(a)
                self.assertIsNot(a, a2)
                self.assertIsNot(a.position, a2.position)
                self.assertIsNot(a.force, a2.force)

                # Check attributes match
                self.assertEqual(a2.symbol, "N")
                self.assertEqual(a2.atomic_number, 7)
                self.assertEqual(a2.name, "N1")
                self.assertEqual(a2.label, "LABEL1")
                self.assertAlmostEqual(a2.charge, -0.5)
                self.assertEqual(a2.path, "/A/1/N")
                self.assertEqual((a2.position.x, a2.position.y, a2.position.z), (1.0, 2.0, 3.0))
                self.assertEqual((a2.force.x, a2.force.y, a2.force.z), (0.1, 0.2, 0.3))

                # Verify coordinate and force independence
                a2.position.x = 42.0
                self.assertEqual(a.position.x, 1.0)
                self.assertEqual(a2.position.x, 42.0)

                a2.force.y = 84.0
                self.assertEqual(a.force.y, 0.2)
                self.assertEqual(a2.force.y, 84.0)


if __name__ == "__main__":
    unittest.main()
