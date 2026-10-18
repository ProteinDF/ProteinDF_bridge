#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for in-place scalar multiplication on AtomGroup."""

import unittest
import numpy as np

import proteindf_bridge as b


class TestAtomGroupInplaceMultiplication(unittest.TestCase):
    """Test suite for in-place scalar multiplication on AtomGroup."""

    def test_reproduction_basic(self):
        """Verify basic scalar multiplication scales coordinates."""
        ag = b.AtomGroup()
        ag.set_atom("h", b.Atom(symbol="H", position=b.Position(1.0, 2.0, 3.0)))
        ag *= 2.0
        pos = ag.get_atom("h").position
        self.assertAlmostEqual(pos.x, 2.0)
        self.assertAlmostEqual(pos.y, 4.0)
        self.assertAlmostEqual(pos.z, 6.0)

    def test_returns_same_object(self):
        """Verify that ag *= scalar returns the exact same object instance."""
        ag = b.AtomGroup()
        ag.set_atom("c", b.Atom(symbol="C", position=b.Position(1.0, 2.0, 3.0)))
        initial = ag
        ag *= 1.5
        self.assertIs(ag, initial)

    def test_nested_tree_coordinates(self):
        """Verify that coordinates of all atoms in a nested hierarchy are multiplied."""
        root = b.AtomGroup(name="root")
        model = b.AtomGroup(name="model")
        chain = b.AtomGroup(name="chain")
        res = b.AtomGroup(name="res")

        a1 = b.Atom(symbol="C", name="C1", position=b.Position(1.0, 2.0, 3.0))
        a2 = b.Atom(symbol="N", name="N1", position=b.Position(-1.0, 0.5, 2.5))
        a_root = b.Atom(symbol="O", name="O1", position=b.Position(0.0, 1.0, 4.0))

        res.set_atom("C1", a1)
        res.set_atom("N1", a2)
        chain.set_group("res", res)
        model.set_group("chain", chain)
        root.set_group("model", model)
        root.set_atom("O1", a_root)

        root *= 2.0

        c1_pos = root.get_group("model").get_group("chain").get_group("res").get_atom("C1").position
        self.assertAlmostEqual(c1_pos.x, 2.0)
        self.assertAlmostEqual(c1_pos.y, 4.0)
        self.assertAlmostEqual(c1_pos.z, 6.0)

        n1_pos = root.get_group("model").get_group("chain").get_group("res").get_atom("N1").position
        self.assertAlmostEqual(n1_pos.x, -2.0)
        self.assertAlmostEqual(n1_pos.y, 1.0)
        self.assertAlmostEqual(n1_pos.z, 5.0)

        o1_pos = root.get_atom("O1").position
        self.assertAlmostEqual(o1_pos.x, 0.0)
        self.assertAlmostEqual(o1_pos.y, 2.0)
        self.assertAlmostEqual(o1_pos.z, 8.0)

    def test_scalar_types(self):
        """Verify multiplication accepts int, float, and numpy scalar types."""
        # int
        ag_int = b.AtomGroup()
        ag_int.set_atom("c", b.Atom(symbol="C", position=b.Position(1.0, 2.0, 3.0)))
        ag_int *= 3
        pos_int = ag_int.get_atom("c").position
        self.assertAlmostEqual(pos_int.x, 3.0)
        self.assertAlmostEqual(pos_int.y, 6.0)
        self.assertAlmostEqual(pos_int.z, 9.0)

        # numpy.float64
        ag_np = b.AtomGroup()
        ag_np.set_atom("c", b.Atom(symbol="C", position=b.Position(1.0, 2.0, 3.0)))
        ag_np *= np.float64(0.5)
        pos_np = ag_np.get_atom("c").position
        self.assertAlmostEqual(pos_np.x, 0.5)
        self.assertAlmostEqual(pos_np.y, 1.0)
        self.assertAlmostEqual(pos_np.z, 1.5)

    def test_reflected_in_pre_retrieved_references(self):
        """Verify that mutations are reflected in previously extracted child groups, atoms, and positions."""
        root = b.AtomGroup(name="root")
        child = b.AtomGroup(name="child")
        atom = b.Atom(symbol="C", position=b.Position(2.0, 4.0, 6.0))
        child.set_atom("c", atom)
        root.set_group("child", child)

        # Retrieve references beforehand
        retrieved_child = root.get_group("child")
        retrieved_atom = retrieved_child.get_atom("c")
        retrieved_pos = retrieved_atom.position

        root *= 0.5

        # Reflected in retrieved references
        self.assertAlmostEqual(retrieved_pos.x, 1.0)
        self.assertAlmostEqual(retrieved_pos.y, 2.0)
        self.assertAlmostEqual(retrieved_pos.z, 3.0)

        self.assertAlmostEqual(retrieved_atom.position.x, 1.0)
        self.assertAlmostEqual(retrieved_atom.position.y, 2.0)
        self.assertAlmostEqual(retrieved_atom.position.z, 3.0)

        self.assertAlmostEqual(retrieved_child.get_atom("c").position.x, 1.0)
        self.assertAlmostEqual(retrieved_child.get_atom("c").position.y, 2.0)
        self.assertAlmostEqual(retrieved_child.get_atom("c").position.z, 3.0)

    def test_force_is_unchanged(self):
        """Verify that force vector remains unchanged when multiplying AtomGroup."""
        ag = b.AtomGroup()
        atom = b.Atom(
            symbol="H",
            position=b.Position(1.0, 2.0, 3.0),
            force=b.Position(10.0, 20.0, 30.0),
        )
        ag.set_atom("h", atom)

        ag *= 2.0

        atm = ag.get_atom("h")
        self.assertAlmostEqual(atm.position.x, 2.0)
        self.assertAlmostEqual(atm.position.y, 4.0)
        self.assertAlmostEqual(atm.position.z, 6.0)

        self.assertAlmostEqual(atm.force.x, 10.0)
        self.assertAlmostEqual(atm.force.y, 20.0)
        self.assertAlmostEqual(atm.force.z, 30.0)

    def test_invalid_type_raises_type_error(self):
        """Verify that multiplying by a non-numeric type raises TypeError."""
        ag = b.AtomGroup()
        ag.set_atom("h", b.Atom(symbol="H", position=b.Position(1.0, 2.0, 3.0)))

        with self.assertRaises(TypeError):
            ag *= "a"

        with self.assertRaises(TypeError):
            ag *= [1.0, 2.0]


if __name__ == "__main__":
    unittest.main()
