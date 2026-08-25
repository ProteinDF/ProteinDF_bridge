#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from proteindf_bridge.atomgroup import AtomGroup
from proteindf_bridge.atom import Atom
from proteindf_bridge.position import Position
from proteindf_bridge.superposer import Superposer


class TestSuperposer(unittest.TestCase):
    def setUp(self):
        self.ag1 = AtomGroup("mol1")
        self.ag1.set_atom("A1", Atom(name="A1", xyz=Position([1.0, 1.0, 1.0])))
        self.ag1.set_atom("A2", Atom(name="A2", xyz=Position([1.0, -1.0, -1.0])))
        self.ag1.set_atom("A3", Atom(name="A3", xyz=Position([-1.0, 1.0, -1.0])))
        self.ag1.set_atom("A4", Atom(name="A4", xyz=Position([-1.0, -1.0, 1.0])))

        # 平行移動したもの (x+1.0, y+2.0, z+3.0)
        self.ag2 = AtomGroup("mol2")
        self.ag2.set_atom("A1", Atom(name="A1", xyz=Position([2.0, 3.0, 4.0])))
        self.ag2.set_atom("A2", Atom(name="A2", xyz=Position([2.0, 1.0, 2.0])))
        self.ag2.set_atom("A3", Atom(name="A3", xyz=Position([0.0, 3.0, 2.0])))
        self.ag2.set_atom("A4", Atom(name="A4", xyz=Position([0.0, 1.0, 4.0])))

    def test_rmsd(self):
        sp = Superposer(self.ag1, self.ag2)
        self.assertAlmostEqual(sp.rmsd, 0.0, places=5)

    def test_superimpose(self):
        sp = Superposer(self.ag1, self.ag2)
        superimposed = sp.superimpose(self.ag1)
        self.assertEqual(superimposed.get_number_of_atoms(), 4)
        for key in ["A1", "A2", "A3", "A4"]:
            p1 = superimposed.get_atom(key).xyz
            p2 = self.ag2.get_atom(key).xyz
            dist = p1.distance_from(p2)
            self.assertAlmostEqual(dist, 0.0, places=4)


if __name__ == "__main__":
    unittest.main()
