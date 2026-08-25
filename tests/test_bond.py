#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from proteindf_bridge.bond import Bond
from proteindf_bridge.atomgroup import AtomGroup
from proteindf_bridge.atom import Atom
from proteindf_bridge.position import Position


class TestBond(unittest.TestCase):
    def test_bond_setup(self):
        ag = AtomGroup("mol")
        # 共有結合距離（〜1.5Å）にある2つの炭素原子
        c1 = Atom(symbol="C", xyz=Position([0.0, 0.0, 0.0]))
        c2 = Atom(symbol="C", xyz=Position([1.5, 0.0, 0.0]))
        # 遠く離れた炭素原子（10.0Å）
        c3 = Atom(symbol="C", xyz=Position([10.0, 0.0, 0.0]))

        ag.set_atom("1", c1)
        ag.set_atom("2", c2)
        ag.set_atom("3", c3)

        bond = Bond()
        bond.setup(ag)

        self.assertIsNotNone(bond._distmat)
        self.assertIsNotNone(bond._bondmat)
        # c1 と c2 の距離は約 1.5
        self.assertAlmostEqual(bond._distmat.get(1, 0), 1.5)
        # c1 と c2 は結合あり(1)、c1/c2 と c3 は結合なし(0)
        self.assertEqual(bond._bondmat.get(1, 0), 1)
        self.assertEqual(bond._bondmat.get(2, 0), 0)
        self.assertEqual(bond._bondmat.get(2, 1), 0)


if __name__ == "__main__":
    unittest.main()
