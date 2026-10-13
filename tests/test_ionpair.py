#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from proteindf_bridge.ionpair import IonPair
from proteindf_bridge.atomgroup import AtomGroup
from proteindf_bridge.atom import Atom
from proteindf_bridge.position import Position
from proteindf_bridge.error import BrValueError


class TestIonPair(unittest.TestCase):
    def test_ion_pair_detection(self):
        protein = AtomGroup("protein")
        chainA = AtomGroup("A")

        # GLU (陰イオン性残基: CD, OE1, OE2)
        glu = AtomGroup("1")
        glu.name = "GLU"
        glu.set_atom("CD", Atom(name="CD", xyz=Position([0.0, 0.0, 0.0])))
        glu.set_atom("OE1", Atom(name="OE1", xyz=Position([0.0, 1.0, 0.0])))
        glu.set_atom("OE2", Atom(name="OE2", xyz=Position([1.0, 0.0, 0.0])))
        chainA.set_group("1", glu)

        # LYS (陽イオン性残基: NZ) - 距離 < 4.0Å
        lys = AtomGroup("2")
        lys.name = "LYS"
        lys.set_atom("NZ", Atom(name="NZ", xyz=Position([0.5, 2.0, 0.0])))
        chainA.set_group("2", lys)

        protein.set_group("A", chainA)

        ip = IonPair(protein)
        pairs = ip.get_ion_pairs()

        self.assertEqual(len(pairs), 1)
        self.assertEqual(pairs[0][2], "GLU")
        self.assertEqual(pairs[0][3], "LYS")

    def test_lys_missing_nz_raises_br_value_error(self):
        lys = AtomGroup("2")
        lys.name = "LYS"
        ip = IonPair(AtomGroup())
        with self.assertRaises(BrValueError) as cm:
            ip._get_center_LYS(lys)
        self.assertIn("NZ", str(cm.exception))


if __name__ == "__main__":
    unittest.main()
