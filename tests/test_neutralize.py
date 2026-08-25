#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from proteindf_bridge.neutralize import Neutralize
from proteindf_bridge.atomgroup import AtomGroup
from proteindf_bridge.atom import Atom
from proteindf_bridge.position import Position


class TestNeutralize(unittest.TestCase):
    def test_exempt_list(self):
        protein = AtomGroup("protein")
        model = AtomGroup("model_1")
        chainA = AtomGroup("A")

        glu = AtomGroup("1")
        glu.name = "GLU"
        glu.set_atom("CD", Atom(name="CD", xyz=Position([0.0, 0.0, 0.0])))
        glu.set_atom("OE1", Atom(name="OE1", xyz=Position([0.0, 1.0, 0.0])))
        glu.set_atom("OE2", Atom(name="OE2", xyz=Position([1.0, 0.0, 0.0])))
        chainA.set_group("1", glu)

        lys = AtomGroup("2")
        lys.name = "LYS"
        lys.set_atom("NZ", Atom(name="NZ", xyz=Position([0.5, 2.0, 0.0])))
        lys.set_atom("HZ1", Atom(name="HZ1", xyz=Position([0.5, 2.5, 0.0])))
        lys.set_atom("HZ2", Atom(name="HZ2", xyz=Position([0.0, 2.0, 0.5])))
        lys.set_atom("HZ3", Atom(name="HZ3", xyz=Position([1.0, 2.0, 0.5])))
        chainA.set_group("2", lys)

        model.set_group("A", chainA)
        protein.set_group("model_1", model)

        neut = Neutralize(protein)
        exempt = neut._exempt_list(model)
        self.assertEqual(len(exempt), 2)
        self.assertIsNotNone(neut.neutralized)

    def test_non_destructive_input(self):
        protein = AtomGroup("protein")
        model = AtomGroup("model_1")
        chainA = AtomGroup("A")

        glu = AtomGroup("1")
        glu.name = "GLU"
        glu.set_atom("CD", Atom(name="CD", xyz=Position([0.0, 0.0, 0.0])))
        glu.set_atom("OE1", Atom(name="OE1", xyz=Position([0.0, 1.0, 0.0])))
        glu.set_atom("OE2", Atom(name="OE2", xyz=Position([1.0, 0.0, 0.0])))
        chainA.set_group("1", glu)
        model.set_group("A", chainA)
        protein.set_group("model_1", model)

        initial_atom_count = protein.get_number_of_all_atoms()
        neut = Neutralize(protein)
        # Verify original protein is not modified
        self.assertEqual(protein.get_number_of_all_atoms(), initial_atom_count)
        # Neutralized object has additional ion atom
        self.assertGreater(neut.neutralized.get_number_of_all_atoms(), initial_atom_count)


if __name__ == "__main__":
    unittest.main()
