#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import tempfile
import unittest
from proteindf_bridge.mol2 import SimpleMol2
from proteindf_bridge.atomgroup import AtomGroup
from proteindf_bridge.atom import Atom
from proteindf_bridge.position import Position


class TestSimpleMol2(unittest.TestCase):
    def test_save_and_format(self):
        ag = AtomGroup()
        ag.name = "water"
        ag.set_atom("1", Atom(name="O1", symbol="O", xyz=Position([0.0, 0.0, 0.0]), charge=-0.8))
        ag.set_atom("2", Atom(name="H1", symbol="H", xyz=Position([0.0, 1.0, 0.0]), charge=0.4))
        ag.set_atom("3", Atom(name="H2", symbol="H", xyz=Position([1.0, 0.0, 0.0]), charge=0.4))

        mol2 = SimpleMol2(ag)
        text = str(mol2)

        self.assertIn("@<TRIPOS>MOLECULE", text)
        self.assertIn("water", text)
        self.assertIn("@<TRIPOS>ATOM", text)
        self.assertIn("O1", text)
        self.assertIn("H1", text)

        with tempfile.NamedTemporaryFile(suffix=".mol2", delete=False) as f:
            temp_path = f.name

        try:
            mol2.save(temp_path)
            self.assertTrue(os.path.exists(temp_path))
            with open(temp_path, "r") as rf:
                content = rf.read()
                self.assertIn("@<TRIPOS>MOLECULE", content)
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)


if __name__ == "__main__":
    unittest.main()
