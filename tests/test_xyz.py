#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import tempfile
import unittest
from proteindf_bridge.xyz import Xyz
from proteindf_bridge.atomgroup import AtomGroup
from proteindf_bridge.atom import Atom
from proteindf_bridge.position import Position


class TestXyz(unittest.TestCase):
    def test_init_with_atomgroup(self):
        ag = AtomGroup("water")
        ag.set_atom("1", Atom(symbol="O", xyz=Position([0.0, 0.0, 0.0])))
        ag.set_atom("2", Atom(symbol="H", xyz=Position([0.0, 1.0, 0.0])))
        ag.set_atom("3", Atom(symbol="H", xyz=Position([1.0, 0.0, 0.0])))

        xyz = Xyz(ag)
        text = xyz.get_text()
        self.assertTrue(text.startswith("3\n"))
        self.assertIn("O", text)
        self.assertIn("H", text)

    def test_save_and_load(self):
        ag = AtomGroup("test_mol")
        ag.set_atom("1", Atom(symbol="C", xyz=Position([1.2, 3.4, 5.6])))

        with tempfile.NamedTemporaryFile(suffix=".xyz", delete=False) as f:
            temp_path = f.name

        try:
            xyz1 = Xyz(ag)
            xyz1.save(temp_path)

            # ファイルパス文字列を渡して初期化する経路をテスト
            xyz2 = Xyz(temp_path)
            loaded_ag = xyz2.get_atom_group()
            self.assertEqual(loaded_ag.get_number_of_atoms(), 1)
            atom = loaded_ag.get_atom("0")
            self.assertEqual(atom.symbol, "C")
            self.assertAlmostEqual(atom.xyz.x, 1.2, places=4)
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)


if __name__ == "__main__":
    unittest.main()
