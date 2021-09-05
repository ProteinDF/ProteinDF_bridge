#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import doctest

from proteindf_bridge.utils import Utils
from proteindf_bridge.biopdb import Pdb


class UtilsTest(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_get_sequential_residue_id(self):
        pdb = Pdb("./data/3i3zH.pdb")
        models = pdb.get_atomgroup()
        protein = models["model_1"]

        self.assertEqual(Utils.get_sequential_residue_id(protein, "A", "10"), 10)
        self.assertEqual(Utils.get_sequential_residue_id(protein, "A", "21"), 21)
        self.assertEqual(Utils.get_sequential_residue_id(protein, "B", "1"), 22)
        self.assertEqual(Utils.get_sequential_residue_id(protein, "B", "10"), 31)
        self.assertEqual(Utils.get_sequential_residue_id(protein, "C", "10"), 0)


# def load_tests(loader, tests, ignore):
#     from proteindf_bridge import Utils
#     tests.addTests(doctest.DocTestSuite(Utils))
#     return tests


if __name__ == '__main__':
    unittest.main()
