#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import pickle
import doctest

from proteindf_bridge.biopdb import Pdb

class PdbTests(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_init(self):
        pdb = Pdb("./data/2MGO.pdb")

    def test_get_atomgroup(self):
        pdb = Pdb("./data/2MGO.pdb")
        ag = pdb.get_atomgroup()

        # model
        self.assertEqual(ag.get_number_of_groups(), 20)
        # chain
        self.assertEqual(ag["model_1"].get_number_of_groups(), 1)
        # residues
        self.assertEqual(ag["model_1"]["A"].get_number_of_groups(), 9)

        self.assertEqual(ag["model_1"]["A"][1].get_number_of_atoms(), 12) # CYS1
        self.assertEqual(ag["model_1"]["A"][1].name, "CYS")
        self.assertEqual(ag["model_1"]["A"][2].get_number_of_atoms(), 21) # TYR2
        self.assertEqual(ag["model_1"]["A"][2].name, "TYR")
        self.assertEqual(ag["model_1"]["A"][3].get_number_of_atoms(), 19) # ILE3
        self.assertEqual(ag["model_1"]["A"][3].name, "ILE")
        self.assertEqual(ag["model_1"]["A"][4].get_number_of_atoms(), 17) # GLN4
        self.assertEqual(ag["model_1"]["A"][4].name, "GLN")
        self.assertEqual(ag["model_1"]["A"][5].get_number_of_atoms(), 14) # ASN5
        self.assertEqual(ag["model_1"]["A"][5].name, "ASN")
        self.assertEqual(ag["model_1"]["A"][6].get_number_of_atoms(), 10) # CYS6
        self.assertEqual(ag["model_1"]["A"][6].name, "CYS")
        self.assertEqual(ag["model_1"]["A"][7].get_number_of_atoms(), 14) # PRO7
        self.assertEqual(ag["model_1"]["A"][7].name, "PRO")
        self.assertEqual(ag["model_1"]["A"][8].get_number_of_atoms(), 19) # LEU8
        self.assertEqual(ag["model_1"]["A"][8].name, "LEU")
        self.assertEqual(ag["model_1"]["A"][9].get_number_of_atoms(),  8) # GLY9
        self.assertEqual(ag["model_1"]["A"][9].name, "GLY")

    def test_set_by_atomgroup(self):
        pass


def load_tests(loader, tests, ignore):
    from proteindf_bridge import biopdb
    tests.addTests(doctest.DocTestSuite(biopdb))
    return tests


if __name__ == '__main__':
    unittest.main()
