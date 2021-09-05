#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import doctest

from proteindf_bridge.biopdb import Pdb
from proteindf_bridge.atomgroup import AtomGroup
from proteindf_bridge.format import Format


class FormatTest(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_formats1(self):
        pdb = Pdb("./data/2MGO.pdb")

        models = pdb.get_atomgroup()
        self.assertTrue(Format.is_models(models))
        self.assertFalse(Format.is_protein(models))
        self.assertFalse(Format.is_chain(models))
        self.assertFalse(Format.is_residue(models))

        model = models["model_1"]
        self.assertFalse(Format.is_models(model))
        self.assertTrue(Format.is_protein(model))
        self.assertFalse(Format.is_chain(model))
        self.assertFalse(Format.is_residue(model))

        chain = model["A"]
        self.assertFalse(Format.is_models(chain))
        self.assertFalse(Format.is_protein(chain))
        self.assertTrue(Format.is_chain(chain))
        self.assertFalse(Format.is_residue(chain))

        res = chain["3"]
        self.assertFalse(Format.is_models(res))
        self.assertFalse(Format.is_protein(res))
        self.assertFalse(Format.is_chain(res))
        self.assertTrue(Format.is_residue(res))


def load_tests(loader, tests, ignore):
    from proteindf_bridge import format
    tests.addTests(doctest.DocTestSuite(format))
    return tests


if __name__ == '__main__':
    unittest.main()
