#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import pickle
import doctest

from proteindf_bridge.gro import SimpleGro

class GroTests(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_load(self):
        gro = SimpleGro()
        gro.load("./data/sample.gro")
        # print(gro)

        ag = gro.get_atomgroup()
        print(ag)


def load_tests(loader, tests, ignore):
    from proteindf_bridge import gro
    tests.addTests(doctest.DocTestSuite(gro))
    return tests

if __name__ == '__main__':
    unittest.main()


