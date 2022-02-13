#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import doctest

from proteindf_bridge.biopdb import Pdb
from proteindf_bridge.path import Path


class PathTest(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


def load_tests(loader, tests, ignore):
    from proteindf_bridge import path
    tests.addTests(doctest.DocTestSuite(path))
    return tests


if __name__ == '__main__':
    unittest.main()
