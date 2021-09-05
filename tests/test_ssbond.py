#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import doctest

from proteindf_bridge.ssbond import SSBond


class SSBondTest(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    # def test_check(self):
    #     tmp_pdb = Pdb('./data/1hls.pdb')
    #     models = tmp_pdb.get_atomgroup()
    #     model = models.get_group('model_1')
    #     ssbond = SSBond()
    #     ssbond.check(models)


def load_tests(loader, tests, ignore):
    from proteindf_bridge import ssbond
    tests.addTests(doctest.DocTestSuite(ssbond))
    return tests


if __name__ == '__main__':
    unittest.main()
