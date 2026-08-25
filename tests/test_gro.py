import os
import unittest
import pickle
import doctest

from proteindf_bridge.gro import SimpleGro

DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "proteindf_bridge", "data"))

class GroTests(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_load(self):
        gro = SimpleGro()
        gro.load(os.path.join(DATA_DIR, "sample.gro"))
        ag = gro.get_atomgroup()
        self.assertIsNotNone(ag)


def load_tests(loader, tests, ignore):
    from proteindf_bridge import gro
    tests.addTests(doctest.DocTestSuite(gro))
    return tests

if __name__ == '__main__':
    unittest.main()


