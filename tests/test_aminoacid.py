#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from proteindf_bridge.aminoacid import AminoAcid
from proteindf_bridge.atomgroup import AtomGroup


class TestAminoAcid(unittest.TestCase):
    def test_is_aminoacid(self):
        ag_ala = AtomGroup()
        ag_ala.name = "ALA"
        self.assertTrue(AminoAcid.is_aminoacid(ag_ala))

        ag_other = AtomGroup()
        ag_other.name = "HOH"
        self.assertFalse(AminoAcid.is_aminoacid(ag_other))


if __name__ == "__main__":
    unittest.main()
