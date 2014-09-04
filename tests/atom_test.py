#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2002-2014 The ProteinDF project
# see also AUTHORS and README.
# 
# This file is part of ProteinDF.
# 
# ProteinDF is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# ProteinDF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

import unittest
import pickle
import bridge

class AtomTests(unittest.TestCase):
    def test_init(self):
        atom = bridge.Atom()
        atom.symbol = 'Fe'
        self.assertEqual(atom.atomic_number, 26)

    def test_pickle(self):
        atom1 = bridge.Atom()
        atom1.symbol = 'Na'
        b = pickle.dumps(atom1)
        atom2 = pickle.loads(b)

        self.assertIsInstance(atom2, bridge.Atom)
        self.assertEqual(atom2.symbol, 'Na')

def test_suite():
    """
    builds the test suite.
    """
    def _suite(test_class):
        return unittest.makeSuite(test_class)

    suite = unittest.TestSuite()
    suite.addTests(_suite(AtomTests))
    return suite

if __name__ == '__main__':
    unittest.main(defaultTest = 'test_suite')
