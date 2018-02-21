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

from proteindf_bridge.atom import Atom

class AtomTests(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_init(self):
        atom = Atom()
        atom.symbol = 'Fe'
        self.assertEqual(atom.atomic_number, 26)

    def test_init2(self):
        atom = Atom(symbol='Ni', position=[0.0, 1.0, 2.0])
        self.assertEqual(atom.atomic_number, 28)
        self.assertAlmostEqual(atom.xyz[0], 0.0)
        self.assertAlmostEqual(atom.xyz[1], 1.0)
        self.assertAlmostEqual(atom.xyz[2], 2.0)

    def test_init3(self):
        atom = Atom(symbol='C', xyz="0.0 1.0 2.0")
        self.assertEqual(atom.atomic_number, 6)
        self.assertAlmostEqual(atom.xyz[0], 0.0)
        self.assertAlmostEqual(atom.xyz[1], 1.0)
        self.assertAlmostEqual(atom.xyz[2], 2.0)

    def test_pickle(self):
        atom1 = Atom()
        atom1.symbol = 'Na'
        b = pickle.dumps(atom1)
        atom2 = pickle.loads(b)

        self.assertIsInstance(atom2, Atom)
        self.assertEqual(atom2.symbol, 'Na')


if __name__ == '__main__':
    unittest.main()
