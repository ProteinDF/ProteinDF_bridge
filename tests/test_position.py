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

from proteindf_bridge.position import Position

class PositionTests(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_init1(self):
        pos1 = Position()
        self.assertAlmostEqual(pos1.x, 0.0)
        self.assertAlmostEqual(pos1.y, 0.0)
        self.assertAlmostEqual(pos1.z, 0.0)

    def test_init2(self):
        pos1 = Position([1.0, 2.0, -3.0])
        self.assertAlmostEqual(pos1.x, 1.0)
        self.assertAlmostEqual(pos1.y, 2.0)
        self.assertAlmostEqual(pos1.z, -3.0)

    def test_init3(self):
        pos1 = Position("1.0  2.0 -3.0")
        self.assertAlmostEqual(pos1.x, 1.0)
        self.assertAlmostEqual(pos1.y, 2.0)
        self.assertAlmostEqual(pos1.z, -3.0)

        pos2 = Position("1.0,  2.0,-3.0")
        self.assertAlmostEqual(pos2.x, 1.0)
        self.assertAlmostEqual(pos2.y, 2.0)
        self.assertAlmostEqual(pos2.z, -3.0)

    def test_pickle(self):
        pos1 = Position([1.0, 2.0, 3.0])
        b = pickle.dumps(pos1)
        pos2 = pickle.loads(b)

        self.assertIsInstance(pos2, Position)
        self.assertAlmostEqual(pos2.x, 1.0)
        self.assertAlmostEqual(pos2.y, 2.0)
        self.assertAlmostEqual(pos2.z, 3.0)


if __name__ == '__main__':
    unittest.main()
