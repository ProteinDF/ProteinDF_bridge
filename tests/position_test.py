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

class PositionTests(unittest.TestCase):
    def test_init(self):
        pos1 = bridge.Position()
        self.assertAlmostEqual(pos1.x, 0.0)
        self.assertAlmostEqual(pos1.y, 0.0)
        self.assertAlmostEqual(pos1.z, 0.0)

    def test_pickle(self):
        pos1 = bridge.Position([1.0, 2.0, 3.0])
        b = pickle.dumps(pos1)
        pos2 = pickle.loads(b)

        self.assertIsInstance(pos2, bridge.Position)
        self.assertAlmostEqual(pos2.x, 1.0)
        self.assertAlmostEqual(pos2.y, 2.0)
        self.assertAlmostEqual(pos2.z, 3.0)

def test_suite():
    """
    builds the test suite.
    """
    def _suite(test_class):
        return unittest.makeSuite(test_class)

    suite = unittest.TestSuite()
    suite.addTests(_suite(PositionTests))
    return suite

if __name__ == '__main__':
    unittest.main(defaultTest = 'test_suite')
