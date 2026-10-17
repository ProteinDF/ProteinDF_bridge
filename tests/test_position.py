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
import doctest

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

    def test_norm_zero_vector(self):
        from proteindf_bridge.error import BrValueError
        pos = Position([0.0, 0.0, 0.0])
        with self.assertRaises(BrValueError):
            pos.norm()

    def test_setstate(self):
        pos = Position()
        pos.__setstate__([4.0, 5.0, 6.0])
        self.assertEqual(pos.xyz, [4.0, 5.0, 6.0])

        pos.__setstate__((7.0, 8.0, 9.0))
        self.assertEqual(pos.xyz, [7.0, 8.0, 9.0])

        with self.assertRaises(TypeError):
            pos.__setstate__({1.0, 2.0, 3.0})

    def test_indexing_and_iteration(self):
        p = Position([1.0, 2.0, 3.0])
        self.assertEqual(list(p), [1.0, 2.0, 3.0])
        self.assertEqual(len(p), 3)
        self.assertEqual(p[-1], 3.0)
        self.assertEqual(p[-2], 2.0)
        self.assertEqual(p[-3], 1.0)
        self.assertEqual(p[0], 1.0)
        self.assertEqual(p[1], 2.0)
        self.assertEqual(p[2], 3.0)

        with self.assertRaises(IndexError):
            _ = p[3]
        with self.assertRaises(IndexError):
            _ = p[-4]

        with self.assertRaises(IndexError):
            p[3] = 0.0
        with self.assertRaises(IndexError):
            p[-4] = 0.0

        p[-1] = 9.0
        self.assertEqual(p[2], 9.0)
        self.assertEqual(p[-1], 9.0)

    def test_multiplication_scalars_and_dot(self):
        import numpy as np
        from proteindf_bridge.error import BrInputError

        p = Position([1.0, 2.0, 3.0])

        r1 = p * 2
        self.assertEqual(list(r1), [2.0, 4.0, 6.0])
        r2 = 2 * p
        self.assertEqual(list(r2), [2.0, 4.0, 6.0])

        r3 = p * np.int64(2)
        self.assertEqual(list(r3), [2.0, 4.0, 6.0])
        r4 = np.int64(2) * p
        self.assertEqual(list(r4), [2.0, 4.0, 6.0])

        r5 = p * np.float64(2)
        self.assertEqual(list(r5), [2.0, 4.0, 6.0])
        r6 = np.float64(2) * p
        self.assertEqual(list(r6), [2.0, 4.0, 6.0])

        p2 = Position([2.0, 3.0, 4.0])
        self.assertEqual(p * p2, 20.0)

        with self.assertRaises(TypeError):
            _ = p * "a"

    def test_parity_with_rust(self):
        try:
            import proteindf_bridge.rs as rs_br
        except ImportError:
            self.skipTest("proteindf_bridge.rs is not available")

        import numpy as np

        py_p = Position([1.0, 2.0, 3.0])
        rs_p = rs_br.Position(1.0, 2.0, 3.0)

        self.assertEqual(list(py_p), list(rs_p))

        for idx in range(-3, 3):
            self.assertEqual(py_p[idx], rs_p[idx])

        for out_idx in [3, 4, -4, -5]:
            with self.assertRaises(IndexError):
                _ = py_p[out_idx]
            with self.assertRaises(IndexError):
                _ = rs_p[out_idx]

        for out_idx in [3, -4]:
            with self.assertRaises(IndexError):
                py_p[out_idx] = 0.0
            with self.assertRaises(IndexError):
                rs_p[out_idx] = 0.0

        for factor in [2, 2.0, np.int64(2), np.float64(2)]:
            py_mul = py_p * factor
            rs_mul = rs_p * factor
            self.assertEqual(list(py_mul), list(rs_mul))

            py_rmul = factor * py_p
            rs_rmul = factor * rs_p
            self.assertEqual(list(py_rmul), list(rs_rmul))

        py_p2 = Position([2.0, 3.0, 4.0])
        rs_p2 = rs_br.Position(2.0, 3.0, 4.0)
        self.assertEqual(py_p * py_p2, rs_p * rs_p2)


def load_tests(loader, tests, ignore):
    from proteindf_bridge import position
    tests.addTests(doctest.DocTestSuite(position))
    return tests


if __name__ == '__main__':
    unittest.main()
