#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from proteindf_bridge.vector import Vector


class TestVector(unittest.TestCase):
    def test_init_with_size(self):
        v = Vector(5)
        self.assertEqual(len(v), 5)
        for i in range(5):
            self.assertEqual(v[i], 0.0)

    def test_init_with_list(self):
        v = Vector([1.0, 2.0, 3.0])
        self.assertEqual(len(v), 3)
        self.assertEqual(v[0], 1.0)
        self.assertEqual(v[1], 2.0)
        self.assertEqual(v[2], 3.0)

    def test_set_and_get(self):
        v = Vector(3)
        v.set(1, 4.5)
        self.assertEqual(v.get(1), 4.5)
        v[2] = 9.0
        self.assertEqual(v[2], 9.0)

    def test_resize(self):
        v = Vector([1.0, 2.0, 3.0])
        v.resize(5)
        self.assertEqual(len(v), 5)
        self.assertEqual(v[0], 1.0)
        self.assertEqual(v[3], 0.0)

    def test_dot_and_abs(self):
        v1 = Vector([3.0, 4.0])
        v1_abs = v1.abs()
        self.assertEqual(list(v1_abs), [3.0, 4.0])

        v2 = Vector([1.0, 2.0])
        v3 = Vector([3.0, 4.0])
        dot = v2 * v3  # Vector 同士の乗算は内積
        self.assertAlmostEqual(dot, 11.0)

    def test_arithmetic_operations(self):
        v1 = Vector([1.0, 2.0, 3.0])
        v2 = Vector([4.0, 5.0, 6.0])

        v_add = v1 + v2
        self.assertEqual(list(v_add), [5.0, 7.0, 9.0])

        v_sub = v2 - v1
        self.assertEqual(list(v_sub), [3.0, 3.0, 3.0])

        v_mul = v1 * 2.0
        self.assertEqual(list(v_mul), [2.0, 4.0, 6.0])


if __name__ == "__main__":
    unittest.main()
