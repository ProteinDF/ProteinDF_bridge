#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import numpy
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

    def test_buffer_roundtrip(self):
        v = Vector([1.0, 2.5, -3.0])
        buf = v.get_buffer()
        v2 = Vector()
        v2.set_buffer(buf)
        self.assertEqual(list(v), list(v2))

    def test_buffer_resize(self):
        v = Vector(2)
        v_target = Vector([1.0, 2.0, 3.0, 4.0])
        v.set_buffer(v_target.get_buffer())
        self.assertEqual(len(v), 4)
        self.assertEqual(list(v), [1.0, 2.0, 3.0, 4.0])

    def test_buffer_invalid_length(self):
        v = Vector([1.0, 2.0, 3.0])
        with self.assertRaises(ValueError):
            v.set_buffer(b"1234567")
        self.assertEqual(len(v), 3)
        self.assertEqual(list(v), [1.0, 2.0, 3.0])

    def test_buffer_writable(self):
        v = Vector(3)
        buf = v.get_buffer()
        v2 = Vector(3)
        v2.set_buffer(buf)
        v2[0] = 42.0
        self.assertEqual(v2[0], 42.0)

    def test_buffer_memoryview_float(self):
        v = Vector([1.0, 2.5, -3.0])
        mv = memoryview(numpy.array([1.0, 2.5, -3.0], dtype=float))
        v2 = Vector()
        v2.set_buffer(mv)
        self.assertEqual(list(v), list(v2))


    def test_init_with_dict(self):
        v = Vector({"size": 3, "data": [1.0, 2.0, 3.0]})
        self.assertEqual(len(v), 3)
        self.assertEqual(list(v), [1.0, 2.0, 3.0])

        v_empty = Vector({"size": 3})
        self.assertEqual(len(v_empty), 3)
        self.assertEqual(list(v_empty), [0.0, 0.0, 0.0])

        with self.assertRaises(ValueError):
            Vector({"size": 3, "data": [1.0, 2.0]})

    def test_multiplication_scalars_and_types(self):
        v = Vector([1.0, 2.0, 3.0])

        # int
        res_i1 = v * 2
        self.assertIsInstance(res_i1, Vector)
        self.assertEqual(list(res_i1), [2.0, 4.0, 6.0])

        res_i2 = 3 * v
        self.assertIsInstance(res_i2, Vector)
        self.assertEqual(list(res_i2), [3.0, 6.0, 9.0])

        v_inplace = Vector([1.0, 2.0, 3.0])
        v_inplace *= 2
        self.assertIsInstance(v_inplace, Vector)
        self.assertEqual(list(v_inplace), [2.0, 4.0, 6.0])

        # float
        res_f1 = v * 2.0
        self.assertIsInstance(res_f1, Vector)
        self.assertEqual(list(res_f1), [2.0, 4.0, 6.0])

        res_f2 = 3.0 * v
        self.assertIsInstance(res_f2, Vector)
        self.assertEqual(list(res_f2), [3.0, 6.0, 9.0])

        v_inplace_f = Vector([1.0, 2.0, 3.0])
        v_inplace_f *= 2.0
        self.assertIsInstance(v_inplace_f, Vector)
        self.assertEqual(list(v_inplace_f), [2.0, 4.0, 6.0])

        # numpy real numbers
        np_f = numpy.float64(2.5)
        np_i = numpy.int64(3)

        res_np_f1 = v * np_f
        self.assertIsInstance(res_np_f1, Vector)
        self.assertEqual(list(res_np_f1), [2.5, 5.0, 7.5])

        res_np_f2 = np_f * v
        self.assertIsInstance(res_np_f2, Vector)
        self.assertEqual(list(res_np_f2), [2.5, 5.0, 7.5])

        v_inplace_np = Vector([1.0, 2.0, 3.0])
        v_inplace_np *= np_f
        self.assertIsInstance(v_inplace_np, Vector)
        self.assertEqual(list(v_inplace_np), [2.5, 5.0, 7.5])

        res_np_i1 = v * np_i
        self.assertIsInstance(res_np_i1, Vector)
        self.assertEqual(list(res_np_i1), [3.0, 6.0, 9.0])

        res_np_i2 = np_i * v
        self.assertIsInstance(res_np_i2, Vector)
        self.assertEqual(list(res_np_i2), [3.0, 6.0, 9.0])

        v_inplace_np2 = Vector([1.0, 2.0, 3.0])
        v_inplace_np2 *= np_i
        self.assertIsInstance(v_inplace_np2, Vector)
        self.assertEqual(list(v_inplace_np2), [3.0, 6.0, 9.0])

        # unsupported types (TypeError)
        with self.assertRaises(TypeError):
            _ = v * "str"
        with self.assertRaises(TypeError):
            _ = "str" * v
        with self.assertRaises(TypeError):
            v_err = Vector([1.0, 2.0, 3.0])
            v_err *= "str"

        # bool must raise TypeError
        with self.assertRaises(TypeError):
            _ = v * True
        with self.assertRaises(TypeError):
            _ = True * v
        with self.assertRaises(TypeError):
            v_err = Vector([1.0, 2.0, 3.0])
            v_err *= True


if __name__ == "__main__":
    unittest.main()
