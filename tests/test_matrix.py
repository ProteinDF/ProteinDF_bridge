#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import numpy
from proteindf_bridge.matrix import Matrix, SymmetricMatrix
from proteindf_bridge.vector import Vector


class TestMatrix(unittest.TestCase):
    def test_init_and_shape(self):
        m = Matrix(2, 3)
        self.assertEqual(m.rows, 2)
        self.assertEqual(m.cols, 3)

    def test_init_with_nested_list(self):
        data = [[1.0, 2.0], [3.0, 4.0]]
        m = Matrix(data)
        self.assertEqual(m.rows, 2)
        self.assertEqual(m.cols, 2)
        self.assertAlmostEqual(m.get(0, 1), 2.0)
        self.assertAlmostEqual(m.get(1, 0), 3.0)

    def test_set_and_get(self):
        m = Matrix(2, 2)
        m.set(0, 1, 5.5)
        self.assertAlmostEqual(m.get(0, 1), 5.5)
        m.add(0, 1, 2.0)
        self.assertAlmostEqual(m.get(0, 1), 7.5)

    def test_multiplication(self):
        c = Matrix([[7, 4, -1], [3, 0, 5]])
        d = Matrix([[8, 4, 2], [1, 3, -6], [-7, 0, 5]])
        cd = c * d
        expected = Matrix([[67, 40, -15], [-11, 12, 31]])
        self.assertEqual(cd, expected)

    def test_vector_multiplication(self):
        m = Matrix([[1.0, 2.0], [3.0, 4.0]])
        v = Vector([1.0, 1.0])
        res = m * v
        self.assertAlmostEqual(res[0], 3.0)
        self.assertAlmostEqual(res[1], 7.0)

    def test_buffer_roundtrip_square(self):
        m = Matrix([[1.0, 2.0], [3.0, 4.0]])
        buf = m.get_buffer()
        m2 = Matrix(2, 2)
        m2.set_buffer(buf)
        self.assertEqual(m, m2)

    def test_buffer_roundtrip_non_square(self):
        m = Matrix([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        buf = m.get_buffer()
        m2 = Matrix(2, 3)
        m2.set_buffer(buf)
        self.assertEqual(m, m2)

    def test_buffer_invalid_length(self):
        m = Matrix([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        with self.assertRaises(ValueError):
            m.set_buffer(b"too_short")
        self.assertEqual(m.rows, 2)
        self.assertEqual(m.cols, 3)
        self.assertAlmostEqual(m.get(0, 0), 1.0)
        self.assertAlmostEqual(m.get(1, 2), 6.0)

    def test_buffer_writable_after_set(self):
        m = Matrix([[1.0, 2.0], [3.0, 4.0]])
        buf = m.get_buffer()
        m2 = Matrix(2, 2)
        m2.set_buffer(buf)
        m2.set(0, 1, 99.0)
        self.assertAlmostEqual(m2.get(0, 1), 99.0)

    def test_buffer_bytearray(self):
        m = Matrix([[1.0, 2.0], [3.0, 4.0]])
        ba = bytearray(m.get_buffer())
        m2 = Matrix(2, 2)
        m2.set_buffer(ba)
        self.assertEqual(m, m2)

    def test_buffer_memoryview_float(self):
        m = Matrix([[1.0, 2.0], [3.0, 4.0]])
        mv = memoryview(numpy.array([[1.0, 2.0], [3.0, 4.0]], dtype=float))
        m2 = Matrix(2, 2)
        m2.set_buffer(mv)
        self.assertEqual(m, m2)

    def test_get_symmetric_matrix_warning(self):
        m = Matrix([[1.0, 2.0], [3.0, 4.0]])
        with self.assertLogs("proteindf_bridge.matrix", level="WARNING") as cm:
            sm = m.get_symmetric_matrix()
        self.assertIsNotNone(sm)
        self.assertIsInstance(sm, SymmetricMatrix)
        self.assertTrue(any("warning:" in msg for msg in cm.output))


class TestSymmetricMatrix(unittest.TestCase):
    def test_symmetric_indexing(self):
        sm = SymmetricMatrix(3)
        self.assertEqual(sm.dim, 3)
        sm.set(0, 1, 2.5)
        # 対称なので (1, 0) も 2.5
        self.assertAlmostEqual(sm.get(1, 0), 2.5)
        self.assertAlmostEqual(sm.get(0, 1), 2.5)

    def test_symmetric_add(self):
        sm = SymmetricMatrix(3)
        sm.add(0, 1, 3.5)
        self.assertAlmostEqual(sm.get(1, 0), 3.5)
        self.assertAlmostEqual(sm.get(0, 1), 3.5)
        # さらに加算
        sm.add(1, 0, 1.5)
        self.assertAlmostEqual(sm.get(1, 0), 5.0)
        self.assertAlmostEqual(sm.get(0, 1), 5.0)

    def test_get_raw_data(self):
        sm = SymmetricMatrix(3)
        sm.set(0, 0, 1.0)
        sm.set(1, 0, 2.0)
        sm.set(1, 1, 3.0)
        raw = sm.get_raw_data()
        self.assertEqual(raw["row"], 3)
        self.assertEqual(raw["col"], 3)
        self.assertEqual(raw["type"], "SP")
        self.assertEqual(len(raw["data"]), 6)

    def test_eigenvalues(self):
        sm = SymmetricMatrix(2)
        sm.set(0, 0, 2.0)
        sm.set(1, 1, 2.0)
        sm.set(0, 1, 1.0)
        eigvals, eigvecs = sm.eig()
        self.assertEqual(len(eigvals), 2)
        self.assertAlmostEqual(eigvals[0], 1.0)
        self.assertAlmostEqual(eigvals[1], 3.0)

    def test_symmetric_buffer_roundtrip(self):
        sm = SymmetricMatrix(3)
        sm.set(0, 0, 1.0)
        sm.set(0, 1, 2.0)
        sm.set(1, 1, 3.0)
        sm.set(2, 0, 4.0)
        buf = sm.get_buffer()
        sm2 = SymmetricMatrix(3)
        sm2.set_buffer(buf)
        self.assertEqual(sm, sm2)


if __name__ == "__main__":
    unittest.main()
