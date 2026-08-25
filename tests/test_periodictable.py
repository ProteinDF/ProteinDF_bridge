#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from proteindf_bridge.periodictable import PeriodicTable


class TestPeriodicTable(unittest.TestCase):
    def test_get_symbol(self):
        self.assertEqual(PeriodicTable.get_symbol(1), "H")
        self.assertEqual(PeriodicTable.get_symbol(6), "C")
        self.assertEqual(PeriodicTable.get_symbol(8), "O")
        self.assertEqual(PeriodicTable.get_symbol(20), "Ca")

    def test_get_atomic_number(self):
        self.assertEqual(PeriodicTable.get_atomic_number("H"), 1)
        self.assertEqual(PeriodicTable.get_atomic_number("C"), 6)
        self.assertEqual(PeriodicTable.get_atomic_number("N"), 7)
        self.assertEqual(PeriodicTable.get_atomic_number("Cu"), 29)

    def test_get_weight_and_vdw(self):
        self.assertAlmostEqual(PeriodicTable.atomic_weight(6), 12.01, places=2)
        self.assertAlmostEqual(PeriodicTable.vdw(6), 1.70, places=2)


if __name__ == "__main__":
    unittest.main()
