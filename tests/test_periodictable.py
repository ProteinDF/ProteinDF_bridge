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

    def test_covalent_radius(self):
        # Rust版periodic_table.rsの値: H=0.31, C=0.76, N=0.71, O=0.66
        self.assertAlmostEqual(PeriodicTable.covalent_radius("H"), 0.31, places=2)
        self.assertAlmostEqual(PeriodicTable.covalent_radius(1), 0.31, places=2)
        self.assertAlmostEqual(PeriodicTable.covalent_radius("C"), 0.76, places=2)
        self.assertAlmostEqual(PeriodicTable.covalent_radius(6), 0.76, places=2)
        self.assertAlmostEqual(PeriodicTable.covalent_radius("N"), 0.71, places=2)
        self.assertAlmostEqual(PeriodicTable.covalent_radius(7), 0.71, places=2)
        self.assertAlmostEqual(PeriodicTable.covalent_radius("O"), 0.66, places=2)
        self.assertAlmostEqual(PeriodicTable.covalent_radius(8), 0.66, places=2)

        # 原子番号0や表にない原子で例外
        with self.assertRaises((IndexError, KeyError)):
            PeriodicTable.covalent_radius(0)
        with self.assertRaises((IndexError, KeyError)):
            PeriodicTable.covalent_radius(999)
        with self.assertRaises((ValueError, IndexError, KeyError)):
            PeriodicTable.covalent_radius("Unknown")


if __name__ == "__main__":
    unittest.main()
