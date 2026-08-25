#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from proteindf_bridge.str_processing import StrUtils


class TestStrProcessing(unittest.TestCase):
    def test_sort_nicely(self):
        input_list = ["item10", "item2", "item1"]
        sorted_list = StrUtils.sort_nicely(input_list)
        self.assertEqual(sorted_list, ["item1", "item2", "item10"])

    def test_get_common_str(self):
        self.assertEqual(StrUtils.get_common_str("abcdef", "abcxyz"), "abc")
        self.assertEqual(StrUtils.get_common_str("hello", "world"), "")

    def test_str_to_bool(self):
        self.assertTrue(StrUtils.str_to_bool("1"))
        self.assertTrue(StrUtils.str_to_bool("true"))
        self.assertTrue(StrUtils.str_to_bool("YES"))
        self.assertFalse(StrUtils.str_to_bool("0"))
        self.assertFalse(StrUtils.str_to_bool("false"))
        self.assertFalse(StrUtils.str_to_bool("NO"))

    def test_unicode_conversion(self):
        b = b"hello"
        u = StrUtils.to_unicode(b)
        self.assertEqual(u, "hello")
        self.assertIsInstance(u, str)

        b_converted = StrUtils.to_bytes("hello")
        self.assertEqual(b_converted, b"hello")
        self.assertIsInstance(b_converted, bytes)


if __name__ == "__main__":
    unittest.main()
