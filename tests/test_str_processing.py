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

    def test_to_unicode_dict(self):
        d = {
            0: 1.0,
            1: b"val1",
            b"str_key": b"val2",
            "normal_key": [b"elem", {2: b"nested_val"}],
            3: {4: 5.0, b"sub_key": [0, 1]},
        }
        res = StrUtils.to_unicode_dict(d)
        expected = {
            0: 1.0,
            1: "val1",
            "str_key": "val2",
            "normal_key": ["elem", {2: "nested_val"}],
            3: {4: 5.0, "sub_key": [0, 1]},
        }
        self.assertEqual(res, expected)
        self.assertIn(0, res)
        self.assertNotIn("0", res)
        self.assertIn("str_key", res)
        self.assertIn(2, res["normal_key"][1])
        self.assertIn(4, res[3])
        self.assertIn("sub_key", res[3])


if __name__ == "__main__":
    unittest.main()
