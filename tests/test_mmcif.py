#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import tempfile
import unittest
from proteindf_bridge.mmcif import SimpleMmcif


class TestSimpleMmcif(unittest.TestCase):
    def test_load_mmcif_basic(self):
        mmcif_content = """data_test_entry
_entry.id test_entry
_cell.length_a 10.0
_cell.length_b 20.0
_cell.length_c 30.0
"""
        with tempfile.NamedTemporaryFile(suffix=".cif", mode="w", delete=False) as f:
            f.write(mmcif_content)
            temp_path = f.name

        try:
            cif = SimpleMmcif(temp_path)
            self.assertIn("data_test_entry", cif._data)
            data_block, tables = cif._data["data_test_entry"]
            self.assertEqual(data_block["_entry.id"], "test_entry")
            self.assertEqual(data_block["_cell.length_a"], "10.0")
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)


if __name__ == "__main__":
    unittest.main()
