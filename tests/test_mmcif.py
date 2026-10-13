#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import tempfile
import unittest

from proteindf_bridge.mmcif import SimpleMmcif
import proteindf_bridge.rs as rs_br

DATA_DIR = os.path.abspath(
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "rust",
        "crates",
        "proteindf-bridge",
        "tests",
        "data",
    )
)


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

    def test_load_pdb_entries(self):
        """Test that 5 PDB entries load without error and match Rust version data."""
        cif_files = ["1HLS.cif", "1WCT.cif", "2FB4.cif", "2MGO.cif", "3I3Z.cif"]
        for filename in cif_files:
            cif_path = os.path.join(DATA_DIR, filename)
            self.assertTrue(os.path.exists(cif_path), f"File not found: {cif_path}")

            py_cif = SimpleMmcif(cif_path)
            rs_cif = rs_br.SimpleMmcif(cif_path)

            self.assertEqual(
                py_cif.get_molecule_names(),
                rs_cif.get_molecule_names(),
                f"Block names mismatch for {filename}",
            )

            # Compare _data structure with Rust data
            for block_name, (py_kv, py_tables) in py_cif._data.items():
                self.assertIn(block_name, rs_cif.data)
                rs_kv, rs_tables = rs_cif.data[block_name]

                self.assertEqual(py_kv, rs_kv, f"Key-values mismatch in {filename} [{block_name}]")
                self.assertEqual(len(py_tables), len(rs_tables), f"Table count mismatch in {filename} [{block_name}]")
                for i, (py_table, rs_table) in enumerate(zip(py_tables, rs_tables)):
                    self.assertEqual(
                        len(py_table),
                        len(rs_table),
                        f"Row count mismatch in {filename} table {i}",
                    )
                    self.assertEqual(
                        py_table,
                        rs_table,
                        f"Table rows mismatch in {filename} table {i}",
                    )

    def test_load_ala_quotes(self):
        """Test ALA.cif quoted values are parsed correctly without concatenating subsequent fields."""
        cif_path = os.path.join(DATA_DIR, "ALA.cif")
        self.assertTrue(os.path.exists(cif_path), f"File not found: {cif_path}")

        py_cif = SimpleMmcif(cif_path)
        self.assertIn("data_ALA", py_cif._data)
        kv, tables = py_cif._data["data_ALA"]

        # Find _pdbx_chem_comp_identifier table
        identifier_table = None
        for table in tables:
            if table and "_pdbx_chem_comp_identifier.type" in table[0]:
                identifier_table = table
                break

        self.assertIsNotNone(identifier_table, "_pdbx_chem_comp_identifier table not found")
        # Check that "SYSTEMATIC NAME" is read as an individual field value
        types = [row["_pdbx_chem_comp_identifier.type"] for row in identifier_table]
        self.assertIn("SYSTEMATIC NAME", types)

        # Check atomgroup extraction
        ag = py_cif.get_atomgroup("data_ALA")
        self.assertEqual(ag.name, "ALA")
        self.assertEqual(ag.get_number_of_atoms(), 13)

    def test_save_load_msgpack_roundtrip(self):
        """Test save_msgpack and load_msgpack roundtrip."""
        cif_path = os.path.join(DATA_DIR, "ALA.cif")
        py_cif = SimpleMmcif(cif_path)

        with tempfile.NamedTemporaryFile(suffix=".msgpack", delete=False) as f:
            temp_path = f.name

        try:
            py_cif.save_msgpack(temp_path)

            loaded_cif = SimpleMmcif()
            loaded_cif.load_msgpack(temp_path)

            self.assertEqual(loaded_cif.get_molecule_names(), py_cif.get_molecule_names())
            # msgpack serializes tuples as lists, but keys and contents are identical
            for block_name in py_cif._data:
                orig_kv, orig_tables = py_cif._data[block_name]
                loaded_kv, loaded_tables = loaded_cif._data[block_name]
                self.assertEqual(orig_kv, loaded_kv)
                self.assertEqual(orig_tables, loaded_tables)
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)

    def test_nonexistent_file_raises_error(self):
        """Test that attempting to load a nonexistent file raises an error from Rust."""
        nonexistent = "/path/to/nonexistent/file.cif"
        with self.assertRaises(Exception):
            SimpleMmcif(nonexistent)


if __name__ == "__main__":
    unittest.main()
