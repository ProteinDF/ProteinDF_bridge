#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from proteindf_bridge.dbmanager import DbManager


class TestDbManager(unittest.TestCase):
    def setUp(self):
        self.db = DbManager(":memory:")

    def test_create_table_list(self):
        self.db.create_table("users", ["id", "name", "email"], primary_key="id")
        self.assertTrue(self.db.has_table("users"))
        self.assertEqual(self.db.get_table_names(), ["users"])
        self.assertEqual(self.db.get_field_names("users"), ["id", "name", "email"])

    def test_create_table_dict(self):
        self.db.create_table("products", {"id": "INTEGER", "title": "TEXT", "price": "REAL"}, primary_key="id")
        self.assertTrue(self.db.has_table("products"))
        self.assertEqual(self.db.get_field_names("products"), ["id", "title", "price"])

    def test_insert_and_select(self):
        self.db.create_table("users", ["id", "name"], primary_key="id")
        self.db.insert("users", {"id": 1, "name": "Alice"})
        self.db.insert("users", {"id": 2, "name": "Bob"})

        rows = self.db.select("users")
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["name"], "Alice")
        self.assertEqual(rows[1]["name"], "Bob")

        filtered = self.db.select("users", where="name = 'Alice'")
        self.assertEqual(len(filtered), 1)
        self.assertEqual(filtered[0]["id"], 1)

    def test_get_results_multiline(self):
        self.db.create_table("items", ["id", "val"], primary_key="id")
        for i in range(5):
            self.db.insert("items", {"id": i, "val": f"v_{i}"})

        results = self.db.get_results("SELECT * FROM items ORDER BY id ASC")
        self.assertEqual(len(results), 5)
        for i in range(5):
            self.assertEqual(results[i]["id"], i)
            self.assertEqual(results[i]["val"], f"v_{i}")

    def test_get_results_empty(self):
        self.db.create_table("items", ["id", "val"], primary_key="id")
        results = self.db.get_results("SELECT * FROM items WHERE id = 999")
        self.assertEqual(results, [])

    def test_update_and_delete(self):
        self.db.create_table("items", ["id", "val"], primary_key="id")
        self.db.insert("items", {"id": 1, "val": "old"})

        self.db.update("items", {"val": "new"}, where="id = 1")
        res = self.db.select("items", where="id = 1")
        self.assertEqual(res[0]["val"], "new")

        self.db.delete("items", where="id = 1")
        res2 = self.db.select("items", where="id = 1")
        self.assertEqual(res2, [])


if __name__ == "__main__":
    unittest.main()
