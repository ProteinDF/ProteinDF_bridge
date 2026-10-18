#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Validation tests for PyAtomGroup to CoreAtomGroup conversion.
Ensures that single-pass path generation preserves tree structure, paths,
bonds, raw data, serialization, and atom selection for 2FB4 and 2MGO.
"""

import hashlib
import os
import pickle
import unittest

import proteindf_bridge.rs as rs


class TestRsAgToCore(unittest.TestCase):
    def setUp(self):
        self.data_dir = os.path.abspath(
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
        self.pdb_2fb4 = os.path.join(self.data_dir, "2FB4.pdb")
        self.pdb_2mgo = os.path.join(self.data_dir, "2MGO.pdb")

    def test_2fb4_to_core_operations(self):
        """Verify get_bond_list, get_raw_data, str, pickle roundtrip, and select on 2FB4."""
        self.assertTrue(os.path.exists(self.pdb_2fb4), f"2FB4.pdb not found at {self.pdb_2fb4}")
        ag = rs.Pdb(self.pdb_2fb4).get_atomgroup()

        # 1. get_bond_list
        bonds = ag.get_bond_list()
        self.assertEqual(len(bonds), 6)
        expected_bonds = [
            ("/model_1/L/22/153_SG", "/model_1/L/87/648_SG", 1),
            ("/model_1/L/136/1009_SG", "/model_1/L/195/1463_SG", 1),
            ("/model_1/H/22/1760_SG", "/model_1/H/96/2356_SG", 1),
            ("/model_1/H/101D/2419_SG", "/model_1/H/104B/2448_SG", 1),
            ("/model_1/H/145/2740_SG", "/model_1/H/201/3154_SG", 1),
            ("/model_1/L/213/1595_SG", "/model_1/H/221/3312_SG", 1),
        ]
        self.assertEqual(bonds, expected_bonds)

        # 2. get_raw_data
        raw = ag.get_raw_data()
        self.assertIn("groups", raw)
        self.assertIn("model_1", raw["groups"])
        raw_hash = hashlib.sha256(str(raw).encode("utf-8")).hexdigest()
        self.assertEqual(
            raw_hash,
            "816a9b04875a2e1a5097c3e99c269dc5463725e64148368b37fb56519221cff9",
        )

        # 3. str(ag)
        s = str(ag)
        self.assertEqual(len(s), 355369)
        str_hash = hashlib.sha256(s.encode("utf-8")).hexdigest()
        self.assertEqual(
            str_hash,
            "b71c00066f126fa85459e74e3c4697133f8d56f6841ba41abc45340d237494cd",
        )

        # 4. pickle roundtrip
        pickled = pickle.dumps(ag)
        ag_restored = pickle.loads(pickled)
        self.assertEqual(ag_restored.get_number_of_all_atoms(), 3417)
        restored_str = str(ag_restored)
        self.assertEqual(len(restored_str), 355369)
        self.assertEqual(
            hashlib.sha256(restored_str.encode("utf-8")).hexdigest(),
            "48c651479c5e92410ca63c6ba031d223cecad5235921348dd15f1a28d6eaf99a",
        )
        self.assertEqual(set(ag_restored.get_bond_list()), set(expected_bonds))

        # 5. select(Select_Name("CA"))
        sel = rs.Select_Name("CA")
        ag_ca = ag.select(sel)
        ca_atoms = ag_ca.get_atom_list()
        self.assertEqual(len(ca_atoms), 445)
        ca_paths = [a.path for a in ca_atoms]
        self.assertEqual(ca_paths[0], "/model_1/L/1/2_CA")
        self.assertEqual(ca_paths[-1], "/model_1/H/221/3308_CA")
        paths_hash = hashlib.sha256(str(ca_paths).encode("utf-8")).hexdigest()
        self.assertEqual(
            paths_hash,
            "980b3b98b014729e73b5f13e97ef5edb04d19888e43793f137464a4e657cec06",
        )

    def test_2mgo_to_core_operations(self):
        """Verify get_bond_list, get_raw_data, str, pickle roundtrip, and select on 2MGO."""
        self.assertTrue(os.path.exists(self.pdb_2mgo), f"2MGO.pdb not found at {self.pdb_2mgo}")
        ag = rs.Pdb(self.pdb_2mgo).get_atomgroup()

        # 1. get_bond_list
        bonds = ag.get_bond_list()
        self.assertEqual(len(bonds), 20)
        expected_bonds = [
            (f"/model_{i}/A/1/6_SG", f"/model_{i}/A/6/89_SG", 1)
            for i in range(1, 21)
        ]
        self.assertEqual(bonds, expected_bonds)

        # 2. get_raw_data
        raw = ag.get_raw_data()
        self.assertIn("groups", raw)
        self.assertIn("model_1", raw["groups"])
        raw_hash = hashlib.sha256(str(raw).encode("utf-8")).hexdigest()
        self.assertEqual(
            raw_hash,
            "09e991a60c32df13e43ecee8e827ac6fae0838a9cbd803bfbe721970b7595ef1",
        )

        # 3. str(ag)
        s = str(ag)
        self.assertEqual(len(s), 267175)
        str_hash = hashlib.sha256(s.encode("utf-8")).hexdigest()
        self.assertEqual(
            str_hash,
            "4fcec5b42cadfc05c9784707b02a24fe489501f3a4a50ed9c2e0068342bb89ac",
        )

        # 4. pickle roundtrip
        pickled = pickle.dumps(ag)
        ag_restored = pickle.loads(pickled)
        self.assertEqual(ag_restored.get_number_of_all_atoms(), 2680)
        restored_str = str(ag_restored)
        self.assertEqual(len(restored_str), 267175)
        self.assertEqual(
            hashlib.sha256(restored_str.encode("utf-8")).hexdigest(),
            "4fcec5b42cadfc05c9784707b02a24fe489501f3a4a50ed9c2e0068342bb89ac",
        )
        self.assertEqual(ag_restored.get_bond_list(), expected_bonds)

        # 5. select(Select_Name("CA"))
        sel = rs.Select_Name("CA")
        ag_ca = ag.select(sel)
        ca_atoms = ag_ca.get_atom_list()
        self.assertEqual(len(ca_atoms), 180)
        ca_paths = [a.path for a in ca_atoms]
        self.assertEqual(ca_paths[0], "/model_1/A/1/2_CA")
        self.assertEqual(ca_paths[-1], "/model_20/A/9/128_CA")
        paths_hash = hashlib.sha256(str(ca_paths).encode("utf-8")).hexdigest()
        self.assertEqual(
            paths_hash,
            "69d63de3c1a3a636180c396b7a7207af14bc8babcc4b4ce6972a2d454445c934",
        )


if __name__ == "__main__":
    unittest.main()
