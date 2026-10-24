#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit tests validating direct Python tree traversal for get_raw_data,
get_bond_list, and pickle serialization on AtomGroup.
"""

import pickle
import unittest

import proteindf_bridge.rs as rs


class TestRsAgRawDataAndBonds(unittest.TestCase):
    def test_empty_tree(self):
        """Validate get_raw_data, get_bond_list, and pickle roundtrip on empty group."""
        ag = rs.AtomGroup("empty")

        expected_raw = {"name": "empty"}
        expected_bonds = []

        self.assertEqual(ag.get_raw_data(), expected_raw)
        self.assertEqual(ag.get_bond_list(), expected_bonds)

        pickled = pickle.dumps(ag)
        restored = pickle.loads(pickled)

        self.assertEqual(restored.name, "empty")
        self.assertEqual(restored.get_raw_data(), expected_raw)
        self.assertEqual(restored.get_bond_list(), expected_bonds)
        self.assertEqual(restored.get_number_of_all_atoms(), 0)

    def test_single_atom_tree(self):
        """Validate get_raw_data, get_bond_list, and pickle roundtrip on single-atom group."""
        ag = rs.AtomGroup("single")
        atom = rs.Atom(
            name="CA",
            symbol="C",
            charge=0.0,
            xyz=[1.0, 2.0, 3.0],
            force=[0.0, 0.0, 0.0],
        )
        ag.set_atom("1_CA", atom)

        expected_raw = {
            "atoms": {
                "1_CA": {
                    "Z": 6,
                    "name": "CA",
                    "Q": 0.0,
                    "xyz": [1.0, 2.0, 3.0],
                    "force": [0.0, 0.0, 0.0],
                }
            },
            "name": "single",
        }
        expected_bonds = []

        self.assertEqual(ag.get_raw_data(), expected_raw)
        self.assertEqual(ag.get_bond_list(), expected_bonds)

        pickled = pickle.dumps(ag)
        restored = pickle.loads(pickled)

        self.assertEqual(restored.name, "single")
        self.assertEqual(restored.get_raw_data(), expected_raw)
        self.assertEqual(restored.get_bond_list(), expected_bonds)
        self.assertEqual(restored.get_number_of_all_atoms(), 1)

    def test_nested_tree_with_bonds_and_secondary_structure(self):
        """Validate nested tree with relative bonds, secondary structures, and subtree calls."""
        ag_root = rs.AtomGroup("root")

        sub = rs.AtomGroup("sub")
        sub.secondary_structure = "H"
        a1 = rs.Atom(name="N", symbol="N", xyz=[0.0, 0.0, 0.0], charge=0.0, force=[0.0, 0.0, 0.0])
        a2 = rs.Atom(name="CA", symbol="C", xyz=[1.0, 0.0, 0.0], charge=0.0, force=[0.0, 0.0, 0.0])
        sub.set_atom("1_N", a1)
        sub.set_atom("2_CA", a2)

        sub2 = rs.AtomGroup("sub2")
        a3 = rs.Atom(name="C", symbol="C", xyz=[2.0, 0.0, 0.0], charge=0.0, force=[0.0, 0.0, 0.0])
        sub2.set_atom("3_C", a3)

        ag_root.set_group("sub", sub)
        ag_root.set_group("sub2", sub2)

        # Relative bond inside child group
        ag_root["sub"].add_bond(
            ag_root["sub"].get_atom("1_N"),
            ag_root["sub"].get_atom("2_CA"),
            1,
        )

        # Cross-group bond in root
        ag_root.add_bond(
            ag_root["sub"].get_atom("2_CA"),
            ag_root["sub2"].get_atom("3_C"),
            2,
        )

        expected_root_raw = {
            "groups": {
                "sub": {
                    "atoms": {
                        "1_N": {
                            "Z": 7,
                            "name": "N",
                            "Q": 0.0,
                            "xyz": [0.0, 0.0, 0.0],
                            "force": [0.0, 0.0, 0.0],
                        },
                        "2_CA": {
                            "Z": 6,
                            "name": "CA",
                            "Q": 0.0,
                            "xyz": [1.0, 0.0, 0.0],
                            "force": [0.0, 0.0, 0.0],
                        },
                    },
                    "name": "sub",
                    "bonds": [["1_N", "2_CA", 1]],
                },
                "sub2": {
                    "atoms": {
                        "3_C": {
                            "Z": 6,
                            "name": "C",
                            "Q": 0.0,
                            "xyz": [2.0, 0.0, 0.0],
                            "force": [0.0, 0.0, 0.0],
                        },
                    },
                    "name": "sub2",
                },
            },
            "name": "root",
            "bonds": [["sub/2_CA", "sub2/3_C", 2]],
        }

        expected_root_bonds = [
            ("/sub/1_N", "/sub/2_CA", 1),
            ("/sub/2_CA", "/sub2/3_C", 2),
        ]

        expected_sub_raw = {
            "atoms": {
                "1_N": {
                    "Z": 7,
                    "name": "N",
                    "Q": 0.0,
                    "xyz": [0.0, 0.0, 0.0],
                    "force": [0.0, 0.0, 0.0],
                },
                "2_CA": {
                    "Z": 6,
                    "name": "CA",
                    "Q": 0.0,
                    "xyz": [1.0, 0.0, 0.0],
                    "force": [0.0, 0.0, 0.0],
                },
            },
            "name": "sub",
            "bonds": [["1_N", "2_CA", 1]],
        }

        expected_sub_bonds = [
            ("/sub/1_N", "/sub/2_CA", 1),
        ]

        # 1. Root group checks
        self.assertEqual(ag_root.get_raw_data(), expected_root_raw)
        self.assertEqual(ag_root.get_bond_list(), expected_root_bonds)

        # 2. Non-root subtree checks
        self.assertEqual(ag_root["sub"].get_raw_data(), expected_sub_raw)
        self.assertEqual(ag_root["sub"].get_bond_list(), expected_sub_bonds)

        # 3. Pickle roundtrip of root
        pickled_root = pickle.dumps(ag_root)
        restored_root = pickle.loads(pickled_root)

        self.assertEqual(restored_root.name, "root")
        self.assertEqual(restored_root.get_raw_data(), expected_root_raw)
        self.assertEqual(restored_root.get_bond_list(), expected_root_bonds)
        self.assertEqual(restored_root["sub"].secondary_structure, "H")
        self.assertEqual(restored_root["sub"].get_raw_data(), expected_sub_raw)
        self.assertEqual(restored_root["sub"].get_bond_list(), expected_sub_bonds)

        # 4. Pickle roundtrip of non-root subtree
        pickled_sub = pickle.dumps(ag_root["sub"])
        restored_sub = pickle.loads(pickled_sub)

        self.assertEqual(restored_sub.name, "sub")
        self.assertEqual(restored_sub.get_raw_data(), expected_sub_raw)
        self.assertEqual(restored_sub.get_bond_list(), expected_sub_bonds)
        self.assertEqual(restored_sub.secondary_structure, "H")


if __name__ == "__main__":
    unittest.main()
