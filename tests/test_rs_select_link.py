# -*- coding: utf-8 -*-
import os
import unittest

from proteindf_bridge import (
    AtomGroup,
    Select_Path_simple,
    Select_Path_wildcard,
    load_atomgroup,
)


class TestRsSelectLink(unittest.TestCase):
    def _create_dict_tree(self):
        return {
            "name": "root",
            "groups": {
                "model_1": {
                    "name": "model_1",
                    "groups": {
                        "A": {
                            "name": "A",
                            "groups": {
                                "1": {
                                    "name": "res_1",
                                    "atoms": {
                                        "N": {
                                            "name": "N",
                                            "symbol": "N",
                                            "xyz": [1.0, 2.0, 3.0],
                                        },
                                        "CA": {
                                            "name": "CA",
                                            "symbol": "C",
                                            "xyz": [4.0, 5.0, 6.0],
                                        },
                                    },
                                },
                                "2": {
                                    "name": "res_2",
                                    "atoms": {
                                        "N": {
                                            "name": "N",
                                            "symbol": "N",
                                            "xyz": [7.0, 8.0, 9.0],
                                        },
                                        "CA": {
                                            "name": "CA",
                                            "symbol": "C",
                                            "xyz": [10.0, 11.0, 12.0],
                                        },
                                        "C": {
                                            "name": "C",
                                            "symbol": "C",
                                            "xyz": [13.0, 14.0, 15.0],
                                        },
                                    },
                                },
                            },
                        }
                    },
                }
            },
        }

    def test_select_path_wildcard_without_prior_traversal(self):
        data = self._create_dict_tree()

        # Without prior traversal
        ag_no_traverse = AtomGroup(data)
        sel_wildcard = Select_Path_wildcard("/model_1/A/2/")
        res_no_traverse = ag_no_traverse.select(sel_wildcard)

        # With prior traversal
        ag_traverse = AtomGroup(data)
        _ = ag_traverse.get_group("model_1").get_group("A").get_group("2")
        res_traverse = ag_traverse.select(sel_wildcard)

        self.assertEqual(res_no_traverse.get_number_of_all_atoms(), 3)
        self.assertEqual(
            res_no_traverse.get_number_of_all_atoms(),
            res_traverse.get_number_of_all_atoms(),
        )

    def test_select_various_path_queries_consistency(self):
        data = self._create_dict_tree()

        # Prior traversed baseline
        ag_baseline = AtomGroup(data)
        _ = ag_baseline.get_group("model_1").get_group("A").get_group("2")
        expected_atoms = ag_baseline.select(Select_Path_wildcard("/model_1/A/2/")).get_number_of_all_atoms()
        self.assertEqual(expected_atoms, 3)

        # Select_Path_simple without prior traversal
        ag_simple = AtomGroup(data)
        res_simple = ag_simple.select(Select_Path_simple("/model_1/A/2/"))
        self.assertEqual(res_simple.get_number_of_all_atoms(), 3)
        self.assertEqual(res_simple.get_number_of_all_atoms(), expected_atoms)

        # Select_Path_wildcard with wildcard * without prior traversal
        ag_star = AtomGroup(data)
        res_star = ag_star.select(Select_Path_wildcard("/model_1/A/2/*"))
        self.assertEqual(res_star.get_number_of_all_atoms(), 3)
        self.assertEqual(res_star.get_number_of_all_atoms(), expected_atoms)

    def test_selected_result_paths(self):
        data = self._create_dict_tree()

        ag_no_traverse = AtomGroup(data)
        res_no_traverse = ag_no_traverse.select(Select_Path_wildcard("/model_1/A/2/"))

        ag_traverse = AtomGroup(data)
        _ = ag_traverse.get_group("model_1").get_group("A").get_group("2")
        res_traverse = ag_traverse.select(Select_Path_wildcard("/model_1/A/2/"))

        group_no_traverse = res_no_traverse.get_group("model_1").get_group("A").get_group("2")
        group_traverse = res_traverse.get_group("model_1").get_group("A").get_group("2")

        # Explicit expected path
        expected_path = "/model_1/A/2/"
        self.assertEqual(group_no_traverse.path, expected_path)
        self.assertEqual(group_no_traverse.path, group_traverse.path)

        model_no_traverse = res_no_traverse.get_group("model_1")
        model_traverse = res_traverse.get_group("model_1")
        self.assertEqual(model_no_traverse.path, "/model_1/")
        self.assertEqual(model_no_traverse.path, model_traverse.path)

        chain_no_traverse = res_no_traverse.get_group("model_1").get_group("A")
        chain_traverse = res_traverse.get_group("model_1").get_group("A")
        self.assertEqual(chain_no_traverse.path, "/model_1/A/")
        self.assertEqual(chain_no_traverse.path, chain_traverse.path)

    def test_select_brd_file_without_prior_traversal(self):
        brd_path = os.path.join(
            os.path.dirname(__file__), "..", "proteindf_bridge", "data", "ACE_ALA_NME.brd"
        )
        ag_no_traverse = load_atomgroup(brd_path)
        res_no_traverse = ag_no_traverse.select(Select_Path_wildcard("/2/"))

        ag_traverse = load_atomgroup(brd_path)
        _ = ag_traverse.get_group("2")
        res_traverse = ag_traverse.select(Select_Path_wildcard("/2/"))

        self.assertEqual(res_no_traverse.get_number_of_all_atoms(), 10)
        self.assertEqual(
            res_no_traverse.get_number_of_all_atoms(),
            res_traverse.get_number_of_all_atoms(),
        )
        self.assertEqual(res_no_traverse.get_group("2").path, "/2/")


if __name__ == "__main__":
    unittest.main()
