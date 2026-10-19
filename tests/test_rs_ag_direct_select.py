import unittest

import proteindf_bridge.rs as rs


class TestRsAgDirectSelect(unittest.TestCase):
    def _create_nested_tree(self):
        ag = rs.AtomGroup("root")
        m1 = rs.AtomGroup("model_1")
        c1 = rs.AtomGroup("chain_A")
        r1 = rs.AtomGroup("res_1")
        r1.secondary_structure = "H"
        a1 = rs.Atom(name="CA", symbol="C", xyz=rs.Position(1.0, 2.0, 3.0))
        a2 = rs.Atom(name="CB", symbol="C", xyz=rs.Position(4.0, 5.0, 6.0))
        a3 = rs.Atom(name="N", symbol="N", xyz=rs.Position(7.0, 8.0, 9.0))
        r1.set_atom("CA", a1)
        r1.set_atom("CB", a2)
        r1.set_atom("N", a3)
        r1.add_bond(a1, a2, 1)

        r2 = rs.AtomGroup("res_2")
        r2.secondary_structure = "E"
        a4 = rs.Atom(name="CA", symbol="C", xyz=rs.Position(10.0, 11.0, 12.0))
        a5 = rs.Atom(name="O", symbol="O", xyz=rs.Position(13.0, 14.0, 15.0))
        r2.set_atom("CA", a4)
        r2.set_atom("O", a5)

        c1.set_group("1", r1)
        c1.set_group("2", r2)
        m1.set_group("A", c1)
        ag.set_group("model_1", m1)
        return ag

    def test_nested_tree_builtin_and_custom_selectors_raw_data(self):
        ag = self._create_nested_tree()

        # 1. SelectName("CA")
        res_name_ca = ag.select(rs.SelectName("CA"))
        expected_select_name_ca = {
            "groups": {
                "model_1": {
                    "groups": {
                        "A": {
                            "groups": {
                                "1": {
                                    "atoms": {
                                        "CA": {
                                            "Z": 6,
                                            "name": "CA",
                                            "Q": 0.0,
                                            "xyz": [1.0, 2.0, 3.0],
                                            "force": [0.0, 0.0, 0.0],
                                        }
                                    },
                                    "name": "res_1",
                                },
                                "2": {
                                    "atoms": {
                                        "CA": {
                                            "Z": 6,
                                            "name": "CA",
                                            "Q": 0.0,
                                            "xyz": [10.0, 11.0, 12.0],
                                            "force": [0.0, 0.0, 0.0],
                                        }
                                    },
                                    "name": "res_2",
                                },
                            },
                            "name": "chain_A",
                        }
                    },
                    "name": "model_1",
                }
            },
            "name": "root",
        }
        self.assertEqual(res_name_ca.get_raw_data(), expected_select_name_ca)
        self.assertIsNone(
            res_name_ca.get_group("model_1").get_group("A").get_group("1").secondary_structure
        )

        # 2. SelectSymbol("C")
        res_sym_c = ag.select(rs.SelectSymbol("C"))
        expected_select_symbol_c = {
            "groups": {
                "model_1": {
                    "groups": {
                        "A": {
                            "groups": {
                                "1": {
                                    "atoms": {
                                        "CA": {
                                            "Z": 6,
                                            "name": "CA",
                                            "Q": 0.0,
                                            "xyz": [1.0, 2.0, 3.0],
                                            "force": [0.0, 0.0, 0.0],
                                        },
                                        "CB": {
                                            "Z": 6,
                                            "name": "CB",
                                            "Q": 0.0,
                                            "xyz": [4.0, 5.0, 6.0],
                                            "force": [0.0, 0.0, 0.0],
                                        },
                                    },
                                    "name": "res_1",
                                },
                                "2": {
                                    "atoms": {
                                        "CA": {
                                            "Z": 6,
                                            "name": "CA",
                                            "Q": 0.0,
                                            "xyz": [10.0, 11.0, 12.0],
                                            "force": [0.0, 0.0, 0.0],
                                        }
                                    },
                                    "name": "res_2",
                                },
                            },
                            "name": "chain_A",
                        }
                    },
                    "name": "model_1",
                }
            },
            "name": "root",
        }
        self.assertEqual(res_sym_c.get_raw_data(), expected_select_symbol_c)

        # 3. SelectPath("/model_1/A/1/") - matches subtree with bonds and secondary_structure
        res_path_r1 = ag.select(rs.SelectPath("/model_1/A/1/"))
        expected_select_path_r1 = {
            "groups": {
                "model_1": {
                    "groups": {
                        "A": {
                            "groups": {
                                "1": {
                                    "atoms": {
                                        "CA": {
                                            "Z": 6,
                                            "name": "CA",
                                            "Q": 0.0,
                                            "xyz": [1.0, 2.0, 3.0],
                                            "force": [0.0, 0.0, 0.0],
                                        },
                                        "CB": {
                                            "Z": 6,
                                            "name": "CB",
                                            "Q": 0.0,
                                            "xyz": [4.0, 5.0, 6.0],
                                            "force": [0.0, 0.0, 0.0],
                                        },
                                        "N": {
                                            "Z": 7,
                                            "name": "N",
                                            "Q": 0.0,
                                            "xyz": [7.0, 8.0, 9.0],
                                            "force": [0.0, 0.0, 0.0],
                                        },
                                    },
                                    "name": "res_1",
                                    "bonds": [["", "", 1]],
                                }
                            },
                            "name": "chain_A",
                        }
                    },
                    "name": "model_1",
                }
            },
            "name": "root",
        }
        self.assertEqual(res_path_r1.get_raw_data(), expected_select_path_r1)
        self.assertEqual(
            res_path_r1.get_group("model_1").get_group("A").get_group("1").secondary_structure,
            "H",
        )

        # 4. Selection on non-root subtree
        sub_c1 = ag.get_group("model_1").get_group("A")
        res_sub_ca = sub_c1.select(rs.SelectName("CA"))
        expected_subtree_ca = {
            "groups": {
                "1": {
                    "atoms": {
                        "CA": {
                            "Z": 6,
                            "name": "CA",
                            "Q": 0.0,
                            "xyz": [1.0, 2.0, 3.0],
                            "force": [0.0, 0.0, 0.0],
                        }
                    },
                    "name": "res_1",
                },
                "2": {
                    "atoms": {
                        "CA": {
                            "Z": 6,
                            "name": "CA",
                            "Q": 0.0,
                            "xyz": [10.0, 11.0, 12.0],
                            "force": [0.0, 0.0, 0.0],
                        }
                    },
                    "name": "res_2",
                },
            },
            "name": "chain_A",
        }
        self.assertEqual(res_sub_ca.get_raw_data(), expected_subtree_ca)
        self.assertEqual(res_sub_ca.path, "/model_1/A/")

        # 5. Custom selector matching res_1
        class Res1Selector:
            def is_match(self, obj):
                if hasattr(obj, "name"):
                    return obj.name == "res_1"
                return False

        res_custom_r1 = ag.select(Res1Selector())
        self.assertEqual(res_custom_r1.get_raw_data(), expected_select_path_r1)
        self.assertEqual(
            res_custom_r1.get_group("model_1").get_group("A").get_group("1").secondary_structure,
            "H",
        )

        # 6. Custom selector matching N atoms
        class NitrogenSelector:
            def is_match(self, obj):
                if hasattr(obj, "symbol"):
                    return obj.symbol == "N"
                return False

        res_custom_n = ag.select(NitrogenSelector())
        expected_custom_n = {
            "groups": {
                "model_1": {
                    "groups": {
                        "A": {
                            "groups": {
                                "1": {
                                    "atoms": {
                                        "N": {
                                            "Z": 7,
                                            "name": "N",
                                            "Q": 0.0,
                                            "xyz": [7.0, 8.0, 9.0],
                                            "force": [0.0, 0.0, 0.0],
                                        }
                                    },
                                    "name": "res_1",
                                }
                            },
                            "name": "chain_A",
                        }
                    },
                    "name": "model_1",
                }
            },
            "name": "root",
        }
        self.assertEqual(res_custom_n.get_raw_data(), expected_custom_n)

    def test_empty_tree_selection(self):
        empty_ag = rs.AtomGroup("empty_ag")

        # Built-in selector on empty tree (no match)
        res_empty = empty_ag.select(rs.SelectName("CA"))
        self.assertEqual(res_empty.get_raw_data(), {"name": "empty_ag"})
        self.assertEqual(res_empty.path, "/")

        # Built-in selector matching empty tree itself
        res_empty_match = empty_ag.select(rs.SelectName("empty_ag"))
        self.assertEqual(res_empty_match.get_raw_data(), {"name": "empty_ag"})
        self.assertEqual(res_empty_match.path, "/")

        # Custom selector on empty tree
        class AnySelector:
            def is_match(self, obj):
                return False

        res_custom_empty = empty_ag.select(AnySelector())
        self.assertEqual(res_custom_empty.get_raw_data(), {"name": "empty_ag"})

    def test_custom_selector_receives_original_objects(self):
        ag = self._create_nested_tree()
        model_1 = ag.get_group("model_1")
        chain_A = model_1.get_group("A")
        res_1 = chain_A.get_group("1")
        atom_ca = res_1.get_atom("CA")

        received_objects = []

        class ObjectRecorder:
            def is_match(self, obj):
                received_objects.append(obj)
                return False

        ag.select(ObjectRecorder())

        # Verify that original objects were passed to is_match without cloning
        self.assertTrue(any(obj is ag for obj in received_objects))
        self.assertTrue(any(obj is model_1 for obj in received_objects))
        self.assertTrue(any(obj is chain_A for obj in received_objects))
        self.assertTrue(any(obj is res_1 for obj in received_objects))
        self.assertTrue(any(obj is atom_ca for obj in received_objects))

    def test_returned_tree_is_independent_from_original(self):
        ag = self._create_nested_tree()
        selected = ag.select(rs.SelectName("CA"))

        orig_atom = ag.get_group("model_1").get_group("A").get_group("1").get_atom("CA")
        sel_atom = selected.get_group("model_1").get_group("A").get_group("1").get_atom("CA")

        # Verify objects are not identical
        self.assertIsNot(orig_atom, sel_atom)
        self.assertIsNot(ag, selected)

        # Modifying selected tree atom position does not mutate original tree
        sel_atom.xyz.x = 999.0
        self.assertEqual(sel_atom.xyz.x, 999.0)
        self.assertEqual(orig_atom.xyz.x, 1.0)

        # Modifying original tree atom position does not mutate selected tree
        orig_atom.xyz.x = 888.0
        self.assertEqual(orig_atom.xyz.x, 888.0)
        self.assertEqual(sel_atom.xyz.x, 999.0)

        # Modifying group hierarchy in selected tree does not affect original
        selected.get_group("model_1").name = "mutated_model"
        self.assertEqual(selected.get_group("model_1").name, "mutated_model")
        self.assertEqual(ag.get_group("model_1").name, "model_1")


if __name__ == "__main__":
    unittest.main()
