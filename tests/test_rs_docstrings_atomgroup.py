#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests validating that AtomGroup class and its public methods and properties have docstrings."""

import inspect
import unittest

from proteindf_bridge import AtomGroup


class TestAtomGroupDocstrings(unittest.TestCase):
    """Test suite ensuring all public interfaces of AtomGroup have non-empty docstrings."""

    EXPECTED_PUBLIC_MEMBERS = [
        "name",
        "parent",
        "path",
        "path_depth",
        "is_model_level",
        "is_chain_level",
        "is_residue_level",
        "validate_schema",
        "secondary_structure",
        "get_number_of_atoms",
        "get_number_of_all_atoms",
        "get_number_of_groups",
        "get_number_of_bonds",
        "get_atom",
        "set_atom",
        "has_atom",
        "has_atomkey",
        "has_atomname",
        "del_atom",
        "remove_atom",
        "erase_atom",
        "get_group",
        "set_group",
        "has_group",
        "has_groupkey",
        "has_groupname",
        "del_group",
        "remove_group",
        "erase_group",
        "sort_atoms",
        "sort_groups",
        "atoms",
        "groups",
        "get_group_list",
        "get_atom_keys",
        "get_atom_list",
        "get_path_list",
        "pickup_atoms",
        "assign_charges",
        "get_family",
        "divide_path",
        "charge",
        "real_nuclei_charge",
        "nuclei_charge",
        "get_atom_kinds",
        "get_xyz",
        "get_atom_kinds_count",
        "get_formula",
        "formula",
        "weight",
        "center",
        "box",
        "get_box",
        "merge",
        "select",
        "add_bond",
        "bonds",
        "get_bond_list",
        "setup",
        "setup_with_db",
        "add_missing_hydrogens",
        "apply_secondary_structure",
        "shift_by",
        "rotate",
        "sum_of_atomic_number",
        "restructure",
        "get_raw_data",
        "set_by_dict_data",
        "save_csv",
    ]

    def test_class_docstring(self):
        """Verify that the AtomGroup class itself has a non-empty docstring."""
        self.assertIsNotNone(AtomGroup.__doc__)
        self.assertTrue(
            bool(AtomGroup.__doc__.strip()),
            "AtomGroup class docstring must not be empty",
        )

    def test_expected_public_members_docstrings(self):
        """Verify that all expected public methods and properties have non-empty docstrings."""
        for member_name in self.EXPECTED_PUBLIC_MEMBERS:
            with self.subTest(member=member_name):
                self.assertTrue(
                    hasattr(AtomGroup, member_name),
                    f"AtomGroup should have attribute '{member_name}'",
                )
                attr = inspect.getattr_static(AtomGroup, member_name)
                doc = inspect.getdoc(attr)
                self.assertIsNotNone(
                    doc,
                    f"Docstring of '{member_name}' should not be None",
                )
                self.assertTrue(
                    bool(doc.strip()),
                    f"Docstring of '{member_name}' should not be empty",
                )

    def test_all_public_members_docstrings(self):
        """Verify that any public member found on AtomGroup has a non-empty docstring."""
        public_members = [m for m in dir(AtomGroup) if not m.startswith("_")]
        self.assertTrue(len(public_members) > 0)
        for member_name in public_members:
            with self.subTest(member=member_name):
                attr = inspect.getattr_static(AtomGroup, member_name)
                doc = inspect.getdoc(attr)
                self.assertIsNotNone(
                    doc,
                    f"Public attribute '{member_name}' has no docstring",
                )
                self.assertTrue(
                    bool(doc.strip()),
                    f"Public attribute '{member_name}' docstring is empty",
                )


if __name__ == "__main__":
    unittest.main()
