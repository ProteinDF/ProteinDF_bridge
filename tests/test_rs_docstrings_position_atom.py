#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test suite validating docstrings on Position, Atom, and exception classes
(BrError, BrInputError, BrValueError) as well as their public methods and properties.
"""

import unittest

from proteindf_bridge import Atom, Position
from proteindf_bridge.error import BrError, BrInputError, BrValueError


class TestRsDocstringsPositionAtom(unittest.TestCase):
    """Test that Position, Atom, and exception classes have non-empty docstrings."""

    def test_exception_classes_docstrings(self):
        """Verify that BrError, BrInputError, and BrValueError have non-empty docstrings."""
        for exc in (BrError, BrInputError, BrValueError):
            doc = exc.__doc__
            self.assertIsNotNone(doc, f"{exc.__name__}.__doc__ is None")
            self.assertTrue(
                isinstance(doc, str) and len(doc.strip()) > 0,
                f"{exc.__name__}.__doc__ is empty",
            )

    def test_position_class_docstring(self):
        """Verify that Position class has a non-empty docstring."""
        doc = Position.__doc__
        self.assertIsNotNone(doc, "Position.__doc__ is None")
        self.assertTrue(
            isinstance(doc, str) and len(doc.strip()) > 0,
            "Position.__doc__ is empty",
        )

    def test_position_public_members_docstrings(self):
        """Verify that all public methods and properties of Position have non-empty docstrings."""
        expected_properties = ["x", "y", "z", "xyz", "position", "epsilon"]
        expected_methods = [
            "move_to",
            "square_distance_from",
            "distance_from",
            "norm",
            "rotate",
            "dot",
            "cross",
            "get_raw_data",
            "dihedral",
        ]

        for prop_name in expected_properties:
            prop = getattr(Position, prop_name)
            doc = prop.__doc__
            self.assertIsNotNone(
                doc, f"Position.{prop_name} (property) docstring is None"
            )
            self.assertTrue(
                isinstance(doc, str) and len(doc.strip()) > 0,
                f"Position.{prop_name} (property) docstring is empty",
            )

        for method_name in expected_methods:
            method = getattr(Position, method_name)
            doc = method.__doc__
            self.assertIsNotNone(
                doc, f"Position.{method_name} (method) docstring is None"
            )
            self.assertTrue(
                isinstance(doc, str) and len(doc.strip()) > 0,
                f"Position.{method_name} (method) docstring is empty",
            )

        # Ensure every public attribute on Position has a non-empty docstring
        for name in dir(Position):
            if name.startswith("_"):
                continue
            attr = getattr(Position, name)
            doc = attr.__doc__
            self.assertIsNotNone(
                doc, f"Public attribute Position.{name} docstring is None"
            )
            self.assertTrue(
                isinstance(doc, str) and len(doc.strip()) > 0,
                f"Public attribute Position.{name} docstring is empty",
            )

    def test_atom_class_docstring(self):
        """Verify that Atom class has a non-empty docstring."""
        doc = Atom.__doc__
        self.assertIsNotNone(doc, "Atom.__doc__ is None")
        self.assertTrue(
            isinstance(doc, str) and len(doc.strip()) > 0,
            "Atom.__doc__ is empty",
        )

    def test_atom_public_members_docstrings(self):
        """Verify that all public methods and properties of Atom have non-empty docstrings."""
        expected_properties = [
            "atomic_number",
            "symbol",
            "name",
            "label",
            "charge",
            "path",
            "xyz",
            "position",
            "force",
            "is_real",
            "vdw",
        ]
        expected_methods = [
            "weight",
            "move_to",
            "shift_by",
            "rotate",
            "get_raw_data",
            "set_by_raw_data",
        ]

        for prop_name in expected_properties:
            prop = getattr(Atom, prop_name)
            doc = prop.__doc__
            self.assertIsNotNone(
                doc, f"Atom.{prop_name} (property) docstring is None"
            )
            self.assertTrue(
                isinstance(doc, str) and len(doc.strip()) > 0,
                f"Atom.{prop_name} (property) docstring is empty",
            )

        for method_name in expected_methods:
            method = getattr(Atom, method_name)
            doc = method.__doc__
            self.assertIsNotNone(
                doc, f"Atom.{method_name} (method) docstring is None"
            )
            self.assertTrue(
                isinstance(doc, str) and len(doc.strip()) > 0,
                f"Atom.{method_name} (method) docstring is empty",
            )

        # Ensure every public attribute on Atom has a non-empty docstring
        for name in dir(Atom):
            if name.startswith("_"):
                continue
            attr = getattr(Atom, name)
            doc = attr.__doc__
            self.assertIsNotNone(
                doc, f"Public attribute Atom.{name} docstring is None"
            )
            self.assertTrue(
                isinstance(doc, str) and len(doc.strip()) > 0,
                f"Public attribute Atom.{name} docstring is empty",
            )


if __name__ == "__main__":
    unittest.main()
