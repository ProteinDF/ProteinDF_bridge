#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2014 The ProteinDF development team.
# see also AUTHORS and README if provided.
#
# This file is a part of the ProteinDF software package.
#
# The ProteinDF is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The ProteinDF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

"""Unit tests for unified exception classes between Rust (rs) and Python."""

import pickle
import unittest

import proteindf_bridge
import proteindf_bridge.error as br_error
import proteindf_bridge.rs as rs_br


class TestRsErrorUnify(unittest.TestCase):
    """Test unification of BrError, BrInputError, and BrValueError across rs and pure Python."""

    def test_class_identity(self):
        """Verify that proteindf_bridge.error classes are identical to rs classes."""
        self.assertIs(proteindf_bridge.BrError, rs_br.BrError)
        self.assertIs(proteindf_bridge.BrInputError, rs_br.BrInputError)
        self.assertIs(proteindf_bridge.BrValueError, rs_br.BrValueError)

        self.assertIs(br_error.BrError, rs_br.BrError)
        self.assertIs(br_error.BrInputError, rs_br.BrInputError)
        self.assertIs(br_error.BrValueError, rs_br.BrValueError)

    def test_inheritance_hierarchy(self):
        """Verify inheritance hierarchy: BrInputError and BrValueError subclass BrError, which subclasses Exception."""
        self.assertTrue(issubclass(rs_br.BrError, Exception))
        self.assertTrue(issubclass(rs_br.BrInputError, rs_br.BrError))
        self.assertTrue(issubclass(rs_br.BrValueError, rs_br.BrError))

        self.assertTrue(issubclass(br_error.BrInputError, br_error.BrError))
        self.assertTrue(issubclass(br_error.BrValueError, br_error.BrError))

    def test_format_and_errmsg(self):
        """Verify constructor arguments, errmsg attribute, and str() formatting."""
        # BrError() - default
        e1 = rs_br.BrError()
        self.assertEqual(e1.errmsg, "")
        self.assertEqual(str(e1), "Bridge module error: ")

        # BrError(msg)
        e2 = rs_br.BrError("some general error")
        self.assertEqual(e2.errmsg, "some general error")
        self.assertEqual(str(e2), "Bridge module error: some general error")

        # BrInputError(expr, msg)
        e3 = rs_br.BrInputError("f()", "m")
        self.assertEqual(e3.errmsg, "Input Error: m (f())")
        self.assertEqual(str(e3), "Bridge module error: Input Error: m (f())")

        # BrValueError(expr, msg)
        e4 = rs_br.BrValueError("Position.norm()", "cannot normalize zero vector")
        self.assertEqual(e4.errmsg, "Value Error: cannot normalize zero vector (Position.norm())")
        self.assertEqual(
            str(e4),
            "Bridge module error: Value Error: cannot normalize zero vector (Position.norm())",
        )

    def test_catch_pure_python_exception(self):
        """Verify that pure Python-raised exceptions can be caught by rs or error exception types."""
        # Function raising BrValueError as pure Python modules do
        def raise_val_err():
            raise br_error.BrValueError("Position.norm()", "test error")

        # Catch by rs_br.BrValueError
        caught = False
        try:
            raise_val_err()
        except rs_br.BrValueError as e:
            caught = True
            self.assertEqual(e.errmsg, "Value Error: test error (Position.norm())")
            self.assertEqual(
                str(e),
                "Bridge module error: Value Error: test error (Position.norm())",
            )
        self.assertTrue(caught)

        # Catch by br_error.BrError (base class)
        caught = False
        try:
            raise_val_err()
        except br_error.BrError as e:
            caught = True
            self.assertEqual(e.errmsg, "Value Error: test error (Position.norm())")
        self.assertTrue(caught)

    def test_catch_rust_exception(self):
        """Verify that Rust-raised exceptions can be caught by proteindf_bridge.error.BrError and have correct format."""
        # Non-existent file via rs.load_atomgroup
        caught = False
        try:
            rs_br.load_atomgroup("__nonexistent_file_xyz_12345.brd__")
        except br_error.BrError as e:
            caught = True
            self.assertTrue(hasattr(e, "errmsg"))
            self.assertIn("I/O error", e.errmsg)
            self.assertEqual(str(e), f"Bridge module error: {e.errmsg}")
        self.assertTrue(caught)

        # PeriodicTable error via rs
        caught = False
        try:
            rs_br.PeriodicTable.get_symbol(9999)
        except br_error.BrError as e:
            caught = True
            self.assertTrue(hasattr(e, "errmsg"))
            self.assertIn("atomic number 9999 not found", e.errmsg)
            self.assertEqual(str(e), f"Bridge module error: {e.errmsg}")
            self.assertIsInstance(e, rs_br.BrValueError)
        self.assertTrue(caught)

    def test_pickle_roundtrip(self):
        """Verify that exceptions can be pickled and unpickled with class and errmsg preserved."""
        test_exceptions = [
            rs_br.BrError(),
            rs_br.BrError("m"),
            rs_br.BrInputError("f()", "m"),
            rs_br.BrValueError("Position.norm()", "zero vector"),
            br_error.BrError("from_error_module"),
            br_error.BrInputError("expr_test", "msg_test"),
            br_error.BrValueError("val_expr", "val_msg"),
        ]

        for orig in test_exceptions:
            serialized = pickle.dumps(orig)
            restored = pickle.loads(serialized)

            self.assertIs(type(restored), type(orig))
            self.assertEqual(restored.errmsg, orig.errmsg)
            self.assertEqual(str(restored), str(orig))


if __name__ == "__main__":
    unittest.main()
