#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test suite validating subclassing support for rs.Atom in Python.
"""

import unittest

import proteindf_bridge.rs as rs


class TestRsAtomSubclass(unittest.TestCase):
    """Test suite for Python subclassing of rs.Atom."""

    def test_subclass_definition(self):
        """Verify that a Python class can inherit from rs.Atom without TypeError."""
        class SubAtom(rs.Atom):
            pass

        self.assertTrue(issubclass(SubAtom, rs.Atom))

    def test_subclass_instantiation_and_methods(self):
        """Verify that a subclass instance can be initialized with standard Atom arguments and use methods."""
        class SubAtom(rs.Atom):
            pass

        sub = SubAtom("C", name="CA", xyz=[1.0, 2.0, 3.0], charge=0.0)
        self.assertIsInstance(sub, SubAtom)
        self.assertIsInstance(sub, rs.Atom)
        self.assertEqual(sub.symbol, "C")
        self.assertEqual(sub.atomic_number, 6)
        self.assertEqual(sub.name, "CA")
        self.assertEqual(sub.xyz.x, 1.0)
        self.assertEqual(sub.xyz.y, 2.0)
        self.assertEqual(sub.xyz.z, 3.0)
        self.assertEqual(sub.charge, 0.0)

        # Verify attributes and methods on subclass instance
        self.assertTrue(sub.is_real)
        self.assertAlmostEqual(sub.vdw, 1.70, places=2)
        self.assertAlmostEqual(sub.weight(), 12.011, places=2)

        sub.shift_by([1.0, -1.0, 0.5])
        self.assertAlmostEqual(sub.xyz.x, 2.0)
        self.assertAlmostEqual(sub.xyz.y, 1.0)
        self.assertAlmostEqual(sub.xyz.z, 3.5)

    def test_subclass_init_super_without_args(self):
        """Verify that calling super().__init__() without arguments succeeds in a subclass."""
        class SubAtomWithInit(rs.Atom):
            def __init__(self, *args, **kwargs):
                super().__init__()
                self.initialized = True

        sub = SubAtomWithInit("N", name="N", xyz=[0.0, 0.0, 0.0])
        self.assertEqual(sub.symbol, "N")
        self.assertEqual(sub.name, "N")
        self.assertTrue(sub.initialized)

    def test_subclass_with_custom_init_and_new(self):
        """Verify subclass instantiation accepting custom arguments via super().__new__ and default __new__."""
        # Case 1: Custom __new__ delegating base arguments to super().__new__
        class CustomNewAtom(rs.Atom):
            def __new__(cls, *args, tag="default_tag", **kwargs):
                instance = super().__new__(cls, *args, **kwargs)
                return instance

            def __init__(self, *args, tag="default_tag", **kwargs):
                super().__init__()
                self.tag = tag

        obj1 = CustomNewAtom("C", name="C1", xyz=[0.0, 1.0, 2.0], tag="custom_tag")
        self.assertEqual(obj1.symbol, "C")
        self.assertEqual(obj1.name, "C1")
        self.assertEqual(obj1.tag, "custom_tag")
        self.assertEqual(obj1.xyz.y, 1.0)

        # Case 2: Subclass without overriding __new__, receiving standard Atom arguments in __init__
        class CustomInitAtom(rs.Atom):
            def __init__(self, symbol, name=None, **kwargs):
                super().__init__()
                self.label_info = f"{symbol}-{name}"

        obj2 = CustomInitAtom("O", name="OG")
        self.assertEqual(obj2.symbol, "O")
        self.assertEqual(obj2.name, "OG")
        self.assertEqual(obj2.label_info, "O-OG")
