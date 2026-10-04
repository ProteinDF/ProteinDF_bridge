#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test suite validating that all Rust pyclasses and exception types exposed by
proteindf_bridge.rs have their __module__ attribute correctly set to
'proteindf_bridge.rs'.
"""

import unittest
import proteindf_bridge.rs as rs


class TestRsModuleName(unittest.TestCase):
    """Verify that all classes and exceptions in proteindf_bridge.rs have __module__ == 'proteindf_bridge.rs'."""

    def test_all_classes_and_exceptions_module_name(self):
        """Enumerate all type objects in proteindf_bridge.rs and verify their __module__."""
        types_checked = {}
        for name in dir(rs):
            if name.startswith("_"):
                continue
            obj = getattr(rs, name)
            if isinstance(obj, type):
                types_checked[name] = obj
                self.assertEqual(
                    obj.__module__,
                    "proteindf_bridge.rs",
                    f"Type {name} has incorrect __module__: {obj.__module__!r} (expected 'proteindf_bridge.rs')",
                )

        # Expected classes across various modules (PR#14-16, PR#44-48)
        expected_types = [
            # Exceptions
            "BrError",
            "BrInputError",
            "BrValueError",
            # Foundation & Data model
            "PeriodicTable",
            "Vector",
            "Matrix",
            "SymmetricMatrix",
            "Position",
            "Atom",
            "Bond",
            "AtomGroup",
            # Format I/O
            "Format",
            "Xyz",
            "SimpleGro",
            "SimpleMol2",
            "AmberPrmtop",
            "Pdb",
            "SimpleMmcif",
            "MmcifStructureReport",
            "StructConnPartnerUnresolved",
            "UnresolvedStructConn",
            # Structural operations
            "AminoAcid",
            "SSBond",
            "IonPair",
            "Superposer",
            "SuperposerQuaternion",
            # Selectors
            "SelectAtom",
            "SelectAtomGroup",
            "SelectName",
            "SelectSymbol",
            "SelectPath",
            "SelectPathSimple",
            "SelectPathWildcard",
            "SelectPathRegex",
            "SelectRange",
            # Analysis & Modeling (PR#45-48)
            "CcdAtom",
            "CcdBondTemplate",
            "CcdTemplateDb",
            "SchemaViolation",
            "SecondaryStructure",
            "HydrogenationReport",
            "OverallHydrogenationReport",
            "HydrogenBond",
            "SidechainHydrogenBond",
            "ChPiInteraction",
            "Interaction",
            "InteractionSet",
            "Modeling",
            "Neutralize",
            "RamachandranAngle",
        ]

        for exp in expected_types:
            self.assertIn(
                exp,
                types_checked,
                f"Expected type {exp} was not found in proteindf_bridge.rs",
            )

        # Ensure a comprehensive set of types was checked (currently 60 types)
        self.assertGreaterEqual(
            len(types_checked),
            50,
            f"Expected at least 50 types to be checked, found {len(types_checked)}",
        )


if __name__ == "__main__":
    unittest.main()
