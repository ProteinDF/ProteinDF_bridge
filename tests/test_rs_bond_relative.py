#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for relative bond paths in rs.AtomGroup."""

import os
import tempfile
import unittest

import proteindf_bridge.rs as rs
from proteindf_bridge.functions import load_atomgroup, save_atomgroup

DATA_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "proteindf_bridge", "data")
)


class TestRsBondRelative(unittest.TestCase):
    def test_add_bond_sibling_20_2(self):
        """Criterion 1 (a): Sibling 20 and 2; add_bond on mol['2'] atoms stores in '2' only."""
        mol = rs.AtomGroup()
        grp20 = rs.AtomGroup()
        grp20.set_atom("N", rs.Atom(symbol="N"))
        grp20.set_atom("CA", rs.Atom(symbol="C"))
        grp2 = rs.AtomGroup()
        grp2.set_atom("N", rs.Atom(symbol="N"))
        grp2.set_atom("CA", rs.Atom(symbol="C"))

        mol.set_group("20", grp20)
        mol.set_group("2", grp2)

        mol.add_bond(mol["2"]["N"], mol["2"]["CA"])
        self.assertEqual(len(mol["2"].bonds()), 1)
        self.assertEqual(len(mol["20"].bonds()), 0)
        self.assertEqual(len(mol.bonds()), 0)
        self.assertEqual(len(mol["2"].get_bond_list()), 1)
        self.assertEqual(len(mol["20"].get_bond_list()), 0)
        self.assertEqual(len(mol.get_bond_list()), 1)

    def test_add_bond_sibling_2_20(self):
        """Criterion 1 (b): Sibling 2 and 20; add_bond on mol['2'] atoms stores in '2' only."""
        mol = rs.AtomGroup()
        grp2 = rs.AtomGroup()
        grp2.set_atom("N", rs.Atom(symbol="N"))
        grp2.set_atom("CA", rs.Atom(symbol="C"))
        grp20 = rs.AtomGroup()
        grp20.set_atom("N", rs.Atom(symbol="N"))
        grp20.set_atom("CA", rs.Atom(symbol="C"))

        mol.set_group("2", grp2)
        mol.set_group("20", grp20)

        mol.add_bond(mol["2"]["N"], mol["2"]["CA"])
        self.assertEqual(len(mol["2"].bonds()), 1)
        self.assertEqual(len(mol["20"].bonds()), 0)
        self.assertEqual(len(mol.bonds()), 0)
        self.assertEqual(len(mol["2"].get_bond_list()), 1)
        self.assertEqual(len(mol["20"].get_bond_list()), 0)
        self.assertEqual(len(mol.get_bond_list()), 1)

    def test_unparented_group_add_bond_then_set_group(self):
        """Criterion 2: Adding a bond to an unparented group, then setting it into a parent tree.

        get_bond_list() returns paths with the new parent prefix.
        """
        model = rs.AtomGroup()
        res = rs.AtomGroup()
        res.set_atom("N", rs.Atom(symbol="N"))
        res.set_atom("CA", rs.Atom(symbol="C"))
        res.add_bond(res["N"], res["CA"])

        # Relative bond in res
        self.assertEqual(res.bonds(), [("N", "CA", 1)])
        # get_bond_list before parenting
        self.assertEqual(res.get_bond_list(), [("/N", "/CA", 1)])

        # Insert res into model
        model.set_group("res1", res)

        # get_bond_list from inside the tree
        self.assertEqual(model["res1"].get_bond_list(), [("/res1/N", "/res1/CA", 1)])
        self.assertEqual(model.get_bond_list(), [("/res1/N", "/res1/CA", 1)])

        # Insert model into a root
        root = rs.AtomGroup()
        root.set_group("model_1", model)
        self.assertEqual(
            root["model_1"]["res1"].get_bond_list(),
            [("/model_1/res1/N", "/model_1/res1/CA", 1)],
        )
        self.assertEqual(
            root["model_1"].get_bond_list(),
            [("/model_1/res1/N", "/model_1/res1/CA", 1)],
        )
        self.assertEqual(
            root.get_bond_list(),
            [("/model_1/res1/N", "/model_1/res1/CA", 1)],
        )

    def test_2mgo_pdb_bonds_and_brd_roundtrip(self):
        """Criterion 3: 2MGO.pdb bonds match expectations and survive brd save/load."""
        pdb_path = os.path.join(DATA_DIR, "2MGO.pdb")
        self.assertTrue(os.path.exists(pdb_path), f"File not found: {pdb_path}")

        pdb = rs.Pdb(pdb_path)
        ag = pdb.get_atomgroup()
        bonds_orig = ag.get_bond_list()

        # 2MGO has 20 models, each having 1 SSBOND (between 1 CYS 6_SG and 6 CYS 89_SG) -> 20 bonds
        self.assertEqual(len(bonds_orig), 20)
        self.assertEqual(
            bonds_orig[0],
            ("/model_1/A/1/6_SG", "/model_1/A/6/89_SG", 1),
        )

        from proteindf_bridge.biopdb import Pdb as PyPdb
        py_pdb = PyPdb(pdb_path)
        py_bonds = [tuple(b) for b in py_pdb.get_atomgroup().get_bond_list()]
        self.assertEqual(bonds_orig, py_bonds)

        with tempfile.TemporaryDirectory() as tmpdir:
            brd_path = os.path.join(tmpdir, "2MGO.brd")

            # Save via proteindf_bridge.functions.save_atomgroup
            save_atomgroup(ag, brd_path)
            self.assertTrue(os.path.exists(brd_path))

            # Reload via rs.load_atomgroup
            reloaded_rs = rs.load_atomgroup(brd_path)
            self.assertEqual(reloaded_rs.get_bond_list(), bonds_orig)

            # Reload via proteindf_bridge.load_atomgroup
            reloaded_py = load_atomgroup(brd_path)
            self.assertEqual(
                [tuple(b) for b in reloaded_py.get_bond_list()],
                bonds_orig,
            )


if __name__ == "__main__":
    unittest.main()
