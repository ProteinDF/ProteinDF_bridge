#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Comparison test suite for Phase 5 (PR#16):
Validating proteindf_bridge.rs (PyO3 Rust bindings) structural operations against pure Python proteindf_bridge.
"""

import math
import os
import unittest
import numpy as np

# Pure Python implementation
from proteindf_bridge.aminoacid import AminoAcid as PyAminoAcid
from proteindf_bridge.biopdb import Pdb as PyPdb
from proteindf_bridge.ionpair import IonPair as PyIonPair
from proteindf_bridge.select import (
    Select_Name as PySelect_Name,
    Select_Path_simple as PySelect_Path_simple,
    Select_Path_wildcard as PySelect_Path_wildcard,
    Select_PathRegex as PySelect_PathRegex,
    Select_Range as PySelect_Range,
    Select_Symbol as PySelect_Symbol,
)
from proteindf_bridge.ssbond import SSBond as PySSBond
from proteindf_bridge.superposer import Superposer as PySuperposer

# Rust PyO3 bindings
import proteindf_bridge.rs as rs_br

DATA_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "proteindf_bridge", "data")
)


class TestRsPhase3(unittest.TestCase):

    # --------------------------------------------------------------------------
    # 1. AminoAcid
    # --------------------------------------------------------------------------
    def test_amino_acid(self):
        # Recognizable amino acids
        for aa in ["ALA", "CYS", "CYX", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "PRO", "TYR", "VAL"]:
            ag_rs = rs_br.AtomGroup()
            ag_rs.name = aa
            self.assertTrue(rs_br.AminoAcid.is_aminoacid(ag_rs))
            self.assertEqual(
                rs_br.AminoAcid.is_aminoacid(ag_rs),
                PyAminoAcid.is_aminoacid(ag_rs),
            )
            # Direct name check
            self.assertTrue(rs_br.AminoAcid.is_aminoacid_name(aa))

        # Non-amino acids
        for other in ["HOH", "WAT", "NA", "CL", "LIG", "UNK"]:
            ag_rs = rs_br.AtomGroup()
            ag_rs.name = other
            self.assertFalse(rs_br.AminoAcid.is_aminoacid(ag_rs))
            self.assertFalse(rs_br.AminoAcid.is_aminoacid_name(other))

        # aa_list
        self.assertEqual(rs_br.AminoAcid.aa_list(), PyAminoAcid._AA_list)
        self.assertEqual(rs_br.AminoAcid._AA_list, PyAminoAcid._AA_list)

    # --------------------------------------------------------------------------
    # 2. SSBond
    # --------------------------------------------------------------------------
    def test_ssbond_1hls(self):
        pdb_path = os.path.join(DATA_DIR, "1hls.pdb")
        rs_pdb = rs_br.Pdb(pdb_path)
        py_pdb = PyPdb(pdb_path)

        rs_models = rs_pdb.get_atomgroup()
        py_models = py_pdb.get_atomgroup()

        rs_model = rs_models["model_1"]
        py_model = py_models["model_1"]

        rs_ssb = rs_br.SSBond(rs_model)
        py_ssb = PySSBond(py_model)

        rs_bonds = rs_ssb.get_bonds()
        py_bonds = py_ssb.get_bonds()

        expected = [
            ("/model_1/A/6/", "/model_1/A/11/"),
            ("/model_1/A/7/", "/model_1/B/7/"),
            ("/model_1/A/20/", "/model_1/B/19/"),
        ]
        self.assertEqual(rs_bonds, expected)
        self.assertEqual(rs_bonds, py_bonds)

        # find_bonds static method
        self.assertEqual(rs_br.SSBond.find_bonds(rs_model), expected)

    # --------------------------------------------------------------------------
    # 3. IonPair
    # --------------------------------------------------------------------------
    def test_ion_pair_synthetic(self):
        protein_rs = rs_br.AtomGroup("protein")
        chainA_rs = rs_br.AtomGroup("A")

        # GLU: CD, OE1, OE2
        glu_rs = rs_br.AtomGroup("1")
        glu_rs.name = "GLU"
        glu_rs.set_atom("CD", rs_br.Atom(name="CD", xyz=rs_br.Position(0.0, 0.0, 0.0)))
        glu_rs.set_atom("OE1", rs_br.Atom(name="OE1", xyz=rs_br.Position(0.0, 1.0, 0.0)))
        glu_rs.set_atom("OE2", rs_br.Atom(name="OE2", xyz=rs_br.Position(1.0, 0.0, 0.0)))
        chainA_rs.set_group("1", glu_rs)

        # LYS: NZ
        lys_rs = rs_br.AtomGroup("2")
        lys_rs.name = "LYS"
        lys_rs.set_atom("NZ", rs_br.Atom(name="NZ", xyz=rs_br.Position(0.5, 2.0, 0.0)))
        chainA_rs.set_group("2", lys_rs)

        protein_rs.set_group("A", chainA_rs)

        ip_rs = rs_br.IonPair(protein_rs)
        ip_py = PyIonPair(protein_rs)

        pairs_rs = ip_rs.get_ion_pairs()
        pairs_py = ip_py.get_ion_pairs()

        expected = [("/A/1/", "/A/2/", "GLU", "LYS")]
        self.assertEqual(pairs_rs, expected)
        self.assertEqual(pairs_py, expected)

        # Static find_ion_pairs
        self.assertEqual(rs_br.IonPair.find_ion_pairs(protein_rs), expected)

    def test_ion_pair_1hls(self):
        pdb_path = os.path.join(DATA_DIR, "1hls.pdb")
        rs_pdb = rs_br.Pdb(pdb_path)
        py_pdb = PyPdb(pdb_path)

        rs_model = rs_pdb.get_atomgroup()["model_1"]
        py_model = py_pdb.get_atomgroup()["model_1"]

        pairs_rs = rs_br.IonPair(rs_model).get_ion_pairs()
        pairs_py = PyIonPair(py_model).get_ion_pairs()

        self.assertEqual(len(pairs_rs), len(pairs_py))
        self.assertEqual(set(pairs_rs), set(pairs_py))

    # --------------------------------------------------------------------------
    # 4. Selector
    # --------------------------------------------------------------------------
    def _create_selector_test_group(self):
        # group with sub1 and sub2
        rs_root = rs_br.AtomGroup()
        sub1 = rs_br.AtomGroup()
        sub1.set_atom("C1", rs_br.Atom(name="C1", symbol="C", xyz=rs_br.Position(0.0, 0.0, 0.0)))
        sub1.set_atom("H1", rs_br.Atom(name="H1", symbol="H", xyz=rs_br.Position(1.0, 0.0, 0.0)))
        sub1.set_atom("N1", rs_br.Atom(name="N1", symbol="N", xyz=rs_br.Position(2.0, 0.0, 0.0)))
        rs_root.set_group("sub1", sub1)

        sub2 = rs_br.AtomGroup()
        sub2.set_atom("C1", rs_br.Atom(name="C1", symbol="C", xyz=rs_br.Position(0.0, 1.0, 0.0)))
        sub2.set_atom("H1", rs_br.Atom(name="H1", symbol="H", xyz=rs_br.Position(0.0, 2.0, 0.0)))
        sub2.set_atom("N1", rs_br.Atom(name="N1", symbol="N", xyz=rs_br.Position(0.0, 3.0, 0.0)))
        rs_root.set_group("sub2", sub2)

        return rs_root

    def test_select_range(self):
        rs_root = self._create_selector_test_group()

        sel_rs = rs_br.Select_Range("0.0 0.0 0.0", 0.1)
        sel_py = PySelect_Range("0.0 0.0 0.0", 0.1)

        res_rs = rs_root.select(sel_rs)
        res_py = rs_root.select(sel_py)

        self.assertEqual(res_rs.get_number_of_all_atoms(), 1)
        self.assertEqual(res_py.get_number_of_all_atoms(), 1)
        self.assertEqual(res_rs["sub1"]["C1"].symbol, "C")
        self.assertEqual(res_py["sub1"]["C1"].symbol, "C")

        # Range 1.1
        sel_rs2 = rs_br.Select_Range("0.0 0.0 0.0", 1.1)
        sel_py2 = PySelect_Range("0.0 0.0 0.0", 1.1)
        res_rs2 = rs_root.select(sel_rs2)
        res_py2 = rs_root.select(sel_py2)

        self.assertEqual(res_rs2.get_number_of_all_atoms(), 3)
        self.assertEqual(res_py2.get_number_of_all_atoms(), 3)

    def test_select_symbol(self):
        rs_root = self._create_selector_test_group()

        sel_rs = rs_br.Select_Symbol("c")
        sel_py = PySelect_Symbol("c")

        res_rs = rs_root.select(sel_rs)
        res_py = rs_root.select(sel_py)

        self.assertEqual(res_rs.get_number_of_all_atoms(), 2)
        self.assertEqual(res_py.get_number_of_all_atoms(), 2)
        self.assertEqual(res_rs["sub1"]["C1"].symbol, "C")
        self.assertEqual(res_rs["sub2"]["C1"].symbol, "C")
        self.assertEqual(res_py["sub1"]["C1"].symbol, "C")
        self.assertEqual(res_py["sub2"]["C1"].symbol, "C")

    def test_select_name(self):
        rs_root = self._create_selector_test_group()

        sel_rs = rs_br.Select_Name("C1")
        sel_py = PySelect_Name("C1")

        res_rs = rs_root.select(sel_rs)
        res_py = rs_root.select(sel_py)

        self.assertEqual(res_rs.get_number_of_all_atoms(), 2)
        self.assertEqual(res_py.get_number_of_all_atoms(), 2)

    def test_select_path_simple(self):
        rs_root = self._create_selector_test_group()

        sel_rs = rs_br.Select_Path_simple("/sub1/C1")
        sel_py = PySelect_Path_simple("/sub1/C1")

        res_rs = rs_root.select(sel_rs)
        res_py = rs_root.select(sel_py)

        self.assertEqual(res_rs.get_number_of_all_atoms(), 1)
        self.assertEqual(res_py.get_number_of_all_atoms(), 1)
        self.assertEqual(res_rs["sub1"]["C1"].path, "/sub1/C1")
        self.assertEqual(res_py["sub1"]["C1"].path, "/sub1/C1")

    def test_select_path_wildcard(self):
        rs_root = self._create_selector_test_group()

        sel_rs = rs_br.Select_Path_wildcard("/sub1/*")
        sel_py = PySelect_Path_wildcard("/sub1/*")

        res_rs = rs_root.select(sel_rs)
        res_py = rs_root.select(sel_py)

        self.assertEqual(res_rs.get_number_of_all_atoms(), 3)
        self.assertEqual(res_py.get_number_of_all_atoms(), 3)

    def test_select_path_regex(self):
        rs_root = self._create_selector_test_group()

        sel_rs = rs_br.Select_PathRegex(r"/sub\d/C1")
        sel_py = PySelect_PathRegex(r"/sub\d/C1")

        res_rs = rs_root.select(sel_rs)
        res_py = rs_root.select(sel_py)

        self.assertEqual(res_rs.get_number_of_all_atoms(), 2)
        self.assertEqual(res_py.get_number_of_all_atoms(), 2)

    def test_select_atom_and_atomgroup(self):
        rs_root = self._create_selector_test_group()

        ref_atom = rs_br.Atom(symbol="C", xyz=rs_br.Position(0.0, 0.05, 0.0))
        sel_atom = rs_br.Select_Atom(ref_atom, 0.1)
        res = rs_root.select(sel_atom)
        self.assertEqual(res.get_number_of_all_atoms(), 1)
        self.assertEqual(res["sub1"]["C1"].symbol, "C")

        sel_ag = rs_br.Select_AtomGroup(res, 0.01)
        res2 = rs_root.select(sel_ag)
        self.assertEqual(res2.get_number_of_all_atoms(), 1)

    # --------------------------------------------------------------------------
    # 5. Superposer (Kabsch)
    # --------------------------------------------------------------------------
    def test_superposer(self):
        ag1 = rs_br.AtomGroup("mol1")
        ag1.set_atom("A1", rs_br.Atom(name="A1", xyz=rs_br.Position(1.0, 1.0, 1.0)))
        ag1.set_atom("A2", rs_br.Atom(name="A2", xyz=rs_br.Position(1.0, -1.0, -1.0)))
        ag1.set_atom("A3", rs_br.Atom(name="A3", xyz=rs_br.Position([-1.0, 1.0, -1.0])))
        ag1.set_atom("A4", rs_br.Atom(name="A4", xyz=rs_br.Position([-1.0, -1.0, 1.0])))

        # Shifted by (1.0, 2.0, 3.0)
        ag2 = rs_br.AtomGroup("mol2")
        ag2.set_atom("A1", rs_br.Atom(name="A1", xyz=rs_br.Position(2.0, 3.0, 4.0)))
        ag2.set_atom("A2", rs_br.Atom(name="A2", xyz=rs_br.Position(2.0, 1.0, 2.0)))
        ag2.set_atom("A3", rs_br.Atom(name="A3", xyz=rs_br.Position(0.0, 3.0, 2.0)))
        ag2.set_atom("A4", rs_br.Atom(name="A4", xyz=rs_br.Position(0.0, 1.0, 4.0)))

        sp_rs = rs_br.Superposer(ag1, ag2)
        sp_py = PySuperposer(ag1, ag2)

        self.assertAlmostEqual(sp_rs.rmsd, 0.0, places=5)
        self.assertAlmostEqual(sp_py.rmsd, 0.0, places=5)
        self.assertEqual(sp_rs.num_of_positions, 4)
        self.assertEqual(sp_py.num_of_positions, 4)

        superimposed = sp_rs.superimpose(ag1)
        self.assertEqual(superimposed.get_number_of_atoms(), 4)
        for key in ["A1", "A2", "A3", "A4"]:
            p1 = superimposed.get_atom(key).xyz
            p2 = ag2.get_atom(key).xyz
            self.assertAlmostEqual(p1.distance_from(p2), 0.0, places=4)

    # --------------------------------------------------------------------------
    # 6. Superposer_quaternion
    # --------------------------------------------------------------------------
    def test_superposer_quaternion(self):
        ag1 = rs_br.AtomGroup("mol1")
        ag1.set_atom("A1", rs_br.Atom(name="A1", xyz=rs_br.Position(1.0, 1.0, 1.0)))
        ag1.set_atom("A2", rs_br.Atom(name="A2", xyz=rs_br.Position(1.0, -1.0, -1.0)))
        ag1.set_atom("A3", rs_br.Atom(name="A3", xyz=rs_br.Position(-1.0, 1.0, -1.0)))
        ag1.set_atom("A4", rs_br.Atom(name="A4", xyz=rs_br.Position(-1.0, -1.0, 1.0)))

        # Shifted by (1.0, 2.0, 3.0)
        ag2 = rs_br.AtomGroup("mol2")
        ag2.set_atom("A1", rs_br.Atom(name="A1", xyz=rs_br.Position(2.0, 3.0, 4.0)))
        ag2.set_atom("A2", rs_br.Atom(name="A2", xyz=rs_br.Position(2.0, 1.0, 2.0)))
        ag2.set_atom("A3", rs_br.Atom(name="A3", xyz=rs_br.Position(0.0, 3.0, 2.0)))
        ag2.set_atom("A4", rs_br.Atom(name="A4", xyz=rs_br.Position(0.0, 1.0, 4.0)))

        sq_rs = rs_br.Superposer_quaternion(ag1, ag2)
        rmsd = sq_rs.rmsd
        self.assertAlmostEqual(rmsd, 0.0, places=5)
        self.assertAlmostEqual(sq_rs.calc(), rmsd)

        rot_mat = sq_rs.rotation_mat
        self.assertEqual(rot_mat.rows, 3)
        self.assertEqual(rot_mat.cols, 3)

        superimposed = sq_rs.superimpose(ag1)
        self.assertEqual(superimposed.get_number_of_atoms(), 4)
        for key in ["A1", "A2", "A3", "A4"]:
            p1 = superimposed.get_atom(key).xyz
            p2 = ag2.get_atom(key).xyz
            self.assertAlmostEqual(p1.distance_from(p2), 0.0, places=4)

    def test_quaternion_with_arbitrary_axis_rotation(self):
        # Rotation around axis (1,1,1)/sqrt(3) by 37 degrees + translation (2, 3, 4)
        axis = np.array([1.0, 1.0, 1.0]) / math.sqrt(3.0)
        theta = math.radians(37.0)
        c = math.cos(theta)
        s = math.sin(theta)
        C = 1.0 - c
        ux, uy, uz = axis

        rot = np.array([
            [ux * ux * C + c, ux * uy * C - uz * s, ux * uz * C + uy * s],
            [uy * ux * C + uz * s, uy * uy * C + c, uy * uz * C - ux * s],
            [uz * ux * C - uy * s, uz * uy * C + ux * s, uz * uz * C + c],
        ])
        trans = np.array([2.0, 3.0, 4.0])

        coords1 = [
            np.array([1.0, 1.0, 1.0]),
            np.array([1.0, -1.0, -1.0]),
            np.array([-1.0, 1.0, -1.0]),
            np.array([-1.0, -1.0, 1.0]),
        ]

        ag1 = rs_br.AtomGroup("mol1")
        ag2 = rs_br.AtomGroup("mol2")
        for i, p in enumerate(coords1):
            key = f"A{i + 1}"
            ag1.set_atom(key, rs_br.Atom(name=key, xyz=rs_br.Position(p[0], p[1], p[2])))
            p_rot = rot @ p + trans
            ag2.set_atom(key, rs_br.Atom(name=key, xyz=rs_br.Position(p_rot[0], p_rot[1], p_rot[2])))

        sq = rs_br.Superposer_quaternion(ag1, ag2)
        sp = rs_br.Superposer(ag1, ag2)

        # Both Kabsch and quaternion methods should give near-zero RMSD
        self.assertAlmostEqual(sq.rmsd, 0.0, places=4)
        self.assertAlmostEqual(sp.rmsd, 0.0, places=4)
        self.assertAlmostEqual(sq.rmsd, sp.rmsd, places=4)

        superimposed = sq.superimpose(ag1)
        for key in ["A1", "A2", "A3", "A4"]:
            p1 = superimposed.get_atom(key).xyz
            p2 = ag2.get_atom(key).xyz
            self.assertAlmostEqual(p1.distance_from(p2), 0.0, places=4)

    def test_custom_selector_exception_propagation(self):
        class FaultySelector:
            def is_match(self, obj):
                raise ValueError("custom selector exception test")

        ag = rs_br.AtomGroup("test")
        ag.set_atom("C1", rs_br.Atom(name="C1", symbol="C", xyz=rs_br.Position(0.0, 0.0, 0.0)))

        with self.assertRaises(ValueError) as ctx:
            ag.select(FaultySelector())

        self.assertIn("custom selector exception test", str(ctx.exception))

    def test_custom_selector_success(self):
        class CarbonSelector:
            def is_match(self, obj):
                if hasattr(obj, "symbol"):
                    return obj.symbol == "C"
                return False

        ag = rs_br.AtomGroup("test")
        ag.set_atom("C1", rs_br.Atom(name="C1", symbol="C", xyz=rs_br.Position(0.0, 0.0, 0.0)))
        ag.set_atom("H1", rs_br.Atom(name="H1", symbol="H", xyz=rs_br.Position(1.0, 0.0, 0.0)))

        selected = ag.select(CarbonSelector())
        self.assertEqual(selected.get_number_of_atoms(), 1)
        self.assertEqual(selected.get_atom("C1").symbol, "C")


if __name__ == "__main__":
    unittest.main()
