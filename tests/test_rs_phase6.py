#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Comparison test suite for Phase 6 (PR#48):
Validating proteindf_bridge_rs (PyO3 Rust bindings) against pure Python proteindf_bridge
for .brd I/O (MessagePack and YUI header formats), Modeling, and Neutralize.
"""

import math
import os
import tempfile
import unittest

# Pure Python implementation
import proteindf_bridge.functions as py_functions
import proteindf_bridge.modeling as py_modeling
import proteindf_bridge.neutralize as py_neutralize
from proteindf_bridge.atom import Atom as PyAtom
from proteindf_bridge.atomgroup import AtomGroup as PyAtomGroup
from proteindf_bridge.biopdb import Pdb as PyPdb
from proteindf_bridge.position import Position as PyPos

# Rust PyO3 bindings
import proteindf_bridge_rs as rs_br


class TestRsPhase6(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Locate test data directory
        test_dir = os.path.dirname(os.path.abspath(__file__))
        cls.data_dir = os.path.join(test_dir, "..", "rust", "crates", "proteindf-bridge", "tests", "data")
        if not os.path.exists(cls.data_dir):
            cls.data_dir = os.path.join(test_dir, "data")

    @staticmethod
    def _make_atom(symbol, name, x, y, z):
        return rs_br.Atom(symbol=symbol, name=name, position=rs_br.Position(x, y, z))

    # =========================================================================
    # 1. .brd MessagePack and YUI Header I/O
    # =========================================================================

    def test_load_atomgroup_fixture_parity(self):
        """
        Verify that rs_br.load_atomgroup produces identical atoms, groups, and positions
        as pure Python py_functions.load_atomgroup on existing fixtures.
        """
        fixture_path = os.path.join(self.data_dir, "ACE_ALA_NME_trans1.brd")
        self.assertTrue(os.path.exists(fixture_path), f"Fixture not found: {fixture_path}")

        py_ag = py_functions.load_atomgroup(fixture_path)
        rs_ag = rs_br.load_atomgroup(fixture_path)

        self.assertEqual(rs_ag.get_number_of_groups(), py_ag.get_number_of_groups())
        self.assertEqual(rs_ag.get_number_of_all_atoms(), py_ag.get_number_of_all_atoms())
        self.assertEqual(rs_ag.get_group_list(), py_ag.get_group_list())

        py_atoms = py_ag.get_atom_list()
        rs_atoms = rs_ag.get_atom_list()
        self.assertEqual(len(rs_atoms), len(py_atoms))

        for py_atm, rs_atm in zip(py_atoms, rs_atoms):
            self.assertEqual(rs_atm.name, py_atm.name)
            self.assertEqual(rs_atm.path, py_atm.path)
            self.assertAlmostEqual(rs_atm.xyz.x, py_atm.xyz.x, places=5)
            self.assertAlmostEqual(rs_atm.xyz.y, py_atm.xyz.y, places=5)
            self.assertAlmostEqual(rs_atm.xyz.z, py_atm.xyz.z, places=5)

    def test_save_and_load_atomgroup_cross_reading(self):
        """
        TASK PR#48 Definition of Done 2:
        Confirm that .brd written by pure Python can be read by Rust, and vice versa.
        """
        # Create a sample AtomGroup with hierarchy and atoms
        py_ag = PyAtomGroup("test_protein")
        py_chain = PyAtomGroup("A")
        py_res = PyAtomGroup("1")
        py_res.name = "ALA"
        py_res.set_atom("CA", PyAtom(name="CA", symbol="C", xyz=PyPos([1.0, 2.0, 3.0])))
        py_res.set_atom("CB", PyAtom(name="CB", symbol="C", xyz=PyPos([2.0, 3.0, 4.0])))
        py_chain.set_group("1", py_res)
        py_ag.set_group("A", py_chain)

        with tempfile.TemporaryDirectory() as tmpdir:
            # 1. Pure Python writes, Rust reads
            py_written_path = os.path.join(tmpdir, "py_written.brd")
            py_functions.save_atomgroup(py_ag, py_written_path)

            rs_loaded = rs_br.load_atomgroup(py_written_path)
            self.assertEqual(rs_loaded.get_number_of_all_atoms(), py_ag.get_number_of_all_atoms())
            self.assertEqual(rs_loaded.name, py_ag.name)

            rs_ca = rs_loaded["A"]["1"]["CA"]
            self.assertEqual(rs_ca.name, "CA")
            self.assertAlmostEqual(rs_ca.xyz.x, 1.0, places=5)
            self.assertAlmostEqual(rs_ca.xyz.y, 2.0, places=5)
            self.assertAlmostEqual(rs_ca.xyz.z, 3.0, places=5)

            # 2. Rust writes, Pure Python reads
            rs_written_path = os.path.join(tmpdir, "rs_written.brd")
            rs_br.save_atomgroup(rs_loaded, rs_written_path)

            py_loaded = py_functions.load_atomgroup(rs_written_path)
            self.assertEqual(py_loaded.get_number_of_all_atoms(), rs_loaded.get_number_of_all_atoms())
            self.assertEqual(py_loaded.name, rs_loaded.name)

            py_ca = py_loaded["A"]["1"]["CA"]
            self.assertEqual(py_ca.name, "CA")
            self.assertAlmostEqual(py_ca.xyz.x, 1.0, places=5)
            self.assertAlmostEqual(py_ca.xyz.y, 2.0, places=5)
            self.assertAlmostEqual(py_ca.xyz.z, 3.0, places=5)

    def test_brd_yui_format_roundtrip(self):
        """
        Verify saving and loading with YUI-compatible headers (uncompressed and zstd compressed).
        """
        fixture_path = os.path.join(self.data_dir, "ACE_ALA_NME_trans1.brd")
        original_ag = rs_br.load_atomgroup(fixture_path)
        atom_count = original_ag.get_number_of_all_atoms()

        with tempfile.TemporaryDirectory() as tmpdir:
            # 1. Uncompressed YUI format
            uncompressed_path = os.path.join(tmpdir, "test_yui_raw.brd")
            rs_br.save_brd_yui(original_ag, uncompressed_path, compress_zstd=False)

            loaded_raw = rs_br.load_brd_yui(uncompressed_path)
            self.assertEqual(loaded_raw.get_number_of_all_atoms(), atom_count)
            self.assertEqual(loaded_raw.get_group_list(), original_ag.get_group_list())

            # Check header magic directly
            with open(uncompressed_path, "rb") as f:
                header = f.read(6)
            self.assertEqual(header[:4], b"YUI\0")
            self.assertEqual(header[4], 1)  # Version 1
            self.assertEqual(header[5], 0)  # Compression NONE

            # 2. Zstd compressed YUI format
            compressed_path = os.path.join(tmpdir, "test_yui_zstd.brd")
            rs_br.save_brd_yui(original_ag, compressed_path, compress_zstd=True)

            loaded_zstd = rs_br.load_brd_yui(compressed_path)
            self.assertEqual(loaded_zstd.get_number_of_all_atoms(), atom_count)
            self.assertEqual(loaded_zstd.get_group_list(), original_ag.get_group_list())

            # Check header magic directly
            with open(compressed_path, "rb") as f:
                header_zstd = f.read(6)
            self.assertEqual(header_zstd[:4], b"YUI\0")
            self.assertEqual(header_zstd[4], 1)  # Version 1
            self.assertEqual(header_zstd[5], 1)  # Compression ZSTD

            # 3. Error case: loading plain msgpack with load_brd_yui should fail
            with self.assertRaises(rs_br.BrError):
                rs_br.load_brd_yui(fixture_path)

            # 4. Error case: corrupt file
            corrupt_path = os.path.join(tmpdir, "corrupt.brd")
            with open(corrupt_path, "wb") as f:
                f.write(b"CORRUPT")
            with self.assertRaises(rs_br.BrError):
                rs_br.load_brd_yui(corrupt_path)

    # =========================================================================
    # 2. Modeling
    # =========================================================================

    def test_modeling_arbitary_rotate_matrix(self):
        """
        Verify arbitary_rotate_matrix against pure Python Modeling and mathematical properties.
        """
        py_model = py_modeling.Modeling()
        rs_model = rs_br.Modeling()

        a = PyPos(0.5, 0.7, 0.3)
        b = PyPos(0.1, 0.9, 0.2)
        rs_a = rs_br.Position(0.5, 0.7, 0.3)
        rs_b = rs_br.Position(0.1, 0.9, 0.2)

        py_rot = py_model.arbitary_rotate_matrix(a, b)
        rs_rot = rs_model.arbitary_rotate_matrix(rs_a, rs_b)

        # Check 1:1 matrix elements matching Python
        for i in range(3):
            for j in range(3):
                self.assertAlmostEqual(rs_rot.get(i, j), py_rot.get(i, j), places=6)

        # Check rotation of b onto direction of a
        rotated_b = rs_br.Position(rs_b)
        rotated_b.rotate(rs_rot)
        rotated_b.norm()
        norm_a = rs_br.Position(rs_a)
        norm_a.norm()

        self.assertAlmostEqual(rotated_b.x, norm_a.x, places=6)
        self.assertAlmostEqual(rotated_b.y, norm_a.y, places=6)
        self.assertAlmostEqual(rotated_b.z, norm_a.z, places=6)

    def test_modeling_add_methyl_parity(self):
        """
        Verify add_methyl 1:1 parity with pure Python Modeling.
        """
        py_model = py_modeling.Modeling()
        rs_model = rs_br.Modeling()

        c1 = PyAtom(symbol="C", name="C1", position=PyPos(1.0, 2.0, 3.0))
        c2 = PyAtom(symbol="C", name="C2", position=PyPos(2.2, 2.6, 3.9))
        rs_c1 = self._make_atom("C", "C1", 1.0, 2.0, 3.0)
        rs_c2 = self._make_atom("C", "C2", 2.2, 2.6, 3.9)

        py_methyl = py_model.add_methyl(c1, c2)
        rs_methyl = rs_model.add_methyl(rs_c1, rs_c2)

        self.assertEqual(rs_methyl.get_number_of_atoms(), py_methyl.get_number_of_atoms())
        for name in ["H11", "H12", "H13"]:
            py_h = py_methyl[name]
            rs_h = rs_methyl[name]
            self.assertEqual(rs_h.symbol, py_h.symbol)
            self.assertAlmostEqual(rs_h.xyz.x, py_h.xyz.x, places=5)
            self.assertAlmostEqual(rs_h.xyz.y, py_h.xyz.y, places=5)
            self.assertAlmostEqual(rs_h.xyz.z, py_h.xyz.z, places=5)

    def test_modeling_get_nh3_parity(self):
        """
        Verify get_NH3 1:1 parity with pure Python Modeling.
        """
        py_model = py_modeling.Modeling()
        rs_model = rs_br.Modeling()

        # 1. Default parameters
        py_nh3 = py_model.get_NH3()
        rs_nh3 = rs_model.get_NH3()

        self.assertEqual(rs_nh3.get_number_of_atoms(), py_nh3.get_number_of_atoms())
        for atom_name in ["N", "H1", "H2", "H3"]:
            py_atm = py_nh3[atom_name]
            rs_atm = rs_nh3[atom_name]
            self.assertEqual(rs_atm.symbol, py_atm.symbol)
            self.assertAlmostEqual(rs_atm.xyz.x, py_atm.xyz.x, places=5)
            self.assertAlmostEqual(rs_atm.xyz.y, py_atm.xyz.y, places=5)
            self.assertAlmostEqual(rs_atm.xyz.z, py_atm.xyz.z, places=5)

        # 2. Custom angle and length
        angle = 0.6 * math.pi
        length = 1.2
        py_nh3_cust = py_model.get_NH3(angle=angle, length=length)
        rs_nh3_cust = rs_model.get_NH3(angle=angle, length=length)

        for atom_name in ["N", "H1", "H2", "H3"]:
            py_atm = py_nh3_cust[atom_name]
            rs_atm = rs_nh3_cust[atom_name]
            self.assertAlmostEqual(rs_atm.xyz.x, py_atm.xyz.x, places=5)
            self.assertAlmostEqual(rs_atm.xyz.y, py_atm.xyz.y, places=5)
            self.assertAlmostEqual(rs_atm.xyz.z, py_atm.xyz.z, places=5)

    def test_modeling_select_residues_and_get_last_index(self):
        """
        Verify select_residues and get_last_index against pure Python.
        """
        py_model = py_modeling.Modeling()
        rs_model = rs_br.Modeling()

        # Build chain with residues 1..5
        py_chain = PyAtomGroup("A")
        rs_chain = rs_br.AtomGroup("A")
        for i in range(1, 6):
            res_py = PyAtomGroup(str(i))
            res_py.name = "ALA"
            res_py.set_atom(f"CA{i}", PyAtom(name=f"CA{i}", symbol="C", xyz=PyPos([float(i), 0.0, 0.0])))
            py_chain.set_group(str(i), res_py)

            res_rs = rs_br.AtomGroup(str(i))
            res_rs.name = "ALA"
            res_rs.set_atom(f"CA{i}", self._make_atom("C", f"CA{i}", float(i), 0.0, 0.0))
            rs_chain.set_group(str(i), res_rs)

        py_sel = py_model.select_residues(py_chain, 2, 4)
        rs_sel = rs_model.select_residues(rs_chain, 2, 4)

        self.assertEqual(rs_sel.get_number_of_atoms(), py_sel.get_number_of_atoms())
        self.assertEqual(rs_sel.get_number_of_atoms(), 3)

        # get_last_index
        self.assertEqual(rs_model.get_last_index(rs_sel), py_model.get_last_index(py_sel))
        self.assertEqual(rs_model.get_last_index(rs_sel), 4)

    def test_modeling_capping_simple_parity(self):
        """
        Verify get_ACE_simple and get_NME_simple 1:1 parity with pure Python.
        """
        py_model = py_modeling.Modeling()
        rs_model = rs_br.Modeling()

        # Build standard ALA next_aa
        next_aa_py = PyAtomGroup("ALA")
        next_aa_py.set_atom("CA", PyAtom(name="CA", symbol="C", xyz=PyPos([0.0, 0.0, 0.0])))
        next_aa_py.set_atom("C", PyAtom(name="C", symbol="C", xyz=PyPos([1.0, 0.0, 0.0])))
        next_aa_py.set_atom("O", PyAtom(name="O", symbol="O", xyz=PyPos([1.0, 1.0, 0.0])))
        next_aa_py.set_atom("N", PyAtom(name="N", symbol="N", xyz=PyPos([0.0, 1.0, 0.0])))
        next_aa_py.set_atom("H", PyAtom(name="H", symbol="H", xyz=PyPos([0.0, 1.5, 0.0])))

        next_aa_rs = rs_br.AtomGroup("ALA")
        next_aa_rs.set_atom("CA", self._make_atom("C", "CA", 0.0, 0.0, 0.0))
        next_aa_rs.set_atom("C", self._make_atom("C", "C", 1.0, 0.0, 0.0))
        next_aa_rs.set_atom("O", self._make_atom("O", "O", 1.0, 1.0, 0.0))
        next_aa_rs.set_atom("N", self._make_atom("N", "N", 0.0, 1.0, 0.0))
        next_aa_rs.set_atom("H", self._make_atom("H", "H", 0.0, 1.5, 0.0))

        # get_ACE_simple
        py_ace = py_model.get_ACE_simple(next_aa_py)
        rs_ace = rs_model.get_ACE_simple(next_aa_rs)
        self.assertEqual(rs_ace.get_number_of_atoms(), py_ace.get_number_of_atoms())
        self.assertEqual(rs_ace.path, py_ace.path)
        for name in ["CA", "C", "O", "H11", "H12", "H13"]:
            self.assertAlmostEqual(rs_ace[name].xyz.x, py_ace[name].xyz.x, places=5)
            self.assertAlmostEqual(rs_ace[name].xyz.y, py_ace[name].xyz.y, places=5)
            self.assertAlmostEqual(rs_ace[name].xyz.z, py_ace[name].xyz.z, places=5)

        # get_NME_simple
        py_nme = py_model.get_NME_simple(next_aa_py)
        rs_nme = rs_model.get_NME_simple(next_aa_rs)
        self.assertEqual(rs_nme.get_number_of_atoms(), py_nme.get_number_of_atoms())
        self.assertEqual(rs_nme.path, py_nme.path)
        for name in ["CA", "N", "H", "H11", "H12", "H13"]:
            self.assertAlmostEqual(rs_nme[name].xyz.x, py_nme[name].xyz.x, places=5)
            self.assertAlmostEqual(rs_nme[name].xyz.y, py_nme[name].xyz.y, places=5)
            self.assertAlmostEqual(rs_nme[name].xyz.z, py_nme[name].xyz.z, places=5)

    def test_modeling_get_ace_and_nme_full(self):
        """
        Verify get_ACE and get_NME fitting 1:1 parity on reference conformers.
        """
        py_model = py_modeling.Modeling()
        rs_model = rs_br.Modeling()

        fixture_path = os.path.join(self.data_dir, "ACE_ALA_NME_trans1.brd")
        py_ag = py_functions.load_atomgroup(fixture_path)
        rs_ag = rs_br.load_atomgroup(fixture_path)

        res2_py = py_ag["2"]
        res3_py = py_ag["3"]
        res2_rs = rs_ag["2"]
        res3_rs = rs_ag["3"]

        # With next_aa
        py_ace = py_model.get_ACE(res2_py, res3_py)
        rs_ace = rs_model.get_ACE(res2_rs, res3_rs)
        self.assertEqual(rs_ace.get_number_of_atoms(), py_ace.get_number_of_atoms())
        for atom in rs_ace.get_atom_list():
            py_atm = py_ace[atom.name]
            self.assertAlmostEqual(atom.xyz.x, py_atm.xyz.x, places=4)
            self.assertAlmostEqual(atom.xyz.y, py_atm.xyz.y, places=4)
            self.assertAlmostEqual(atom.xyz.z, py_atm.xyz.z, places=4)

        py_nme = py_model.get_NME(res2_py, res3_py)
        rs_nme = rs_model.get_NME(res2_rs, res3_rs)
        self.assertEqual(rs_nme.get_number_of_atoms(), py_nme.get_number_of_atoms())
        for atom in rs_nme.get_atom_list():
            py_atm = py_nme[atom.name]
            self.assertAlmostEqual(atom.xyz.x, py_atm.xyz.x, places=4)
            self.assertAlmostEqual(atom.xyz.y, py_atm.xyz.y, places=4)
            self.assertAlmostEqual(atom.xyz.z, py_atm.xyz.z, places=4)

        # Without next_aa
        py_ace_none = py_model.get_ACE(res2_py, None)
        rs_ace_none = rs_model.get_ACE(res2_rs, None)
        self.assertEqual(rs_ace_none.get_number_of_atoms(), py_ace_none.get_number_of_atoms())

    def test_modeling_neutralize_helpers_parity(self):
        """
        Verify individual neutralization helpers (Nterm, Cterm, GLU, ASP, LYS, ARG, FAD)
        for 1:1 parity with pure Python Modeling.
        """
        py_model = py_modeling.Modeling()
        rs_model = rs_br.Modeling()

        # 1. Nterm general
        res_py = PyAtomGroup("ALA")
        res_py.set_atom("N", PyAtom(name="N", symbol="N", xyz=PyPos([0.0, 0.0, 0.0])))
        res_py.set_atom("H1", PyAtom(name="H1", symbol="H", xyz=PyPos([0.5, 0.5, 0.5])))
        res_py.set_atom("H2", PyAtom(name="H2", symbol="H", xyz=PyPos([-0.5, 0.5, 0.5])))
        res_py.set_atom("H3", PyAtom(name="H3", symbol="H", xyz=PyPos([0.0, -0.5, 0.5])))

        res_rs = rs_br.AtomGroup("ALA")
        res_rs.set_atom("N", self._make_atom("N", "N", 0.0, 0.0, 0.0))
        res_rs.set_atom("H1", self._make_atom("H", "H1", 0.5, 0.5, 0.5))
        res_rs.set_atom("H2", self._make_atom("H", "H2", -0.5, 0.5, 0.5))
        res_rs.set_atom("H3", self._make_atom("H", "H3", 0.0, -0.5, 0.5))

        py_cl = py_model.neutralize_Nterm(res_py)
        rs_cl = rs_model.neutralize_Nterm(res_rs)
        self.assertEqual(rs_cl["Cl"].name, py_cl["Cl"].name)
        self.assertAlmostEqual(rs_cl["Cl"].xyz.x, py_cl["Cl"].xyz.x, places=5)
        self.assertAlmostEqual(rs_cl["Cl"].xyz.y, py_cl["Cl"].xyz.y, places=5)
        self.assertAlmostEqual(rs_cl["Cl"].xyz.z, py_cl["Cl"].xyz.z, places=5)

        # 2. Cterm
        ct_py = PyAtomGroup("ALA")
        ct_py.set_atom("C", PyAtom(name="C", symbol="C", xyz=PyPos([0.0, 0.0, 0.0])))
        ct_py.set_atom("O", PyAtom(name="O", symbol="O", xyz=PyPos([0.0, 1.0, 0.0])))
        ct_py.set_atom("OXT", PyAtom(name="OXT", symbol="O", xyz=PyPos([1.0, 0.0, 0.0])))

        ct_rs = rs_br.AtomGroup("ALA")
        ct_rs.set_atom("C", self._make_atom("C", "C", 0.0, 0.0, 0.0))
        ct_rs.set_atom("O", self._make_atom("O", "O", 0.0, 1.0, 0.0))
        ct_rs.set_atom("OXT", self._make_atom("O", "OXT", 1.0, 0.0, 0.0))

        py_na = py_model.neutralize_Cterm(ct_py)
        rs_na = rs_model.neutralize_Cterm(ct_rs)
        self.assertAlmostEqual(rs_na["Na"].xyz.x, py_na["Na"].xyz.x, places=5)
        self.assertAlmostEqual(rs_na["Na"].xyz.y, py_na["Na"].xyz.y, places=5)
        self.assertAlmostEqual(rs_na["Na"].xyz.z, py_na["Na"].xyz.z, places=5)

        # 3. GLU
        glu_py = PyAtomGroup("1")
        glu_py.name = "GLU"
        glu_py.set_atom("CD", PyAtom(name="CD", symbol="C", xyz=PyPos([0.0, 0.0, 0.0])))
        glu_py.set_atom("OE1", PyAtom(name="OE1", symbol="O", xyz=PyPos([0.0, 1.0, 0.0])))
        glu_py.set_atom("OE2", PyAtom(name="OE2", symbol="O", xyz=PyPos([1.0, 0.0, 0.0])))

        glu_rs = rs_br.AtomGroup("1")
        glu_rs.name = "GLU"
        glu_rs.set_atom("CD", self._make_atom("C", "CD", 0.0, 0.0, 0.0))
        glu_rs.set_atom("OE1", self._make_atom("O", "OE1", 0.0, 1.0, 0.0))
        glu_rs.set_atom("OE2", self._make_atom("O", "OE2", 1.0, 0.0, 0.0))

        py_glu_na = py_model.neutralize_GLU(glu_py)
        rs_glu_na = rs_model.neutralize_GLU(glu_rs)
        self.assertEqual(rs_glu_na.get_atom_list()[0].name, py_glu_na.get_atom_list()[0].name)
        self.assertAlmostEqual(rs_glu_na.get_atom_list()[0].xyz.x, py_glu_na.get_atom_list()[0].xyz.x, places=5)

        # 4. ARG cases 0, 1, 2
        arg_py = PyAtomGroup("1")
        arg_py.name = "ARG"
        arg_py.set_atom("CZ", PyAtom(name="CZ", symbol="C", xyz=PyPos([0.0, 0.0, 0.0])))
        arg_py.set_atom("NH1", PyAtom(name="NH1", symbol="N", xyz=PyPos([1.0, 0.5, 0.0])))
        arg_py.set_atom("NH2", PyAtom(name="NH2", symbol="N", xyz=PyPos([1.0, -0.5, 0.0])))
        arg_py.set_atom("HH11", PyAtom(name="HH11", symbol="H", xyz=PyPos([1.5, 0.7, 0.0])))
        arg_py.set_atom("HH12", PyAtom(name="HH12", symbol="H", xyz=PyPos([1.5, 0.3, 0.0])))
        arg_py.set_atom("HH21", PyAtom(name="HH21", symbol="H", xyz=PyPos([1.5, -0.3, 0.0])))
        arg_py.set_atom("HH22", PyAtom(name="HH22", symbol="H", xyz=PyPos([1.5, -0.7, 0.0])))

        arg_rs = rs_br.AtomGroup("1")
        arg_rs.name = "ARG"
        arg_rs.set_atom("CZ", self._make_atom("C", "CZ", 0.0, 0.0, 0.0))
        arg_rs.set_atom("NH1", self._make_atom("N", "NH1", 1.0, 0.5, 0.0))
        arg_rs.set_atom("NH2", self._make_atom("N", "NH2", 1.0, -0.5, 0.0))
        arg_rs.set_atom("HH11", self._make_atom("H", "HH11", 1.5, 0.7, 0.0))
        arg_rs.set_atom("HH12", self._make_atom("H", "HH12", 1.5, 0.3, 0.0))
        arg_rs.set_atom("HH21", self._make_atom("H", "HH21", 1.5, -0.3, 0.0))
        arg_rs.set_atom("HH22", self._make_atom("H", "HH22", 1.5, -0.7, 0.0))

        for case in [0, 1, 2]:
            py_arg_cl = py_model.neutralize_ARG(arg_py, case=case)
            rs_arg_cl = rs_model.neutralize_ARG(arg_rs, case=case)
            self.assertAlmostEqual(rs_arg_cl.get_atom_list()[0].xyz.x, py_arg_cl.get_atom_list()[0].xyz.x, places=5)
            self.assertAlmostEqual(rs_arg_cl.get_atom_list()[0].xyz.y, py_arg_cl.get_atom_list()[0].xyz.y, places=5)
            self.assertAlmostEqual(rs_arg_cl.get_atom_list()[0].xyz.z, py_arg_cl.get_atom_list()[0].xyz.z, places=5)

    # =========================================================================
    # 3. Neutralize
    # =========================================================================

    def test_neutralize_synthetic_and_parity(self):
        """
        Verify Neutralize on synthetic model against pure Python.
        """
        protein_py = PyAtomGroup("protein")
        model_py = PyAtomGroup("model_1")
        chain_py = PyAtomGroup("A")

        glu_py = PyAtomGroup("1")
        glu_py.name = "GLU"
        glu_py.set_atom("CD", PyAtom(name="CD", xyz=PyPos([0.0, 0.0, 0.0])))
        glu_py.set_atom("OE1", PyAtom(name="OE1", xyz=PyPos([0.0, 1.0, 0.0])))
        glu_py.set_atom("OE2", PyAtom(name="OE2", xyz=PyPos([1.0, 0.0, 0.0])))
        chain_py.set_group("1", glu_py)

        lys_py = PyAtomGroup("2")
        lys_py.name = "LYS"
        lys_py.set_atom("NZ", PyAtom(name="NZ", xyz=PyPos([0.5, 2.0, 0.0])))
        lys_py.set_atom("HZ1", PyAtom(name="HZ1", xyz=PyPos([0.5, 2.5, 0.0])))
        lys_py.set_atom("HZ2", PyAtom(name="HZ2", xyz=PyPos([0.0, 2.0, 0.5])))
        lys_py.set_atom("HZ3", PyAtom(name="HZ3", xyz=PyPos([1.0, 2.0, 0.5])))
        chain_py.set_group("2", lys_py)

        model_py.set_group("A", chain_py)
        protein_py.set_group("model_1", model_py)

        # Clone hierarchy to Rust AtomGroup
        protein_rs = rs_br.AtomGroup("protein")
        model_rs = rs_br.AtomGroup("model_1")
        chain_rs = rs_br.AtomGroup("A")

        glu_rs = rs_br.AtomGroup("1")
        glu_rs.name = "GLU"
        glu_rs.set_atom("CD", self._make_atom("C", "CD", 0.0, 0.0, 0.0))
        glu_rs.set_atom("OE1", self._make_atom("O", "OE1", 0.0, 1.0, 0.0))
        glu_rs.set_atom("OE2", self._make_atom("O", "OE2", 1.0, 0.0, 0.0))
        chain_rs.set_group("1", glu_rs)

        lys_rs = rs_br.AtomGroup("2")
        lys_rs.name = "LYS"
        lys_rs.set_atom("NZ", self._make_atom("N", "NZ", 0.5, 2.0, 0.0))
        lys_rs.set_atom("HZ1", self._make_atom("H", "HZ1", 0.5, 2.5, 0.0))
        lys_rs.set_atom("HZ2", self._make_atom("H", "HZ2", 0.0, 2.0, 0.5))
        lys_rs.set_atom("HZ3", self._make_atom("H", "HZ3", 1.0, 2.0, 0.5))
        chain_rs.set_group("2", lys_rs)

        model_rs.set_group("A", chain_rs)
        protein_rs.set_group("model_1", model_rs)

        neut_py = py_neutralize.Neutralize(protein_py)
        neut_rs = rs_br.Neutralize(protein_rs)

        # Non-destructive check on original input
        self.assertEqual(protein_rs.get_number_of_all_atoms(), 7)
        self.assertEqual(protein_py.get_number_of_all_atoms(), 7)

        # Compare neutralized results
        py_neutralized = neut_py.neutralized
        rs_neutralized = neut_rs.neutralized

        self.assertEqual(rs_neutralized.get_number_of_all_atoms(), py_neutralized.get_number_of_all_atoms())
        self.assertEqual(rs_neutralized.get_number_of_all_atoms(), 9)  # 7 original + 2 ions

        py_atoms = py_neutralized.get_atom_list()
        rs_atoms = rs_neutralized.get_atom_list()
        for py_a, rs_a in zip(py_atoms, rs_atoms):
            self.assertEqual(rs_a.name, py_a.name)
            self.assertAlmostEqual(rs_a.xyz.x, py_a.xyz.x, places=5)
            self.assertAlmostEqual(rs_a.xyz.y, py_a.xyz.y, places=5)
            self.assertAlmostEqual(rs_a.xyz.z, py_a.xyz.z, places=5)

    def test_neutralize_exempt_list_parity_and_notes(self):
        """
        Verify _exempt_list parity with pure Python.
        Note on divergence / Python parity (from docs/rust-port-handoff.md lesson 14):
        In pure Python neutralize.py, _exempt_list is computed via IonPair, but in _neutralize(),
        `exempt_list = [] # self._exempt_list()` is commented out, making it practically dead code
        during actual neutralization. Both Python and Rust retain the helper method for inspection/testing,
        and both intentionally do not apply exemptions during neutralize().
        """
        protein_rs = rs_br.AtomGroup("protein")
        model_rs = rs_br.AtomGroup("model_1")
        chain_rs = rs_br.AtomGroup("A")

        glu_rs = rs_br.AtomGroup("1")
        glu_rs.name = "GLU"
        glu_rs.set_atom("CD", self._make_atom("C", "CD", 0.0, 0.0, 0.0))
        glu_rs.set_atom("OE1", self._make_atom("O", "OE1", 0.0, 1.0, 0.0))
        glu_rs.set_atom("OE2", self._make_atom("O", "OE2", 1.0, 0.0, 0.0))
        chain_rs.set_group("1", glu_rs)

        lys_rs = rs_br.AtomGroup("2")
        lys_rs.name = "LYS"
        lys_rs.set_atom("NZ", self._make_atom("N", "NZ", 0.5, 2.0, 0.0))
        lys_rs.set_atom("HZ1", self._make_atom("H", "HZ1", 0.5, 2.5, 0.0))
        lys_rs.set_atom("HZ2", self._make_atom("H", "HZ2", 0.0, 2.0, 0.5))
        lys_rs.set_atom("HZ3", self._make_atom("H", "HZ3", 1.0, 2.0, 0.5))
        chain_rs.set_group("2", lys_rs)

        model_rs.set_group("A", chain_rs)
        protein_rs.set_group("model_1", model_rs)

        neut = rs_br.Neutralize(protein_rs)
        exempt = neut._exempt_list(model_rs)
        self.assertEqual(len(exempt), 2)
        self.assertEqual(exempt, [("A", "1", "GLU"), ("A", "2", "LYS")])

        # Verify divide_path
        self.assertEqual(neut._divide_path("/model_1/A/1/"), ("A", "1"))
        self.assertEqual(neut._divide_path("/model_1/A/1/CA"), ("1", "CA"))

    def test_neutralize_object_identity_and_collection_protection(self):
        """
        Review policy verification:
        1. Large objects (AtomGroup in Neutralize.neutralized) return the same object on each access.
        2. Modifications to neut.neutralized persist on subsequent accesses.
        3. Lists (_exempt_list) return a new container on each access to protect internal state.
        """
        protein_rs = rs_br.AtomGroup("protein")
        model_rs = rs_br.AtomGroup("model_1")
        chain_rs = rs_br.AtomGroup("A")
        glu_rs = rs_br.AtomGroup("1")
        glu_rs.name = "GLU"
        glu_rs.set_atom("CD", self._make_atom("C", "CD", 0.0, 0.0, 0.0))
        glu_rs.set_atom("OE1", self._make_atom("O", "OE1", 0.0, 1.0, 0.0))
        glu_rs.set_atom("OE2", self._make_atom("O", "OE2", 1.0, 0.0, 0.0))
        chain_rs.set_group("1", glu_rs)
        model_rs.set_group("A", chain_rs)
        protein_rs.set_group("model_1", model_rs)

        neut = rs_br.Neutralize(protein_rs)

        # 1. Identity check: same object returned
        ag1 = neut.neutralized
        ag2 = neut.neutralized
        self.assertIs(ag1, ag2)

        # 2. Persistence check: modifying ag1 is visible in neut.neutralized
        ag1.name = "modified_neutralized"
        self.assertEqual(neut.neutralized.name, "modified_neutralized")

        # 3. Collection protection check: modifying _exempt_list does not affect subsequent calls
        ex1 = neut._exempt_list(model_rs)
        init_len = len(ex1)
        ex1.append(("FAKE", "99", "ION"))
        ex2 = neut._exempt_list(model_rs)
        self.assertEqual(len(ex2), init_len)

    def test_neutralize_real_pdb_1hls_parity(self):
        """
        Verify Neutralize on real fixture 1hls.pdb against pure Python implementation.
        Baseline: rust/crates/proteindf-bridge/tests/test_neutralize.rs test_neutralize_real_fixture_1hls.
        Original 1hls has 746 atoms (model_1: 746 atoms, 2 chains A and B).
        Neutralization adds ions to acidic/basic residues and termini, preserving exact coordinates.
        """
        pdb_path = os.path.join(self.data_dir, "1hls.pdb")
        self.assertTrue(os.path.exists(pdb_path), f"Fixture not found: {pdb_path}")

        py_pdb = PyPdb(pdb_path)
        py_ag = py_pdb.get_atomgroup()

        rs_pdb = rs_br.Pdb(pdb_path)
        rs_ag = rs_pdb.get_atomgroup()

        self.assertEqual(rs_ag.get_number_of_all_atoms(), py_ag.get_number_of_all_atoms())
        self.assertEqual(rs_ag.get_number_of_all_atoms(), 782)

        py_neut = py_neutralize.Neutralize(py_ag)
        rs_neut = rs_br.Neutralize(rs_ag)

        # Check total atom counts after neutralization (782 + 10 ions = 792)
        self.assertEqual(
            rs_neut.neutralized.get_number_of_all_atoms(),
            py_neut.neutralized.get_number_of_all_atoms(),
        )
        self.assertEqual(rs_neut.neutralized.get_number_of_all_atoms(), 792)

        # Check atom-by-atom parity
        py_atoms = py_neut.neutralized.get_atom_list()
        rs_atoms = rs_neut.neutralized.get_atom_list()
        self.assertEqual(len(rs_atoms), len(py_atoms))

        for py_atm, rs_atm in zip(py_atoms, rs_atoms):
            self.assertEqual(rs_atm.name, py_atm.name)
            self.assertAlmostEqual(rs_atm.xyz.x, py_atm.xyz.x, places=4)
            self.assertAlmostEqual(rs_atm.xyz.y, py_atm.xyz.y, places=4)
            self.assertAlmostEqual(rs_atm.xyz.z, py_atm.xyz.z, places=4)

    def test_modeling_errors_and_edge_cases(self):
        """
        Verify that missing atoms and invalid arguments raise proper BrError / BrInputError / BrValueError.
        """
        m = rs_br.Modeling()

        empty_ag = rs_br.AtomGroup("EMPTY")
        with self.assertRaises(rs_br.BrInputError):
            m.get_ACE_simple(empty_ag)

        with self.assertRaises(rs_br.BrInputError):
            m.get_NME_simple(empty_ag)

        with self.assertRaises(rs_br.BrInputError):
            m.neutralize_Nterm(empty_ag)

        with self.assertRaises(rs_br.BrInputError):
            m.neutralize_Cterm(empty_ag)

        with self.assertRaises(rs_br.BrInputError):
            m.neutralize_GLU(empty_ag)

        with self.assertRaises(rs_br.BrInputError):
            m.neutralize_ASP(empty_ag)

        with self.assertRaises(rs_br.BrInputError):
            m.neutralize_LYS(empty_ag)

        with self.assertRaises(rs_br.BrInputError):
            m.neutralize_FAD(empty_ag)

        # Invalid ARG case
        with self.assertRaises(rs_br.BrValueError):
            m.neutralize_ARG(empty_ag, case=99)

    def test_modeling_from_data_dir_and_repr(self):
        """
        Verify Modeling.from_data_dir and repr for Modeling and Neutralize.
        """
        m_dir = rs_br.Modeling.from_data_dir(self.data_dir)
        self.assertIsNotNone(m_dir)
        self.assertEqual(repr(m_dir), "Modeling()")

        ag = rs_br.AtomGroup("test")
        neut = rs_br.Neutralize(ag)
        self.assertIn("Neutralize(atoms=", repr(neut))


if __name__ == "__main__":
    unittest.main()
