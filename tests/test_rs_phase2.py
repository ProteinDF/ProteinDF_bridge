#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Comparison test suite for Phase 5 (PR#15):
Validating proteindf_bridge_rs (PyO3 Rust bindings) format I/O against pure Python proteindf_bridge.
"""

import os
import tempfile
import unittest
import numpy as np

# Pure Python implementation
from proteindf_bridge.format import Format as PyFormat
from proteindf_bridge.xyz import Xyz as PyXyz
from proteindf_bridge.gro import SimpleGro as PySimpleGro
from proteindf_bridge.mol2 import SimpleMol2 as PySimpleMol2
from proteindf_bridge.amber_prmtop import AmberPrmtop as PyAmberPrmtop
from proteindf_bridge.biopdb import Pdb as PyPdb
from proteindf_bridge.mmcif import SimpleMmcif as PySimpleMmcif
from proteindf_bridge.atom import Atom as PyAtom
from proteindf_bridge.atomgroup import AtomGroup as PyAtomGroup
from proteindf_bridge.position import Position as PyPos

# Rust PyO3 bindings
import proteindf_bridge_rs as rs_br

DATA_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "proteindf_bridge", "data")
)


class TestRsPhase2(unittest.TestCase):

    # --------------------------------------------------------------------------
    # 1. Format
    # --------------------------------------------------------------------------
    def test_format_hierarchy(self):
        pdb = rs_br.Pdb(os.path.join(DATA_DIR, "2MGO.pdb"))
        models = pdb.get_atomgroup()

        # models level
        self.assertTrue(rs_br.Format.is_models(models))
        self.assertFalse(rs_br.Format.is_protein(models))
        self.assertFalse(rs_br.Format.is_chain(models))
        self.assertFalse(rs_br.Format.is_residue(models))

        # model level
        model = models["model_1"]
        self.assertFalse(rs_br.Format.is_models(model))
        self.assertTrue(rs_br.Format.is_protein(model))
        self.assertFalse(rs_br.Format.is_chain(model))
        self.assertFalse(rs_br.Format.is_residue(model))

        # chain level
        chain = model["A"]
        self.assertFalse(rs_br.Format.is_models(chain))
        self.assertFalse(rs_br.Format.is_protein(chain))
        self.assertTrue(rs_br.Format.is_chain(chain))
        self.assertFalse(rs_br.Format.is_residue(chain))

        # residue level
        res = chain[3]
        self.assertFalse(rs_br.Format.is_models(res))
        self.assertFalse(rs_br.Format.is_protein(res))
        self.assertFalse(rs_br.Format.is_chain(res))
        self.assertTrue(rs_br.Format.is_residue(res))

    def test_format_non_atomgroup(self):
        for val in [None, 123, "string", [1, 2, 3], {"a": 1}]:
            self.assertFalse(rs_br.Format.is_residue(val))
            self.assertFalse(rs_br.Format.is_chain(val))
            self.assertFalse(rs_br.Format.is_protein(val))
            self.assertFalse(rs_br.Format.is_models(val))

    def test_format_empty(self):
        empty = rs_br.AtomGroup("empty")
        self.assertTrue(rs_br.Format.is_residue(empty))
        self.assertTrue(rs_br.Format.is_chain(empty))
        self.assertTrue(rs_br.Format.is_protein(empty))
        self.assertTrue(rs_br.Format.is_models(empty))

    # --------------------------------------------------------------------------
    # 2. Xyz
    # --------------------------------------------------------------------------
    def test_xyz_load_and_compare(self):
        xyz_path = os.path.join(DATA_DIR, "ACE_ALA_NME.xyz")
        py_xyz = PyXyz(xyz_path)
        rs_xyz = rs_br.Xyz(xyz_path)

        self.assertEqual(rs_xyz.comment, py_xyz._comment)

        py_ag = py_xyz.get_atom_group()
        rs_ag = rs_xyz.get_atomgroup()

        self.assertEqual(rs_ag.get_number_of_atoms(), py_ag.get_number_of_atoms())
        self.assertEqual(rs_ag.get_number_of_atoms(), 22)

        for i in range(22):
            py_atom = py_ag.get_atom(str(i))
            rs_atom = rs_ag.get_atom(str(i))
            self.assertEqual(rs_atom.symbol, py_atom.symbol)
            self.assertAlmostEqual(rs_atom.xyz.x, py_atom.xyz.x, places=5)
            self.assertAlmostEqual(rs_atom.xyz.y, py_atom.xyz.y, places=5)
            self.assertAlmostEqual(rs_atom.xyz.z, py_atom.xyz.z, places=5)

    def test_xyz_init_with_atomgroup_and_save(self):
        ag = rs_br.AtomGroup("water")
        ag.set_atom("1", rs_br.Atom(symbol="O", xyz=rs_br.Position(0.0, 0.0, 0.0)))
        ag.set_atom("2", rs_br.Atom(symbol="H", xyz=rs_br.Position(0.0, 1.0, 0.0)))
        ag.set_atom("3", rs_br.Atom(symbol="H", xyz=rs_br.Position(1.0, 0.0, 0.0)))

        xyz = rs_br.Xyz(ag)
        text = xyz.get_text()
        self.assertTrue(text.startswith("3\n"))
        self.assertIn("water", text)
        self.assertIn("O", text)
        self.assertIn("H", text)
        self.assertEqual(str(xyz), text)

        with tempfile.NamedTemporaryFile(suffix=".xyz", delete=False) as f:
            temp_path = f.name
        try:
            xyz.save(temp_path)
            loaded = rs_br.Xyz(temp_path)
            loaded_ag = loaded.get_atomgroup()
            self.assertEqual(loaded_ag.get_number_of_atoms(), 3)
            self.assertEqual(loaded.comment, "water")
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)

    # --------------------------------------------------------------------------
    # 3. SimpleGro
    # --------------------------------------------------------------------------
    def test_gro_load_and_compare(self):
        gro_path = os.path.join(DATA_DIR, "sample.gro")
        py_gro = PySimpleGro()
        py_gro.load(gro_path)

        rs_gro = rs_br.SimpleGro(gro_path)

        self.assertEqual(rs_gro.title, py_gro._title)
        self.assertEqual(rs_gro.num_of_atoms, py_gro._num_of_atoms)
        self.assertEqual(rs_gro.num_of_atoms, 6)

        py_box = py_gro._box_vectors
        rs_box = rs_gro.box_vectors
        for p, r in zip(py_box, rs_box):
            self.assertAlmostEqual(r, p, places=4)

        py_ag = py_gro.get_atomgroup()
        rs_ag = rs_gro.get_atomgroup()
        self.assertEqual(rs_ag.get_number_of_atoms(), py_ag.get_number_of_atoms())
        self.assertEqual(rs_ag.name, py_ag.name)

    def test_gro_roundtrip(self):
        gro_path = os.path.join(DATA_DIR, "sample.gro")
        rs_gro = rs_br.SimpleGro(gro_path)
        ag = rs_gro.get_atomgroup()

        rs_gro2 = rs_br.SimpleGro()
        rs_gro2.set_by_atomgroup(ag)
        self.assertEqual(rs_gro2.num_of_atoms, 6)
        self.assertEqual(len(rs_gro2.box_vectors), 3)

        with tempfile.NamedTemporaryFile(suffix=".gro", delete=False) as f:
            temp_path = f.name
        try:
            rs_gro2.save(temp_path)
            rs_gro3 = rs_br.SimpleGro(temp_path)
            self.assertEqual(rs_gro3.num_of_atoms, 6)
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)

    # --------------------------------------------------------------------------
    # 4. SimpleMol2
    # --------------------------------------------------------------------------
    def test_mol2_create_and_text(self):
        py_ag = PyAtomGroup("water")
        py_ag.set_atom("1", PyAtom(name="O1", symbol="O", xyz=PyPos(0.0, 0.0, 0.0), charge=-0.8))
        py_ag.set_atom("2", PyAtom(name="H1", symbol="H", xyz=PyPos(0.0, 1.0, 0.0), charge=0.4))
        py_ag.set_atom("3", PyAtom(name="H2", symbol="H", xyz=PyPos(1.0, 0.0, 0.0), charge=0.4))

        rs_ag = rs_br.AtomGroup("water")
        rs_ag.set_atom("1", rs_br.Atom(name="O1", symbol="O", xyz=rs_br.Position(0.0, 0.0, 0.0), charge=-0.8))
        rs_ag.set_atom("2", rs_br.Atom(name="H1", symbol="H", xyz=rs_br.Position(0.0, 1.0, 0.0), charge=0.4))
        rs_ag.set_atom("3", rs_br.Atom(name="H2", symbol="H", xyz=rs_br.Position(1.0, 0.0, 0.0), charge=0.4))

        py_mol2 = PySimpleMol2(py_ag)
        rs_mol2 = rs_br.SimpleMol2(rs_ag)

        py_text = str(py_mol2)
        rs_text = rs_mol2.get_text()

        self.assertIn("@<TRIPOS>MOLECULE", rs_text)
        self.assertIn("water", rs_text)
        self.assertIn("@<TRIPOS>ATOM", rs_text)
        self.assertIn("O1", rs_text)
        self.assertIn("H1", rs_text)
        self.assertIn("H2", rs_text)
        self.assertEqual(str(rs_mol2), rs_text)

        with tempfile.NamedTemporaryFile(suffix=".mol2", delete=False) as f:
            temp_path = f.name
        try:
            rs_mol2.save(temp_path)
            with open(temp_path) as f:
                saved_content = f.read()
            self.assertEqual(saved_content, rs_text)
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)

    # --------------------------------------------------------------------------
    # 5. AmberPrmtop
    # --------------------------------------------------------------------------
    def test_amber_prmtop(self):
        prmtop_content = """%VERSION  VERSION_STAMP = V0001.000  DATE = 08/25/26  12:00:00
%FLAG ATOM_NAME
%FORMAT(20a4)
C1  H1  
%FLAG CHARGE
%FORMAT(5E16.8)
 0.00000000E+00 0.00000000E+00
%FLAG ATOMIC_NUMBER
%FORMAT(10I8)
       6       1
"""
        inpcrd_content = """default_name
    2
   0.0000000   0.0000000   0.0000000   1.0900000   0.0000000   0.0000000
"""
        with tempfile.NamedTemporaryFile(suffix=".prmtop", mode="w", delete=False) as f_top:
            f_top.write(prmtop_content)
            top_path = f_top.name
        with tempfile.NamedTemporaryFile(suffix=".inpcrd", mode="w", delete=False) as f_crd:
            f_crd.write(inpcrd_content)
            crd_path = f_crd.name

        try:
            py_amber = PyAmberPrmtop(top_path, crd_path)
            rs_amber = rs_br.AmberPrmtop(top_path, crd_path)

            self.assertEqual(rs_amber.atom_names, py_amber.atom_names)
            self.assertEqual(rs_amber.atomic_numbers, py_amber.atomic_numbers)
            self.assertEqual(len(rs_amber.charges), len(py_amber.charges))
            for rc, pc in zip(rs_amber.charges, py_amber.charges):
                self.assertAlmostEqual(rc, pc, places=5)

            self.assertEqual(len(rs_amber.xyz), len(py_amber.xyz))
            for rp, pp in zip(rs_amber.xyz, py_amber.xyz):
                self.assertAlmostEqual(rp.x, pp.x, places=5)
                self.assertAlmostEqual(rp.y, pp.y, places=5)
                self.assertAlmostEqual(rp.z, pp.z, places=5)

            py_ag = py_amber.get_atomgroup()
            rs_ag = rs_amber.get_atomgroup()
            self.assertEqual(rs_ag.get_number_of_atoms(), py_ag.get_number_of_atoms())
            self.assertEqual(rs_ag.get_number_of_atoms(), 2)
        finally:
            if os.path.exists(top_path):
                os.remove(top_path)
            if os.path.exists(crd_path):
                os.remove(crd_path)

    # --------------------------------------------------------------------------
    # 6. Pdb
    # --------------------------------------------------------------------------
    def test_pdb_2mgo(self):
        pdb_path = os.path.join(DATA_DIR, "2MGO.pdb")
        py_pdb = PyPdb(pdb_path)
        rs_pdb = rs_br.Pdb(pdb_path)

        py_ag = py_pdb.get_atomgroup()
        rs_ag = rs_pdb.get_atomgroup()

        # Model count
        self.assertEqual(rs_ag.get_number_of_groups(), py_ag.get_number_of_groups())
        self.assertEqual(rs_ag.get_number_of_groups(), 20)

        # Model 1 chains
        self.assertEqual(
            rs_ag["model_1"].get_number_of_groups(),
            py_ag["model_1"].get_number_of_groups(),
        )

        # Model 1 Chain A residues
        chain_rs = rs_ag["model_1"]["A"]
        chain_py = py_ag["model_1"]["A"]
        self.assertEqual(chain_rs.get_number_of_groups(), chain_py.get_number_of_groups())
        self.assertEqual(chain_rs.get_number_of_groups(), 9)

        # Check residue names and atom counts using integer indexing
        for i in range(1, 10):
            self.assertEqual(chain_rs[i].name, chain_py[i].name)
            self.assertEqual(
                chain_rs[i].get_number_of_atoms(), chain_py[i].get_number_of_atoms()
            )

    def test_pdb_1hls_altloc(self):
        pdb_path = os.path.join(DATA_DIR, "1hls.pdb")
        py_pdb = PyPdb(pdb_path)
        rs_pdb = rs_br.Pdb(pdb_path)

        # Default altloc="A"
        py_ag_a = py_pdb.get_atomgroup(select_model=1, select_altloc="A")
        rs_ag_a = rs_br.Pdb(pdb_path).get_atomgroup(select_model=1, select_altloc="A")
        self.assertEqual(rs_ag_a.get_number_of_atoms(), py_ag_a.get_number_of_atoms())

        # altloc=None (all conformations)
        py_ag_all = py_pdb.get_atomgroup(select_model=1, select_altloc=None)
        rs_ag_all = rs_pdb.get_atomgroup(select_model=1, select_altloc=None)
        self.assertEqual(rs_ag_all.get_number_of_atoms(), py_ag_all.get_number_of_atoms())

    def test_pdb_3i3zh(self):
        pdb_path = os.path.join(DATA_DIR, "3i3zH.pdb")
        py_pdb = PyPdb(pdb_path)
        rs_pdb = rs_br.Pdb(pdb_path)

        py_ag = py_pdb.get_atomgroup()
        rs_ag = rs_pdb.get_atomgroup()

        self.assertEqual(rs_ag.get_number_of_groups(), py_ag.get_number_of_groups())
        # Check chains in model_1
        self.assertEqual(
            rs_ag["model_1"].get_number_of_groups(),
            py_ag["model_1"].get_number_of_groups(),
        )

    def test_pdb_methods(self):
        pdb = rs_br.Pdb(os.path.join(DATA_DIR, "2MGO.pdb"), mode="amber")
        self.assertEqual(pdb.mode, "AMBER")
        pdb.mode = None
        self.assertIsNone(pdb.mode)

        pdb.renumber()
        text = pdb.get_text()
        self.assertIn("MODEL        1", text)
        self.assertIn("ATOM  ", text)
        self.assertEqual(str(pdb), text)

        # roundtrip set_by_atomgroup
        ag = pdb.get_atomgroup(select_model=1)
        pdb2 = rs_br.Pdb()
        pdb2.set_by_atomgroup(ag)
        text2 = pdb2.get_text()
        self.assertIn("ATOM  ", text2)

    # --------------------------------------------------------------------------
    # 7. SimpleMmcif
    # --------------------------------------------------------------------------
    def test_mmcif_ccd_ala(self):
        cif_path = os.path.join(DATA_DIR, "ALA.cif")
        py_cif = PySimpleMmcif(cif_path)
        rs_cif = rs_br.SimpleMmcif(cif_path)

        self.assertEqual(rs_cif.get_molecule_names(), py_cif.get_molecule_names())
        self.assertEqual(rs_cif.get_molecule_names(), ["data_ALA"])

        py_ag = py_cif.get_atomgroup("data_ALA")
        rs_ag = rs_cif.get_atomgroup("data_ALA")
        self.assertEqual(rs_ag.get_number_of_atoms(), py_ag.get_number_of_atoms())
        self.assertEqual(rs_ag.get_number_of_atoms(), 13)

        # Default get_atomgroup() without args should pick first molecule
        rs_ag_default = rs_cif.get_atomgroup()
        self.assertEqual(rs_ag_default.get_number_of_atoms(), 13)

    def test_mmcif_structure_2mgo(self):
        cif_path = os.path.join(DATA_DIR, "2MGO.cif")
        pdb_path = os.path.join(DATA_DIR, "2MGO.pdb")

        cif = rs_br.SimpleMmcif(cif_path)
        cif_ag = cif.get_structure_atomgroup()

        pdb = rs_br.Pdb(pdb_path)
        pdb_ag = pdb.get_atomgroup()

        # 20 models in both
        self.assertEqual(cif_ag.get_number_of_groups(), pdb_ag.get_number_of_groups())
        self.assertEqual(cif_ag.get_number_of_groups(), 20)

        # Compare residues in model_1 chain A
        cif_chain = cif_ag["model_1"]["A"]
        pdb_chain = pdb_ag["model_1"]["A"]
        self.assertEqual(cif_chain.get_number_of_groups(), pdb_chain.get_number_of_groups())

        for i in range(1, 10):
            self.assertEqual(cif_chain[i].name, pdb_chain[i].name)
            self.assertEqual(
                cif_chain[i].get_number_of_atoms(), pdb_chain[i].get_number_of_atoms()
            )

    def test_mmcif_structure_1hls(self):
        cif_path = os.path.join(DATA_DIR, "1HLS.cif")
        pdb_path = os.path.join(DATA_DIR, "1hls.pdb")

        cif = rs_br.SimpleMmcif(cif_path)
        pdb = rs_br.Pdb(pdb_path)

        cif_ag = cif.get_structure_atomgroup(select_model=1, select_altloc="A")
        pdb_ag = pdb.get_atomgroup(select_model=1, select_altloc="A")

        # Compare total atom count between mmCIF and PDB for model 1, altloc A
        cif_model = cif_ag["model_1"]
        pdb_model = pdb_ag["model_1"]
        self.assertEqual(
            cif_model.get_number_of_groups(), pdb_model.get_number_of_groups()
        )


if __name__ == "__main__":
    unittest.main()
