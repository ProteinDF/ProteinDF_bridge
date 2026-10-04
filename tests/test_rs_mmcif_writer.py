#!/usr/bin/env python
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: The ProteinDF development team
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Test suite for PR#44:
Validating proteindf_bridge.rs (PyO3 Rust bindings) PDBx/mmCIF structure writing
methods on PySimpleMmcif, including write_structure, save_structure, save,
set_by_atomgroup, keyword options, roundtrip fidelity, and error propagation.
"""

import os
import tempfile
import unittest

# Rust PyO3 bindings
import proteindf_bridge.rs as rs_br

DATA_DIR = os.path.abspath(
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "rust",
        "crates",
        "proteindf-bridge",
        "tests",
        "data",
    )
)


class TestRsMmcifWriter(unittest.TestCase):
    """Tests for mmCIF structure writing in PySimpleMmcif (PR#44)."""

    def setUp(self):
        self.tmp_files = []

    def tearDown(self):
        for f in self.tmp_files:
            if os.path.exists(f):
                try:
                    os.remove(f)
                except OSError:
                    pass

    def _get_tmp_path(self, suffix=".cif"):
        fd, path = tempfile.mkstemp(suffix=suffix)
        os.close(fd)
        self.tmp_files.append(path)
        return path

    def _assert_atomgroups_match_roundtrip(self, original, reloaded, tol=1e-3):
        """Helper to assert that original and reloaded AtomGroups match in hierarchy and atom properties."""
        self.assertEqual(
            original.get_number_of_groups(),
            reloaded.get_number_of_groups(),
            "model group count mismatch",
        )
        self.assertEqual(
            original.get_number_of_atoms(),
            reloaded.get_number_of_atoms(),
            "total atom count mismatch",
        )

        for model_key, orig_model in original.groups():
            self.assertTrue(reloaded.has_group(model_key), f"missing {model_key}")
            re_model = reloaded.get_group(model_key)
            self.assertEqual(
                orig_model.get_number_of_groups(),
                re_model.get_number_of_groups(),
                f"chain count mismatch in {model_key}",
            )

            for chain_key, orig_chain in orig_model.groups():
                self.assertTrue(re_model.has_group(chain_key), f"missing {chain_key}")
                re_chain = re_model.get_group(chain_key)
                self.assertEqual(
                    orig_chain.get_number_of_groups(),
                    re_chain.get_number_of_groups(),
                    f"residue count mismatch in {model_key}/{chain_key}",
                )

                for res_key, orig_res in orig_chain.groups():
                    self.assertTrue(re_chain.has_group(res_key), f"missing {res_key}")
                    re_res = re_chain.get_group(res_key)
                    self.assertEqual(orig_res.name, re_res.name)
                    self.assertEqual(
                        orig_res.get_number_of_atoms(),
                        re_res.get_number_of_atoms(),
                        f"atom count mismatch in residue {res_key}",
                    )

                    orig_atoms = {
                        atom.name: atom for _, atom in orig_res.atoms()
                    }
                    re_atoms = {
                        atom.name: atom for _, atom in re_res.atoms()
                    }
                    self.assertEqual(
                        len(orig_atoms),
                        orig_res.get_number_of_atoms(),
                        f"atom names not unique in residue {res_key}",
                    )

                    for atom_name, o_atom in orig_atoms.items():
                        self.assertIn(
                            atom_name,
                            re_atoms,
                            f"missing atom {atom_name} in {model_key}/{chain_key}/{res_key}",
                        )
                        r_atom = re_atoms[atom_name]

                        self.assertEqual(
                            o_atom.atomic_number,
                            r_atom.atomic_number,
                            f"atomic number mismatch for {atom_name}",
                        )
                        self.assertAlmostEqual(
                            o_atom.charge,
                            r_atom.charge,
                            places=3,
                            msg=f"charge mismatch for {atom_name}",
                        )

                        o_pos = o_atom.position
                        r_pos = r_atom.position
                        self.assertAlmostEqual(
                            o_pos.x,
                            r_pos.x,
                            delta=tol,
                            msg=f"x mismatch for {atom_name}",
                        )
                        self.assertAlmostEqual(
                            o_pos.y,
                            r_pos.y,
                            delta=tol,
                            msg=f"y mismatch for {atom_name}",
                        )
                        self.assertAlmostEqual(
                            o_pos.z,
                            r_pos.z,
                            delta=tol,
                            msg=f"z mismatch for {atom_name}",
                        )

    # --------------------------------------------------------------------------
    # 1. Roundtrip tests (TASK PR#44 Completion Definition #1)
    # --------------------------------------------------------------------------
    def test_roundtrip_1hls_real_data(self):
        """Test roundtrip of 1HLS.cif (20-model NMR structure) via save_structure."""
        cif_path = os.path.join(DATA_DIR, "1HLS.cif")
        cif = rs_br.SimpleMmcif(cif_path)
        orig_ag = cif.get_structure_atomgroup()
        self.assertEqual(orig_ag.get_number_of_groups(), 20)

        out_path = self._get_tmp_path()
        rs_br.SimpleMmcif.save_structure(
            orig_ag, out_path, data_block_name="1HLS"
        )

        reloaded_cif = rs_br.SimpleMmcif(out_path)
        reloaded_ag = reloaded_cif.get_structure_atomgroup()

        self._assert_atomgroups_match_roundtrip(orig_ag, reloaded_ag)

    def test_roundtrip_2fb4_insertion_codes(self):
        """Test roundtrip of 2FB4.cif (with insertion codes) via set_by_atomgroup and save."""
        cif_path = os.path.join(DATA_DIR, "2FB4.cif")
        cif = rs_br.SimpleMmcif(cif_path)
        orig_ag = cif.get_structure_atomgroup()

        writer = rs_br.SimpleMmcif()
        writer.set_by_atomgroup(orig_ag, data_block_name="2FB4")

        out_path = self._get_tmp_path()
        writer.save(out_path)

        reloaded_cif = rs_br.SimpleMmcif(out_path)
        reloaded_ag = reloaded_cif.get_structure_atomgroup()

        self._assert_atomgroups_match_roundtrip(orig_ag, reloaded_ag)

    def test_roundtrip_1wct_struct_conn(self):
        """Test roundtrip of 1WCT.cif (with covale and disulf struct_conn) via save(atomgroup=...)."""
        cif_path = os.path.join(DATA_DIR, "1WCT.cif")
        cif = rs_br.SimpleMmcif(cif_path)
        orig_ag = cif.get_structure_atomgroup()

        out_path = self._get_tmp_path()
        writer = rs_br.SimpleMmcif()
        writer.save(out_path, atomgroup=orig_ag, data_block_name="1WCT")

        reloaded_cif = rs_br.SimpleMmcif(out_path)
        reloaded_ag = reloaded_cif.get_structure_atomgroup()

        self._assert_atomgroups_match_roundtrip(orig_ag, reloaded_ag)

        # Check that inter-residue bonds (_struct_conn) are preserved
        orig_bonds = orig_ag.get_bond_list()
        reloaded_bonds = reloaded_ag.get_bond_list()
        self.assertEqual(len(orig_bonds), len(reloaded_bonds))

    # --------------------------------------------------------------------------
    # 2. Text generation & options (write_structure, get_text, options)
    # --------------------------------------------------------------------------
    def test_write_structure_and_get_text(self):
        """Test string generation via write_structure, get_text, and __str__."""
        cif_path = os.path.join(DATA_DIR, "ALA.cif")
        cif = rs_br.SimpleMmcif(cif_path)
        # Build a valid hierarchy for ALA CCD ligand
        ag_ala = cif.get_atomgroup("data_ALA")

        root = rs_br.AtomGroup("root")
        model = rs_br.AtomGroup("1")
        chain = rs_br.AtomGroup("A")
        res = rs_br.AtomGroup("1")
        res.name = "ALA"
        for k, a in ag_ala.atoms():
            res.set_atom(k, a)
        chain.set_group("1", res)
        model.set_group("A", chain)
        root.set_group("model_1", model)

        # Static write_structure
        text_static = rs_br.SimpleMmcif.write_structure(
            root, data_block_name="ALA_BLOCK"
        )
        self.assertIn("data_ALA_BLOCK", text_static)
        self.assertIn("_atom_site.group_PDB", text_static)
        self.assertIn("ATOM", text_static)
        self.assertIn("ALA", text_static)

        # Stateful get_text and __str__
        cif_obj = rs_br.SimpleMmcif()
        cif_obj.set_by_atomgroup(root, data_block_name="ALA_BLOCK")
        text_stateful = cif_obj.get_text()
        self.assertEqual(text_static, text_stateful)
        self.assertEqual(str(cif_obj), text_stateful)

    def test_charge_to_b_factor_options(self):
        """Test charge_to_b_factor keyword option."""
        # Create hierarchy with partial charge
        root = rs_br.AtomGroup("root")
        model = rs_br.AtomGroup("1")
        chain = rs_br.AtomGroup("A")
        res = rs_br.AtomGroup("1")
        res.name = "ALA"

        atom = rs_br.Atom()
        atom.name = "CA"
        atom.atomic_number = 6
        atom.charge = -0.4285
        atom.position = rs_br.Position(1.0, 2.0, 3.0)

        res.set_atom("CA", atom)
        chain.set_group("1", res)
        model.set_group("A", chain)
        root.set_group("model_1", model)

        # Non-integer charge without charge_to_b_factor should raise BrInputError
        with self.assertRaises((rs_br.BrInputError, rs_br.BrError)):
            rs_br.SimpleMmcif.write_structure(root, charge_to_b_factor=False)

        # With charge_to_b_factor=True, it should succeed and write charge to B_iso_or_equiv
        text = rs_br.SimpleMmcif.write_structure(root, charge_to_b_factor=True)
        self.assertIn("-0.4285", text)
        self.assertIn("?", text)  # pdbx_formal_charge is ?

        # In set_by_atomgroup
        cif = rs_br.SimpleMmcif()
        cif.set_by_atomgroup(root, charge_to_b_factor=True)
        self.assertIn("-0.4285", cif.get_text())

    def test_unknown_kwarg_is_charge2tempfactor_raises_type_error(self):
        """Test that passing removed is_charge2tempfactor keyword raises TypeError."""
        root = rs_br.AtomGroup("root")
        model = rs_br.AtomGroup("1")
        chain = rs_br.AtomGroup("A")
        res = rs_br.AtomGroup("1")
        res.name = "ALA"

        atom = rs_br.Atom()
        atom.name = "CA"
        atom.atomic_number = 6
        atom.position = rs_br.Position(1.0, 2.0, 3.0)

        res.set_atom("CA", atom)
        chain.set_group("1", res)
        model.set_group("A", chain)
        root.set_group("model_1", model)

        cif = rs_br.SimpleMmcif()
        out_path = self._get_tmp_path()

        # write_structure should raise TypeError on is_charge2tempfactor
        with self.assertRaises(TypeError):
            rs_br.SimpleMmcif.write_structure(root, is_charge2tempfactor=True)

        # save_structure should raise TypeError on is_charge2tempfactor
        with self.assertRaises(TypeError):
            rs_br.SimpleMmcif.save_structure(root, out_path, is_charge2tempfactor=True)

        # set_by_atomgroup should raise TypeError on is_charge2tempfactor
        with self.assertRaises(TypeError):
            cif.set_by_atomgroup(root, is_charge2tempfactor=True)

        # save should raise TypeError on is_charge2tempfactor
        with self.assertRaises(TypeError):
            cif.save(out_path, atomgroup=root, is_charge2tempfactor=True)

    # --------------------------------------------------------------------------
    # 3. Error propagation (TASK PR#44 Completion Definition #2)
    # --------------------------------------------------------------------------
    def test_error_schema_violation_raises_br_error(self):
        """Test that protein schema violation raises BrInputError/BrError."""
        # Invalid hierarchy: direct atoms in chain (depth 2)
        root = rs_br.AtomGroup("root")
        model = rs_br.AtomGroup("1")
        chain = rs_br.AtomGroup("A")

        atom = rs_br.Atom()
        atom.name = "CA"
        atom.atomic_number = 6
        atom.position = rs_br.Position(0.0, 0.0, 0.0)

        chain.set_atom("CA", atom)
        model.set_group("A", chain)
        root.set_group("model_1", model)

        out_path = self._get_tmp_path()

        # save_structure
        with self.assertRaises((rs_br.BrInputError, rs_br.BrError)):
            rs_br.SimpleMmcif.save_structure(root, out_path)

        # set_by_atomgroup
        cif = rs_br.SimpleMmcif()
        with self.assertRaises((rs_br.BrInputError, rs_br.BrError)):
            cif.set_by_atomgroup(root)

    def test_error_unparseable_residue_key(self):
        """Test that unparseable residue key raises BrInputError/BrError."""
        root = rs_br.AtomGroup("root")
        model = rs_br.AtomGroup("1")
        chain = rs_br.AtomGroup("A")
        # "INVALID" cannot be parsed into seq_id + ins_code
        res = rs_br.AtomGroup("INVALID")
        res.name = "ALA"

        atom = rs_br.Atom()
        atom.name = "CA"
        atom.atomic_number = 6
        atom.position = rs_br.Position(0.0, 0.0, 0.0)

        res.set_atom("CA", atom)
        chain.set_group("INVALID", res)
        model.set_group("A", chain)
        root.set_group("model_1", model)

        out_path = self._get_tmp_path()
        with self.assertRaises((rs_br.BrInputError, rs_br.BrError)):
            rs_br.SimpleMmcif.save_structure(root, out_path)

    def test_error_no_atomgroup_specified(self):
        """Test that calling save() or get_text() without AtomGroup raises BrInputError/BrError."""
        cif = rs_br.SimpleMmcif()
        out_path = self._get_tmp_path()

        with self.assertRaises((rs_br.BrInputError, rs_br.BrError)):
            cif.save(out_path)

        with self.assertRaises((rs_br.BrInputError, rs_br.BrError)):
            cif.get_text()


if __name__ == "__main__":
    unittest.main()
