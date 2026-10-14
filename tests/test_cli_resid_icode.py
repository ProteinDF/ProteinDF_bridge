#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
import os
import shutil
import sys
import tempfile
import unittest

import proteindf_bridge.rs as rs
from proteindf_bridge.cli import brd_renumber_resid, pdb2brd, reorder


class TestCliResidIcode(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.pdb_2fb4 = os.path.abspath(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                "rust",
                "crates",
                "proteindf-bridge",
                "tests",
                "data",
                "2FB4.pdb",
            )
        )
        self.assertTrue(
            os.path.exists(self.pdb_2fb4), f"2FB4.pdb not found at {self.pdb_2fb4}"
        )

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def _run_cli(self, func, args):
        """Execute CLI entry point in-process with captured stdout/stderr."""
        old_argv = sys.argv
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        sys.argv = args
        sys.stdout = io.StringIO()
        sys.stderr = io.StringIO()
        exc = None
        try:
            func()
        except Exception as e:
            exc = e
        finally:
            sys.argv = old_argv
            out = sys.stdout.getvalue()
            err = sys.stderr.getvalue()
            sys.stdout = old_stdout
            sys.stderr = old_stderr
        return out, err, exc

    def test_reorder_get_max_resid_with_2fb4(self):
        pdb = rs.Pdb(self.pdb_2fb4)
        ag = pdb.get_atomgroup()
        model = ag["model_1"]
        max_resid = reorder.get_max_resid(model)
        self.assertEqual(max_resid, 289)

    def test_cli_brd_renumber_resid_with_2fb4(self):
        # 1. Convert 2FB4.pdb to .brd via pdb2brd CLI
        brd_path = os.path.join(self.temp_dir, "2fb4.brd")
        out, err, exc = self._run_cli(
            pdb2brd.main, ["pdb2brd.py", self.pdb_2fb4, brd_path]
        )
        self.assertIsNone(exc)

        orig_ag = rs.load_atomgroup(brd_path)
        orig_h = orig_ag["model_1"]["H"]
        self.assertTrue(orig_h.has_group("101A"))
        orig_h_residue_count = len(list(orig_h.groups()))
        orig_total_residue_count = sum(
            len(list(chain.groups()))
            for _, model in orig_ag.groups()
            for _, chain in model.groups()
        )

        # 2. Renumber with increment 10 via brd_renumber_resid CLI
        renumbered_brd = os.path.join(self.temp_dir, "2fb4_renum.brd")
        out, err, exc = self._run_cli(
            brd_renumber_resid.main,
            ["brd-renumber-resid.py", brd_path, "10", "-o", renumbered_brd],
        )
        self.assertIsNone(exc)

        # 3. Check renumbered atomgroup
        renum_ag = rs.load_atomgroup(renumbered_brd)
        renum_h = renum_ag["model_1"]["H"]

        # '101A' must be renumbered to '111A'
        self.assertTrue(renum_h.has_group("111A"))
        self.assertFalse(renum_h.has_group("101A"))
        self.assertEqual(renum_h["111A"].name, orig_h["101A"].name)

        # Number of residues must not change
        renum_h_residue_count = len(list(renum_h.groups()))
        self.assertEqual(renum_h_residue_count, orig_h_residue_count)

        renum_total_residue_count = sum(
            len(list(chain.groups()))
            for _, model in renum_ag.groups()
            for _, chain in model.groups()
        )
        self.assertEqual(renum_total_residue_count, orig_total_residue_count)

    def test_cli_reorder_with_2fb4(self):
        # Verify that reorder.py CLI runs successfully on 2FB4.brd
        brd_path = os.path.join(self.temp_dir, "2fb4.brd")
        out, err, exc = self._run_cli(
            pdb2brd.main, ["pdb2brd.py", self.pdb_2fb4, brd_path]
        )
        self.assertIsNone(exc)

        reordered_brd = os.path.join(self.temp_dir, "2fb4_reordered.brd")
        out, err, exc = self._run_cli(
            reorder.main, ["reorder.py", brd_path, reordered_brd]
        )
        self.assertIsNone(exc)
        self.assertTrue(os.path.exists(reordered_brd))


if __name__ == "__main__":
    unittest.main()
