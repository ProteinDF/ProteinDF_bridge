#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Golden test suite for CLI migration (PR#52).
Freezes the output (written files, standard outputs, exceptions) of all 26
structure-processing CLI tools running under the pure Python implementation.

Target CLIs (26 tools):
- Format conversion (10):
    pdb2brd.py, brd2pdb.py, brd2gro.py, gro2brd.py, brd2xyz.py, xyz2brd.py,
    brd2txt.py, mmcif2txt.py, mmcif2mol2.py, read_amber_prmtop.py
- Structure operation & display (14):
    brd-box.py, brd-density.py, brd-formula.py, brd-show-res.py, brd-show-bonds.py,
    brd-setup-bond.py, brd-select.py, brd-select-path.py, brd-divide.py,
    brd-divide-mainchain.py, brd-renumber-resid.py, brd-restructure.py,
    reorder.py, crystallize.py
- Modeling (2):
    neutralize.py, superposer.py

Excluded tools (6 tools):
    db2txt.py, doctest_runner.py, module_inspect.py, mpac2txt.py, mpac2yml.py, yml2mpac.py
    Reason: General-purpose utilities that do not handle molecular structure data
    (per RUST_PORT_SPEC.md §4.5 scope decision).
"""

import io
import json
import os
import re
import shutil
import sys
import tempfile
import unittest

import proteindf_bridge.functions as py_fn
from proteindf_bridge.cli import (
    brd2gro,
    brd2pdb,
    brd2txt,
    brd2xyz,
    brd_box,
    brd_density,
    brd_divide,
    brd_divide_mainchain,
    brd_formula,
    brd_renumber_resid,
    brd_restructure,
    brd_select,
    brd_select_path,
    brd_setup_bond,
    brd_show_bonds,
    brd_show_res,
    crystallize,
    gro2brd,
    mmcif2mol2,
    mmcif2txt,
    neutralize,
    pdb2brd,
    read_amber_prmtop,
    reorder,
    superposer,
    xyz2brd,
)


class TestCliGolden(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        test_dir = os.path.dirname(os.path.abspath(__file__))
        cls.repo_root = os.path.abspath(os.path.join(test_dir, ".."))
        cls.data_dir = os.path.join(cls.repo_root, "proteindf_bridge", "data")
        cls.golden_dir = os.path.join(test_dir, "cli_golden")

        # Standard fixtures
        cls.p_1hls_pdb = os.path.join(cls.data_dir, "1hls.pdb")
        cls.p_2mgo_pdb = os.path.join(cls.data_dir, "2MGO.pdb")
        cls.p_gro = os.path.join(cls.data_dir, "sample.gro")
        cls.p_xyz = os.path.join(cls.data_dir, "ACE_ALA_NME.xyz")
        cls.p_cif = os.path.join(cls.data_dir, "ALA.cif")
        cls.p_cis1 = os.path.join(cls.data_dir, "ACE_ALA_NME_cis1.brd")
        cls.p_cis2 = os.path.join(cls.data_dir, "ACE_ALA_NME_cis2.brd")

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        # Generate standard protein .brd inputs for tests:
        # - 1hls.brd (782 atoms, standard protein model/chain/residue hierarchy, charged residues)
        # - 2mgo.brd (136 atoms, NMR model, has explicit SSBOND disulfide bonds)
        self.brd_1hls = os.path.join(self.temp_dir, "1hls_input.brd")
        self._run_cli(pdb2brd.main, ["pdb2brd.py", self.p_1hls_pdb, self.brd_1hls])

        self.brd_2mgo = os.path.join(self.temp_dir, "2mgo_input.brd")
        self._run_cli(pdb2brd.main, ["pdb2brd.py", self.p_2mgo_pdb, self.brd_2mgo])

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def _run_cli(self, func, args):
        """
        Execute CLI entry point in-process with captured stdout/stderr.
        """
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

    @staticmethod
    def _canonicalize_brd(obj):
        """
        Recursively normalize a loaded .brd data dictionary.
        - Sort dictionary keys
        - Round float values to 6 decimal places to prevent floating-point representation drift
        """
        if isinstance(obj, dict):
            return {str(k): TestCliGolden._canonicalize_brd(v) for k, v in sorted(obj.items())}
        elif isinstance(obj, list):
            return [TestCliGolden._canonicalize_brd(x) for x in obj]
        elif isinstance(obj, float):
            return round(obj, 6)
        else:
            return obj

    @staticmethod
    def _normalize_text(content):
        """
        Normalize text:
        - Convert CRLF to LF
        - Strip trailing whitespace on each line
        - Strip outer leading/trailing empty lines, ensure single final newline
        """
        lines = [line.rstrip() for line in content.strip().splitlines()]
        return "\n".join(lines) + "\n"

    def _assert_brd_matches_golden(self, actual_path, golden_filename):
        self.assertTrue(os.path.exists(actual_path), f"Output file does not exist: {actual_path}")
        actual_data = py_fn.load_msgpack(actual_path)
        actual_canon = self._canonicalize_brd(actual_data)

        golden_path = os.path.join(self.golden_dir, golden_filename)
        self.assertTrue(os.path.exists(golden_path), f"Golden file does not exist: {golden_path}")
        with open(golden_path, "r", encoding="utf-8") as f:
            expected_canon = json.load(f)

        self.assertEqual(actual_canon, expected_canon)

    def _assert_text_matches_golden(self, actual_text, golden_filename):
        norm_actual = self._normalize_text(actual_text)
        golden_path = os.path.join(self.golden_dir, golden_filename)
        self.assertTrue(os.path.exists(golden_path), f"Golden file does not exist: {golden_path}")
        with open(golden_path, "r", encoding="utf-8") as f:
            expected_text = f.read()

        self.assertEqual(norm_actual, expected_text)

    # =========================================================================
    # Group 1: Format Conversion (10 tools)
    # =========================================================================

    def test_cli_pdb2brd(self):
        out_path = os.path.join(self.temp_dir, "pdb2brd.brd")
        out, err, exc = self._run_cli(pdb2brd.main, ["pdb2brd.py", self.p_1hls_pdb, out_path])
        self.assertIsNone(exc)
        self._assert_brd_matches_golden(out_path, "pdb2brd.brd.json")
        self._assert_text_matches_golden(out, "pdb2brd.stdout.txt")

    def test_cli_brd2pdb(self):
        out_path = os.path.join(self.temp_dir, "brd2pdb.pdb")
        out, err, exc = self._run_cli(brd2pdb.main, ["brd2pdb.py", self.brd_1hls, "-o", out_path])
        self.assertIsNone(exc)
        with open(out_path, "r", encoding="utf-8") as f:
            pdb_content = f.read()
        self._assert_text_matches_golden(pdb_content, "brd2pdb.pdb.txt")
        self._assert_text_matches_golden(out, "brd2pdb.stdout.txt")

    def test_cli_brd2gro(self):
        out, err, exc = self._run_cli(brd2gro.main, ["brd2gro.py", self.brd_1hls])
        self.assertIsNone(exc)
        self._assert_text_matches_golden(out, "brd2gro.stdout.txt")

    def test_cli_gro2brd(self):
        out_path = os.path.join(self.temp_dir, "gro2brd.brd")
        out, err, exc = self._run_cli(gro2brd.main, ["gro2brd.py", self.p_gro, out_path])
        self.assertIsNone(exc)
        self._assert_brd_matches_golden(out_path, "gro2brd.brd.json")
        self._assert_text_matches_golden(out, "gro2brd.stdout.txt")

    def test_cli_brd2xyz(self):
        out_path = os.path.join(self.temp_dir, "brd2xyz.xyz")
        out, err, exc = self._run_cli(brd2xyz.main, ["brd2xyz.py", self.brd_1hls, "-o", out_path])
        self.assertIsNone(exc)
        with open(out_path, "r", encoding="utf-8") as f:
            xyz_content = f.read()
        self._assert_text_matches_golden(xyz_content, "brd2xyz.xyz.txt")
        self._assert_text_matches_golden(out, "brd2xyz.stdout.txt")

    def test_cli_xyz2brd(self):
        out_path = os.path.join(self.temp_dir, "xyz2brd.brd")
        out, err, exc = self._run_cli(xyz2brd.main, ["xyz2brd.py", self.p_xyz, out_path])
        self.assertIsNone(exc)
        self._assert_brd_matches_golden(out_path, "xyz2brd.brd.json")
        self._assert_text_matches_golden(out, "xyz2brd.stdout.txt")

    def test_cli_brd2txt(self):
        out_path = os.path.join(self.temp_dir, "brd2txt.csv")
        out, err, exc = self._run_cli(brd2txt.main, ["brd2txt.py", self.brd_1hls, "-c", out_path])
        self.assertIsNone(exc)
        with open(out_path, "r", encoding="utf-8") as f:
            csv_content = f.read()
        self._assert_text_matches_golden(csv_content, "brd2txt.csv.txt")
        self._assert_text_matches_golden(out, "brd2txt.stdout.txt")

    def test_cli_mmcif2txt(self):
        out_path = os.path.join(self.temp_dir, "mmcif2txt.mpac")
        out, err, exc = self._run_cli(mmcif2txt.main, ["mmcif2txt.py", self.p_cif, "-w", out_path])
        self.assertIsNone(exc)
        self._assert_brd_matches_golden(out_path, "mmcif2txt.mpac.json")
        self._assert_text_matches_golden(out, "mmcif2txt.stdout.txt")

    def test_cli_mmcif2mol2(self):
        out, err, exc = self._run_cli(mmcif2mol2.main, ["mmcif2mol2.py", self.p_cif, "-o", self.temp_dir])
        self.assertIsNone(exc)
        mol2_path = os.path.join(self.temp_dir, "data_ALA.mol2")
        self.assertTrue(os.path.exists(mol2_path))
        with open(mol2_path, "r", encoding="utf-8") as f:
            mol2_content = f.read()
        self._assert_text_matches_golden(mol2_content, "mmcif2mol2.mol2.txt")
        self._assert_text_matches_golden(out, "mmcif2mol2.stdout.txt")

    def test_cli_read_amber_prmtop(self):
        top_path = os.path.join(self.temp_dir, "test.prmtop")
        crd_path = os.path.join(self.temp_dir, "test.inpcrd")
        with open(top_path, "w", encoding="utf-8") as f:
            f.write(
                "%VERSION  VERSION_STAMP = V0001.000  DATE = 08/25/26  12:00:00\n"
                "%FLAG ATOM_NAME\n"
                "%FORMAT(20a4)\n"
                "C1  H1  \n"
                "%FLAG CHARGE\n"
                "%FORMAT(5E16.8)\n"
                " 0.00000000E+00 0.00000000E+00\n"
                "%FLAG ATOMIC_NUMBER\n"
                "%FORMAT(10I8)\n"
                "       6       1\n"
            )
        with open(crd_path, "w", encoding="utf-8") as f:
            f.write(
                "default_name\n"
                "    2\n"
                "   0.0000000   0.0000000   0.0000000   1.0900000   0.0000000   0.0000000\n"
            )
        orig_cwd = os.getcwd()
        os.chdir(self.temp_dir)
        try:
            out, err, exc = self._run_cli(read_amber_prmtop.main, ["read_amber_prmtop.py", "test.prmtop", "test.inpcrd"])
            self.assertIsNone(exc)
            self._assert_brd_matches_golden("output.brd", "read_amber_prmtop.brd.json")
            self._assert_text_matches_golden(out, "read_amber_prmtop.stdout.txt")
        finally:
            os.chdir(orig_cwd)

    # =========================================================================
    # Group 2: Structure Operation & Display (14 tools)
    # =========================================================================

    def test_cli_brd_box(self):
        out, err, exc = self._run_cli(brd_box.main, ["brd-box.py", self.brd_1hls, "1", "1.0"])
        self.assertIsNone(exc)
        self._assert_text_matches_golden(out, "brd_box.stdout.txt")

    def test_cli_brd_density(self):
        out, err, exc = self._run_cli(brd_density.main, ["brd-density.py", self.brd_1hls, "1", "--box", "30.0", "30.0", "30.0"])
        self.assertIsNone(exc)
        self._assert_text_matches_golden(out, "brd_density.stdout.txt")

    def test_cli_brd_formula(self):
        out, err, exc = self._run_cli(brd_formula.main, ["brd-formula.py", self.brd_1hls])
        self.assertIsNone(exc)
        self._assert_text_matches_golden(out, "brd_formula.stdout.txt")

    def test_cli_brd_show_res(self):
        out, err, exc = self._run_cli(brd_show_res.main, ["brd-show-res.py", self.brd_1hls])
        self.assertIsNone(exc)
        self._assert_text_matches_golden(out, "brd_show_res.stdout.txt")

    def test_cli_brd_show_bonds(self):
        # 2MGO has SSBOND records and yields explicit bonds in brd
        out, err, exc = self._run_cli(brd_show_bonds.main, ["brd-show-bonds.py", self.brd_2mgo])
        self.assertIsNone(exc)
        self._assert_text_matches_golden(out, "brd_show_bonds.stdout.txt")

    def test_cli_brd_setup_bond(self):
        # Pure Python implementation has AttributeError: module 'proteindf_bridge' has no attribute 'Bond'
        # (Bond import is commented out in proteindf_bridge/__init__.py).
        # We record and test this expected exception so that the behavior is pinned until PR#54 migration.
        out_path = os.path.join(self.temp_dir, "bonds.brd")
        out, err, exc = self._run_cli(brd_setup_bond.main, ["brd-setup-bond.py", self.brd_1hls, "-o", out_path])
        self.assertIsNotNone(exc)
        self.assertIsInstance(exc, AttributeError)
        exc_str = f"{type(exc).__name__}: {exc}"
        self._assert_text_matches_golden(exc_str, "brd_setup_bond.exception.txt")

    def test_cli_brd_select(self):
        out_path = os.path.join(self.temp_dir, "select.brd")
        out, err, exc = self._run_cli(brd_select.main, ["brd-select.py", self.brd_1hls, "-q", "*/A/1/*", "-o", out_path])
        self.assertIsNone(exc)
        self._assert_brd_matches_golden(out_path, "brd_select.brd.json")
        self._assert_text_matches_golden(out, "brd_select.stdout.txt")

    def test_cli_brd_select_path(self):
        out_path = os.path.join(self.temp_dir, "select_path.brd")
        out, err, exc = self._run_cli(brd_select_path.main, ["brd-select-path.py", self.brd_1hls, "-q", "/model_1/A/1/", "-o", out_path])
        self.assertIsNone(exc)
        self._assert_brd_matches_golden(out_path, "brd_select_path.brd.json")
        self._assert_text_matches_golden(out, "brd_select_path.stdout.txt")

    def test_cli_brd_divide(self):
        out_path = os.path.join(self.temp_dir, "divide.brd")
        out, err, exc = self._run_cli(brd_divide.main, ["brd-divide.py", self.brd_1hls, "-o", out_path])
        self.assertIsNone(exc)
        self._assert_brd_matches_golden(out_path, "brd_divide.brd.json")
        self._assert_text_matches_golden(out, "brd_divide.stdout.txt")

    def test_cli_brd_divide_mainchain(self):
        out_path = os.path.join(self.temp_dir, "divide_mc.brd")
        out, err, exc = self._run_cli(brd_divide_mainchain.main, ["brd-divide-mainchain.py", self.brd_1hls, "-o", out_path])
        self.assertIsNone(exc)
        self._assert_brd_matches_golden(out_path, "brd_divide_mainchain.brd.json")
        self._assert_text_matches_golden(out, "brd_divide_mainchain.stdout.txt")

    def test_cli_brd_renumber_resid(self):
        out_path = os.path.join(self.temp_dir, "renumber_resid.brd")
        out, err, exc = self._run_cli(brd_renumber_resid.main, ["brd-renumber-resid.py", self.brd_1hls, "10", "-o", out_path])
        self.assertIsNone(exc)
        self._assert_brd_matches_golden(out_path, "brd_renumber_resid.brd.json")
        self._assert_text_matches_golden(out, "brd_renumber_resid.stdout.txt")

    def test_cli_brd_restructure(self):
        out_path = os.path.join(self.temp_dir, "restructure.brd")
        out, err, exc = self._run_cli(brd_restructure.main, ["brd-restructure.py", self.brd_1hls, self.brd_1hls, "-o", out_path])
        self.assertIsNone(exc)
        self._assert_brd_matches_golden(out_path, "brd_restructure.brd.json")
        self._assert_text_matches_golden(out, "brd_restructure.stdout.txt")

    def test_cli_reorder(self):
        out_path = os.path.join(self.temp_dir, "reorder.brd")
        out, err, exc = self._run_cli(reorder.main, ["reorder.py", self.brd_1hls, out_path])
        self.assertIsNone(exc)
        self._assert_brd_matches_golden(out_path, "reorder.brd.json")
        self._assert_text_matches_golden(out, "reorder.stdout.txt")

    def test_cli_crystallize(self):
        out_path = os.path.join(self.temp_dir, "crystallize.brd")
        out, err, exc = self._run_cli(crystallize.main, ["crystallize.py", self.brd_1hls, out_path, "--num_x", "1", "--num_y", "1", "--num_z", "1"])
        self.assertIsNone(exc)
        self._assert_brd_matches_golden(out_path, "crystallize.brd.json")
        # crystallize.py prints args Namespace containing dynamic input/output file paths.
        # Normalize paths so the test is deterministic across runs and machines.
        norm_out = re.sub(r"INPUT_BRD_PATH=\['[^']*'\]", "INPUT_BRD_PATH=['<INPUT_BRD_PATH>']", out)
        norm_out = re.sub(r"OUTPUT_BRD_PATH=\['[^']*'\]", "OUTPUT_BRD_PATH=['<OUTPUT_BRD_PATH>']", norm_out)
        self._assert_text_matches_golden(norm_out, "crystallize.stdout.txt")

    # =========================================================================
    # Group 3: Modeling (2 tools)
    # =========================================================================

    def test_cli_neutralize(self):
        out_path = os.path.join(self.temp_dir, "neutralize.brd")
        out, err, exc = self._run_cli(neutralize.main, ["neutralize.py", self.brd_1hls, out_path])
        self.assertIsNone(exc)
        self._assert_brd_matches_golden(out_path, "neutralize.brd.json")
        self._assert_text_matches_golden(out, "neutralize.stdout.txt")

    def test_cli_superposer(self):
        out, err, exc = self._run_cli(superposer.main, ["superposer.py", self.p_cis1, self.p_cis2])
        self.assertIsNone(exc)
        self._assert_text_matches_golden(out, "superposer.stdout.txt")


if __name__ == "__main__":
    unittest.main()
