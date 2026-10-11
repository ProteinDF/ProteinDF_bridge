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
        - Sort and canonicalize bond lists (endpoints order-independent) since heuristic bond
          search uses spatial cell hashing which may iterate cells non-deterministically.
        """
        if isinstance(obj, dict):
            res = {}
            for k, v in sorted(obj.items()):
                if k == "bonds" and isinstance(v, list):
                    norm_b = []
                    for b in v:
                        b_can = TestCliGolden._canonicalize_brd(b)
                        if isinstance(b_can, list) and len(b_can) >= 2:
                            p1, p2 = str(b_can[0]), str(b_can[1])
                            if p1 > p2:
                                norm_b.append([p2, p1] + list(b_can[2:]))
                            else:
                                norm_b.append([p1, p2] + list(b_can[2:]))
                        else:
                            norm_b.append(b_can)
                    norm_b.sort(key=lambda x: (str(x[0]), str(x[1])))
                    res[str(k)] = norm_b
                else:
                    res[str(k)] = TestCliGolden._canonicalize_brd(v)
            return res
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

    def test_cli_mmcif2mol2_scratch_roundtrip(self):
        scratch_path = os.path.join(self.temp_dir, "scratch.mpac")
        out_dir_read = os.path.join(self.temp_dir, "from_scratch")
        os.makedirs(out_dir_read, exist_ok=True)

        # 1. Run with -w to write scratch file
        out_w, err_w, exc_w = self._run_cli(
            mmcif2mol2.main,
            ["mmcif2mol2.py", self.p_cif, "-w", scratch_path],
        )
        self.assertIsNone(exc_w)
        self.assertTrue(os.path.exists(scratch_path))
        self.assertGreater(os.path.getsize(scratch_path), 0)

        # 2. Run with -r to read scratch file and output mol2
        out_r, err_r, exc_r = self._run_cli(
            mmcif2mol2.main,
            ["mmcif2mol2.py", "dummy.cif", "-r", scratch_path, "-o", out_dir_read],
        )
        self.assertIsNone(exc_r)
        mol2_path = os.path.join(out_dir_read, "data_ALA.mol2")
        self.assertTrue(os.path.exists(mol2_path))
        with open(mol2_path, "r", encoding="utf-8") as f:
            mol2_content = f.read()
        self._assert_text_matches_golden(mol2_content, "mmcif2mol2.mol2.txt")

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
        out_path = os.path.join(self.temp_dir, "bonds.brd")
        out, err, exc = self._run_cli(brd_setup_bond.main, ["brd-setup-bond.py", self.brd_1hls, "-o", out_path])
        self.assertIsNone(exc)
        self._assert_brd_matches_golden(out_path, "brd_setup_bond.brd.json")

        # Correctness criteria per TASK_cli-rust-migration.md / docs/cli-migration-notes.md §3.5:
        # 1. Real PDB 1hls: GLU 4 backbone C=O bond is unique and has order 2 (from CCD template DB)
        # 2. Chain B ARG 22 has >= 2 double bonds (C=O and CZ=NH2)
        # 3. Inter-residue peptide bonds (Chain A 1-2, Chain B 1-2) detected with order 1
        # 4. Disulfide bond (Chain A CYS6 SG - CYS11 SG) detected with order 1
        # 5. No duplicate bonds exist across all detected bonds
        import proteindf_bridge.rs as rs
        ag = rs.load_atomgroup(out_path)
        bonds = ag.get_bond_list()

        # 1. GLU 4 C=O bond
        glu4_co = [
            b for b in bonds
            if (('/A/4/' in b[0] and (b[0].endswith('_C') or b[0].endswith('/C'))) and ('/A/4/' in b[1] and (b[1].endswith('_O') or b[1].endswith('/O'))))
            or (('/A/4/' in b[1] and (b[1].endswith('_C') or b[1].endswith('/C'))) and ('/A/4/' in b[0] and (b[0].endswith('_O') or b[0].endswith('/O'))))
        ]
        self.assertEqual(len(glu4_co), 1, "GLU 4 C=O bond should appear exactly once")
        self.assertEqual(glu4_co[0][2], 2, "GLU 4 C=O bond should have order 2 from CCD template")

        # 2. ARG 22 in Chain B double bonds
        arg22_double = [b for b in bonds if ('/B/22/' in b[0] or '/B/22/' in b[1]) and b[2] == 2]
        self.assertGreaterEqual(len(arg22_double), 2, "ARG 22 should have at least 2 double bonds")

        # 3. Peptide bonds
        pep_a = [
            b for b in bonds
            if (('/A/1/' in b[0] and (b[0].endswith('_C') or b[0].endswith('/C')) and '/A/2/' in b[1] and (b[1].endswith('_N') or b[1].endswith('/N'))))
            or (('/A/2/' in b[0] and (b[0].endswith('_N') or b[0].endswith('/N')) and '/A/1/' in b[1] and (b[1].endswith('_C') or b[1].endswith('/C'))))
        ]
        self.assertEqual(len(pep_a), 1, "Peptide bond between Chain A residue 1 C and residue 2 N should be detected")
        self.assertEqual(pep_a[0][2], 1, "Peptide bond order should be 1")

        pep_b = [
            b for b in bonds
            if (('/B/1/' in b[0] and (b[0].endswith('_C') or b[0].endswith('/C')) and '/B/2/' in b[1] and (b[1].endswith('_N') or b[1].endswith('/N'))))
            or (('/B/2/' in b[0] and (b[0].endswith('_N') or b[0].endswith('/N')) and '/B/1/' in b[1] and (b[1].endswith('_C') or b[1].endswith('/C'))))
        ]
        self.assertEqual(len(pep_b), 1, "Peptide bond between Chain B residue 1 C and residue 2 N should be detected")
        self.assertEqual(pep_b[0][2], 1, "Peptide bond order should be 1")

        # 4. Disulfide bond
        ssbond = [
            b for b in bonds
            if (('/A/6/' in b[0] and (b[0].endswith('_SG') or b[0].endswith('/SG')) and '/A/11/' in b[1] and (b[1].endswith('_SG') or b[1].endswith('/SG'))))
            or (('/A/11/' in b[0] and (b[0].endswith('_SG') or b[0].endswith('/SG')) and '/A/6/' in b[1] and (b[1].endswith('_SG') or b[1].endswith('/SG'))))
        ]
        self.assertEqual(len(ssbond), 1, "Disulfide bond between Chain A CYS6 SG and CYS11 SG should be detected")
        self.assertEqual(ssbond[0][2], 1, "Disulfide bond order should be 1")

        # 5. No duplicates
        seen = set()
        for b in bonds:
            pair = tuple(sorted([b[0], b[1]]))
            self.assertNotIn(pair, seen, f"Duplicate bond found: {pair}")
            seen.add(pair)

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

        # Verify physical and chemical properties of neutralized output (782 -> 792 atoms, 10 ions)
        ag = py_fn.load_msgpack(out_path)
        all_atoms = []

        def collect(g):
            for a in g.get("atoms", {}).values():
                all_atoms.append(a)
            for sub in g.get("groups", {}).values():
                collect(sub)

        collect(ag)
        self.assertEqual(len(all_atoms), 792, "1hls (782 atoms) should have 792 atoms after neutralizing 10 charged sites")

        # Ion counts and elements
        na_ions = [a for a in all_atoms if a.get("symbol") == "Na" or a.get("name") == "Na"]
        cl_ions = [a for a in all_atoms if a.get("symbol") == "Cl" or a.get("name") == "Cl"]
        self.assertEqual(len(na_ions), 6, "Expected 6 Na+ counter-ions")
        self.assertEqual(len(cl_ions), 4, "Expected 4 Cl- counter-ions")

        # Total charge
        total_charge = sum(float(a.get("charge", 0.0)) for a in all_atoms)
        self.assertAlmostEqual(total_charge, 0.0, places=4, msg="Overall system charge should be 0.0")

        # Non-bonded distance checks (ions should not clash with other atoms or each other)
        import math
        all_ions = na_ions + cl_ions
        non_ions = [a for a in all_atoms if a not in all_ions]
        for ion in all_ions:
            ix, iy, iz = ion["xyz"]
            # Distance to protein atoms (>= 1.5A to avoid clash)
            for patm in non_ions:
                px, py, pz = patm["xyz"]
                d = math.sqrt((ix - px)**2 + (iy - py)**2 + (iz - pz)**2)
                self.assertGreater(d, 1.5, f"Ion {ion.get('name')} too close to atom {patm.get('name')}: {d:.3f}A")
            # Distance to other ions (>= 5.0A)
            for ion2 in all_ions:
                if ion2 is not ion:
                    i2x, i2y, i2z = ion2["xyz"]
                    d = math.sqrt((ix - i2x)**2 + (iy - i2y)**2 + (iz - i2z)**2)
                    self.assertGreater(d, 5.0, f"Ion {ion.get('name')} too close to ion {ion2.get('name')}: {d:.3f}A")

        # Determinism check across runs
        repeat_path = os.path.join(self.temp_dir, "neutralize_repeat.brd")
        out2, err2, exc2 = self._run_cli(neutralize.main, ["neutralize.py", self.brd_1hls, repeat_path])
        self.assertIsNone(exc2)
        with open(out_path, "rb") as f1, open(repeat_path, "rb") as f2:
            self.assertEqual(f1.read(), f2.read(), "Neutralize output must be strictly deterministic and bit-identical across runs")

    def test_cli_superposer(self):
        out, err, exc = self._run_cli(superposer.main, ["superposer.py", self.p_cis1, self.p_cis2])
        self.assertIsNone(exc)

        # Parse RMSD line: pure Python baseline is 1.6536679544491055, Rust (nalgebra Kabsch SVD) is 1.6536679544491049.
        # Difference is ~6e-16 (1-2 ULPs of 64-bit float representation).
        # We enforce delta=1e-14 so that any meaningful drift is caught while allowing floating-point rounding variance.
        lines = out.strip().splitlines()
        self.assertTrue(lines[0].startswith("rmsd: "), f"First line must start with 'rmsd: ', got {lines[0]}")
        actual_rmsd = float(lines[0].split()[1])
        expected_py_rmsd = 1.6536679544491055
        self.assertAlmostEqual(
            actual_rmsd,
            expected_py_rmsd,
            delta=1e-14,
            msg="Kabsch RMSD must match pure Python baseline within 64-bit float rounding precision (1e-14)"
        )

        # Structure representation matches Rust Display golden output
        self._assert_text_matches_golden(out, "superposer.stdout.txt")

    def test_cli_superposer_quaternion(self):
        out, err, exc = self._run_cli(superposer.main, ["superposer.py", "-q", self.p_cis1, self.p_cis2])
        self.assertIsNone(exc)

        # Parse RMSD line with quaternion: pure Python baseline is 1.7371258580366855, Rust is 1.7371258580366848.
        # Difference is ~6.66e-16 (1-2 ULPs of 64-bit float representation).
        lines = out.strip().splitlines()
        self.assertTrue(lines[0].startswith("rmsd: "), f"First line must start with 'rmsd: ', got {lines[0]}")
        actual_rmsd = float(lines[0].split()[1])
        expected_py_rmsd = 1.7371258580366855
        self.assertAlmostEqual(
            actual_rmsd,
            expected_py_rmsd,
            delta=1e-14,
            msg="Quaternion RMSD must match pure Python baseline within 64-bit float rounding precision (1e-14)"
        )
        self.assertIn("# group name=", out)
        self.assertIn("1_CH3 C (CH3 )", out)


if __name__ == "__main__":
    unittest.main()
