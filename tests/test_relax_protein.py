#!/usr/bin/env python
# -*- coding: utf-8 -*-

import importlib.util
import os
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch

import proteindf_bridge as bridge


def load_relax_module():
    script_path = os.path.abspath(
        os.path.join(
            os.path.dirname(__file__),
            "..",
            "scripts",
            "relax_protein.py",
        )
    )
    spec = importlib.util.spec_from_file_location("relax_protein", script_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


relax_module = load_relax_module()


class TestRelaxExecIon(unittest.TestCase):
    def setUp(self):
        # Relax instance without running full __init__ pipeline
        self.relax = relax_module.Relax.__new__(relax_module.Relax)
        self.relax._logger = MagicMock()

    def test_exec_cmd_success(self):
        # Exit code 0 returns True
        cmdline = f'{sys.executable} -c "import sys; sys.exit(0)"'
        self.assertTrue(self.relax._exec_cmd(cmdline))

    def test_exec_cmd_failures(self):
        # Exit codes 1, 2, and negative (signal) should all raise bridge.BrError
        for exit_code in (1, 2):
            cmdline = f'{sys.executable} -c "import sys; sys.exit({exit_code})"'
            with self.assertRaises(bridge.BrError) as cm:
                self.relax._exec_cmd(cmdline)
            self.assertIn(f"exit code {exit_code}", str(cm.exception))

        # Negative exit code via mock Popen (e.g., killed by SIGKILL = -9)
        mock_proc = MagicMock()
        mock_proc.communicate.return_value = (b"", b"Killed")
        mock_proc.returncode = -9
        with patch("subprocess.Popen", return_value=mock_proc):
            with self.assertRaises(bridge.BrError) as cm:
                self.relax._exec_cmd("dummy_cmd")
            self.assertIn("exit code -9", str(cm.exception))

    def test_get_ion_residues(self):
        # PDB containing Na+ and Cl- with insertion codes
        pdb_content = (
            "ATOM      1  NA  Na+ A  10      10.000  10.000  10.000  1.00  0.00          NA\n"
            "ATOM      2  NA  Na+ A  20A     12.000  12.000  12.000  1.00  0.00          NA\n"
            "ATOM      3  CL  Cl- A  15      15.000  15.000  15.000  1.00  0.00          CL\n"
            "ATOM      4  CL  Cl- A  35B     20.000  20.000  20.000  1.00  0.00          CL\n"
            "END\n"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            orig_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                with open("md3_before.pdb", "w") as f:
                    f.write(pdb_content)

                start, end = self.relax._get_ion_residues()
                self.assertIsInstance(start, int)
                self.assertIsInstance(end, int)
                self.assertEqual(start, 10)
                self.assertEqual(end, 35)
            finally:
                os.chdir(orig_cwd)


if __name__ == "__main__":
    unittest.main()
