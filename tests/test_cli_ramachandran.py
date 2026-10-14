#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for brd-ramachandran.py CLI."""

import builtins
import csv
import io
import math
import os
import shutil
import sys
import tempfile
import unittest
from unittest import mock

import proteindf_bridge.rs as rs
from proteindf_bridge.cli import brd_ramachandran


class TestBrdRamachandranCli(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.brd_path = os.path.join(self.temp_dir, "helix.brd")
        self._create_alpha_helix_brd(self.brd_path)

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def _create_alpha_helix_brd(self, path):
        """Construct a minimal 3-residue alpha-helix backbone."""
        def nerf(a, b, c, length, theta_deg, chi_deg):
            theta = math.radians(theta_deg)
            chi = math.radians(chi_deg)

            def sub(x, y):
                return [x[i] - y[i] for i in range(3)]

            def cross(x, y):
                return [
                    x[1] * y[2] - x[2] * y[1],
                    x[2] * y[0] - x[0] * y[2],
                    x[0] * y[1] - x[1] * y[0],
                ]

            def dot(x, y):
                return sum(x[i] * y[i] for i in range(3))

            def norm(x):
                return math.sqrt(dot(x, x))

            bc = [x / norm(sub(c, b)) for x in sub(c, b)]
            ab = sub(b, a)
            n = cross(ab, bc)
            n = [x / norm(n) for x in n]
            nbc = cross(n, bc)
            d = [
                -length * math.cos(theta),
                length * math.sin(theta) * math.cos(chi),
                length * math.sin(theta) * math.sin(chi),
            ]
            return [
                c[i] + d[0] * bc[i] + d[1] * nbc[i] + d[2] * n[i]
                for i in range(3)
            ]

        n1 = [0.0, 0.0, 0.0]
        ca1 = [1.46, 0.0, 0.0]
        c1 = [
            1.46 - 1.52 * math.cos(math.radians(111.0)),
            1.52 * math.sin(math.radians(111.0)),
            0.0,
        ]

        n2 = nerf(n1, ca1, c1, 1.33, 116.0, -47.0)
        ca2 = nerf(ca1, c1, n2, 1.46, 122.0, 180.0)
        c2 = nerf(c1, n2, ca2, 1.52, 111.0, -57.0)

        n3 = nerf(n2, ca2, c2, 1.33, 116.0, -47.0)
        ca3 = nerf(ca2, c2, n3, 1.46, 122.0, 180.0)
        c3 = nerf(c2, n3, ca3, 1.52, 111.0, -57.0)

        chain = rs.AtomGroup("A")
        chain.name = "A"

        for resid, resname, (n, ca, c) in [
            ("1", "ALA", (n1, ca1, c1)),
            ("2", "ALA", (n2, ca2, c2)),
            ("3", "ALA", (n3, ca3, c3)),
        ]:
            res = rs.AtomGroup(resid)
            res.name = resname
            res.set_atom(
                "N", rs.Atom(symbol="N", name="N", position=rs.Position(*n))
            )
            res.set_atom(
                "CA", rs.Atom(symbol="C", name="CA", position=rs.Position(*ca))
            )
            res.set_atom(
                "C", rs.Atom(symbol="C", name="C", position=rs.Position(*c))
            )
            chain.set_group(resid, res)

        model = rs.AtomGroup("model_1")
        model.name = "model_1"
        model.set_group("A", chain)

        models = rs.AtomGroup("models")
        models.set_group("model_1", model)

        rs.save_atomgroup(models, path)

    def test_csv_output_and_plot(self):
        """Test CSV columns and values, and optional plot generation."""
        csv_path = os.path.join(self.temp_dir, "output.csv")
        brd_ramachandran.main([self.brd_path, "-o", csv_path])

        self.assertTrue(os.path.exists(csv_path))

        with open(csv_path, newline="", encoding="utf-8") as f:
            reader = csv.reader(f)
            rows = list(reader)

        # Header check
        self.assertEqual(
            rows[0],
            ["model", "chain", "residue_key", "residue_name", "phi", "psi"],
        )
        self.assertEqual(len(rows), 4)  # header + 3 residues

        # Residue 1: phi is None (empty string), psi is present
        self.assertEqual(rows[1][0], "model_1")
        self.assertEqual(rows[1][1], "A")
        self.assertEqual(rows[1][2], "1")
        self.assertEqual(rows[1][3], "ALA")
        self.assertEqual(rows[1][4], "")
        self.assertNotEqual(rows[1][5], "")

        # Residue 2 (center): phi and psi are both negative
        self.assertEqual(rows[2][2], "2")
        self.assertEqual(rows[2][3], "ALA")
        phi_val = float(rows[2][4])
        psi_val = float(rows[2][5])
        self.assertLess(phi_val, 0.0)
        self.assertLess(psi_val, 0.0)
        self.assertAlmostEqual(phi_val, -57.0, places=1)
        self.assertAlmostEqual(psi_val, -47.0, places=1)

        # Residue 3: phi is present, psi is None (empty string)
        self.assertEqual(rows[3][2], "3")
        self.assertNotEqual(rows[3][4], "")
        self.assertEqual(rows[3][5], "")

        # Test --plot if matplotlib is available
        try:
            import matplotlib  # noqa: F401
        except ImportError:
            self.skipTest("matplotlib not available for --plot test")

        plot_path = os.path.join(self.temp_dir, "ramachandran.png")
        stdout_buf = io.StringIO()
        with mock.patch("sys.stdout", stdout_buf):
            brd_ramachandran.main([self.brd_path, "--plot", plot_path])
        self.assertTrue(os.path.exists(plot_path))
        self.assertGreater(os.path.getsize(plot_path), 0)

    def test_missing_matplotlib_error(self):
        """Test error message and exit code 1 when matplotlib is missing for --plot."""
        plot_path = os.path.join(self.temp_dir, "ramachandran.png")
        stderr_buf = io.StringIO()
        stdout_buf = io.StringIO()

        with mock.patch.dict(
            sys.modules, {"matplotlib": None, "matplotlib.pyplot": None}
        ):
            with mock.patch("sys.stdout", stdout_buf):
                with mock.patch("sys.stderr", stderr_buf):
                    with self.assertRaises(SystemExit) as cm:
                        brd_ramachandran.main(
                            [self.brd_path, "--plot", plot_path]
                        )

        self.assertEqual(cm.exception.code, 1)
        stderr_output = stderr_buf.getvalue()
        self.assertIn("pip install proteindf_bridge[plot]", stderr_output)

    def test_csv_does_not_import_matplotlib(self):
        """Test that running CSV output does not import matplotlib."""
        csv_path = os.path.join(self.temp_dir, "output.csv")
        orig_import = builtins.__import__
        imported = []

        def spy_import(name, *args, **kwargs):
            if "matplotlib" in name:
                imported.append(name)
            return orig_import(name, *args, **kwargs)

        with mock.patch("builtins.__import__", side_effect=spy_import):
            brd_ramachandran.main([self.brd_path, "-o", csv_path])

        self.assertEqual(
            imported, [], f"matplotlib was imported during CSV run: {imported}"
        )

    def test_generate_plot_with_gly_and_pro(self):
        """Test plot generation with GLY and PRO residues."""
        try:
            import matplotlib  # noqa: F401
        except ImportError:
            self.skipTest("matplotlib not available")

        records = [
            {
                "model": "model_1",
                "chain": "A",
                "residue_key": "1",
                "residue_name": "ALA",
                "phi": -57.0,
                "psi": -47.0,
            },
            {
                "model": "model_1",
                "chain": "A",
                "residue_key": "2",
                "residue_name": "GLY",
                "phi": 60.0,
                "psi": 40.0,
            },
            {
                "model": "model_1",
                "chain": "A",
                "residue_key": "3",
                "residue_name": "PRO",
                "phi": -65.0,
                "psi": -30.0,
            },
        ]
        plot_path = os.path.join(self.temp_dir, "test_plot.png")
        brd_ramachandran.generate_plot(records, plot_path)
        self.assertTrue(os.path.exists(plot_path))
        self.assertGreater(os.path.getsize(plot_path), 0)


if __name__ == "__main__":
    unittest.main()
