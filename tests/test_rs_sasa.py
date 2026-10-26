#!/usr/bin/env python
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: The ProteinDF development team
# SPDX-License-Identifier: GPL-3.0-or-later

"""Test suite for Solvent-Accessible Surface Area (rs.AtomGroup.sasa).

Validates:
1. Single isolated atom: total_area is close to 4 * pi * (r_vdw + probe_radius)^2.
2. Two overlapping atoms: individual atom area is smaller than isolated atom area.
3. Parameter options: probe_radius=0 and varying n_points produce different SASA values.
4. Hydrogen handling: hydrogens are ignored in SASA calculation, atoms count matches heavy atoms.
5. Real structure 1hls.pdb: 403 heavy atoms, 51 residues, total_area close to 3992.15 Å².
6. Hierarchy level consistency: sasa() called from root, model, chain, or residue yields consistent
   results, and tree is not modified (get_raw_data() unchanged).
7. Serialization and formatting: to_dict() produces JSON-serializable structure, __str__ and __repr__.
"""

import json
import math
import os
import timeit
import unittest

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


class TestRsSasa(unittest.TestCase):
    """Test suite for rs.AtomGroup.sasa calculation and data structures."""

    def setUp(self):
        self.pdb_1hls = os.path.join(DATA_DIR, "1hls.pdb")

    def _load_1hls(self) -> rs_br.AtomGroup:
        pdb = rs_br.Pdb()
        pdb.load(self.pdb_1hls)
        return pdb.get_atomgroup()

    # -------------------------------------------------------------------------
    # 1. Single atom: total_area is close to 4 * pi * (r + 1.4)^2
    # -------------------------------------------------------------------------
    def test_single_atom_sasa(self):
        res = rs_br.AtomGroup(name="ALA")
        res.path = "/model_1/A/1/"
        res.set_atom("CA", rs_br.Atom(symbol="C", name="CA", xyz=[0.0, 0.0, 0.0]))

        report = res.sasa()
        self.assertEqual(len(report.atoms), 1)
        self.assertEqual(len(report.residues), 1)

        # C Bondi radius is 1.70, probe is 1.4 -> r_eff = 3.10 Å
        r_eff = 1.70 + 1.4
        expected_area = 4.0 * math.pi * (r_eff**2)

        self.assertAlmostEqual(report.total_area, expected_area, places=6)
        self.assertAlmostEqual(report.atoms[0].area, expected_area, places=6)
        self.assertAlmostEqual(report.residues[0].area, expected_area, places=6)
        self.assertEqual(report.atoms[0].name, "CA")
        self.assertEqual(report.atoms[0].path, "/model_1/A/1/CA")
        self.assertEqual(report.residues[0].comp_id, "ALA")
        self.assertEqual(report.residues[0].path, "/model_1/A/1/")

    # -------------------------------------------------------------------------
    # 2. Two overlapping atoms: individual atom area is smaller than isolated
    # -------------------------------------------------------------------------
    def test_two_overlapping_atoms(self):
        r_eff = 1.70 + 1.4
        single_area = 4.0 * math.pi * (r_eff**2)

        # Two C atoms separated by 2.0 Å (< 2 * 3.10 = 6.20 Å)
        res = rs_br.AtomGroup(name="GLY")
        res.path = "/model_1/A/1/"
        res.set_atom("C1", rs_br.Atom(symbol="C", name="C1", xyz=[0.0, 0.0, 0.0]))
        res.set_atom("C2", rs_br.Atom(symbol="C", name="C2", xyz=[2.0, 0.0, 0.0]))

        report = res.sasa()
        self.assertEqual(len(report.atoms), 2)
        for atom_sasa in report.atoms:
            self.assertGreater(atom_sasa.area, 0.0)
            self.assertLess(atom_sasa.area, single_area)

        self.assertLess(report.total_area, 2.0 * single_area)

        # Analytical exposed area: 4πR² - 2πRh where h = R - d/2
        d = 2.0
        h = r_eff - d / 2.0
        analytical_each = 4.0 * math.pi * (r_eff**2) - 2.0 * math.pi * r_eff * h
        point_step_area = single_area / 92.0
        for atom_sasa in report.atoms:
            self.assertLessEqual(
                abs(atom_sasa.area - analytical_each), 2.5 * point_step_area
            )

    # -------------------------------------------------------------------------
    # 3. Parameter options: probe_radius=0 and varying n_points change values
    # -------------------------------------------------------------------------
    def test_parameter_options(self):
        res = rs_br.AtomGroup(name="GLY")
        res.path = "/model_1/A/1/"
        res.set_atom("C1", rs_br.Atom(symbol="C", name="C1", xyz=[0.0, 0.0, 0.0]))
        res.set_atom("C2", rs_br.Atom(symbol="C", name="C2", xyz=[2.0, 0.0, 0.0]))

        rep_default = res.sasa()
        rep_probe0 = res.sasa(probe_radius=0.0)
        rep_points = res.sasa(n_points=200)

        # probe_radius=0 has smaller sphere radii (1.70 vs 3.10), so total_area must differ significantly
        self.assertNotEqual(rep_default.total_area, rep_probe0.total_area)
        self.assertLess(rep_probe0.total_area, rep_default.total_area)

        # Changing n_points slightly adjusts numerical integration of overlapping spheres
        self.assertNotEqual(rep_default.total_area, rep_points.total_area)

    # -------------------------------------------------------------------------
    # 4. Hydrogen handling: hydrogens ignored, atoms count matches heavy atoms
    # -------------------------------------------------------------------------
    def test_hydrogens_ignored(self):
        heavy = rs_br.AtomGroup(name="GLY")
        heavy.path = "/model_1/A/1/"
        heavy.set_atom("N", rs_br.Atom(symbol="N", name="N", xyz=[0.0, 0.0, 0.0]))
        heavy.set_atom("CA", rs_br.Atom(symbol="C", name="CA", xyz=[1.46, 0.0, 0.0]))
        heavy.set_atom("C", rs_br.Atom(symbol="C", name="C", xyz=[2.0, 1.4, 0.0]))
        heavy.set_atom("O", rs_br.Atom(symbol="O", name="O", xyz=[1.3, 2.4, 0.0]))

        with_h = rs_br.AtomGroup(name="GLY")
        with_h.path = "/model_1/A/1/"
        with_h.set_atom("N", rs_br.Atom(symbol="N", name="N", xyz=[0.0, 0.0, 0.0]))
        with_h.set_atom("CA", rs_br.Atom(symbol="C", name="CA", xyz=[1.46, 0.0, 0.0]))
        with_h.set_atom("C", rs_br.Atom(symbol="C", name="C", xyz=[2.0, 1.4, 0.0]))
        with_h.set_atom("O", rs_br.Atom(symbol="O", name="O", xyz=[1.3, 2.4, 0.0]))
        with_h.set_atom("H", rs_br.Atom(symbol="H", name="H", xyz=[-0.5, -0.8, 0.0]))
        with_h.set_atom(
            "HA2", rs_br.Atom(symbol="H", name="HA2", xyz=[1.8, -0.5, 0.9])
        )
        with_h.set_atom(
            "HA3", rs_br.Atom(symbol="H", name="HA3", xyz=[1.8, -0.5, -0.9])
        )

        rep_heavy = heavy.sasa()
        rep_with_h = with_h.sasa()

        self.assertEqual(len(rep_with_h.atoms), 4)
        self.assertAlmostEqual(rep_with_h.total_area, rep_heavy.total_area, places=6)
        for a_heavy, a_with_h in zip(rep_heavy.atoms, rep_with_h.atoms):
            self.assertEqual(a_heavy.name, a_with_h.name)
            self.assertAlmostEqual(a_heavy.area, a_with_h.area, places=6)

    # -------------------------------------------------------------------------
    # 5. Real structure 1hls: 403 heavy atoms, 51 residues, total_area ~ 3992.15 Å²
    # -------------------------------------------------------------------------
    def test_real_structure_1hls(self):
        ag = self._load_1hls()
        report = ag.sasa()

        self.assertEqual(len(report.atoms), 403)
        self.assertEqual(len(report.residues), 51)
        self.assertAlmostEqual(report.total_area, 3992.15, delta=0.1)

        # Sum of atom areas equals sum of residue areas equals total_area
        atom_sum = sum(a.area for a in report.atoms)
        res_sum = sum(r.area for r in report.residues)
        self.assertAlmostEqual(report.total_area, atom_sum, places=6)
        self.assertAlmostEqual(report.total_area, res_sum, places=6)

        # Benchmark: 1hls sasa() average time over 20 runs
        runs = 20
        total_time = timeit.timeit(lambda: ag.sasa(), number=runs)
        avg_ms = (total_time / runs) * 1000.0
        print(f"\n[BENCHMARK] 1hls sasa() 20-run average: {avg_ms:.2f} ms")

    # -------------------------------------------------------------------------
    # 6. Hierarchy levels & immutability: call from root/model/chain/residue
    # -------------------------------------------------------------------------
    def test_hierarchy_levels_and_immutability(self):
        res = rs_br.AtomGroup(name="ALA")
        res.set_atom("N", rs_br.Atom(symbol="N", name="N", xyz=[0.0, 0.0, 0.0]))
        res.set_atom("CA", rs_br.Atom(symbol="C", name="CA", xyz=[1.45, 0.0, 0.0]))
        res.set_atom("C", rs_br.Atom(symbol="C", name="C", xyz=[2.0, 1.42, 0.0]))
        res.set_atom("O", rs_br.Atom(symbol="O", name="O", xyz=[1.25, 2.39, 0.0]))
        res.set_atom("CB", rs_br.Atom(symbol="C", name="CB", xyz=[1.98, -0.77, 1.21]))

        chain = rs_br.AtomGroup(name="A")
        chain.set_group("1", res)

        model = rs_br.AtomGroup(name="model_1")
        model.set_group("A", chain)

        root = rs_br.AtomGroup()
        root.set_group("model_1", model)

        raw_before = root.get_raw_data()

        rep_root = root.sasa()
        rep_model = root["model_1"].sasa()
        rep_chain = root["model_1"]["A"].sasa()
        rep_res = root["model_1"]["A"]["1"].sasa()

        raw_after = root.get_raw_data()
        self.assertEqual(
            raw_before,
            raw_after,
            "Tree get_raw_data() must not change after sasa()",
        )

        self.assertEqual(len(rep_root.atoms), 5)
        self.assertEqual(len(rep_root.residues), 1)

        self.assertAlmostEqual(rep_root.total_area, rep_model.total_area, places=6)
        self.assertAlmostEqual(rep_root.total_area, rep_chain.total_area, places=6)
        self.assertAlmostEqual(rep_root.total_area, rep_res.total_area, places=6)

        for i in range(len(rep_root.atoms)):
            self.assertEqual(rep_root.atoms[i].name, rep_model.atoms[i].name)
            self.assertAlmostEqual(
                rep_root.atoms[i].area, rep_model.atoms[i].area, places=6
            )

    # -------------------------------------------------------------------------
    # 7. Serialization and string formatting: to_dict(), json.dumps, __str__, __repr__
    # -------------------------------------------------------------------------
    def test_serialization_and_formatting(self):
        res = rs_br.AtomGroup(name="ALA")
        res.path = "/model_1/A/1/"
        res.set_atom("CA", rs_br.Atom(symbol="C", name="CA", xyz=[0.0, 0.0, 0.0]))

        report = res.sasa()

        # to_dict & json.dumps
        d = report.to_dict()
        self.assertIn("total_area", d)
        self.assertIn("atoms", d)
        self.assertIn("residues", d)
        self.assertEqual(len(d["atoms"]), 1)
        self.assertEqual(len(d["residues"]), 1)
        self.assertEqual(d["atoms"][0]["name"], "CA")
        self.assertEqual(d["atoms"][0]["path"], "/model_1/A/1/CA")
        self.assertEqual(d["residues"][0]["comp_id"], "ALA")
        self.assertEqual(d["residues"][0]["path"], "/model_1/A/1/")

        json_str = json.dumps(d)
        loaded = json.loads(json_str)
        self.assertAlmostEqual(loaded["total_area"], report.total_area, places=6)

        # AtomSasa and ResidueSasa to_dict
        atom_dict = report.atoms[0].to_dict()
        res_dict = report.residues[0].to_dict()
        self.assertIsInstance(json.dumps(atom_dict), str)
        self.assertIsInstance(json.dumps(res_dict), str)

        # __str__ and __repr__
        s = str(report)
        r = repr(report)
        self.assertIn("residue", s.lower())
        self.assertIn("atom", s.lower())
        self.assertIn("total_area", s)
        self.assertIn("SasaReport", r)

        atom_repr = repr(report.atoms[0])
        self.assertIn("AtomSasa", atom_repr)
        self.assertIn("CA", atom_repr)

        res_repr = repr(report.residues[0])
        self.assertIn("ResidueSasa", res_repr)
        self.assertIn("ALA", res_repr)


if __name__ == "__main__":
    unittest.main()
