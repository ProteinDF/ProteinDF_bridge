#!/usr/bin/env python
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: The ProteinDF development team
# SPDX-License-Identifier: GPL-3.0-or-later

"""Test suite for geometry inspection (rs.AtomGroup.check_geometry).

Validates:
1. Violations beyond the thresholds are reported.
2. Deviations within the thresholds are not reported.
3. Negative or NaN thresholds raise ValueError.
4. Thresholds are keyword-only arguments with default values (10% / 10% / 0.10 A).
5. 1hls with the default CCD: large thresholds give zero violations.
6. include_hydrogens defaults to False and excludes hydrogen-containing angles.
"""

import math
import os
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

HUGE = 1.0e6


def build_ala(shift_cb: float = 0.0) -> rs_br.AtomGroup:
    """Builds an ALA residue from the builtin CCD ideal coordinates, moving CB along x."""
    tmpl = rs_br.CcdTemplateDb.builtin().lookup("ALA")
    assert tmpl is not None
    group = rs_br.AtomGroup("ALA")
    for idx, a in enumerate(tmpl.atoms, start=1):
        if a.is_hydrogen:
            continue
        x, y, z = a.ideal_xyz
        if a.name == "CB":
            x += shift_cb
        group.set_atom(f"{idx}_{a.name}", rs_br.Atom(symbol=a.element, name=a.name, xyz=rs_br.Position(x, y, z)))
    return group


class TestRsGeometry(unittest.TestCase):
    def check(self, group, length=0.05, angle=0.05, planarity=0.1, include_hydrogens=True):
        return group.check_geometry(
            max_bond_length_ratio=length,
            max_bond_angle_ratio=angle,
            max_planarity_deviation=planarity,
            include_hydrogens=include_hydrogens,
        )

    def test_violation_reported(self):
        """A displaced CB is reported as a bond-length violation on CA-CB."""
        rep = self.check(build_ala(shift_cb=1.0))
        self.assertTrue(rep)
        self.assertGreater(len(rep), 0)
        kinds = {v.kind for v in rep.violations}
        self.assertIn("bond_length", kinds)
        v = next(v for v in rep.violations if v.kind == "bond_length" and set(v.atoms) == {"CA", "CB"})
        self.assertEqual(v.comp_id, "ALA")
        self.assertGreater(v.deviation, 0.1)
        self.assertAlmostEqual(v.deviation, abs(v.observed - v.ideal))

    def test_within_threshold_not_reported(self):
        """Ideal coordinates produce no violations, and a small shift stays within a loose threshold."""
        rep = self.check(build_ala())
        self.assertFalse(rep)
        self.assertEqual(len(rep), 0)
        self.assertEqual(len(rep.violations), 0)

        rep = self.check(build_ala(shift_cb=0.05), length=0.5, angle=0.3, planarity=0.5)
        self.assertEqual(len(rep), 0)

    def test_invalid_threshold_raises_value_error(self):
        """Negative or NaN thresholds raise ValueError for each threshold."""
        ala = build_ala()
        for bad in (-0.1, math.nan):
            with self.assertRaises(ValueError):
                self.check(ala, length=bad)
            with self.assertRaises(ValueError):
                self.check(ala, angle=bad)
            with self.assertRaises(ValueError):
                self.check(ala, planarity=bad)

    def test_thresholds_keyword_only(self):
        """Passing thresholds positionally raises TypeError."""
        ala = build_ala()
        with self.assertRaises(TypeError):
            ala.check_geometry(None, 0.1, 0.1, 0.1, True)

    def test_default_thresholds(self):
        """Omitted thresholds use 10% / 10% / 0.10 A."""
        # CA-CB stretched by 1.0 A (about 65%): reported with the defaults.
        rep = build_ala(shift_cb=1.0).check_geometry(include_hydrogens=False)
        self.assertTrue(any(v.kind == "bond_length" and set(v.atoms) == {"CA", "CB"} for v in rep.violations))
        # A small shift (about 0.05 A, about 3%) stays within the defaults.
        rep = build_ala(shift_cb=0.05).check_geometry(include_hydrogens=False)
        self.assertEqual(len(rep), 0)
        # The defaults are the same as passing them explicitly.
        ala = build_ala(shift_cb=0.3)
        explicit = ala.check_geometry(
            max_bond_length_ratio=0.10,
            max_bond_angle_ratio=0.10,
            max_planarity_deviation=0.10,
            include_hydrogens=False,
        )
        default = ala.check_geometry(include_hydrogens=False)
        self.assertEqual(len(default), len(explicit))

    def test_1hls_default_ccd(self):
        """1hls with the default CCD: tight thresholds report violations; huge ones report none."""
        ag = rs_br.SimpleMmcif(os.path.join(DATA_DIR, "1HLS.cif")).get_structure_atomgroup()
        tight = self.check(ag, length=0.0, angle=0.0, planarity=0.0)
        self.assertGreater(len(tight), 0)
        loose = self.check(ag, length=HUGE, angle=HUGE, planarity=HUGE)
        self.assertEqual(len(loose), 0)
        self.assertFalse(loose)

    def test_1hls_exclude_hydrogens(self):
        """With include_hydrogens=False, no angle violation contains a hydrogen atom."""
        ag = rs_br.SimpleMmcif(os.path.join(DATA_DIR, "1HLS.cif")).get_structure_atomgroup()

        def angles(include_hydrogens):
            rep = self.check(ag, length=HUGE, angle=0.05, planarity=HUGE, include_hydrogens=include_hydrogens)
            return [v for v in rep.violations if v.kind == "bond_angle"]

        def has_h(v):
            return any(n.startswith(("H", "D")) for n in v.atoms)

        with_h = angles(True)
        without_h = angles(False)
        self.assertTrue(any(has_h(v) for v in with_h))
        self.assertFalse(any(has_h(v) for v in without_h))
        self.assertEqual(len([v for v in with_h if not has_h(v)]), len(without_h))

    def test_include_hydrogens_default_is_false(self):
        """Omitting include_hydrogens gives the same result as include_hydrogens=False."""
        ag = rs_br.SimpleMmcif(os.path.join(DATA_DIR, "1HLS.cif")).get_structure_atomgroup()
        rep_omitted = ag.check_geometry()
        rep_explicit = ag.check_geometry(include_hydrogens=False)
        self.assertEqual(len(rep_omitted), len(rep_explicit))
        self.assertEqual(
            [(v.kind, v.comp_id, v.atoms) for v in rep_omitted.violations],
            [(v.kind, v.comp_id, v.atoms) for v in rep_explicit.violations],
        )

    def test_include_hydrogens_true_reports_hydrogen_violations(self):
        """include_hydrogens=True reports violations involving hydrogen atoms on 1HLS."""
        ag = rs_br.SimpleMmcif(os.path.join(DATA_DIR, "1HLS.cif")).get_structure_atomgroup()
        rep = ag.check_geometry(include_hydrogens=True)
        h_violations = [v for v in rep.violations if any(n.startswith(("H", "D")) for n in v.atoms)]
        self.assertGreater(len(h_violations), 0)


if __name__ == "__main__":
    unittest.main()
