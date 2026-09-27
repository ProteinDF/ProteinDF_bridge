#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2002-2014 The ProteinDF project
# see also AUTHORS and README.
#
# This file is part of ProteinDF.
#
# ProteinDF is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ProteinDF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

import math
import unittest

from proteindf_bridge.modeling import Modeling
from proteindf_bridge.position import Position


class ModelingTests(unittest.TestCase):
    def setUp(self):
        self.modeling = Modeling()

    def test_arbitary_rotate_matrix_is_orthonormal_for_generic_vectors(self):
        """Regression test for a bug where the (1, 2) matrix entry used
        `nx * nz` instead of `ny * nz`, producing a non-orthonormal matrix
        (i.e. not a valid rotation) for general, non axis-aligned input
        vectors. Axis-aligned vectors (e.g. one of the two along Z) leave
        nz == 0, masking the bug, so this test deliberately avoids that.
        """
        a = Position(0.5, 0.7, 0.3)
        b = Position(0.1, 0.9, 0.2)
        rot = self.modeling.arbitary_rotate_matrix(a, b)

        # R^T * R must equal the identity matrix for a proper rotation.
        max_off = 0.0
        for i in range(3):
            for j in range(3):
                s = sum(rot.get(k, i) * rot.get(k, j) for k in range(3))
                expected = 1.0 if i == j else 0.0
                max_off = max(max_off, abs(s - expected))
        self.assertLess(max_off, 1.0e-10)

    def test_arbitary_rotate_matrix_maps_in_b_onto_in_a_direction(self):
        """`arbitary_rotate_matrix(in_a, in_b)` returns R such that rotating
        `in_b` by R yields a vector parallel to `in_a` (this is the
        convention actually relied on by add_methyl/get_NH3 callers, which
        rotate a template's own reference axis onto a real bond direction).
        """
        a = Position(0.5, 0.7, 0.3)
        b = Position(0.1, 0.9, 0.2)
        rot = self.modeling.arbitary_rotate_matrix(a, b)

        rotated_b = Position(b)
        rotated_b.rotate(rot)

        la = math.sqrt(a.x * a.x + a.y * a.y + a.z * a.z)
        lb = math.sqrt(b.x * b.x + b.y * b.y + b.z * b.z)
        a_unit = Position(a.x / la, a.y / la, a.z / la)
        rotated_b_unit = Position(
            rotated_b.x / lb, rotated_b.y / lb, rotated_b.z / lb
        )

        err = math.sqrt(
            (rotated_b_unit.x - a_unit.x) ** 2
            + (rotated_b_unit.y - a_unit.y) ** 2
            + (rotated_b_unit.z - a_unit.z) ** 2
        )
        self.assertLess(err, 1.0e-10)

    def test_add_methyl_places_three_hydrogens_at_tetrahedral_bond_length(self):
        """Existing production use of arbitary_rotate_matrix (add_methyl,
        used by neutralize.py / modeling capping helpers): sanity-check
        that the rotated+shifted template still produces the same C-H
        bond length as the hardcoded ethane template itself, regardless
        of the real C1->C2 direction chosen (including a non
        axis-aligned one) -- a rotation+translation must preserve
        distances.
        """
        from proteindf_bridge.atom import Atom

        # Same as modeling.py's own hardcoded ethane template H11 position
        # (C1 at the origin), so this is the expected bond length after any
        # rotation/translation (which preserve distances).
        template_c1_h11_dist = Position(-0.85617, -0.58901, -0.35051).distance_from()

        c1 = Atom(symbol="C", name="C1", position=Position(1.0, 2.0, 3.0))
        c2 = Atom(symbol="C", name="C2", position=Position(2.2, 2.6, 3.9))

        methyl = self.modeling.add_methyl(c1, c2)
        for name in ("H11", "H12", "H13"):
            h = methyl[name]
            dist = (h.xyz - c1.xyz).distance_from()
            self.assertAlmostEqual(dist, template_c1_h11_dist, delta=1.0e-4)


if __name__ == "__main__":
    unittest.main()
