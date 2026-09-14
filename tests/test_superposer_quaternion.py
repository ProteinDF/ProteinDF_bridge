#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from proteindf_bridge.atomgroup import AtomGroup
from proteindf_bridge.atom import Atom
from proteindf_bridge.position import Position
from proteindf_bridge.superposer_quaternion import Superposer_quaternion


class TestSuperposerQuaternion(unittest.TestCase):
    def setUp(self):
        self.ag1 = AtomGroup("mol1")
        self.ag1.set_atom("A1", Atom(name="A1", xyz=Position([1.0, 1.0, 1.0])))
        self.ag1.set_atom("A2", Atom(name="A2", xyz=Position([1.0, -1.0, -1.0])))
        self.ag1.set_atom("A3", Atom(name="A3", xyz=Position([-1.0, 1.0, -1.0])))
        self.ag1.set_atom("A4", Atom(name="A4", xyz=Position([-1.0, -1.0, 1.0])))

        # 平行移動したもの (x+1.0, y+2.0, z+3.0)
        self.ag2 = AtomGroup("mol2")
        self.ag2.set_atom("A1", Atom(name="A1", xyz=Position([2.0, 3.0, 4.0])))
        self.ag2.set_atom("A2", Atom(name="A2", xyz=Position([2.0, 1.0, 2.0])))
        self.ag2.set_atom("A3", Atom(name="A3", xyz=Position([0.0, 3.0, 2.0])))
        self.ag2.set_atom("A4", Atom(name="A4", xyz=Position([0.0, 1.0, 4.0])))

    def test_rmsd_and_calc(self):
        sq = Superposer_quaternion(self.ag1, self.ag2)
        rmsd = sq.rmsd
        self.assertAlmostEqual(rmsd, 0.0, places=5)
        self.assertAlmostEqual(sq.calc(), rmsd)

    def test_rotation_mat(self):
        sq = Superposer_quaternion(self.ag1, self.ag2)
        rot_mat = sq.rotation_mat
        self.assertEqual(rot_mat.rows, 3)
        self.assertEqual(rot_mat.cols, 3)

    def test_superimpose(self):
        sq = Superposer_quaternion(self.ag1, self.ag2)
        superimposed = sq.superimpose(self.ag1)
        self.assertEqual(superimposed.get_number_of_atoms(), 4)
        for key in ["A1", "A2", "A3", "A4"]:
            p1 = superimposed.get_atom(key).xyz
            p2 = self.ag2.get_atom(key).xyz
            dist = p1.distance_from(p2)
            self.assertAlmostEqual(dist, 0.0, places=4)

    def test_arbitrary_axis_rotation(self):
        # 軸 (1,1,1)/sqrt(3) まわりに 37度 回転 + 並進 (2.0, 3.0, 4.0)
        import math
        import numpy as np

        axis = np.array([1.0, 1.0, 1.0]) / math.sqrt(3.0)
        theta = math.radians(37.0)
        c = math.cos(theta)
        s = math.sin(theta)
        C = 1.0 - c
        ux, uy, uz = axis

        rot = np.array([
            [c + ux * ux * C, ux * uy * C - uz * s, ux * uz * C + uy * s],
            [uy * ux * C + uz * s, c + uy * uy * C, uy * uz * C - ux * s],
            [uz * ux * C - uy * s, uz * uy * C + ux * s, c + uz * uz * C],
        ])

        pts1 = [
            [1.2, 2.3, 3.4],
            [4.5, 1.1, 0.2],
            [0.1, 5.6, 2.7],
            [3.3, 0.4, 6.1],
        ]

        ag1 = AtomGroup("mol1")
        ag2 = AtomGroup("mol2")
        trans = np.array([2.0, 3.0, 4.0])

        for i, pt in enumerate(pts1):
            key = f"A{i+1}"
            ag1.set_atom(key, Atom(name=key, xyz=Position(pt)))
            rotated_pt = rot.dot(np.array(pt)) + trans
            ag2.set_atom(key, Atom(name=key, xyz=Position(rotated_pt.tolist())))

        sq = Superposer_quaternion(ag1, ag2)
        rmsd = sq.rmsd
        self.assertLess(rmsd, 1e-10)


if __name__ == "__main__":
    unittest.main()
