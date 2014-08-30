#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2014 The ProteinDF development team.
# see also AUTHORS and README if provided.
# 
# This file is a part of the ProteinDF software package.
# 
# The ProteinDF is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# The ProteinDF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

import math
import copy

import pdfbridge

class Superposer(object):
    """
    """

    def __init__(self, atom_group1, atom_group2):
        """
        @param atom_group1[in] RMSD・回転行列を求めたい分子団
        @param atom_group2[in] 基準となる分子団
        """
        self._atomgroup1 = pdfbridge.AtomGroup(atom_group1)
        self._atomgroup2 = pdfbridge.AtomGroup(atom_group2)
        (positions1, positions2) = self._match_positions(atom_group1, atom_group2)
        self._num_of_positions = len(positions1)
        assert(self._num_of_positions == len(positions2))
        self._positions1 = positions1
        self._positions2 = positions2

        self._center1 = None
        self._center2 = None
        self._shift_positions1 = None
        self._shift_positions2 = None
        self._rotation_mat = None
        self._update_positions1 = None
        self._update_positions2 = None
        self._rmsd = None

        #self._calc(self._atomgroup1,
        #           self._atomgroup2)

    #
    def _get_num_of_positions(self):
      return self._num_of_positions

    num_of_positions = property(_get_num_of_positions)
    # -----------------------------------------------------------------
    def _get_positions1(self):
        return self._positions1

    positions1 = property(_get_positions1)
    # -----------------------------------------------------------------
    def _get_positions2(self):
        return self._positions2

    positions2 = property(_get_positions2)
    # -----------------------------------------------------------------
    def _get_center1(self):
        if self._center1 == None:
            self._center1 = self.get_center(self.positions1)
            #print("center1: {}".format(self._center1))
        return self._center1

    center1 = property(_get_center1)
    # -----------------------------------------------------------------
    def _get_center2(self):
        if self._center2 == None:
            self._center2 = self.get_center(self.positions2)
            #print("center2: {}".format(self._center2))
        return self._center2

    center2 = property(_get_center2)
    # -----------------------------------------------------------------
    def _get_shift_positions1(self):
        if self._shift_positions1 == None:
            self._shift_positions1 = self._shift_positions(self.positions1, self.center1)
            #print(self._shift_positions1)
        return self._shift_positions1

    shift_positions1 = property(_get_shift_positions1)
    # -----------------------------------------------------------------
    def _get_shift_positions2(self):
        if self._shift_positions2 == None:
            self._shift_positions2 = self._shift_positions(self.positions2, self.center2)
            #print(self._shift_positions2)
        return self._shift_positions2

    shift_positions2 = property(_get_shift_positions2)
    # -----------------------------------------------------------------
    def _get_rotation_mat(self):
        if self._rotation_mat == None:
            self._rotation_mat = self._get_rotation_matrix(self.num_of_positions,
                                                           self.shift_positions1,
                                                           self.shift_positions2)
        return self._rotation_mat

    rotation_mat = property(_get_rotation_mat)
    # -----------------------------------------------------------------
    def _get_update_positions1(self):
        # positions1 を回転
        if self._update_positions1 == None:
            self._update_positions1 = copy.deepcopy(self.shift_positions1)
            for i in range(len(self.positions1)):
                self._update_positions1[i].rotate(self.rotation_mat)
        return self._update_positions1

    update_positions1 = property(_get_update_positions1)
    # -----------------------------------------------------------------
    def _get_update_positions2(self):
        # positions2 はそのまま
        if self._update_positions2 == None:
            self._update_positions2 = copy.deepcopy(self.shift_positions2)
        return self._update_positions2

    update_positions2 = property(_get_update_positions2)
    # -----------------------------------------------------------------
    def _get_rmsd(self):
        if self._rmsd == None:
            self._rmsd = self._calc_rmsd(self.update_positions1,
                                         self.update_positions2)
        return self._rmsd

    rmsd = property(_get_rmsd)
    # -----------------------------------------------------------------
    def superimpose(self, atomgroup):
        # positions1 を移動
        answer = pdfbridge.AtomGroup(atomgroup)
        answer.shift_by(- self.center1)

        answer.rotate(self.rotation_mat)

        # positions2 のcenterに移動
        answer.shift_by(self.center2)
        
        return answer
        
    # -----------------------------------------------------------------
    def _calc(self, atom_group1, atom_group2):
        (positions1, positions2) = self._match_positions(atom_group1, atom_group2)
        num_of_positions = len(positions1)
        assert(num_of_positions == len(positions2))

        (translation_vct1, self._translation_vct2) = self._fix_positions(positions1, positions2)

        self._rotation_mat = self._get_rotation_matrix(num_of_positions,
                                                       positions1, positions2)

        self._update_positions(positions1, positions2,
                               self._rotation_mat,
                               self._translation_vct2)

        self._rmsd = self._calc_rmsd(positions1, positions2)
        #print('<<<< superposer')

    def _match_positions(self, atom_group1, atom_group2):
        """
        compare two atom_group objects, and select common points.
        return the list of list corresponding two points.
        """
        positions1 = []
        positions2 = []
        for key, ag1 in atom_group1.groups():
            if atom_group2.has_group(key):
                (p1, p2) = self._match_positions(ag1,
                                                 atom_group2.get_group(key))
                positions1 += p1
                positions2 += p2

        for key, atom1 in atom_group1.atoms():
            if atom_group2.has_atom(key):
                #print(str(atom1), str(atom_group2.get_atom(key)))
                positions1 += [atom1.xyz]
                positions2 += [atom_group2.get_atom(key).xyz]

        return (positions1, positions2)

    def _fix_positions(self, positions1, positions2):
        center1 = self.get_center(positions1)
        center2 = self.get_center(positions2)

        translation_vct1 = - center1
        translation_vct2 =   center2

        for i in range(len(positions1)):
            positions1[i] -= center1
        for i in range(len(positions2)):
            positions2[i] -= center2

        return (translation_vct1, translation_vct2)


    def get_center(self, positions):
        """
        return the center position of input positions
        """
        c = pdfbridge.Position()

        num_of_positions = len(positions)
        for i in range(num_of_positions):
            c += positions[i]
        c /= num_of_positions

        return c

    def _shift_positions(self, positions, center):
        answer = [ p - center for p in positions ]

        # check
        sum_of_positions = pdfbridge.Position()
        for p in answer:
            sum_of_positions += p
        assert(sum_of_positions.distance_from() < 1.0E-5)

        return answer

    def _get_rotation_matrix(self, num_of_points,
                             positions1, positions2):
        r = pdfbridge.Matrix(3, 3)

        # r_ij = Sum_over_k{p2(k, i) * p1(k, j)}
        for k in range(num_of_points):
            x1 = positions1[k].x
            y1 = positions1[k].y
            z1 = positions1[k].z
            x2 = positions2[k].x
            y2 = positions2[k].y
            z2 = positions2[k].z

            r.add(0, 0, x2 * x1)
            r.add(0, 1, x2 * y1)
            r.add(0, 2, x2 * z1)
            r.add(1, 0, y2 * x1)
            r.add(1, 1, y2 * y1)
            r.add(1, 2, y2 * z1)
            r.add(2, 0, z2 * x1)
            r.add(2, 1, z2 * y1)
            r.add(2, 2, z2 * z1)

        #print(' > rot mat: r')
        #print(r)
        tr = r.copy()
        tr.transpose()
        #print(' > rot mat: tr')
        #print(tr)
        trr = tr * r
        #print(' > rot mat: trr')
        #print(trr)
        trr = trr.get_symmetric_matrix()
        #print(trr)

        eigval, eigvec = trr.eig()
        #print('eigval')
        #print(eigval)
        #print('eigvec')
        #print(eigvec)
        a = self._make_right_handed(eigvec)
        #print('a')
        #print(a)

        b = pdfbridge.Matrix(3, 3)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    v = r.get(j , k) * a.get(i, k)
                    b.add(i, j, v)

            # normalize b[i]
            w = 0.0
            for j in range(3):
                w += b.get(i, j) * b.get(i, j)

            t = math.sqrt(1.0 / w)
            for j in range(3):
                v = b.get(i, j)
                b.set(i, j, v * t)
        #print('b')
        #print(b)

        # b[2] = b[0] x b[1]
        tmp_vct = self._calc_vector_product(b.get_row_vector(0),
                                            b.get_row_vector(1))
        for i in range(3):
            b.set(2, i, tmp_vct[i])

        # rotation matrix r_ij = b_ki * a_kj
        mat = self._set_rotation(a, b)
        return mat

    def _make_right_handed(self, mat):
        assert(mat.rows == 3)
        assert(mat.cols == 3)

        v1 = pdfbridge.Vector(3)
        v2 = pdfbridge.Vector(3)
        for i in range(3):
            v1[i] = mat.get(0, i)
            v2[i] = mat.get(1, i)

        v3 = self._calc_vector_product(v1, v2)

        answer = pdfbridge.Matrix(3, 3)
        for i in range(3):
            answer.set(0, i, v1[i])
            answer.set(1, i, v2[i])
            answer.set(2, i, v3[i])

        return answer

    def _calc_vector_product(self, v1, v2):
        assert(len(v1) == 3)
        assert(len(v2) == 3)

        v3 = pdfbridge.Vector(3)
        v3[0] = v1[1] * v2[2] - v1[2] * v2[1]
        v3[1] = v1[2] * v2[0] - v1[0] * v2[2]
        v3[2] = v1[0] * v2[1] - v1[1] * v2[0]

        return v3

    def _set_rotation(self, a, b):
        assert(a.rows == 3)
        assert(a.cols == 3)
        assert(b.rows == 3)
        assert(b.cols == 3)

        r = pdfbridge.Matrix(3, 3)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    v = b.get(k, i) * a.get(k, j)
                    r.add(i, j, v)
        return r

    def _update_positions(self, positions1, positions2,
                          rotation_mat, translation_vct2):
        for i in range(len(positions1)):
            positions1[i].rotate(rotation_mat)

        for i in range(len(positions1)):
            positions1[i] += translation_vct2
        for i in range(len(positions2)):
            positions2[i] += translation_vct2

    def _calc_rmsd(self, positions1, positions2):
        """
        calc rmsd.
        store the value to self._rmsd
        """
        num_of_positions = min(len(self.update_positions1),
                               len(self.update_positions2))
        msd = 0.0
        for i in range(num_of_positions):
            msd += positions1[i].square_distance_from(positions2[i])
        msd /= float(num_of_positions)
        rmsd = math.sqrt(msd)

        return rmsd

    def _debug_positions(self, positions):
        output = ""
        for i, p in enumerate(positions):
            output += "[{}] {}\n".format(i, str(p))
        return output
