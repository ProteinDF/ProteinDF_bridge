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

import pdfbridge

class Superposer_quaternion(object):
    def __init__(self, atomgroup1, atomgroup2):
        self._atomgroup1 = pdfbridge.AtomGroup(atomgroup1)
        self._atomgroup2 = pdfbridge.AtomGroup(atomgroup2)

        (self._positions1, self._positions2) = self._match_positions(self._atomgroup1, self._atomgroup2)

        self._center1 = None
        self._center2 = None
        self._r_A = None
        self._r_B = None
        self._va = None
        self._vb = None
        self._matB = None
        self._eigval = None
        self._matR = None
        self._rmsd = None

    # -----------------------------------------------------------------
    def _get_center1(self):
        if self._center1 == None:
            self._center1 = self._calc_center(self._positions1)
        return self._center1

    center1 = property(_get_center1)
    # -----------------------------------------------------------------
    def _get_center2(self):
        if self._center2 == None:
            self._center2 = self._calc_center(self._positions2)
        return self._center2

    center2 = property(_get_center2)
    # -----------------------------------------------------------------
    def _get_r_A(self):
        if self._r_A == None:
            self._r_A = self._shift_positions(self._positions1, self.center1)
            print('>>>> r_A')
            for v in self._r_A:
                print(v)
        return self._r_A

    r_A = property(_get_r_A)
    # -----------------------------------------------------------------
    def _get_r_B(self):
        if self._r_B == None:
            self._r_B = self._shift_positions(self._positions2, self.center2)
            print('>>>> r_B')
            for v in self._r_B:
                print(v)
        return self._r_B

    r_B = property(_get_r_B)
    # -----------------------------------------------------------------
    def _get_va(self):
        if self._va == None:
            r_A = self.r_A
            r_B = self.r_B
            self._va = self._make_va(r_A, r_B)
        return self._va

    va = property(_get_va)
    # -----------------------------------------------------------------
    def _get_vb(self):
        if self._vb == None:
            r_A = self.r_A
            r_B = self.r_B
            self._vb = self._make_vb(r_A, r_B)
        return self._vb

    vb = property(_get_vb)
    # -----------------------------------------------------------------
    def _get_matB(self):
        if self._matB == None:
            va = self.va
            vb = self.vb
            self._matB = self._make_B(va, vb)
            print('>>>> matB')
            print(self._matB)
        return self._matB

    matB = property(_get_matB)
    # -----------------------------------------------------------------
    def _get_eigval(self):
        if self._eigval == None:
            matB = self.matB
            (eigval, eigvec) = matB.eig()
            print('>>>> eigval')
            print(eigval)
            self._eigval = eigval

        return self._eigval

    eigval = property(_get_eigval)
    # -----------------------------------------------------------------
    def _get_matR(self):
        if self._matR == None:
            q = self.eigval
            self._matR = self._make_R(q)
            print('>>>> matR')
            print(self._matR)
        return self._matR

    matR = property(_get_matR)
    # -----------------------------------------------------------------
    def _get_rmsd(self):
        if self._rmsd == None:
            self._rmsd = self._calc_rmsd(
                self.r_A,
                self.r_B,
                self.matR)

        return self._rmsd

    rmsd = property(_get_rmsd)
    # -----------------------------------------------------------------
    def calc(self):
        (positions1, positions2) = self.match_positions(atom_group1, atom_group2)
        num_of_positions1 = len(positions1)
        num_of_positions2 = len(positions2)
        assert(num_of_positions1 == num_of_positions2)

        center1 = self.calc_center(position1)
        center2 = self.calc_center(position2)

        # shift
        r_A = self._shift_positions(position1, center1)
        r_B = self._shift_positions(position2, center2)

        # make a, b
        va = self._make_va(r_A, r_B)
        vb = self._make_vb(r_A, r_B)

        matB = self.make_B(va, vb)
        (eigval, eigvec) = matB.eig()

        matR = self.make_R()


    def _match_positions(self, atomgroup1, atomgroup2):
        """
        compare two atom_group objects, and select common points.
        return the list of list corresponding two points.
        """
        positions1 = []
        positions2 = []
        for key, ag1 in atomgroup1.groups():
            if atomgroup2.has_group(key):
                (p1, p2) = self._match_positions(ag1,
                                                 atomgroup2.get_group(key))
                positions1 += p1
                positions2 += p2

        for key, atom1 in atomgroup1.atoms():
            if atomgroup2.has_atom(key):
                #print(str(atom1), str(atom_group2.get_atom(key)))
                positions1 += [atom1.xyz]
                positions2 += [atomgroup2.get_atom(key).xyz]

        return (positions1, positions2)


    def _calc_center(self, positions):
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
        print('>>>> sum_of_positions:')
        print(sum_of_positions)

        return answer

    def _make_va(self, r_A, r_B):
        num_of_positions = len(r_A)
        assert(num_of_positions == len(r_B))

        va = [ pdfbridge.Position() for x in range(num_of_positions) ]
        for i in range(num_of_positions):
            va[i] = r_B[i] + r_A[i]

        print('>>>> va:')
        for v in va:
            print(v)
        return va

    def _make_vb(self, r_A, r_B):
        num_of_positions = len(r_A)
        assert(num_of_positions == len(r_B))

        vb = [ pdfbridge.Position() for x in range(num_of_positions) ]
        for i in range(num_of_positions):
            vb[i] = r_B[i] - r_A[i]

        print('>>>> vb:')
        for v in vb:
            print(v)
        return vb

    def _make_B(self, va, vb):
        num_of_positions = len(va)
        assert(num_of_positions == len(vb))

        B = pdfbridge.SymmetricMatrix(4)
        for i in range(num_of_positions):
            a = va[i]
            b = vb[i]
            ax = a.x
            ay = a.y
            az = a.z
            bx = b.x
            by = b.y
            bz = b.z
            B.add(0, 0,  bx*bx + by*by + bz*bz)
            B.add(0, 1,  az*by - ay*bz)
            B.add(0, 2, -az*bx + ax*bz)
            B.add(0, 3,  ay*bx - ax*by)
            B.add(1, 1,  bx*bx + ay*ay + az*az)
            B.add(1, 2,  bx*by - ax*ay)
            B.add(1, 3,  bx*bz - ax*az)
            B.add(2, 2,  ax*ax + by*by + az*az)
            B.add(2, 3,  by*bz - ay*az)
            B.add(3, 3,  ax*ax + ay*ay + bz*bz)

        print(B)
        B *= 1.0 / float(num_of_positions * num_of_positions)
        return B

    def _make_R(self, eigvec):
        assert(len(eigvec) == 4)
        q0 = eigvec[0]
        q1 = eigvec[1]
        q2 = eigvec[2]
        q3 = eigvec[3]

        R = Matrix(3, 3)
        R.set(0, 0, 2.0*q0*q0 +2.0*q1*q1 -1.0)
        R.set(0, 1, 2.0*q1*q2 -2.0*q0*q3)
        R.set(0, 2, 2.0*q1*q3 +2.0*q0*q2)

        R.set(1, 0, 2.0*q1*q2 +2.0*q0*q3)
        R.set(1, 1, 2.0*q0*q0 +2.0*q2*q2 -1.0)
        R.set(1, 2, 2.0*q2*q3 -2.0*q0*q1)

        R.set(2, 0, 2.0*q1*q3 -2.0*q0*q2)
        R.set(2, 1, 2.0*q2*q3 +2.0*q0*q1)
        R.set(2, 2, 2.0*q0*q0 +2.0*q3*q3 -1.0)

#        R.set(0, 0, 1.0 - 2.0*q2*q2 -2.0*q3*q3)
#        R.set(0, 1, 2.0*q1*q2 -2.0*q0*q3)
#        R.set(0, 2, 2.0*q1*q3 +2.0*q0*q2)

#        R.set(1, 0, 2.0*q2*q1 +2.0*q0*q3)
#        R.set(1, 1, 1.0 - 2.0*q3*q3 -2.0*q1*q1)
#        R.set(1, 2, 2.0*q2*q3 -2.0*q0*q1)

#        R.set(2, 0, 2.0*q3*q1 -2.0*q0*q2)
#        R.set(2, 1, 2.0*q3*q2 +2.0*q0*q1)
#        R.set(2, 2, 1.0 - 2.0*q1*q1 -2.0*q2*q2)

        return R

    def _calc_rmsd(self, r_A, r_B, R):
        rmsd = 0.0
        num_of_positions = len(r_A)
        assert(num_of_positions == len(r_B))

        for i in range(num_of_positions):
            R_r_A = R * pdfbridge.Vector(r_A[i].xyz)
            r = r_B[i] - pdfbridge.Position(R_r_A)
            rmsd += r.square_distance_from()

        rmsd *= 1.0 / float(num_of_positions)
        rmsd = math.sqrt(rmsd)
        return rmsd
