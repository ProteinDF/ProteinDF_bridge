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

class Bond(object):
    """
    """
    def __init__(self):
        """
        """
        self._atoms = []
        self._distmat = None
        self._bondmat = None

    def setup(self, mol):
        """
        """
        self._list_atoms(mol)
        self._make_distance_matrix()
        self._make_bond_matrix()

        num_of_atoms = len(self._atoms)
        for p in range(num_of_atoms):
            for q in range(p):
                b = self._bondmat.get(p, q)
                if b > 0:
                    mol.add_bond(self._atoms[p], self._atoms[q], b)
        
        
    def _make_bond_matrix(self):
        """
        結合行列を作成します
        """
        num_of_atoms = len(self._atoms)
        self._bondmat = SymmetricMatrix(num_of_atoms)
        for p in range(num_of_atoms):
            vdw_p = self._atoms[p].vdw
            for q in range(p):
                vdw_q = self._atoms[q].vdw
                r = self._distmat.get(p, q)
                if r <= (vdw_p + vdw_q) + 0.4:
                    self._bondmat.set(p, q, 1)
                else:
                    self._bondmat.set(p, q, 0)
        
    def _make_distance_matrix(self):
        """
        距離行列を作成します
        """
        num_of_atoms = len(self._atoms)
        self._distmat = SymmetricMatrix(num_of_atoms)
        for p in range(num_of_atoms):
            for q in range(p):
                d = self._atoms[p].xyz.distance_from(self._atoms[q].xyz)
                self._distmat.set(p, q, d)            
        
    def _list_atoms(self, mol):
        assert(isinstance(mol, AtomGroup))
        for key, atomgroup in mol.groups():
            self._list_atoms(atomgroup)       
        for key, atom in mol.atoms():
            self._atoms.append(atom)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
