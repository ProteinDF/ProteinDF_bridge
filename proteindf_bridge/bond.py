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

from .matrix import SymmetricMatrix
from .atomgroup import AtomGroup
from .periodictable import PeriodicTable

COVALENT_BOND_TOLERANCE = 0.45


def _extract_model(path):
    parts = [p for p in str(path).split("/") if p]
    if not parts:
        return None
    first = parts[0]
    if first.startswith("model_"):
        num_str = first[6:]
        if num_str.isdigit():
            val = int(num_str)
            if val >= 1:
                return val
    return None


def _is_cross_model(m1, m2):
    return m1 is not None and m2 is not None and m1 != m2


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
        if hasattr(mol, "_update_path"):
            mol._update_path(force=True)
        self._list_atoms(mol)
        self._make_distance_matrix()
        self._make_bond_matrix()

        existing_bonds = set()
        for b in mol.get_bond_list():
            p1, p2 = b[0], b[1]
            if p1 <= p2:
                existing_bonds.add((p1, p2))
            else:
                existing_bonds.add((p2, p1))

        num_of_atoms = len(self._atoms)
        for p in range(num_of_atoms):
            for q in range(p):
                b = self._bondmat.get(p, q)
                if b > 0:
                    p1 = self._atoms[p].path
                    p2 = self._atoms[q].path
                    key = (p1, p2) if p1 <= p2 else (p2, p1)
                    if key not in existing_bonds:
                        existing_bonds.add(key)
                        mol.add_bond(self._atoms[p], self._atoms[q], 1)

    def _make_bond_matrix(self):
        """
        Create the bond matrix.
        """
        num_of_atoms = len(self._atoms)
        self._bondmat = SymmetricMatrix(num_of_atoms)
        cov_radii = [
            PeriodicTable.covalent_radius(atom.atomic_number)
            for atom in self._atoms
        ]
        model_ids = [_extract_model(atom.path) for atom in self._atoms]
        for p in range(num_of_atoms):
            cov_p = cov_radii[p]
            model_p = model_ids[p]
            for q in range(p):
                cov_q = cov_radii[q]
                r = self._distmat.get(p, q)
                cross_model = _is_cross_model(model_p, model_ids[q])
                if not cross_model and r <= (cov_p + cov_q) + COVALENT_BOND_TOLERANCE:
                    self._bondmat.set(p, q, 1)
                else:
                    self._bondmat.set(p, q, 0)
        
    def _make_distance_matrix(self):
        """
        Create the distance matrix.
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
