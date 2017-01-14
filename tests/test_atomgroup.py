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

import unittest
import pickle

from pdfbridge.atom import Atom
from pdfbridge.atomgroup import AtomGroup

class AtomGroupTests(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass
        
    def test_init(self):
        group1 = AtomGroup()
        atom1 = Atom(symbol='C')
        atom2 = Atom(symbol='H')
        atom3 = Atom(symbol='N')
        subgrp = AtomGroup()
        subgrp.set_atom('C1', atom1)
        subgrp.set_atom('H1', atom2)
        subgrp.set_atom('N1', atom3)
        group1.set_group('grp', subgrp)

        self.assertEqual(group1['grp']['C1'].symbol, 'C')
        self.assertEqual(group1.get_number_of_atoms(), 0)
        self.assertEqual(group1.get_number_of_groups(), 1)
        self.assertEqual(group1.get_number_of_all_atoms(), 3)
        self.assertAlmostEqual(group1.sum_of_atomic_number(), 14.0)

    def test_pickle(self):
        group1 = AtomGroup()
        atom1 = Atom(symbol='C')
        atom2 = Atom(symbol='H')
        atom3 = Atom(symbol='N')
        subgrp = AtomGroup()
        subgrp.set_atom('C1', atom1)
        subgrp.set_atom('H1', atom2)
        subgrp.set_atom('N1', atom3)
        group1.set_group('grp', subgrp)

        b = pickle.dumps(group1)
        group2 = pickle.loads(b)

        self.assertIsInstance(group2, AtomGroup)
        self.assertEqual(group2.get_number_of_atoms(), 0)
        self.assertEqual(group2.get_number_of_groups(), 1)
        self.assertEqual(group2.get_number_of_all_atoms(), 3)
        self.assertAlmostEqual(group2.sum_of_atomic_number(), 14.0)

    def test_set_atom2(self):
        atom1 = Atom(symbol='C', xyz="1.0 2.0 3.0")
        atom2 = Atom(symbol='C', xyz="1.1 2.1 3.1")

        atomgroup = AtomGroup()
        atomgroup.set_atom('/qwe/C1', atom1)
        atomgroup.set_atom('/qwe/C2', atom2)
        print(atomgroup)
        
        
    def test_get_path_list(self):
        group1 = AtomGroup()
        atom1 = Atom(symbol='C')
        atom2 = Atom(symbol='H')
        atom3 = Atom(symbol='N')
        subgrp = AtomGroup()
        subgrp.set_atom('C1', atom1)
        subgrp.set_atom('H1', atom2)
        subgrp.set_atom('N1', atom3)
        group1.set_group('grp', subgrp)

        path_list = group1.get_path_list()
        self.assertIsInstance(path_list, list)
        self.assertEqual(len(path_list), 3)
        self.assertEqual(path_list[0], '/grp/C1')
        self.assertEqual(path_list[1], '/grp/H1')
        self.assertEqual(path_list[2], '/grp/N1')

        
if __name__ == '__main__':
    unittest.main()
