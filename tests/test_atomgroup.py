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

from pdfbridge.position import Position
from pdfbridge.atom import Atom
from pdfbridge.atomgroup import AtomGroup
from pdfbridge.select import Select_Range

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
        self.assertEqual(group1['grp']['C1'].path, '/grp/C1')
        self.assertEqual(group1['grp']['H1'].symbol, 'H')
        self.assertEqual(group1['grp']['H1'].path, '/grp/H1')
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

        
    def test_set_atom_by_path(self):
        atom1 = Atom(symbol='C', xyz="1.1 2.1 3.1")
        atom2 = Atom(symbol='C', xyz="1.2 2.2 3.2")
        atom3 = Atom(symbol='C', xyz="1.3 2.3 3.3")

        atomgroup = AtomGroup()
        atomgroup.set_atom('/C1', atom1)
        atomgroup.set_atom('/group_A/C2', atom2)
        atomgroup.set_atom('/group_A/group_B/C3', atom3)
        self.assertIsInstance(atomgroup, AtomGroup)
        self.assertEqual(atomgroup.get_number_of_all_atoms(), 3)
        self.assertEqual(atomgroup.get_number_of_atoms(), 1)
        self.assertEqual(atomgroup.has_groupkey('group_A'), True)
        self.assertEqual(atomgroup['group_A'].get_number_of_all_atoms(), 2)
        self.assertEqual(atomgroup['group_A'].get_number_of_atoms(), 1)
        self.assertEqual(atomgroup['group_A'].has_groupkey('group_B'), True)
        self.assertEqual(atomgroup['group_A']['group_B'].get_number_of_atoms(), 1)

        
    def test_path(self):
        atom10 = Atom(symbol="C")
        atom11 = Atom(symbol="H")
        atom12 = Atom(symbol="H")
        atom13 = Atom(symbol="H")
        atomgroup1 = AtomGroup()
        atomgroup1.set_atom("C", atom10)
        atomgroup1.set_atom("H1", atom11)
        atomgroup1.set_atom("H2", atom12)
        atomgroup1.set_atom("H3", atom13)

        self.assertEqual(atomgroup1.path, "/")
        self.assertEqual(atomgroup1["C"].path, "/C")
        self.assertEqual(atomgroup1["H1"].path, "/H1")
        self.assertEqual(atomgroup1["H2"].path, "/H2")
        self.assertEqual(atomgroup1["H3"].path, "/H3")

        atomgroup2 = AtomGroup()
        atomgroup2.set_group("Me", atomgroup1)
        self.assertEqual(atomgroup2.path, "/")
        self.assertEqual(atomgroup2["Me"]["C"].path, "/Me/C")

        atomgroup3 = AtomGroup()
        atomgroup3.set_group("grp3", atomgroup2)
        self.assertEqual(atomgroup3["grp3"]["Me"]["H3"].path, "/grp3/Me/H3")

    def test_path_copy(self):
        atom10 = Atom(symbol="C")
        atomgroup1 = AtomGroup()
        atomgroup1.set_atom("C", atom10)

        self.assertEqual(atomgroup1.path, "/")
        self.assertEqual(atomgroup1["C"].path, "/C")

        grp_cp = AtomGroup(atomgroup1)
        self.assertEqual(grp_cp.path, "/")
        self.assertEqual(grp_cp["C"].path, "/C")

        atomgroup2 = AtomGroup()
        atomgroup2.set_group("Me", atomgroup1)
        grp_cp2 = AtomGroup(atomgroup2)
        self.assertEqual(grp_cp2.path, "/")
        self.assertEqual(grp_cp2["Me"]["C"].path, "/Me/C")

        
    def test_select_range(self):
        atom10 = Atom(symbol='H', xyz="1.0 0.0 0.0")
        atom11 = Atom(symbol='H', xyz="1.1 0.0 0.0")
        atom12 = Atom(symbol='H', xyz="1.2 0.0 0.0")
        atom13 = Atom(symbol='H', xyz="1.3 0.0 0.0")

        grp1 = AtomGroup()
        grp1.set_atom("H0", atom10)
        grp1.set_atom("H1", atom11)
        grp1.set_atom("H2", atom12)
        grp1.set_atom("H3", atom13)

        grp10 = AtomGroup()
        grp10.set_group("g1", grp1)
        grp100 = AtomGroup()
        grp100.set_group("g10", grp10)
        self.assertEqual(grp100["g10"]["g1"]["H0"].path, "/g10/g1/H0")

        selecter1 = Select_Range(Position(0.0, 0.0, 0.0), 1.01)
        part1 = grp100.select(selecter1)
        self.assertEqual(part1.get_number_of_all_atoms(), 1)
        self.assertEqual(part1["g10"]["g1"]["H0"].path, "/g10/g1/H0")

        selecter2 = Select_Range(Position(0.0, 0.0, 0.0), 2.00)
        part2 = grp100.select(selecter2)
        self.assertEqual(part2.get_number_of_all_atoms(), 4)
        self.assertEqual(part2["g10"]["g1"]["H0"].path, "/g10/g1/H0")
        self.assertEqual(part2["g10"]["g1"]["H1"].path, "/g10/g1/H1")
        self.assertEqual(part2["g10"]["g1"]["H2"].path, "/g10/g1/H2")
        self.assertEqual(part2["g10"]["g1"]["H3"].path, "/g10/g1/H3")
        
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


    def test_divide_path(self):
        parts = AtomGroup.divide_path("atom0")
        self.assertEqual(len(parts), 1)
        self.assertEqual(parts[0], "atom0")

        parts = AtomGroup.divide_path("/res1/atom2")
        self.assertEqual(parts[0], "res1")
        self.assertEqual(parts[1], "atom2")
        
        
if __name__ == '__main__':
    unittest.main()
