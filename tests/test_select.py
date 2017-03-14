#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import pickle

from pdfbridge.atom import Atom
from pdfbridge.atomgroup import AtomGroup
from pdfbridge.select import *

class Select_Tests(unittest.TestCase):
    def setUp(self):
        self.group1 = AtomGroup()
        atom11 = Atom(symbol='C', xyz=" 0.0 0.0 0.0")
        atom12 = Atom(symbol='H', xyz=" 1.0 0.0 0.0")
        atom13 = Atom(symbol='N', xyz=" 2.0 0.0 0.0")
        subgrp1 = AtomGroup()
        subgrp1.set_atom('C1', atom11)
        subgrp1.set_atom('H1', atom12)
        subgrp1.set_atom('N1', atom13)
        self.group1.set_group('sub1', subgrp1)

        atom21 = Atom(symbol='C', xyz=" 0.0 1.0 0.0")
        atom22 = Atom(symbol='H', xyz=" 0.0 2.0 0.0")
        atom23 = Atom(symbol='N', xyz=" 0.0 3.0 0.0")
        subgrp2 = AtomGroup()
        subgrp2.set_atom('C1', atom21)
        subgrp2.set_atom('H1', atom22)
        subgrp2.set_atom('N1', atom23)
        
        self.group1.set_group('sub2', subgrp2)

    def tearDown(self):
        pass

    def test_range1(self):
        sel_range = Select_Range("0.0 0.0 0.0", 0.1)
        sel = self.group1.select(sel_range)
        self.assertEqual(isinstance(sel, AtomGroup), True)
        self.assertEqual(sel.get_number_of_all_atoms(), 1)
        self.assertEqual(sel['sub1']['C1'].symbol, 'C')
        self.assertEqual(sel['sub1']['C1'].path, '/sub1/C1')

    def test_range2(self):
        sel_range = Select_Range("0.0 0.0 0.0", 1.1)
        sel = self.group1.select(sel_range)
        self.assertEqual(isinstance(sel, AtomGroup), True)
        self.assertEqual(sel.get_number_of_all_atoms(), 3)
        self.assertEqual(sel['sub1']['C1'].symbol, 'C')
        self.assertEqual(sel['sub1']['C1'].path, '/sub1/C1')
        self.assertEqual(sel['sub1']['H1'].path, '/sub1/H1')
        self.assertEqual(sel['sub2']['C1'].path, '/sub2/C1')
        
if __name__ == '__main__':
    unittest.main()
