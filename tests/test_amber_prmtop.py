#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import tempfile
import unittest
from proteindf_bridge.amber_prmtop import AmberPrmtop


class TestAmberPrmtop(unittest.TestCase):
    def test_load_prmtop_and_inpcrd(self):
        prmtop_content = """%VERSION  VERSION_STAMP = V0001.000  DATE = 08/25/26  12:00:00
%FLAG ATOM_NAME
%FORMAT(20a4)
C1  H1  
%FLAG CHARGE
%FORMAT(5E16.8)
 0.00000000E+00 0.00000000E+00
%FLAG ATOMIC_NUMBER
%FORMAT(10I8)
       6       1
"""
        inpcrd_content = """default_name
    2
   0.0000000   0.0000000   0.0000000   1.0900000   0.0000000   0.0000000
"""
        with tempfile.NamedTemporaryFile(suffix=".prmtop", mode="w", delete=False) as f_top:
            f_top.write(prmtop_content)
            top_path = f_top.name

        with tempfile.NamedTemporaryFile(suffix=".inpcrd", mode="w", delete=False) as f_crd:
            f_crd.write(inpcrd_content)
            crd_path = f_crd.name

        try:
            amber = AmberPrmtop(top_path, crd_path)
            self.assertEqual(len(amber.atom_names), 2)
            self.assertEqual(amber.atom_names[0], "C1")
            self.assertEqual(amber.atom_names[1], "H1")
            self.assertEqual(amber.atomic_numbers, [6, 1])
            self.assertEqual(len(amber.xyz), 2)

            ag = amber.get_atomgroup()
            self.assertEqual(ag.get_number_of_atoms(), 2)
        finally:
            if os.path.exists(top_path):
                os.remove(top_path)
            if os.path.exists(crd_path):
                os.remove(crd_path)


if __name__ == "__main__":
    unittest.main()
