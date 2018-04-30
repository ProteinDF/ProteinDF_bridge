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

class PeriodicTable(object):
    """
    >>> PeriodicTable.get_symbol(1)
    'H'
    >>> PeriodicTable.get_symbol(20)
    'Ca'
    >>> PeriodicTable.get_atomic_number('C')
    6
    >>> PeriodicTable.get_atomic_number('Cu')
    29
    """
    __table = [
        'X',
        'H', 'He',
        'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
        'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
        'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
        'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
        'Cs', 'Ba',
        'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
        'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
        'Fr', 'Ra',
        'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
        'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Fl', 'Lv', 'Uup', 'Uuh', 'Uus', 'Uuo'
        ]

    __atomic_weights = [
        0.0,
        1.008, 4.003,
        6.941, 9.012, 10.81, 12.01, 14.01, 16.00, 19.00, 20.18,
        22.00, 24.31, 26.98, 28.09, 30.97, 32.07, 35.45, 39.95,
        39.10, 40.08, 44.96, 47.87, 50.94, 52.00, 54.94, 55.85, 58.93, 58.69, 63.55, 65.38, 69.72, 72.63, 74.92, 78.97, 79.90, 83.80,
        85.47, 87.62, 88.91, 91.22, 92.91, 95.95, 99.00, 101.1, 102.9, 106.4, 107.9, 112.4, 114.8, 118.7, 121.8, 127.6, 169.9, 131.3,
        132.9, 137.3,
        138.9, 140.1, 140.9, 144.2, 145.0, 150.4, 152.0, 157.3, 158.9, 162.5, 164.9, 167.3, 168.9, 173.1, 175.0,
        178.5, 180.9, 183.8, 186.2, 190.2, 192.2, 195.1, 197.0, 200.6, 204.4, 207.2, 209.0, 210, 210, 222,
        223, 226,
        227, 232.0, 231.0, 238.0, 237, 239, 243, 247, 252, 252, 257, 258, 259, 262, 267, 268, 271, 272, 277, 276, 281, 280, 285, 289, 293
        ]

    __vdw = [
        0.0,
        1.2, 1.4,
        1.82, 0.0, 0.0, 1.70, 1.55, 1.52, 1.47, 1.54,
        2.27, 1.73, 0.0, 2.10, 1.80, 1.80, 1.75, 1.88,
        2.75, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        1.63, 1.40, 1.39, 1.87, 0.0, 1.85, 2.02,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        1.63, 1.72, 1.58, 1.93, 2.17, 0.0, 2.06, 1.98, 2.16
        ]

    @staticmethod
    def get_num_of_atoms():
        return len(PeriodicTable.__table)
    
    @staticmethod
    def get_symbol(atomic_number):
        atomic_number = int(atomic_number)
        try:
            answer = PeriodicTable.__table[atomic_number]
            return answer
        except:
            print("ERROR @PeriodicTable::get_symbol(): not found input:%s." % (atomic_number))
            raise


    @staticmethod
    def get_atomic_number(symbol):
        try:
            symbol = str(symbol)
        except:
            raise

        symbol = symbol.lower()
        symbol = symbol.capitalize()
        
        try:
            answer = PeriodicTable.__table.index(symbol)
            return answer
        except ValueError:
            print("ERROR @PeriodicTable::get_atomic_number(): not found symbol:%s." % (symbol))
            raise
        except:
            raise

    @staticmethod
    def vdw(atom):
        if isinstance(atom, str):
            atom = PeriodicTable.get_atomic_number(atom)
        answer = 0.0
        try:
            answer = PeriodicTable.__vdw[atom]
        except:
            raise
        return answer

    @staticmethod
    def atomic_weight(atom):
        if isinstance(atom, str):
            atom = PeriodicTable.get_atomic_number(atom)
        answer = 0.0
        try:
            answer = PeriodicTable.__atomic_weights[atom]
        except:
            raise
        return answer
    
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
