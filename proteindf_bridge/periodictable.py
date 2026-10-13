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

import logging
logger = logging.getLogger(__name__)


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
        except IndexError:
            logger.error("PeriodicTable.get_symbol(): atomic number %s not found.", atomic_number)
            raise

    @staticmethod
    def get_atomic_number(symbol):
        symbol = str(symbol).lower().capitalize()
        try:
            answer = PeriodicTable.__table.index(symbol)
            return answer
        except ValueError:
            logger.error("PeriodicTable.get_atomic_number(): symbol '%s' not found.", symbol)
            raise

    def __contains__(self, symbol):
        symbol = str(symbol)
        symbol = symbol.lower()
        symbol = symbol.capitalize()

        answer = False
        for i in PeriodicTable.__table:
            if symbol == i:
                answer = True
                break
        return answer

    __covalent_radius = [
        0.0,  # X (dummy)
        0.31, # 1: H
        0.28, # 2: He
        1.28, # 3: Li
        0.96, # 4: Be
        0.84, # 5: B
        0.76, # 6: C
        0.71, # 7: N
        0.66, # 8: O
        0.57, # 9: F
        0.58, # 10: Ne
        1.66, # 11: Na
        1.41, # 12: Mg
        1.21, # 13: Al
        1.11, # 14: Si
        1.07, # 15: P
        1.05, # 16: S
        1.02, # 17: Cl
        1.06, # 18: Ar
        2.03, # 19: K
        1.76, # 20: Ca
        1.70, # 21: Sc
        1.60, # 22: Ti
        1.53, # 23: V
        1.39, # 24: Cr
        1.39, # 25: Mn
        1.32, # 26: Fe
        1.26, # 27: Co
        1.24, # 28: Ni
        1.32, # 29: Cu
        1.22, # 30: Zn
        1.22, # 31: Ga
        1.20, # 32: Ge
        1.19, # 33: As
        1.20, # 34: Se
        1.20, # 35: Br
        1.16, # 36: Kr
        2.20, # 37: Rb
        1.95, # 38: Sr
        1.90, # 39: Y
        1.75, # 40: Zr
        1.64, # 41: Nb
        1.54, # 42: Mo
        1.47, # 43: Tc
        1.46, # 44: Ru
        1.42, # 45: Rh
        1.39, # 46: Pd
        1.45, # 47: Ag
        1.44, # 48: Cd
        1.42, # 49: In
        1.39, # 50: Sn
        1.39, # 51: Sb
        1.38, # 52: Te
        1.39, # 53: I
        1.40, # 54: Xe
        2.44, # 55: Cs
        2.15, # 56: Ba
        2.07, # 57: La
        2.04, # 58: Ce
        2.03, # 59: Pr
        2.01, # 60: Nd
        1.99, # 61: Pm
        1.98, # 62: Sm
        1.98, # 63: Eu
        1.96, # 64: Gd
        1.94, # 65: Tb
        1.92, # 66: Dy
        1.92, # 67: Ho
        1.89, # 68: Er
        1.90, # 69: Tm
        1.87, # 70: Yb
        1.87, # 71: Lu
        1.75, # 72: Hf
        1.70, # 73: Ta
        1.62, # 74: W
        1.51, # 75: Re
        1.44, # 76: Os
        1.41, # 77: Ir
        1.36, # 78: Pt
        1.36, # 79: Au
        1.32, # 80: Hg
        1.45, # 81: Tl
        1.46, # 82: Pb
        1.48, # 83: Bi
        1.40, # 84: Po
        1.50, # 85: At
        1.50, # 86: Rn
        2.60, # 87: Fr
        2.21, # 88: Ra
        2.15, # 89: Ac
        2.06, # 90: Th
        2.00, # 91: Pa
        1.96, # 92: U
        1.90, # 93: Np
        1.87, # 94: Pu
        1.80, # 95: Am
        1.69, # 96: Cm
    ]

    @staticmethod
    def vdw(atom):
        if isinstance(atom, str):
            atom = PeriodicTable.get_atomic_number(atom)
        try:
            return PeriodicTable.__vdw[atom]
        except (IndexError, KeyError):
            logger.error("PeriodicTable.vdw(): no VDW radius for atom %s", atom)
            raise

    @staticmethod
    def covalent_radius(atom):
        if isinstance(atom, str):
            atom = PeriodicTable.get_atomic_number(atom)
        if atom == 0:
            logger.error("PeriodicTable.covalent_radius(): atomic number 0 has no covalent radius")
            raise IndexError("PeriodicTable.covalent_radius(): atomic number 0 has no covalent radius")
        try:
            return PeriodicTable.__covalent_radius[atom]
        except (IndexError, KeyError):
            logger.error("PeriodicTable.covalent_radius(): no covalent radius for atom %s", atom)
            raise

    @staticmethod
    def atomic_weight(atom):
        if isinstance(atom, str):
            atom = PeriodicTable.get_atomic_number(atom)
        try:
            return PeriodicTable.__atomic_weights[atom]
        except (IndexError, KeyError):
            logger.error("PeriodicTable.atomic_weight(): no atomic weight for atom %s", atom)
            raise
    
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
