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
        'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Uut', 'Uuq', 'Uup', 'Uuh', 'Uus', 'Uuo'
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
    def get_symbol(atomic_number):
        try:
            answer = PeriodicTable.__table[atomic_number]
            return answer
        except:
            print("ERROR @PeriodicTable::get_symbol(): not found input:%s." % (atomicNumber))
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
            
if __name__ == "__main__":
    import doctest
    doctest.testmod()
