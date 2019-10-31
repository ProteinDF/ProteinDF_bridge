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

import sys
import os
import argparse
import re

from .error import BrInputError
from .position import Position
from .atom import Atom
from .atomgroup import AtomGroup

class Xyz(object):
    """
    """
    def __init__(self, *args, **kwargs):
        """
        create empty XYZ object
        """
        self._comment = ''
        self._atoms = []

        if len(args) > 0:
            if len(args) == 1:
                rhs = args[0]
                if isinstance(rhs, str):
                    self.load(file_path)
                elif isinstance(rhs, AtomGroup):
                    self.set_by_atomgroup(rhs)
                else:
                    raise BrInputError('Xyz.__init__', 'illegal object type')
            else:
                raise BrInputError('Xyz.__init__', 'illegal the number of args')

    def load(self, file_path):
        if (os.path.isfile(file_path) != True):
            return

        fin = open(file_path, 'r')
        num_of_atoms = int(fin.readline())
        self._comment = fin.readline().rstrip()
        for i in range(num_of_atoms):
            line = fin.readline()
            words = line.split()
            symbol = words[0]
            position = Position([float(words[1]),
                                 float(words[2]),
                                 float(words[3])])
            atom_data = {'symbol': symbol,
                         'position': position}
            self._atoms.append(atom_data)

    def save(self, file_path):
        f = open(file_path, 'w')
        f.write(self.get_text())
        f.close()

    def get_atom_group(self):
        """
        return AtomGroup object
        """
        root = AtomGroup()
        root.name = self._comment
        for i in range(len(self._atoms)):
            atom = Atom(symbol = self._atoms[i]['symbol'],
                        position = self._atoms[i]['position'])
            root.set_atom(str(i), atom)
        return root

    def set_by_atomgroup(self, atomgroup):
        assert(isinstance(atomgroup, AtomGroup))

        for model_key, model in atomgroup.groups():
            self.set_by_atomgroup(model)
        for serial, atom in atomgroup.atoms():
            symbol = atom.symbol
            position = atom.xyz
            atom_data = {'symbol': symbol,
                         'position': position}
            self._atoms.append(atom_data)

    def get_text(self):
        output = ''
        output += '%d\n' % (len(self._atoms))
        output += '%s\n' % (self._comment)
        for i in range(len(self._atoms)):
            symbol = self._atoms[i]['symbol']
            position = self._atoms[i]['position']
            output += '%s % 10.6f % 10.6f % 10.6f\n' % (symbol,
                                                        position.x,
                                                        position.y,
                                                        position.z)
        return output

    def __str__(self):
        return self.get_text()


def main():
    parser = argparse.ArgumentParser(description='read XYZ file')
    parser.add_argument('FILE',
                        nargs = 1,
                        help = 'XYZ file')
    args = parser.parse_args()

    file_path = args.FILE[0]

    xyz_obj = Xyz(file_path)
    print(xyz_obj)


if __name__ == '__main__':
    main()
