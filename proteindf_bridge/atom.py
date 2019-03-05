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

import copy
import math
import logging
logger = logging.getLogger(__name__)

from .error import BrInputError
from .utils import Utils
from .periodictable import PeriodicTable
from .position import Position

class Atom(object):
    """
    >>> a = Atom()
    >>> a.symbol = 'Fe'
    >>> a.atomic_number
    26
    >>> a.symbol
    'Fe'
    >>> a.charge = -0.2
    >>> math.fabs(a.charge - (-0.2)) < 1.0E-10
    True
    >>> b = Atom(a)
    >>> b.symbol = 'Na'
    >>> a.symbol
    'Fe'
    >>> b.symbol
    'Na'
    """
    def __init__(self, *args, **kwargs):
        self._atomic_number = PeriodicTable.get_atomic_number('X')
        self._xyz = Position()
        self._force = Position()
        self.name = ''
        self._label = ''
        self._charge = 0.0
        self._path = ''
        self._parent = None

        if len(args) > 0:
            if len(args) == 1:
                rhs = args[0]
                if (isinstance(rhs, Atom) == True):
                    self._atomic_number = rhs._atomic_number
                    self._xyz = Position(rhs._xyz)
                    self._force = Position(rhs._force)
                    self.name = rhs._name
                    self._label = Utils.to_unicode(rhs.label)
                    self._charge = float(rhs._charge)
                    self._path = Utils.to_unicode(rhs._path)
                elif (isinstance(rhs, dict) == True):
                    self.set_by_raw_data(rhs)
            else:
                raise BrInputError('atom.__init__', 'illegal the number of args')

        if 'symbol' in kwargs:
            self._atomic_number = PeriodicTable.get_atomic_number(kwargs.get('symbol'))
        self._xyz = Position(kwargs.get('position', self._xyz))
        self._xyz = Position(kwargs.get('xyz', self._xyz)) # alias
        self._force = kwargs.get('force', self._force)
        if 'name' in kwargs:
            self.name = kwargs.get('name')
        if 'label' in kwargs:
            self._label = kwargs.get('label')
        if 'charge' in kwargs:
            self._charge = kwargs.get('charge')
        if 'path' in kwargs:
            self._path = kwargs.get('path')
        if 'parent' in kwargs:
            self._parent = kwargs.get('parent')
            #from .atomgroup import AtomGroup
            #assert(isinstance(self._parent, AtomGroup))

    # move ---------------------------------------------------------------------
    def move_to(self, position):
        self.xyz.move_to(position)
        return self

    def shift_by(self, direction):
        direction = Position(direction)
        self.xyz += direction
        return self

    def rotate(self, rotmat):
        self.xyz.rotate(rotmat)

    def __imul__(self, rhs):
        """
        implementation of '*=' operator
        """
        self.xyz *= rhs
        return self

    # --------------------------------------------------------------------------
    def _get_xyz(self):
        return self._xyz

    def _set_xyz(self, p):
        self._xyz = Position(p)

    xyz = property(_get_xyz, _set_xyz)
    position = property(_get_xyz, _set_xyz)

    # --------------------------------------------------------------------------
    def _get_force(self):
        return self._force

    def _set_force(self, f):
        self._force = Position(f)

    force = property(_get_force, _set_force)

    # --------------------------------------------------------------------------
    def _get_atomic_number(self):
        return self._atomic_number

    def _set_atomic_number(self, an):
        an = int(an)
        self._atomic_number = an

    atomic_number = property(_get_atomic_number, _set_atomic_number)

    # --------------------------------------------------------------------------
    def __get_symbol(self):
        answer = PeriodicTable.get_symbol(self.atomic_number)
        return answer

    def __set_symbol(self, symbol):
        self._atomic_number = PeriodicTable.get_atomic_number(symbol)

    symbol = property(__get_symbol, __set_symbol)

    # --------------------------------------------------------------------------
    def __is_real(self):
        atomnum = self.atomic_number
        return atomnum > 0

    is_real = property(__is_real)
    # --------------------------------------------------------------------------
    def _get_name(self):
        return self._name

    def _set_name(self, name):
        name = str(name)
        self._name = Utils.to_unicode(name)

    name = property(_get_name, _set_name)

    # --------------------------------------------------------------------------
    def _get_label(self):
        return self._label

    def _set_label(self, label):
        self._label = Utils.to_unicode(label)

    label = property(_get_label, _set_label)

    # --------------------------------------------------------------------------
    def _get_charge(self):
        return self._charge

    def _set_charge(self, charge):
        self._charge = float(charge)

    charge = property(_get_charge, _set_charge)

    # --------------------------------------------------------------------------
    def __get_vdw(self):
        return PeriodicTable.vdw(self._atomic_number)

    vdw = property(__get_vdw)

    # --------------------------------------------------------------------------
    def _get_path(self):
        return self._path

    def _set_path(self, path):
        self._path = Utils.to_unicode(path)

    path = property(_get_path, _set_path)

    # ==================================================================
    # operator
    # ==================================================================
    def __eq__(self, rhs):
        answer = False
        if ((isinstance(rhs, Atom) == True) and
            (self.atomic_number == rhs.atomic_number) and
            (math.fabs(self.charge - rhs.charge) < 1.0E-10) and
            (self.xyz == rhs.xyz)):
            answer = True
        return answer

    def __ne__(self, rhs):
        return not(self.__eq__(rhs))

    # ==================================================================
    # raw data
    # ==================================================================
    def set_by_raw_data(self, data):
        for key, value in data.items():
            if isinstance(key, bytes):
                key = key.decode('utf-8')

            if key == 'Z':
                self.atomic_number = value
            elif key == 'name':
                self.name = value
            elif key == 'Q':
                self.charge = value
            elif key == 'xyz':
                self.xyz = Position(value)
            elif key == 'force':
                self.force = Position(value)
            else:
                logger.debug("bridge::Atom > unknown key: {}".format(key))
        return self

    def get_raw_data(self):
        data = {}
        data['Z'] = self._atomic_number
        data['name'] = self.name
        data['Q'] = self._charge
        data['xyz'] = self._xyz.get_raw_data()
        data['force'] = self._force.get_raw_data()

        return data

    # ==================================================================
    # debug
    # ==================================================================
    def __str__(self):
        answer = "%2s(name=\"%s\")(% 8.3f, % 8.3f, % 8.3f) Z=% .2f <(% 8.3f, % 8.3f, % 8.3f)>" % (
            self.symbol,
            self.name,
            self.xyz.x,
            self.xyz.y,
            self.xyz.z,
            self.charge,
            self.force.x,
            self.force.y,
            self.force.z)
        return answer

    # ------------------------------------------------------------------
    # serialize
    # ------------------------------------------------------------------
    def __getstate__(self):
        return self.get_raw_data()

    def __setstate__(self, state):
        self.set_by_raw_data(state)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
