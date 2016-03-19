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
import optparse
import math
import copy
import numpy
import logging

import pdfbridge

class Position(object):
    """
    >>> p = Position([5.0, 3.0, -1.2])
    >>> p.x
    5.0
    >>> p.y
    3.0
    >>> p.z
    -1.2

    >>> p.x = 0.0
    >>> p.y = 1.0
    >>> p.z = 2.0
    >>> p == Position([0, 1, 2])
    True

    >>> abs(abs(p) - 2.23606) < 1.0E-5
    True

    >>> p.norm()
    >>> p == Position([0, 1/math.sqrt(5), 2/math.sqrt(5)])
    True

    >>> p.move_to([3, 4, 5])
    >>> p == Position([3.0, 4.0, 5.0])
    True

    >>> tmp = p * 2.0
    >>> tmp == Position([6.0, 8.0, 10.0])
    True

    >>> tmp = (-1.0 * p)
    >>> tmp == Position([-3.0, -4.0, -5.0])
    True

    >>> a = Position([1, 2, 3])
    >>> b = Position([2, 3, 4])
    >>> a * b
    20.0

    >>> p = a + b
    >>> p == Position([3, 5, 7])
    True

    >>> a += b
    >>> a == Position([3, 5, 7])
    True

    >>> p = a - b
    >>> p == Position([1, 2, 3])
    True

    >>> a -= b
    >>> a == Position([1, 2, 3])
    True

    >>> a.dot(b)
    20.0

    >>> n = a.cross(b)
    >>> n == Position([-1, 2, -1])
    True

    """
    def __init__(self, *args, **kwds):
        self._logger = logging.getLogger(__name__)

        self.epsilon = 1.0E-5
        self._position = [0.0, 0.0, 0.0]

        len_args = len(args)
        if len_args > 0:
            if len_args == 3:
                self._position[0] = float(args[0])
                self._position[1] = float(args[1])
                self._position[2] = float(args[2])
            elif len_args == 1:
                if (isinstance(args[0], Position) == True):
                    self._position[0] = args[0]._position[0]
                    self._position[1] = args[0]._position[1]
                    self._position[2] = args[0]._position[2]
                elif ((isinstance(args[0], (list, tuple, numpy.ndarray)) == True) and (len(args[0]) == 3)):
                    self._position[0] = float(args[0][0])
                    self._position[1] = float(args[0][1])
                    self._position[2] = float(args[0][2])
            else:
                raise pdfbridge.InputError("position::__init__", "illegal input")

    # --------------------------------------------------------------------------
    @property
    def xyz(self):
        return self._position

    def __get_x(self):
        return self._position[0]

    def __set_x(self, new_x):
        self._position[0] = float(new_x)

    x = property(__get_x, __set_x)

    def __get_y(self):
        return self._position[1]

    def __set_y(self, new_y):
        self._position[1] = float(new_y)

    y = property(__get_y, __set_y)

    def __get_z(self):
        return self._position[2]

    def __set_z(self, new_z):
        self._position[2] = float(new_z)

    z = property(__get_z, __set_z)

    # --------------------------------------------------------------------------
    def get_raw_data(self):
        return self._position

    def move_to(self, position):
        tmp = Position(position)
        self._position = tmp._position

    def square_distance_from(self, other = None):
        if other == None:
            other = Position()
        other = Position(other)

        d2 = 0.0
        for i in range(3):
            tmp = self._position[i] - other._position[i]
            d2 += tmp * tmp
        return d2

    def distance_from(self, other = None):
        d2 = self.square_distance_from(other)
        return math.sqrt(d2)

    def norm(self):
        n = self.__abs__()
        self._position = [x / n for x in self._position]

    def rotate(self, mat):
        assert(mat.rows == 3)
        assert(mat.cols == 3)
        v1 = pdfbridge.Vector(self._position)
        v2 = mat * v1
        self._position = v2.to_list()
        return self

    def dot(self, rhs):
        """
        内積を求める
        """
        return numpy.dot(self._position, rhs._position)

    def cross(self, rhs):
        """
        外積を求める
        """
        n = numpy.cross(self._position, rhs._position)
        answer = Position(n)
        return answer

    def __str__(self):
        answer = "(% 10.6f, % 10.6f, % 10.6f)" % (self.x,
                                                  self.y,
                                                  self.z)
        return answer

    def __eq__(self, rhs):
        answer = False
        if (isinstance(rhs, Position) == True):
            if (self.distance_from(rhs) < self.epsilon):
                answer = True
        return answer

    def __ne__(self, rhs):
        return not(self.__eq__(rhs))

    def __neg__(self):
        return Position([-x for x in self._position])

    def __abs__(self):
        return math.sqrt(sum([x * x for x in self._position]))

    def __add__(self, rhs):
        if (isinstance(rhs, Position) == True):
            return Position([ x + y for x, y in zip(self._position, rhs._position)])
        else:
            raise pdfbridge.InputError("position.__add__", "illegal input: Position is required.")

    def __sub__(self, rhs):
        return self.__add__(-rhs)

    def __sub__(self, rhs):
        return self.__add__(-rhs)

    def __mul__(self, rhs2):
        rhs1 = Position(self)
        if (isinstance(rhs2, Position) == True):
            return sum([x * y for x, y in zip(rhs1._position, rhs2._position)])
        elif (isinstance(rhs2, float) == True):
            rhs1._position = [x * rhs2 for x in rhs1._position]
            return rhs1
        else:
            raise pdfbridge.InputError("position.__add__", "illegal input: Position is required.")

    __rmul__ = __mul__

    def __imul__(self, rhs):
        v = float(rhs)
        self._position[0] *= v
        self._position[1] *= v
        self._position[2] *= v
        return self

    def __div__(self, rhs):
        v = float(rhs)
        return Position([x / v for x in self._position])

    def __idiv__(self, rhs):
        return self.__imul__(1.0 / float(rhs))

    # ------------------------------------------------------------------
    # serialize
    # ------------------------------------------------------------------
    def __getstate__(self):
        return self.get_raw_data()

    def __setstate__(self, state):
        assert(isinstance(state, (set, list)))
        assert(len(state) == 3)
        self._position[0] = float(state[0])
        self._position[1] = float(state[1])
        self._position[2] = float(state[2])

        
if __name__ == "__main__":
    import doctest
    doctest.testmod()
