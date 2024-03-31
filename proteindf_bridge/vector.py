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

import os
import struct
import copy
import numpy
from types import *


class Vector(object):
    """
    >>> a = Vector(10)
    >>> len(a)
    10
    >>> a.resize(20)
    >>> len(a)
    20
    >>> a.resize(5)
    >>> len(a)
    5
    >>> a.set(1, 3.14)
    >>> (a.get(0) - 0.00 < 1.0E-10)
    True
    >>> (a.get(1) - 3.14 < 1.0E-10)
    True
    >>> (a[1] - 3.14 < 1.0E-10)
    True
    >>> a[2] = 1.73
    >>> (a[2] - 1.73 < 1.0E-10)
    True
    """

    __header_struct_little_endian = "<i"
    __header_struct_big_endian = ">i"
    __body_struct_little_endian = "<d"
    __body_struct_big_endian = ">d"

    def __init__(self, obj=[]):
        """
        初期化

        内部変数self._dataはnumpy.array(float)型
        """
        if isinstance(obj, int):
            size = obj
            self._data = numpy.array([0.0 for x in range(size)])
        elif isinstance(obj, (list, numpy.ndarray)):
            self._data = numpy.array(obj)
        elif isinstance(obj, dict):
            size = raw_data.get("size", 0)
            buf = raw_data.get("data", None)
            if buf != None:
                self._data = numpy.array([0.0 for x in range(size)])
        elif isinstance(obj, Vector):
            self._data = copy.copy(obj._data)
        else:
            print(type(obj))
            raise TypeError

    # --------------------------------------------------------------------------
    @property
    def max(self):
        return self._data.max()

    @property
    def min(self):
        return self._data.min()

    def abs(self):
        answer = Vector(self)
        answer._data = numpy.abs(answer._data)
        return answer

    @property
    def data(self):
        """
        return numpy array
        """
        return self._data

    # --------------------------------------------------------------------------
    def size(self):
        return self.__len__()

    def resize(self, new_size):
        new_data = numpy.array([0.0 for x in range(new_size)])
        for i in range(min(self.size(), new_size)):
            new_data[i] = self._data[i]
        self._data = new_data

    def get(self, index):
        return self._data[index]

    def set(self, index, value):
        self._data[index] = value

    def to_list(self):
        return self._data.tolist()

    def get_buffer(self):
        # return buffer(self._data.tostring())
        return self._data.tostring()

    def set_buffer(self, buf):
        self._data = numpy.fromstring(buf)

    def get_ndarray(self):
        return copy.deepcopy(self._data)

    # --------------------------------------------------------------------------
    def argsort(self):
        return Vector(self._data.argsort())

    def flip(self):
        return Vector(numpy.flip(self._data))

    # --------------------------------------------------------------------------
    def __get_header_struct(self, is_little_endian):
        if is_little_endian:
            return self.__header_struct_little_endian
        else:
            return self.__header_struct_big_endian

    def __get_body_struct(self, is_little_endian):
        if is_little_endian:
            return self.__body_struct_little_endian
        else:
            return self.__body_struct_big_endian

    def __str__(self):
        output = ""
        for order in range(0, len(self), 10):
            output += "\n"
            for j in range(order, min(order + 10, len(self))):
                output += "   %5d th" % (j + 1)
            output += "\n"
            for j in range(order, min(order + 10, len(self))):
                output += "-----------"
            output += "----\n\n"
            for j in range(order, min(order + 10, len(self))):
                output += " %10.6lf" % (self[j])
            output += "\n"
        return output

    def __add__(self, other):
        assert isinstance(other, Vector)
        assert len(self) == len(other)

        answer = Vector(self)
        answer += other
        return answer

    def __iadd__(self, other):
        assert isinstance(other, Vector)
        assert len(self) == len(other)

        self._data += other._data
        return self

    def __sub__(self, other):
        assert isinstance(other, Vector)

        answer = Vector(self)
        answer -= other
        return answer

    def __isub__(self, other):
        assert isinstance(other, Vector)

        self._data -= other._data
        return self

    def __mul__(self, other):
        answer = None
        if isinstance(other, float):
            answer = Vector(self)
            answer *= other
        elif isinstance(other, Vector):
            answer = float(self._data.dot(other._data))
        else:
            raise
        return answer

    def __rmul__(self, other):
        answer = None
        if isinstance(other, float):
            answer = Vector(self)
            answer *= other
        else:
            raise
        return answer

    def __imul__(self, other):
        if isinstance(other, float):
            self._data *= other
        else:
            raise
        return self

    def __neg__(self):
        answer = Vector(self)
        answer *= -1.0
        return answer

    def __len__(self):
        return len(self._data)

    def __getitem__(self, key):
        return self._data[key]

    def __setitem__(self, key, value):
        self._data[key] = value


if __name__ == "__main__":
    import doctest

    doctest.testmod()
