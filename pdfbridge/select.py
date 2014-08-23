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

import re

from .position import Position
from .atom import Atom

class Select(object):
    """
    Selecterインターフェースクラス
    """
    def is_match(self, obj):
        """
        条件に適合した場合、Trueを返す

        サブクラスはこのメソッドを実装すること
        """
        return False


class Select_Name(Select):
    def __init__(self, query):
        assert(isinstance(query, str) == True)
        self.query = query

    def is_match(self, obj):
        answer = False
        name = obj.name.strip().rstrip()
        if (name == self.query):
            answer = True
        return answer

class Select_Path(Select):
    """
    """
    def __init__(self, query):
        assert(isinstance(query, str))
        self._query = query

    def is_match(self, obj):
        answer = False
        path = obj.get_path()
        if (self._query == path):
            answer = True
        return answer


class Select_PathRegex(Select):
    """
    pathに対する正規表現で選択する
    """
    def __init__(self, query):
        assert(isinstance(query, str))
        self._query = query
        self._regex = re.compile(query)

    def is_match(self, obj):
        answer = False
        path = obj.get_path()
        if (self._regex.search(path) != None):
            #print("path=[%s] regex=[%s]" % (path, self._query))
            answer = True
        return answer


class Select_Atom(Select):
    """
    原子記号で選択する
    """
    def __init__(self, atom_symbol):
        assert(isinstance(atom_symbol, str) == True)
        self._atom_symbol = atom_symbol.upper()

    def is_match(self, obj):
        answer = False
        if isinstance(obj, Atom):
            symbol = obj.symbol.upper()
            if symbol == self._atom_symbol:
                answer = True
        return answer

class Select_Range(Select):
    '''
    半径で選択する
    '''
    def __init__(self, pos, d):
        self._pos = Position(pos)
        self._d = float(d)

    def is_match(self, obj):
        answer = False
        if isinstance(obj, Atom):
            d = self._pos.distance_from(obj.xyz)
            if d < self._d:
                answer = True
        return answer
