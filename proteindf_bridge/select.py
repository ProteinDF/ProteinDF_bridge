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

import logging
logger = logging.getLogger(__name__)

from .utils import Utils
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


class Select_Symbol(Select):
    """
    原子記号で選択する
    """
    def __init__(self, atom_symbol):
        self._atom_symbol = atom_symbol.upper()

    def is_match(self, obj):
        answer = False
        if isinstance(obj, Atom):
            symbol = obj.symbol.upper()
            if symbol == self._atom_symbol:
                answer = True
        return answer

class Select_Name(Select):
    def __init__(self, query):
        self.query = Utils.to_unicode(query)

    def is_match(self, obj):
        answer = False
        name = obj.name.strip().rstrip()
        if (name == self.query):
            answer = True
        return answer

class Select_Path(Select):
    """
    """
    def __init__(self, query, use_wildcard = True):
        self._query = Utils.to_unicode(query)
        self._is_used_wildcard = use_wildcard

        if use_wildcard:
            logger.warning("Select_Path() is obsolete. please use Select_Path_wildcard().")
            self._regex_selecter = self._prepare(query)

    def _prepare(self, query):
        query = re.sub("(?<!\\\)\*", '.*', query)
        query = re.sub("(?<!\\\)\?", '?', query)
        query = '^' + query + '$'

        return Select_PathRegex(query)

    def is_match(self, obj):
        if self._is_used_wildcard:
            return self._regex_selecter.is_match(obj)
        else:
            return self._is_match_nowildcard(obj)

    def _is_match_nowildcard(self, obj):
        answer = False
        path = obj.path
        if (self._query == path):
            answer = True
        return answer


class Select_Path_simple(Select):
    """
    """
    def __init__(self, query):
        self._query = Utils.to_unicode(query)

    def is_match(self, obj):
        answer = False
        path = obj.path
        if (self._query == path):
            answer = True
        return answer


class Select_Path_wildcard(Select):
    """selector using path with wildcard
    """
    def __init__(self, query):
        self._query = Utils.to_unicode(query)
        self._regex_selecter = self._prepare(query)

    def _prepare(self, query):
        query = re.sub("(?<!\\\)\*", '.*', query)
        query = re.sub("(?<!\\\)\?", '?', query)
        query = '^' + query + '$'

        return Select_PathRegex(query)

    def is_match(self, obj):
        return self._regex_selecter.is_match(obj)


class Select_PathRegex(Select):
    """
    pathに対する正規表現で選択する
    """
    def __init__(self, query):
        self._query = Utils.to_unicode(query)
        self._regex = re.compile(query)

    def is_match(self, obj):
        answer = False
        path = obj.path
        if (self._regex.search(path) != None):
            #print("path=[%s] regex=[%s]" % (path, self._query))
            answer = True
        return answer

class Select_Range(Select):
    '''
    半径で選択する
    '''
    def __init__(self, pos, d):
        self._pos = Position(pos)
        d = float(d)
        self._d = d
        self._d2 = d * d

    def is_match(self, obj):
        answer = False
        if isinstance(obj, Atom):
            # d = self._pos.distance_from(obj.xyz)
            d2 = self._pos.square_distance_from(obj.xyz)
            if d2 < self._d2:
                answer = True
        return answer

class Select_Atom(Select):
    """
    """
    def __init__(self, atom):
        from .atom import Atom
        assert(isinstance(atom, Atom))
        self._atom = Atom(atom)

    def is_match(self, obj):
        answer = False
        if ((isinstance(obj, Atom)) and
            (self._atom.atomic_number == obj.atomic_number) and
            (self._atom.xyz == obj.xyz)):
                answer = True
        return answer


class Select_AtomGroup(Select):
    """reference atomgroupと同じ原子が存在しているものを返す
    """
    def __init__(self, ref_atomgroup):
        from .atomgroup import AtomGroup
        assert(isinstance(ref_atomgroup, AtomGroup))
        self._ref_atoms = ref_atomgroup.get_atom_list()

    def is_match(self, obj):
        answer = False
        if isinstance(obj, Atom):
            for ref_atom in self._ref_atoms:
                if ref_atom == obj:
                    answer = True
                    break

        return answer
