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
import sys
import re
import copy
import pickle

class Utils(object):
    @classmethod
    def sort_nicely(cls, l):
        """
        Sort the given list in the way that humans expect.
        ref: http://www.codinghorror.com/blog/2007/12/sorting-for-humans-natural-sort-order.html
        """
        l = [ str(x) for x in l ]
        convert = lambda text: int(text) if text.isdigit() else text
        alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
        l.sort(key=alphanum_key)
        return l

    @classmethod
    def add_spaces(cls, s, num_add):
        """
        行頭にnum_add分のスペースを追加する
        """
        spc = ' ' * num_add
        return spc + spc.join(s.splitlines(True))

    @classmethod
    def num_spaces(cls, s):
        """
        行頭のスペース数を返す
        """
        return [len(line) - len(line.lstrip()) for line in s.splitlines()]

    @classmethod
    def del_spaces(cls, s, num_del):
        """
        """
        if num_del < min(cls.num_spaces(s)):
            raise ValueError("removing more spaces than there are!")
        return '\n'.join([ line[num_del:] for line in s.splitlines()])

    @classmethod
    def unindent_block(cls, s):
        return cls.del_spaces(s, min(cls.num_spaces(s)))
        
    @classmethod
    def get_common_str(cls, str1, str2):
        """
        (先頭から)共通文字列を返す

        >>> a = 'abcdef'
        >>> b = 'abcdefg'
        >>> get_common_str(a, b)
        'abcdef'
        """
        answer = ""
        len1 = len(str1)
        len2 = len(str2)
        l = min(len1, len2)
        for i in range(l):
            c1 = str1[i]
            c2 = str2[i]
            if c1 == c2:
                answer += c1
            else:
                break
        return answer

    @classmethod
    def to_unicode_dict(cls, d):
        """
        byteをキーとして保存してある辞書に対して、
        str(utf-8)をキーとする辞書に変換する。
        """
        assert isinstance(d, dict)
        answer = {}
        if isinstance(d, dict):
            for k, v in d.items():
                new_key = cls.to_unicode(k)
                answer[new_key] = v
        return answer

    @classmethod
    def to_unicode(cls, unicode_or_str):
        """
        byteをstr(utf-8)に変換する
        """
        try:
            assert isinstance(unicode_or_str, (unicode, str, bytes))
        except:
            print(type(unicode_or_str))
            print(unicode_or_str)
            raise

        value = unicode_or_str
        if sys.version_info[0] >= 3:
            # Python3
            if isinstance(unicode_or_str, bytes):
                value = unicode_or_str.decode('utf-8')
        else:
            # Python2
            if isinstance(unicode_or_str, str):
                try:
                    value = unicode_or_str.decode('utf-8')
                except:
                    print(type(unicode_or_str))
                    print(unicode_or_str)
                    raise

        return value

    @classmethod
    def to_bytes(cls, unicode_or_str):
        assert isinstance(unicode_or_str, (str, bytes))

        value = unicode_or_str
        if sys.version_info[0] >= 3:
            # Python3
            if isinstance(unicode_or_str, str):
                value = unicode_or_str.encode('utf-8')
        else:
            # Python2
            if isinstance(unicode_or_str, unicode):
                value = unicode_or_str.encode('utf-8')

        return value

    @classmethod
    def str_to_bool(cls, input_str):
        if isinstance(input_str, bool):
            return input_str

        answer = False
        try:
            tmp = int(input_str)
            if tmp != 0:
                answer = True
        except:
            if len(input_str) > 0:
                tmp = input_str.upper()
                if (tmp[0] == 'Y' or
                    tmp[0] == 'T'):
                    answer = True
        return answer
    
    @classmethod
    def check_pickled(cls, data, level=0):
        if isinstance(data, dict):
            for k, v in data.items():
                cls.check_pickled(v)
        else:
            #print('>' * level, data)
            try:
                pickle.dumps(data)
            except:
                print(type(data), data)
                raise

if __name__ == "__main__":
    import doctest
    doctest.testmod()
