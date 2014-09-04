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
        return answer


    #@classmethod
    #def byte2str(cls, s):
    #    """
    #    byteをstr(utf-8)に変換する
    #    """
    #    answer = str(s)
    #    if isinstance(s, bytes):
    #        answer = str(s.decode('utf-8'))
    #    return answer

    @classmethod
    def byte2str_dict_keys(cls, d):
        """
        byteをキーとして保存してある辞書に対して、
        str(utf-8)をキーとする辞書に変換する。
        """
        answer = {}
        if isinstance(d, dict):
            for k, v in d.items():
                if isinstance(k, bytes):
                    old_key = bytes(k)
                    new_key = Utils.byte2str(k)
                    answer[new_key] = v
        return answer

    @classmethod
    def byte2str(cls, obj):
        """
        byteをstr(utf-8)に変換する
        """
        if isinstance(obj, bytes):
            obj = str(obj.decode('utf-8'))
        elif isinstance(obj, list):
            for i, item in enumerate(obj):
                obj[i] = cls.byte2str(item)
        elif isinstance(obj, dict):
            for k, v in obj.items():
                new_key = cls.byte2str(k)
                new_val = cls.byte2str(v)
                del obj[k]
                obj[new_key] = new_val
        return obj

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
