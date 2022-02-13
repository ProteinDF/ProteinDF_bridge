#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
logger = logging.getLogger(__name__)


class Path(object):
    @classmethod
    def get_chain_id(cls, path):
        '''
        >>> Path.get_chain_id("/B/17/561_CA")
        'B'
        '''
        items = cls.split_path(path)
        answer = None
        if len(items) > 1:
            answer = items[0]
        return answer

    @classmethod
    def get_res_id(cls, path):
        '''
        >>> Path.get_res_id("/B/17/561_CA")
        '17'
        '''
        items = cls.split_path(path)
        answer = None
        if len(items) > 2:
            answer = items[1]
        return answer

    @classmethod
    def split_path(cls, path):
        '''
        >>> Path.split_path("/model_1/B/17/561_CA")
        ['model_1', 'B', '17', '561_CA']
        '''
        items = path.split("/")
        if (len(items) > 0) and (items[0] == ''):
            items.pop(0)
        return items
