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

import logging
logger = logging.getLogger(__name__)


class AminoAcid(object):
    _AA_list = [
        'ALA',
        'ASX',  # ASP or ASN
        'ASN',
        'ASP',
        'CYS',
        'CYX',
        'GLU',
        'PHE',
        'GLY',
        'HIS',
        'HIE',
        'HIP',
        'ILE',
        'LYS',
        'LEU',
        'MET',
        'PRO',
        'GLN',
        'ARG',
        'SER',
        'THR',
        'SEC',
        'VAL',
        'TRP',
        'XAA',  # unspecified amino acid
        'TYR',
        'GLX'  # GLN, GLU or GLA, GLP
    ]

    def __init__(self, *args, **kwargs):
        pass

    @classmethod
    def is_aminoacid(self, atomgroup):
        name = atomgroup.name
        answer = False
        if name in self._AA_list:
            answer = True

        return answer
