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

from __future__ import absolute_import


#from .common import NullHandler
from .utils import Utils, str, bytes, basestring, unicode
from .error import BrError, BrInputError, BrValueError

from .vector import Vector
from .matrix import Matrix, SymmetricMatrix, identity_matrix
from .position import Position

from .periodictable import PeriodicTable
from .atom import Atom
from .atomgroup import AtomGroup
#from bond import Bond

from .select import Select, Select_Name, Select_Path, Select_Path_simple, Select_Path_wildcard, Select_PathRegex, Select_Atom, Select_Range

from .biopdb import Pdb
from .mmcif import SimpleMmcif
from .mol2 import SimpleMol2
from .xyz import Xyz

from .aminoacid import AminoAcid

from .amber_prmtop import AmberPrmtop

from .modeling import Modeling
from .ionpair import IonPair
from .neutralize import Neutralize

from .superposer import Superposer
from .superposer_quaternion import Superposer_quaternion

from .dbmanager import DbManager
from .mail import Mail

# __all__ = [
#     'Utils',
#     'BrError', 'BrInputError', 'BrValueError',
#     'Vector',
#     'Matrix', 'SymmetricMatrix',
#     'Position',
#     'PeriodicTable',
#     'Atom',
#     'AtomGroup',
#     'Select_Name', 'Select_Path', 'Select_PathRegex', 'Select_Atom', 'Select_Range',
#     'Pdb',
#     'Xyz',
#     'Modeling',
#     'Ionpair',
#     'Superposer',
#     'Superposer_quaternion',
#     'DbManager',
#     'Mail'
# ]
