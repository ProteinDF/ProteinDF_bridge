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
from collections import OrderedDict

from .error import BrInputError
from .vector import Vector
from .select import Select_Atom, Select_AtomGroup
from .select import Select
from .atom import Atom
from .position import Position
from .periodictable import PeriodicTable
from .str_processing import StrUtils

import logging
logger = logging.getLogger(__name__)


class AtomGroup(object):
    """
    >>> group1 = AtomGroup()
    >>> atom1 = Atom(symbol='C')
    >>> atom2 = Atom(symbol='H')
    >>> atom3 = Atom(symbol='N')
    >>> subgrp = AtomGroup()
    >>> subgrp.set_atom('C1', atom1)
    >>> subgrp.set_atom('H1', atom2)
    >>> subgrp.set_atom('N1', atom3)
    >>> group1.set_group('grp', subgrp)
    >>> group1['grp']['C1'].symbol
    'C'
    >>> group1.get_number_of_atoms()
    0
    >>> group1.get_number_of_groups()
    1
    >>> group1.get_number_of_all_atoms()
    3
    >>> group1.sum_of_atomic_number()
    14.0

    >>> group2 = AtomGroup()
    >>> atom_with_path1 = Atom(symbol='C')
    >>> atom_with_path2 = Atom(symbol='N')
    >>> group2.set_atom('/group1/subgroup1/C1', atom_with_path1)
    >>> group2.set_atom('/group1/subgroup2/N2', atom_with_path2)
    >>> print(group2)
    # group key= name=
      # group key=group1 name=   parent=
        # group key=subgroup1 name=     parent=
        C (    )    0.000    0.000    0.000,  0.00,    0.000    0.000    0.000 /group1/subgroup1/C1
        # group key=subgroup2 name=     parent=
        N (    )    0.000    0.000    0.000,  0.00,    0.000    0.000    0.000 /group1/subgroup2/N2
    <BLANKLINE>
    >>> group2['group1']['subgroup1'].has_atom('C1')
    True
    """

    def __init__(self, *args, **kwargs):
        self._initialize()

        if len(args) > 0:
            if len(args) == 1:
                rhs = args[0]
                if isinstance(rhs, AtomGroup):
                    for k, v in rhs.atoms():
                        self.set_atom(k, v)
                    for k, v in rhs.groups():
                        self.set_group(k, v)
                    self._bonds = copy.copy(rhs._bonds)
                    self.name = rhs.name
                    self._path = StrUtils.to_unicode(rhs._path)
                    self._update_path()
                    self.parent = rhs._parent
                    self._sort_atoms = rhs._sort_atoms
                    self._sort_groups = rhs._sort_groups
                elif isinstance(rhs, dict):
                    self.set_by_dict_data(rhs)
            else:
                raise BrInputError('atomgroup.__init__',
                                   'illegal the number of args')

        if 'name' in kwargs:
            self.name = kwargs.get('name')
        if 'parent' in kwargs:
            self._parent = kwargs.get('parent')

    def _initialize(self):
        self._atoms = OrderedDict()
        self._groups = OrderedDict()
        self._bonds = []
        self.name = ''
        self._path = '/'
        self._parent = None

        # 'nice'を指定すると数字順にアクセスできる
        self._sort_atoms = None
        self._sort_groups = None

    # property -----------------------------------------------------------------
    def _get_sort_atoms(self):
        return self._sort_atoms

    def _set_sort_atoms(self, v):
        self._sort_atoms = v
    sort_atoms = property(_get_sort_atoms, _set_sort_atoms)

    def _get_sort_groups(self):
        return self._sort_groups

    def _set_sort_groups(self, v):
        self._sort_groups = v
    sort_groups = property(_get_sort_groups, _set_sort_groups)

    # --------------------------------------------------------------------------
    def get_number_of_groups(self):
        return len(self._groups)

    def get_number_of_atoms(self):
        return len(self._atoms)

    def get_number_of_all_atoms(self):
        answer = 0
        for key, grp in self.groups():
            answer += grp.get_number_of_all_atoms()
        answer += len(self._atoms)
        return answer

    # parent -------------------------------------------------------------------
    def _get_parent(self):
        return self._parent

    def _set_parent(self, parent):
        self._parent = parent

    parent = property(_get_parent, _set_parent)

    # family -------------------------------------------------------------------
    def get_family(self, query_path):
        if query_path == self.path:
            return self

        common_path = StrUtils.get_common_str(self.path, query_path)
        if len(common_path) == 0:
            return None

        if len(common_path) < len(self.path):
            return self.parent.get_family(query_path)
        else:
            for key, grp in self.groups():
                answer = grp.get_family(query_path)
                if answer is not None:
                    return answer

        return None

    # move ---------------------------------------------------------------------
    def shift_by(self, direction):
        for key, grp in self.groups():
            grp.shift_by(direction)
        for key, atm in self.atoms():
            atm.shift_by(direction)

    def rotate(self, rotmat):
        for key, grp in self.groups():
            grp.rotate(rotmat)
        for key, atm in self.atoms():
            atm.rotate(rotmat)
        return self

    # --------------------------------------------------------------------------
    def sum_of_atomic_number(self):
        """
        原子数の総和を返す
        """
        answer = 0.0
        for key, grp in self.groups():
            answer += grp.sum_of_atomic_number()
        for key, atm in self.atoms():
            answer += atm.atomic_number
        return answer

    def get_atom_kinds(self):
        """
        原子種(シンボル)のリストを返す
        """
        answer = set()
        for key, group in self.groups():
            tmp = group.get_atom_kinds()
            answer.update(tmp)
        for key, atom in self.atoms():
            answer.add(atom.symbol)
        return answer

    def get_atom_kinds_count(self):
        """
        原子種(シンボル)とその数を格納した辞書を返す
        """
        kinds = {}

        for k, subgrp in self.groups():
            subgrp_kinds = subgrp.get_atom_kinds_count()
            for symbol, count in subgrp_kinds.items():
                kinds.setdefault(symbol, 0)
                kinds[symbol] += count

        for k, atom in self.atoms():
            symbol = atom.symbol
            kinds.setdefault(symbol, 0)
            kinds[symbol] += 1

        return kinds

    # --------------------------------------------------------------------------
    def get_num_of_groups(self):
        '''Return the number of groups
        '''
        return len(self._groups)

    def groups(self):
        """
        原子団のリストを返す
        """
        if self._sort_groups == 'nice':
            keys = list(self._groups.keys())
            keys = StrUtils.sort_nicely(keys)
            for k in keys:
                yield(k, self._groups[k])
        else:
            for k, v in self._groups.items():
                yield(k, v)

    def get_group(self, key_or_name):
        """
        入力されたkeyもしくは名前の原子団が含まれている場合、その原子を返す。
        無い場合はNoneを返す。
        """
        key_or_name = StrUtils.to_unicode(key_or_name)
        if key_or_name in self._groups:
            return self._groups.get(key_or_name, None)
        else:
            for k, grp in self.groups():
                if grp.name == key_or_name:
                    return grp
        return None

    def set_group(self, key, value):
        key = str(key)
        key = StrUtils.to_unicode(key)
        assert(isinstance(value, AtomGroup))
        if '_groups' not in self.__dict__:
            self._groups = {}
        self._groups[key] = AtomGroup(value, parent=self)
        self._update_path()

    def has_groupkey(self, key):
        """
        入力されたkeyのグループが含まれている場合、Trueを返す。
        無い場合はFalseを返す。
        """
        answer = False
        key = StrUtils.to_unicode(key)
        if key in self._groups:
            answer = True
        return answer

    def has_groupname(self, name):
        """
        入力されたnameのグループが含まれている場合、Trueを返す。
        無い場合はFalseを返す。
        """
        answer = False
        name = StrUtils.to_unicode(name)
        for k, grp in self.groups():
            if grp.name == name:
                answer = True
                break
        return answer

    def has_group(self, key_or_name):
        """
        入力されたkeyもしくは名前のグループが含まれている場合、Trueを返す。
        無い場合はFalseを返す。
        """
        return (self.has_groupkey(key_or_name)) or (self.has_groupname(key_or_name))

    def erase_group(self, key):
        """remove group

        absolete function.
        """
        logger.info("absolete function: erase_group")
        self.remove_group(key)

    def remove_group(self, key):
        """remove group
        """
        key = StrUtils.to_unicode(key)
        self._groups.pop(key, None)

    def get_group_list(self):
        return [k for k, v in self.groups()]

    # --------------------------------------------------------------------------
    def get_num_of_atoms(self):
        '''Return the number of atoms
        '''
        return len(self._atoms)

    def atoms(self):
        """
        原子のリストを返す
        """
        if self._sort_atoms == 'nice':
            keys = list(self._atoms.keys())
            keys = StrUtils.sort_nicely(keys)
            for k in keys:
                yield(k, self._atoms[k])
        else:
            for k, v in self._atoms.items():
                yield(k, v)

    def get_atom_keys(self):
        return [k for k, v in self.atoms()]

    def get_atom(self, key_or_name):
        """
        入力されたkeyもしくは名前の原子が含まれている場合、その原子を返す。
        無い場合はNoneを返す。
        """
        key_or_name = StrUtils.to_unicode(key_or_name)
        if key_or_name in self._atoms:
            return self._atoms.get(key_or_name, None)
        else:
            for k, atm in self.atoms():
                if atm.name == key_or_name:
                    return atm
        return None

    def _set_atom(self, key, value):
        key = str(key)
        key = StrUtils.to_unicode(key)
        assert(isinstance(value, Atom))
        self._atoms[key] = Atom(value,
                                parent=self,
                                path='{}{}'.format(self.path, key))

    def set_atom(self, key, value):
        assert(isinstance(value, Atom))
        key = str(key)
        keys = key.split('/', 1)
        while (len(keys) > 0) and (len(keys[0]) == 0):
            key = keys[1]
            keys = key.split('/', 1)

        if len(keys) == 1:
            self._set_atom(keys[0], value)
        else:
            grp_key = keys[0]
            rest = keys[1]
            if self.has_groupkey(grp_key) is False:
                self.set_group(grp_key, AtomGroup())
            self._groups[grp_key].set_atom(rest, value)

    def has_atomkey(self, key):
        """
        入力されたkeyの原子が含まれている場合、Trueを返す。
        無い場合はFalseを返す。
        """
        answer = False
        key = StrUtils.to_unicode(key)
        if key in self._atoms:
            answer = True
        return answer

    def has_atomname(self, name):
        """
        入力された名前の原子が含まれている場合、Trueを返す。
        無い場合はFalseを返す。
        """
        answer = False
        name = StrUtils.to_unicode(name)
        name = name.strip().lstrip()
        for k, atm in self.atoms():
            atm_name = atm.name.strip().lstrip()
            if atm_name == name:
                answer = True
                break
        return answer

    def has_atom(self, key_or_name):
        """
        入力されたkeyもしくは名前の原子が含まれている場合、Trueを返す。
        無い場合はFalseを返す。
        """
        return (self.has_atomkey(key_or_name)) or (self.has_atomname(key_or_name))

    def erase_atom(self, key):
        """remove atom

        absolete function.
        """
        logger.info("absolete function: erase_atom()")
        self.remove_atom(key)

    def remove_atom(self, key):
        """remove atom
        """
        key = StrUtils.to_unicode(key)
        self._atoms.pop(key, None)

    def pickup_atoms(self, key_or_name):
        '''
        key または nameが一致した原子の配列を返す
        '''
        answer = []
        for subgrp_key, subgrp in self.groups():
            atomlist = subgrp.pickup_atoms(key_or_name)
            if len(atomlist) > 0:
                answer.extend(atomlist)
        for atm_key, atm in self.atoms():
            if (atm_key == key_or_name) or (atm.name == key_or_name):
                answer.append(atm)
        return answer

    def get_atom_list(self):
        '''
        サブグループ内の原子をリスト型に格納して返す
        '''
        atom_list = []
        for subgrp_key, subgrp in self.groups():
            subgrp_list = subgrp.get_atom_list()
            atom_list.extend(subgrp_list)
        for atom_key, atom in self.atoms():
            atom_list.append(atom)

        return atom_list

    def get_path_list(self):
        '''
        グループ内の原子のパスのリストを返す
        '''
        answer = []
        for subgrp_key, subgrp in self.groups():
            answer.extend(subgrp.get_path_list())
        for atom_key, atom in self.atoms():
            answer.append(atom.path)
        return answer

    def get_formula(self):
        """分子式(組成式; composition formula)を返す
        """
        kinds = self.get_atom_kinds_count()

        formula = ""
        max_atom_id = PeriodicTable.get_num_of_atoms()
        for atom_id in range(1, max_atom_id):
            symbol = PeriodicTable.get_symbol(atom_id)
            if symbol in kinds:
                formula += "{}{}".format(symbol, kinds[symbol])
        if "X" in kinds:
            formula += "X{}".format(kinds["X"])

        return formula

    # name ---------------------------------------------------------------------
    def _get_name(self):
        return self._name

    def _set_name(self, name):
        self._name = StrUtils.to_unicode(name)

    name = property(_get_name, _set_name)

    # charge -------------------------------------------------------------------
    def _get_charge(self):
        charge = 0.0
        for key, ag in self.groups():
            charge += ag.charge
        for key, atom in self.atoms():
            charge += atom.charge
        return charge

    charge = property(_get_charge)
    # nuclei charge ------------------------------------------------------------

    def _get_real_nuclei_charge(self):
        charge = 0.0
        for key, ag in self.groups():
            charge += ag.real_nuclei_charge
        for key, atom in self.atoms():
            if atom.is_real:
                charge += atom.atomic_number
        return charge

    # nuclei charge without dummy atoms
    real_nuclei_charge = property(_get_real_nuclei_charge)
    # nuclei charge ------------------------------------------------------------

    def _get_nuclei_charge(self):
        charge = 0.0
        for key, ag in self.groups():
            charge += ag.nuclei_charge
        for key, atom in self.atoms():
            if atom.is_real:
                charge += atom.atomic_number
            else:
                charge += atom.charge
        return charge

    # nuclei charge including dummy atoms
    nuclei_charge = property(_get_nuclei_charge)
    # --------------------------------------------------------------------------

    def _get_path(self):
        return self._path

    def _set_path(self, value):
        value = StrUtils.to_unicode(value)
        if (len(value) == 0) or (value[-1] != '/'):
            value += '/'

        if (self._path != value):
            self._path = value
            self._update_path()

    path = property(_get_path, _set_path)
    # --------------------------------------------------------------------------

    def merge(self, rhs):
        """
        原子団を結合する
        """
        assert(isinstance(rhs, AtomGroup))
        for key, group in rhs.groups():
            self._merge_group(key, group)
        for key, atom in rhs.atoms():
            self.set_atom(key, Atom(atom))
    # --------------------------------------------------------------------------

    def assign_charges(self, charges):
        assert(isinstance(charges, Vector))
        index = AtomGroup._assign_charges(self, charges, 0)
        print(index, len(charges))
        assert(index == len(charges))

    @staticmethod
    def _assign_charges(atomgroup, charges, charge_index):
        for key, subgrp in atomgroup.groups():
            charge_index = AtomGroup._assign_charges(
                subgrp, charges, charge_index)
        for key, atom in atomgroup.atoms():
            atom.charge = charges.get(charge_index)
            charge_index += 1

        return charge_index

    # formula ------------------------------------------------------------------
    def formula(self):
        atom_list = self.get_atom_list()
        atom_counts = {}
        for atom in atom_list:
            atom_counts.setdefault(atom.symbol, 0)
            atom_counts[atom.symbol] += 1

        answer = ""
        for symbol, counts in atom_counts.items():
            answer += "{}{}".format(symbol, counts)

        return answer

    # weight -------------------------------------------------------------------
    def _get_weight(self):
        atom_list = self.get_atom_list()
        weight = 0.0
        for atom in atom_list:
            weight += atom.weight()

        return weight

    weight = property(_get_weight)

    # --------------------------------------------------------------------------

    def select(self, selector):
        """
        selectorにSelectorオブジェクトを渡すことで
        対応する原子団を返します
        """
        assert(isinstance(selector, Select))
        # self._update_path()

        answer = None
        if (selector.is_match(self)):
            answer = AtomGroup(self)
        else:
            answer = AtomGroup()
            answer.name = self.name
            for key, group in self.groups():
                tmp = group.select(selector)
                if ((tmp.get_number_of_groups() != 0) or (tmp.get_number_of_atoms() != 0)):
                    answer.set_group(key, tmp)
            for key, atom in self.atoms():
                if selector.is_match(atom):
                    answer.set_atom(key, atom)
            answer.path = self.path
            # print("path:{} ::{}".format(self.path, str(answer)))
        return answer

    # --------------------------------------------------------------------------
    def restructure(self, reference, range=1.0E-5):
        """referenceの構造を参照して、データ構造を再構築する。

        フラットな原子リストをPDBデータ構造にビルドアップするときに便利。
        """
        assert(isinstance(reference, AtomGroup))

        # matching
        target_selector = Select_AtomGroup(self, range)
        restructured = reference.select(target_selector)

        # copy attributes: charges
        AtomGroup._copy_attributes(restructured, self)

        # calc the rest
        rest_of_target = AtomGroup._get_rest_of_frame_molecule(
            self, restructured)
        AtomGroup._assign_rest_molecule(rest_of_target, restructured)

        return restructured

    @staticmethod
    def _get_rest_of_frame_molecule(frame_molecule, selected_molecule):
        # calc the rest
        selector = Select_AtomGroup(selected_molecule)
        selected = frame_molecule.select(selector)
        rest_molecule = frame_molecule ^ selected

        return rest_molecule

    @staticmethod
    def _copy_attributes(target, reference):
        for key, subgrp in target.groups():
            AtomGroup._copy_attributes(subgrp, reference)
        for key, atom in target.atoms():
            atom_selector = Select_Atom(atom)
            ref_atoms = reference.select(atom_selector)
            if ref_atoms.get_number_of_all_atoms() > 0:
                ref_atoms = ref_atoms.get_atom_list()
                ref_atom = ref_atoms[0]
                atom.charge = ref_atom.charge

    @staticmethod
    def _assign_rest_molecule(rest_molecule, output_atom_group,
                              model_id="model_1", chain_id="Z", res_name="UNK"):
        chain = AtomGroup()
        res = AtomGroup()
        res.name = res_name
        atom_id = 1
        for atom in rest_molecule.get_atom_list():
            res.set_atom(atom_id, atom)
            atom_id += 1
        chain.set_group(1, res)

        output_atom_group[model_id].set_group(chain_id, chain)
    # --------------------------------------------------------------------------

    def get_number_of_bonds(self):
        return len(self._bonds)

    def get_bond_list(self, bond_list=None):
        """
        タプル('atom1のpath', 'atom2のpath', 結合次数)のリストを返す
        """
        # self._update_path()

        if bond_list is None:
            bond_list = []

        for key, subgrp in self.groups():
            subgrp.get_bond_list(bond_list)

        for b in self._bonds:
            # print("get_bond_list> ", self.path, b[0], b[1])
            bond_info = [None] * 3
            bond_info[0] = '{}{}'.format(self.path, b[0])
            bond_info[1] = '{}{}'.format(self.path, b[1])
            bond_info[2] = b[2]
            bond_list.append(bond_info)

        return bond_list

    def add_bond(self, atom1, atom2, order=1):
        """
        結合情報を追加する
        order = 結合次数
        """
        assert(isinstance(atom1, Atom))
        assert(isinstance(atom2, Atom))
        assert(isinstance(order, int))
        bond_info = (atom1, atom2, order)
        self._add_bond_normalize(bond_info)

    def _add_bond_normalize(self, bond_info):
        """
        結合情報を(正規化しながら)追加する
        """
        assert(len(bond_info) == 3)
        (atom1, atom2, order) = bond_info
        assert(isinstance(atom1, Atom))
        assert(isinstance(atom2, Atom))
        assert(isinstance(order, int))

        common_path = self._get_common_path(atom1.path, atom2.path)
        logger.debug("_add_bond_norm > {} {} {} {}".format(self.path, common_path, atom1.path, atom2.path))
        family = self.get_family(common_path)
        if family is not None:
            family._add_bond(bond_info)
        else:
            self._add_bond(bond_info)

    def _add_bond(self, bond_info):
        (atom1, atom2, order) = bond_info
        atom1_path = atom1.path
        atom2_path = atom2.path
        common_path1 = self._get_common_path(self.path, atom1_path)
        common_path2 = self._get_common_path(self.path, atom2_path)
        if len(common_path1) > 0:
            atom1_path = atom1_path[len(common_path1):]
        if len(common_path2) > 0:
            atom2_path = atom2_path[len(common_path2):]
        logger.debug("_add_bond> {} {} {}".format(self.path, atom1_path, atom2_path))
        self._bonds.append((atom1_path, atom2_path, order))

    def _get_common_path(self, path1, path2):
        common_path = "/"
        common_path_raw = StrUtils.get_common_str(path1, path2)
        if common_path_raw[-1] == '/':
            common_path = common_path_raw
        else:
            last_slash_index = common_path_raw.rfind('/')
            if last_slash_index != -1:
                common_path = common_path_raw[0:last_slash_index + 1]
        return common_path

    # --------------------------------------------------------------------------

    def box(self):
        """
        """
        box_min = self.center()
        box_max = copy.deepcopy(box_min)
        for grpkey, grp in self.groups():
            (grpbox_min, grpbox_max) = grp.box()
            box_min.x = min(box_min.x, grpbox_min.x)
            box_min.y = min(box_min.y, grpbox_min.y)
            box_min.z = min(box_min.z, grpbox_min.z)
            box_max.x = max(box_max.x, grpbox_max.x)
            box_max.y = max(box_max.y, grpbox_max.y)
            box_max.z = max(box_max.z, grpbox_max.z)
        for atmkey, atm in self.atoms():
            box_min.x = min(box_min.x, atm.xyz.x)
            box_min.y = min(box_min.y, atm.xyz.y)
            box_min.z = min(box_min.z, atm.xyz.z)
            box_max.x = max(box_max.x, atm.xyz.x)
            box_max.y = max(box_max.y, atm.xyz.y)
            box_max.z = max(box_max.z, atm.xyz.z)
        return (box_min, box_max)

    def center(self):
        """
        return Position value of the center
        """
        center = Position(0.0, 0.0, 0.0)
        for grpkey, grp in self.groups():
            num_of_grp_atoms = grp.get_number_of_all_atoms()
            if num_of_grp_atoms > 0:
                center += grp.center() * float(num_of_grp_atoms)
        for atmkey, atm in self.atoms():
            center += atm.xyz

        num_of_atoms = self.get_number_of_all_atoms()
        if num_of_atoms > 0:
            center *= (1.0 / num_of_atoms)

        return center

    # file format --------------------------------------------------------------

    def get_xyz(self):
        """
        XYZフォーマット文字列を返す
        """
        output = '%d\n' % (self.get_number_of_all_atoms())
        output += '# \n'
        output += self._get_xyz_recursive()
        return output

    def _get_xyz_recursive(self):
        """
        get_xyz()メソッド内部で再帰的に呼ばれる関数
        """
        output = ""
        for key, grp in self.groups():
            output += grp.get_xyz()
        for key, atm in self.atoms():
            p = atm.xyz
            output += '%2s %8.3f %8.3f %8.3f\n' % (atm.symbol,
                                                   p.x, p.y, p.z)

        return output

    # private method -----------------------------------------------------------
    def _merge_group(self, key, group):
        key = StrUtils.to_unicode(key)
        assert(isinstance(group, AtomGroup))
        if (self.has_group(key)):
            self._groups[key].merge(group)
        else:
            self.set_group(key, group)

    def _update_path(self, force=False):
        for key, group in self.groups():
            group.path = "%s%s/" % (self._path, key)
            if force:
                group._update_path(force)
        for key, atom in self.atoms():
            atom.path = "%s%s" % (self._path, key)

    # --------------------------------------------------------------------------
    def __and__(self, other):
        assert(isinstance(other, AtomGroup))
        answer = AtomGroup()
        for key, group in self.groups():
            if other.has_group(key):
                answer.set_group(key, self.get_group(key) & other.get_group(key))
                if (answer.get_group(key).get_number_of_all_atoms() == 0):
                    answer.remove_group(key)

        for key, atom in self.atoms():
            if other.has_atom(key):
                answer.set_atom(key, atom)

        return answer

    def __iand__(self, rhs):
        """
        implement of '&=' operator
        """
        assert(isinstance(rhs, AtomGroup))
        # self.update_path(self.get_path())
        # rhs.update_path(rhs.get_path())

        for key, group in self.groups():
            if rhs.has_group(key):
                self._groups[key] &= rhs._groups[key]
                if ((self._groups[key].get_number_of_groups() == 0) and (self._groups[key].get_number_of_atoms() == 0)):
                    self.remove_group(key)
            else:
                self.remove_group(key)

        for key, atom in self.atoms():
            if rhs.has_atom(key) is False:
                self.remove_atom(key)

        return self

    def __or__(self, other):
        answer = AtomGroup(self)
        answer.merge(other)

        return answer

    def __ior__(self, rhs):
        """
        implement of '|=' operator
        """
        assert(isinstance(rhs, AtomGroup))
        # self.update_path(self.get_path())
        # rhs.update_path(rhs.get_path())

        self.merge(rhs)
        return self

    def __ixor__(self, rhs):
        """
        implement of '^=' operator
        """
        self = self ^ rhs

        return self

    def __xor__(self, other):
        assert(isinstance(other, AtomGroup))

        answer = AtomGroup()
        subgrp_list = self.get_group_list()
        subgrp_list.extend(other.get_group_list())
        subgrp_list = set(subgrp_list)
        for key in subgrp_list:
            subgrp = AtomGroup()
            if self.has_group(key):
                if other.has_group(key):
                    subgrp = self.get_group(key) ^ other.get_group(key)
                else:
                    subgrp = self.get_group(key)
            else:
                subgrp = other.get_group(key)

            if subgrp.get_number_of_all_atoms() > 0:
                answer.set_group(key, subgrp)

        atom_keys = self.get_atom_keys()
        atom_keys.extend(other.get_atom_keys())
        atom_keys = set(atom_keys)
        for key in atom_keys:
            atom = None
            if self.has_atom(key):
                if other.has_atom(key):
                    pass
                else:
                    atom = self.get_atom(key)
            else:
                atom = other.get_atom(key)

            if atom is not None:
                answer.set_atom(key, atom)

        return answer

    # --------------------------------------------------------------------------
    def __imul__(self, rhs):
        """
        implement of '*=' operator
        """
        for k, subgrp in self.groups():
            subgrp *= rhs
        for k, atom in self.atoms():
            atom *= rhs

        return self

    # --------------------------------------------------------------------------
    def set_by_dict_data(self, data):
        assert(isinstance(data, dict))
        data = StrUtils.to_unicode_dict(data)

        tmp_groups = {}
        tmp_atoms = {}
        for key, value in data.items():
            if isinstance(key, bytes):
                key = key.decode('utf-8')

            if key == 'groups':
                for grp_key, grp_data in value.items():
                    # atomgroup = AtomGroup(grp_data)
                    tmp_groups[grp_key] = AtomGroup(grp_data)
            elif key == 'atoms':
                for atm_key, atm_data in value.items():
                    atom = Atom(atm_data)
                    tmp_atoms[atm_key] = atom
            elif key == 'name':
                self.name = value
            elif key == 'bonds':
                self._bonds = value
            else:
                print('AtomGroup::set_by_dict_dat(): unknown key: {}={}'.format(
                    key, str(value)))

        # store groups and atoms in order
        grp_keys = tmp_groups.keys()
        grp_keys = StrUtils.sort_nicely(grp_keys)
        for grp_key in grp_keys:
            self.set_group(grp_key, tmp_groups[grp_key])
        atom_keys = tmp_atoms.keys()
        atom_keys = StrUtils.sort_nicely(atom_keys)
        for atom_key in atom_keys:
            self.set_atom(atom_key, tmp_atoms[atom_key])

        if "sort_atoms" in data:
            self.sort_atoms = data["sort_atoms"]
        if "sort_groups" in data:
            self.sort_groups = data["sort_groups"]

        # self._update_path()
        return self

    def get_raw_data(self):
        # self._update_path()
        data = {}
        if (len(self._groups) > 0):
            groups = {}
            for key, grp in self.groups():
                groups[key] = grp.get_raw_data()
            data['groups'] = groups
        if (len(self._atoms) > 0):
            atoms = {}
            for key, atm in self.atoms():
                atoms[key] = atm.get_raw_data()
            data['atoms'] = atoms

        data['name'] = self._name
        if len(self._bonds) > 0:
            data['bonds'] = self._bonds

        if self.sort_atoms is not None:
            data['sort_atoms'] = self.sort_atoms
        if self.sort_groups is not None:
            data['sort_groups'] = self.sort_groups

        return data

    def __str__(self):
        # self._update_path()
        return self._get_str()

    def _get_str(self, key='', indent_level=0):
        indent = '  ' * indent_level

        answer = '{indent}# group key={key} name={name}'.format(indent=indent,
                                                                key=key,
                                                                name=self.name)
        if self.parent is not None:
            answer += '{indent} parent={parent}\n'.format(indent=indent,
                                                          parent=self.parent.name)
        else:
            answer += '\n'

        for key, atomgroup in self.groups():
            answer += atomgroup._get_str(key, indent_level + 1)
        for key, atom in self.atoms():
            answer += "{indent}{atom} {atom_path}\n".format(indent=indent,
                                                            atom_path=atom.path,
                                                            atom=str(atom))
        for bond in self._bonds:
            answer += "{indent}bond {atom_path1} <-{order}-> {atom_path2}\n".format(indent=indent,
                                                                                    atom_path1=bond[0],
                                                                                    atom_path2=bond[1],
                                                                                    order=bond[2])
        # answer += '{indent}>\n'.format(indent=indent)

        return answer

    def __getitem__(self, key):
        """
        operator[] for getter

        keyが一致した原子団、原子を返す。
        もしkeyが一致しなければ、名前から検索する。
        """
        key = StrUtils.to_unicode(str(key))
        if (self.has_group(key)):
            return self._groups[key]
        elif key in self._atoms:
            return self._atoms[key]
        else:
            for k, grp in self.groups():
                if grp.name == key:
                    return grp
            for k, atm in self.atoms():
                if atm.name == key:
                    return atm
        raise KeyError(key)

    def __setitem__(self, key, value):
        """operator[] for setter"""
        key = StrUtils.to_unicode(key)
        if isinstance(value, AtomGroup):
            self.set_group(key, value)
        elif isinstance(value, Atom):
            self.set_atom(key, value)
        else:
            raise ValueError(value)

    # ------------------------------------------------------------------
    # serialize
    # ------------------------------------------------------------------

    def __getstate__(self):
        return self.get_raw_data()

    def __setstate__(self, state):
        self._initialize()
        self.set_by_dict_data(state)

    # ------------------------------------------------------------------
    # utilities
    # ------------------------------------------------------------------
    @staticmethod
    def divide_path(path):
        path = str(path)
        parts = path.split("/")

        while '' in parts:
            parts.remove('')

        # num_of_parts = len(parts)
        return parts


if __name__ == "__main__":
    import doctest
    doctest.testmod()
