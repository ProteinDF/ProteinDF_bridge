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
import copy
try:
    from collections import OrderedDict
except ImportError:
    # python 2.6 or earlier, use backport
    from ordereddict import OrderedDict

import pdfbridge

class AtomGroup(object):
    """
    >>> group1 = AtomGroup()
    >>> atom1 = pdfbridge.Atom(symbol='C')
    >>> atom2 = pdfbridge.Atom(symbol='H')
    >>> atom3 = pdfbridge.Atom(symbol='N')
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
    """

    def __init__(self, *args, **kwargs):
        self._atoms = OrderedDict()
        self._groups = OrderedDict()
        self._bonds = []
        self._name = ''
        self._path = '/'
        self._parent = None

        # 'nice'を指定すると数字順にアクセスできる
        self._sort_atoms = None
        self._sort_groups = None

        if len(args) > 0:
            if len(args) == 1:
                rhs = args[0]
                if (isinstance(rhs, AtomGroup) == True):
                    for k, v in rhs.atoms():
                        self.set_atom(k, v)
                    for k, v in rhs.groups():
                        self.set_group(k, v)
                    self._bonds = copy.copy(rhs._bonds)
                    self._name = pdfbridge.Utils.byte2str(rhs._name)
                    self._path = pdfbridge.Utils.byte2str(rhs._path)
                    self.parent = rhs._parent
                    self._sort_atoms = rhs._sort_atoms
                    self._sort_groups = rhs._sort_groups
                elif (isinstance(rhs, dict) == True):
                    self.set_by_dict_data(rhs)
            else:
                raise pdfbridge.InputError('atomgroup.__init__', 'illegal the number of args')

        if 'name' in kwargs:
            self._name = kwargs.get('name')
        if 'parent' in kwargs:
            self._parent = kwargs.get('parent')

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
        
        common_path = pdfbridge.Utils.get_common_str(self.path, query_path)
        if len(common_path) == 0:
            return None
            
        if len(common_path) < len(self.path):
            return self.parent.get_family(query_path)
        else:
            for key, grp in self.groups():
                answer = grp.get_family(query_path)
                if answer != None:
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

    # --------------------------------------------------------------------------
    def groups(self):
        """
        原子団のリストを返す
        """
        if self._sort_atoms == 'nice':
            keys = list(self._groups.keys())
            keys = pdfbridge.Utils.sort_nicely(keys)
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
        key_or_name = pdfbridge.Utils.byte2str(key_or_name)
        if key_or_name in self._groups:
            return self._groups.get(key_or_name, None)
        else:
            for k, grp in self.groups():
                if grp.name == key_or_name:
                    return grp
        return None

    def set_group(self, key, value):
        key = str(key)
        key = pdfbridge.Utils.byte2str(key)
        assert(isinstance(value, AtomGroup))
        if '_groups' not in self.__dict__:
            self._groups = {}
        self._groups[key] = AtomGroup(value, parent = self)
        self._update_path()

    def has_groupkey(self, key):
        """
        入力されたkeyのグループが含まれている場合、Trueを返す。
        無い場合はFalseを返す。
        """
        answer = False
        key = pdfbridge.Utils.byte2str(key)
        if key in self._groups:
            answer = True
        return answer

    def has_groupname(self, name):
        """
        入力されたnameのグループが含まれている場合、Trueを返す。
        無い場合はFalseを返す。
        """
        answer = False
        name = pdfbridge.Utils.byte2str(name)
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
        key = pdfbridge.Utils.byte2str(key)
        self._groups.pop(key, None)

    #def get_group_list(self):
    #    return self.data['groups'].keys()

    # --------------------------------------------------------------------------
    def atoms(self):
        """
        原子のリストを返す
        """
        if self._sort_atoms == 'nice':
            keys = list(self._atoms.keys())
            keys = pdfbridge.Utils.sort_nicely(keys)
            for k in keys:
                yield(k, self._atoms[k])
        else:
            for k, v in self._atoms.items():
                yield(k, v)

    def get_atom(self, key_or_name):
        """
        入力されたkeyもしくは名前の原子が含まれている場合、その原子を返す。
        無い場合はNoneを返す。
        """
        key_or_name = pdfbridge.Utils.byte2str(key_or_name)
        if key_or_name in self._atoms:
            return self._atoms.get(key_or_name, None)
        else:
            for k, atm in self.atoms():
                if atm.name == key_or_name:
                    return atm
        return None

    def set_atom(self, key, value):
        key = str(key)
        key = pdfbridge.Utils.byte2str(key)
        assert(isinstance(value, pdfbridge.Atom))
        self._atoms[key] = pdfbridge.Atom(value,
                                          parent=self,
                                          path='{}{}'.format(self.path, key))
        
    def has_atomkey(self, key):
        """
        入力されたkeyの原子が含まれている場合、Trueを返す。
        無い場合はFalseを返す。
        """
        answer = False
        key = pdfbridge.Utils.byte2str(key)
        if key in self._atoms:
            answer = True
        return answer

    def has_atomname(self, name):
        """
        入力された名前の原子が含まれている場合、Trueを返す。
        無い場合はFalseを返す。
        """
        answer = False
        name = pdfbridge.Utils.byte2str(name)
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
        key = pdfbridge.Utils.byte2str(key)
        self._atoms.pop(key, None)

    #def get_atom_list(self):
    #    return self.data['atoms'].keys()

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
    
    # name ---------------------------------------------------------------------
    def _get_name(self):
        return self._name

    def _set_name(self, name):
        name = str(name)
        self._name = pdfbridge.Utils.byte2str(name)

    name = property(_get_name, _set_name)
    # --------------------------------------------------------------------------
    def _get_charge(self):
        charge = 0.0
        for key, ag in self.groups():
            charge += ag.charge
        for key, atom in self.atoms():
            charge += atom.charge
        return charge

    charge = property(_get_charge)
    # --------------------------------------------------------------------------
    def _get_path(self):
        return self._path

    def _set_path(self, value):
        value = str(value)
        if (self._path != value):
            self._path = value
            if (self._path[-1] != '/'):
                self._path += '/'
            self._update_path()

    path = property(_get_path, _set_path)
    # --------------------------------------------------------------------------
    def merge(self, rhs):
        """
        原子団を結合する
        """
        assert(isinstance(rhs, AtomGroup) == True)
        for key, group in rhs.groups():
            self._merge_group(key, group)
        for key, atom in rhs.atoms():
            self.set_atom(key, pdfbridge.Atom(atom))

    # --------------------------------------------------------------------------
    def select(self, selecter):
        """
        selecterにSelecterオブジェクトを渡すことで
        対応する原子団を返します
        """
        assert(isinstance(selecter, pdfbridge.Select) == True)
        self._update_path()

        answer = None
        if (selecter.is_match(self) == True):
            answer = AtomGroup(self)
        else:
            answer = AtomGroup()
            answer.name = self.name
            for key, group in self.groups():
                tmp = group.select(selecter)
                if ((tmp.get_number_of_groups() != 0) or
                    (tmp.get_number_of_atoms() != 0)):
                    answer.set_group(key, tmp)
            for key, atom in self.atoms():
                if (selecter.is_match(atom) == True):
                    answer.set_atom(key, atom)
            answer.path = self.path
        return answer

    # --------------------------------------------------------------------------
    def get_bond_list(self, bond_list = None):
        """
        タプル('atom1のpath', 'atom2のpath', 結合次数)のリストを返す
        """
        self._update_path()
        
        if bond_list == None:
            bond_list = []

        for key, subgrp in self.groups():
            subgrp.get_bond_list(bond_list)

        for b in self._bonds:
            #print("get_bond_list> ", self.path, b[0], b[1])
            b[0] = '{}{}'.format(self.path, b[0])
            b[1] = '{}{}'.format(self.path, b[1])
            bond_list.append(b)

        return bond_list

    def add_bond(self, atom1, atom2, order =1):
        """
        結合情報を追加する
        order = 結合次数
        """
        bond_info = (atom1, atom2, order)
        self._add_bond_normalize(bond_info)

    def _add_bond_normalize(self, bond_info):
        """
        結合情報を(正規化しながら)追加する
        """
        assert(len(bond_info) == 3)
        (atom1, atom2, order) = bond_info
        assert(isinstance(atom1, pdfbridge.Atom))
        assert(isinstance(atom2, pdfbridge.Atom))
        assert(isinstance(order, int))

        common_path = pdfbridge.Utils.get_common_str(atom1.path, atom2.path)
        family = self.get_family(common_path)
        #print("_add_bond_norm> ", self.path, common_path, atom1.path, atom2.path)
        if family != None:
            family._add_bond(bond_info)
        else:
            self._add_bond(bond_info)

    def _add_bond(self, bond_info):
        (atom1, atom2, order) = bond_info
        atom1_path = atom1.path
        atom2_path = atom2.path
        common_path1 = pdfbridge.Utils.get_common_str(self.path, atom1_path)
        common_path2 = pdfbridge.Utils.get_common_str(self.path, atom2_path)
        if len(common_path1) > 0:
            atom1_path = atom1_path[len(common_path1):]
        if len(common_path2) > 0:
            atom2_path = atom2_path[len(common_path2):]
        #print("_add_bond> ", self.path, atom1_path, atom2_path)
        self._bonds.append((atom1_path, atom2_path, order))

    # --------------------------------------------------------------------------
    def box(self):
        """
        """
        box_min = box_max = self.center()
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
        center = pdfbridge.Position(0.0, 0.0, 0.0)
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
        output  = '%d\n' % (self.get_number_of_all_atoms())
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
        key = pdfbridge.Utils.byte2str(key)
        assert(isinstance(group, AtomGroup) == True)
        if (self.has_group(key) == True):
            self._groups[key].merge(group)
        else:
            self.set_group(key, group)

    def _update_path(self):
        for key, group in self.groups():
            group.path = "%s%s/" % (self._path, key)
        for key, atom in self.atoms():
            atom.path = "%s%s" % (self._path, key)

    # --------------------------------------------------------------------------
    def __iand__(self, rhs):
        """
        implement of '&=' operator
        """
        assert(isinstance(rhs, AtomGroup) == True)
        #self.update_path(self.get_path())
        #rhs.update_path(rhs.get_path())

        for key, group in self.groups():
            if rhs.has_group(key):
                self._groups[key] &= rhs._groups[key]
                if ((self._groups[key].get_number_of_groups() == 0) and
                    (self._groups[key].get_number_of_atoms() == 0)):
                    self.erase_group(key)
            else:
                self.erase_group(key)

        for key, atom in self.atoms():
            if rhs.has_atom(key) != True:
                self.erase_atom(key)

        return self

    def __ior__(self, rhs):
        """
        implement of '|=' operator
        """
        assert(isinstance(rhs, AtomGroup) == True)
        #self.update_path(self.get_path())
        #rhs.update_path(rhs.get_path())

        self.merge(rhs)
        return self

    def __ixor__(self, rhs):
        """
        implement of '^=' operator
        """
        assert(isinstance(rhs, AtomGroup) == True)
        #self.update_path(self.get_path())
        #rhs.update_path(rhs.get_path())

        for key, group in rhs.groups():
            if (self.has_group(key) == True):
                self._groups[key].__ixor__(group)
                if ((self._groups[key].get_number_of_groups() == 0) and
                    (self._groups[key].get_number_of_atoms() == 0)):
                    self.erase_group(key)
            else:
                self.set_group(key, group)

        for key, atom in rhs.atoms():
            if key in self._atoms:
                self.erase_atom(key)
            else:
                self.set_atom(key, atom)

        return self

    # --------------------------------------------------------------------------
    def set_by_dict_data(self, data):
        assert(isinstance(data, dict) == True)
        data = pdfbridge.Utils.byte2str(data)

        tmp_groups = {}
        tmp_atoms = {}
        for key, value in data.items():
            if isinstance(key, bytes):
                key = key.decode('utf-8')

            if key == 'groups':
                for grp_key, grp_data in value.items():
                    atomgroup = AtomGroup(grp_data)
                    tmp_groups[grp_key] = AtomGroup(grp_data)
            elif key == 'atoms':
                for atm_key, atm_data in value.items():
                    atom = pdfbridge.Atom(atm_data)
                    tmp_atoms[atm_key] = atom
            elif key == 'name':
                self.name = value
            elif key == 'bonds':
                self._bonds = value
            else:
                print('unknown key: {}'.format(key))

        # store groups and atoms in order
        grp_keys = tmp_groups.keys()
        grp_keys = pdfbridge.Utils.sort_nicely(grp_keys)
        for grp_key in grp_keys:
            self.set_group(grp_key, tmp_groups[grp_key])
        atom_keys = tmp_atoms.keys()
        atom_keys = pdfbridge.Utils.sort_nicely(atom_keys)
        for atom_key in atom_keys:
            self.set_atom(atom_key, tmp_atoms[atom_key])
                
        self._update_path()
        return self

    def get_raw_data(self):
        self._update_path()
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
        return data

    def __str__(self):
        self._update_path()

        return self._get_str()

    def _get_str(self, indent_level = 0):
        indent = '  ' * indent_level

        answer = '{indent}<grp name={name}'.format(indent=indent,
                                                   name=self.name)
        if self.parent is not None:
            answer += '{indent} parent={parent}'.format(indent=indent,
                                                        parent=self.parent.name)
        answer += '\n'
        for key, atomgroup in self.groups():
            answer += atomgroup._get_str(indent_level +1)
        for key, atom in self.atoms():
            answer += "{indent}{atom_path}:{atom}\n".format(indent=indent,
                                                            atom_path=atom.path,
                                                            atom=str(atom))
        for bond in self._bonds:
            answer += "{indent}bond {atom_path1} <-{order}-> {atom_path2}\n".format(indent=indent,
                                                                                    atom_path1=bond[0],
                                                                                    atom_path2=bond[1],
                                                                                    order=bond[2])
        answer += '{indent}>\n'.format(indent=indent)

        return answer
        
        
    def __getitem__(self, key):
        """
        operator[] for getter

        keyが一致した原子団、原子を返す。
        もしkeyが一致しなければ、名前から検索する。
        """
        key = pdfbridge.Utils.byte2str(str(key))
        if (self.has_group(key) == True):
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
        key = pdfbridge.Utils.byte2str(key)
        if (isinstance(value, AtomGroup) == True):
            self.set_group(key, value)
        elif (isinstance(value, pdfbridge.Atom) == True):
            self.set_atom(key, value)
        else:
            raise ValueError(value)


    # ------------------------------------------------------------------
    # serialize
    # ------------------------------------------------------------------
    def __getstate__(self):
        return self.get_dict_data()

    def __setstate__(self, state):
        self._atoms = {}
        self._groups = {}
        self._bonds = []
        self._name = ''
        self._path = '/'
        self._parent = None
        self.set_by_dict_data(state)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
