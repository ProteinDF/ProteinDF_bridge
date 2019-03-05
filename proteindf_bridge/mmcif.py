#!/usr/bin/env python

# Copyright (C) 2019 The ProteinDF development team.
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
import pprint
import msgpack

import logging
logger = logging.getLogger(__name__)

from .error import BrInputError, BrValueError
from .utils import Utils
from .position import Position
from .atom import Atom
from .atomgroup import AtomGroup

class SimpleMmcif(object):
    """
    """
    _re_data_block = re.compile("(data_\S+)")
    _re_keyvalue = re.compile("^(_\S+)\s+(\S+|\".*\"|\'.*\')\s*$")
    _re_loop_block = re.compile("^\s*loop_\s*$")
    _re_loop_header = re.compile("^(_\S+)\s*$")

    def __init__(self, file_path = None):
        self._data = {}

        if file_path:
            self.load(file_path)


    def load(self, file_path):
        with open(file_path) as file_obj:
            while True:
                line = self._get_line(file_obj)
                if len(line) == 0:
                    break

                match_obj = self._re_data_block.match(line)
                if match_obj:
                    name = match_obj.group(1)
                    # print("data: {}".format(name))
                    self._data[name] = self._load_data_block(file_obj)


    def _get_line(self, file_obj):
        line = ""
        start_semicolon = False # flag of start ";" block
        while True:
            failback_pos = file_obj.tell()
            current_line = file_obj.readline()
            #print("{}> '{}'".format(str(start_semicolon)[0], current_line))
            if (len(current_line) == 0) and (start_semicolon == False):
                break

            if current_line[0] == "#":
                continue

            current_line = current_line.strip()
            if len(line) == 0:
                # if "line" is empty, then this reading time is the first content.
                line = current_line
                continue

            if len(current_line) == 0:
                continue

            #print("{}>> '{}'".format(str(start_semicolon)[0], current_line))
            if current_line[0] == ";":
                if start_semicolon == False:
                    line += ' "'
                    start_semicolon = True
                    line += current_line[1:]
                    continue
                else:
                    line += '"'
                    start_semicolon = False
                    break
            elif current_line[0] == '"':
                line += " " + current_line
                continue
            else:
                #if start_semicolon:
                #    line += '"'
                #    start_semicolon = False
                if start_semicolon == True:
                    line += current_line
                    continue
                #print("getline cont.: [{}]".format(current_line))
                file_obj.seek(failback_pos)
                break

        #print("getline: [{}]".format(line))
        return line


    def _get_value(self, value):
        if (value[0] == '"' and value[-1] == '"'):
            value = value[1:-1]

        return value

    def _load_data_block(self, file_obj):
        kv = {}
        tables = []
        while True:
            failback_pos = file_obj.tell()
            line = self._get_line(file_obj)
            if len(line) == 0:
                # print("break data block: {}".format(line))
                break
            #print("check 109: {}".format(line))

            kv_match_obj = self._re_keyvalue.match(line)
            if kv_match_obj:
                # read key-value component
                key = kv_match_obj.group(1)
                value = kv_match_obj.group(2)
                kv[key] = self._get_value(value)
                continue
            else:
                loop_block_match_obj = self._re_loop_block.match(line)
                if loop_block_match_obj:
                    # read "_loop" component
                    tables.append(self._load_loop_block(file_obj))
                else:
                    # check end of "_data"
                    match_obj = self._re_data_block.match(line)
                    if match_obj:
                        logger.debug("end of data: {}".format(line))
                        # print("end of data: {}".format(line))
                        file_obj.seek(failback_pos)
                        break
                    else:
                        logger.warning("illegal end of data: {}".format(line))
                    continue
                # break

        return (kv, tables)


    def _load_loop_block(self, file_obj):
        data = []

        re_contents = None
        is_reading_loop_header = True
        last_contents_filepointer = None
        header = []
        # print(">>>> begin loop!")
        while True:
            line = self._get_line(file_obj)
            # print(line)

            if is_reading_loop_header:
                loop_header_match_obj = self._re_loop_header.match(line)
                if loop_header_match_obj:
                    name = loop_header_match_obj.group(1)
                    header.append(name)
                    continue
                else:
                    is_reading_header = False

                    # make contents re_contents
                    num_of_columns = len(header)

                    re_str = "^\s*" + "(\S+|\".*\")\s+" * num_of_columns
                    re_str = re_str[:-1]
                    re_str += "*$"
                    # print("re_str:> ", re_str)
                    re_contents = re.compile(re_str)

            data_match_obj = self._re_data_block.match(line)
            if data_match_obj:
                # logger.debug("end of loop by data_: {}".format(line))
                # print("<<<< end of loop by data_: {}".format(line))
                file_obj.seek(last_contents_filepointer)
                last_contents_filepointer = None
                break

            contents_match_obj = re_contents.match(line)
            if contents_match_obj:
                row = {}
                num_of_columns = len(header)
                for i in range(num_of_columns):
                    value = contents_match_obj.group(i +1)
                    row[header[i]] = self._get_value(value)
                data.append(row)
                last_contents_filepointer = file_obj.tell()
            else:
                # logger.debug("end of loop: {}".format(line))
                # print("<<<< end of loop: {}".format(line))
                file_obj.seek(last_contents_filepointer)
                last_contents_filepointer = None
                break

        return data


    def load_msgpack(self, file_path):
        packed = None
        with open(file_path, "rb") as file_obj:
            packed = file_obj.read()
        self._data = msgpack.unpackb(packed)
        if isinstance(self._data, list):
            self._data = Utils.to_unicode_list(self._data)
        elif isinstance(self._data, dict):
            self._data = Utils.to_unicode_dict(self._data)


    def save_msgpack(self, file_path):
        packed = msgpack.packb(self._data)
        with open(file_path, "wb") as file_obj:
            file_obj.write(packed)


    def get_molecule_names(self):
        return list(self._data.keys())


    def get_atomgroup(self, name):
        ag = AtomGroup()

        mol_data = self._data[name]
        try:
            self._get_atomgroup_list(mol_data, ag)
            self._get_atomgroup_bond_list(mol_data, ag)
        except BrInputError as e:
            raise BrInputError(name, "Invalid mmcif data: name={}".format(name))

        return ag


    def _get_atomgroup_list(self, list_item, ag):
        assert(isinstance(ag, AtomGroup))
        try:
            for item in list_item:
                if isinstance(item, list):
                    self._get_atomgroup_list(item, ag)
                elif isinstance(item, dict):
                    self._get_atomgroup_dict(item, ag)
        except BrInputError as e:
            raise e


    def _get_atomgroup_dict(self, dict_item, output_atomgroup):
        """output_atomgroupに出力する
        """
        assert(isinstance(output_atomgroup, AtomGroup))
        assert(isinstance(dict_item, dict))

        re_numbers = re.compile("^[\+\-]?\d*(\.\d+)?$")
        # check input coordinates(x, y, z)
        def get_coordinates(xyz, dict_item):
            assert((xyz == "x") or (xyz == "y") or (xyz == "z"))
            answer = None

            ideal_key = "_chem_comp_atom.pdbx_model_Cartn_{}_ideal".format(xyz)
            model_key = "_chem_comp_atom.model_Cartn_{}".format(xyz)
            xyz_ideal = dict_item.get(ideal_key, "N/A")
            xyz_model = dict_item.get(model_key, "N/A")

            if re_numbers.match(xyz_ideal):
                answer = float(xyz_ideal)
            elif re_numbers.match(xyz_model):
                answer = float(xyz_model)
            else:
                logger.warning("Invalid coordinate({}) [{}, {}]".format(xyz, xyz_ideal, xyz_model))
                # raise BrValueError([xyz_ideal, xyz_model], "Invalid coordinate({})".format(xyz))

            return answer

        # output key-value
        #print(">" * 10)
        #for key, value in dict_item.items():
        #    print("{}> {}".format(key, value))
        #print("<" * 10)

        if "_chem_comp.id" in dict_item:
            output_atomgroup.name = dict_item["_chem_comp.id"]

        if "_chem_comp_atom.atom_id" in dict_item:
            atom = Atom()
            id = dict_item["_chem_comp_atom.atom_id"]
            atom.name = id

            symbol = dict_item["_chem_comp_atom.type_symbol"]
            if symbol == "D":
                symbol = "H"
            atom.symbol = symbol

            x = get_coordinates("x", dict_item)
            y = get_coordinates("y", dict_item)
            z = get_coordinates("z", dict_item)
            if (x != None) and (y != None) and (z != None):
                atom.position = (x, y, z)
            else:
                update = False

            charge = dict_item["_chem_comp_atom.charge"]
            if charge != "?":
                atom.charge = charge
            output_atomgroup.set_atom(id, atom)


        # if "_chem_comp_bond.comp_id" in dict_item:
        #     atom1_name = dict_item["_chem_comp_bond.atom_id_1"]
        #     atom2_name = dict_item["_chem_comp_bond.atom_id_2"]
        #     bond_order_str = dict_item["_chem_comp_bond.value_order"]
        #     bond_order = 0
        #     if bond_order_str == "SING":
        #         bond_order = 1
        #     elif bond_order_str == "DOUB":
        #         bond_order = 2
        #     elif bond_order_str == "TRIP":
        #         bond_order = 3
        #     else:
        #         print("illegal input: {}".format(bond_order_str))
        #     atom1 = output_atomgroup.get_atom(atom1_name)
        #     atom2 = output_atomgroup.get_atom(atom2_name)
        #     output_atomgroup.add_bond(atom1, atom2, bond_order)


    def _get_atomgroup_bond_list(self, list_item, ag):
        """ bond 専用
        """
        assert(isinstance(ag, AtomGroup))
        try:
            for item in list_item:
                if isinstance(item, list):
                    self._get_atomgroup_bond_list(item, ag)
                elif isinstance(item, dict):
                    self._get_atomgroup_bond_dict(item, ag)
        except BrInputError as e:
            raise e


    def _get_atomgroup_bond_dict(self, dict_item, output_atomgroup):
        """ bond 専用
        """
        assert(isinstance(output_atomgroup, AtomGroup))
        assert(isinstance(dict_item, dict))

        if "_chem_comp_bond.comp_id" in dict_item:
            atom1_name = dict_item["_chem_comp_bond.atom_id_1"]
            atom2_name = dict_item["_chem_comp_bond.atom_id_2"]
            bond_order_str = dict_item["_chem_comp_bond.value_order"]
            bond_order = 0
            if bond_order_str == "SING":
                bond_order = 1
            elif bond_order_str == "DOUB":
                bond_order = 2
            elif bond_order_str == "TRIP":
                bond_order = 3
            else:
                logger.warning("illegal input: {}".format(bond_order_str))
            atom1 = output_atomgroup.get_atom(atom1_name)
            atom2 = output_atomgroup.get_atom(atom2_name)
            output_atomgroup.add_bond(atom1, atom2, bond_order)


    def __str__(self):
        answer = ""
        for name in self._data.keys():
            answer += "data: {}\n".format(name)
            (kv, tables) = self._data[name]
            for key, value in kv.items():
                answer += "  {key}: {value}\n".format(key=key, value=value)
            for table in tables:
                answer += "  ----\n"
                for row in range(len(table)):
                    for key, value in table[row].items():
                        answer += "  [{id}]{key}: {value}\n".format(id=row, key=key, value=value)

        return answer
