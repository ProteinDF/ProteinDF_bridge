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

from .atomgroup import AtomGroup
from .atom import Atom
from .error import BrInputError
from .functions import load_msgpack, save_msgpack
from . import rs

import re
import logging
logger = logging.getLogger(__name__)


class SimpleMmcif(object):
    """
    """

    def __init__(self, file_path=None):
        self._data = {}

        if file_path:
            self.load(file_path)

    def load(self, file_path):
        rs_cif = rs.SimpleMmcif(str(file_path))
        for block_name, (kv, tables) in rs_cif.data.items():
            self._data[block_name] = (kv, tables)

    def load_msgpack(self, file_path):
        self._data = load_msgpack(file_path)

    def save_msgpack(self, file_path):
        save_msgpack(self._data, file_path)

    def get_molecule_names(self):
        return list(self._data.keys())

    def get_atomgroup(self, name):
        ag = AtomGroup()

        mol_data = self._data[name]
        try:
            self._get_atomgroup_list(mol_data, ag)
            self._get_atomgroup_bond_list(mol_data, ag)
        except BrInputError as e:
            raise BrInputError(
                name, "Invalid mmcif data: name={}".format(name))

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
        """Write into output_atomgroup.
        """
        assert(isinstance(output_atomgroup, AtomGroup))
        assert(isinstance(dict_item, dict))

        re_numbers = re.compile(r"^[\+\-]?\d*(\.\d+)?$")
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
                logger.warning("Invalid coordinate({}) [{}, {}]".format(
                    xyz, xyz_ideal, xyz_model))
                # raise BrValueError([xyz_ideal, xyz_model], "Invalid coordinate({})".format(xyz))

            return answer

        # output key-value
        #print(">" * 10)
        # for key, value in dict_item.items():
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
        """ For bonds only.
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
        """ For bonds only.
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
                        answer += "  [{id}]{key}: {value}\n".format(
                            id=row, key=key, value=value)

        return answer
