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

from .atomgroup import AtomGroup
from .atom import Atom
from .position import Position
import sys
import os
import optparse
import re
import copy
import logging

logger = logging.getLogger(__name__)


class Pdb(object):
    """ """

    def __init__(self, file_path=None, mode=None):
        """
        create empty PDB object

        mode: None or amber
        """
        self._data = {}
        self._ssbonds = []
        if file_path:
            self.load(file_path)

        self._mode = mode
        if isinstance(self._mode, str):
            self._mode = self._mode.upper()

        # modpdb table
        self._modpdb_amber_atm_tbl = [
            {
                "name": "NA",
                "symbol": "Na",
                "rename": "Na+",
            },
            {
                "name": "CL",
                "symbol": "Cl",
                "rename": "Cl-",
            },
        ]
        self._modpdb_formal_atm_tbl = [
            {
                "name": "NA",
                "symbol": "Na",
                "rename": "NA",
            },
            {
                "name": "CL",
                "symbol": "Cl",
                "rename": "CL",
            },
        ]

        self._modpdb_amber_res_tbl = {"NA": "Na+", "CL": "Cl-"}
        self._modpdb_formal_res_tbl = {"NA": "NA ", "CL": "CL "}

        # leapではNH2はHとみなされず、NH2を重ねて付加してしまう
        self._modpdb_amber_resatom_table = {
            "NME": {
                "HN2": "H",
                "H1": "HH31",
                "H2": "HH32",
                "H3": "HH33",
            },
        }
        # reduceではNH2が無いとNH2を付加してしまう
        self._modpdb_formal_resatom_table = {
            "NME": {
                "H": "HN2",
                "HH31": "H1",
                "HH32": "H2",
                "HH33": "H3",
            },
        }

    def __get_debug(self):
        if not "_debug" in self.__dict__:
            self._debug = False
        return self._debug

    def __set_debug(self, yn):
        self._debug = yn

    debug = property(__get_debug, __set_debug)

    # def get_number_of_items(self):
    #    return len(self._data)

    # def get_name(self, index):
    #    answer = None
    #    if ((0 <= index) and (index < self.get_number_of_items())):
    #        answer = self._data[index].get('name', None)
    #    return answer

    # def set_name(self, index, name):
    #    if ((0 <= index) and (index < self.get_number_of_items())):
    #        self._data[index]['name'] = name

    # def get_element(self, index):
    #    answer = None
    #    if ((0 <= index) and (index < self.get_number_of_items())):
    #        answer = self._data[index].get('element', None)
    #    return answer

    # def set_element(self, index, element):
    #    if ((0 <= index) and (index < self.get_number_of_items())):
    #        self._data[index]['element'] = element

    # def get_posision(self, index):
    #    answer = None
    #    if ((0 <= index) and (index < self.get_number_of_items())):
    #        answer = self._data[index].get('coord', None)
    #    return answer

    # def set_position(self, index, position):
    #    if ((0 <= index) and (index < self.get_number_of_items())):
    #        self._data[index]['coord'] = position

    # def get_charge(self, index):
    #    answer = 0
    #    if ((0 <= index) and (index < self.get_number_of_items())):
    #        answer = self._data[index].get('charge', 0)
    #    return answer

    # def set_charge(self, index, charge):
    #    if ((0 <= index) and (index < self.get_number_of_items())):
    #        self._data[index]['charge'] = charge

    # def get_occupancy(self, index):
    #    answer = None
    #    if ((0 <= index) and (index < self.get_number_of_items())):
    #        answer = self._data[index].get('occupancy', None)
    #    return answer

    # def get_temp_factor(self, index):
    #    answer = None
    #    if ((0 <= index) and (index < self.get_number_of_items())):
    #        answer = self._data[index].get('temp_factor', None)
    #    return answer

    # def set_temp_factor(self, serial, temp_factor):
    #    serial = int(serial)
    #    self._data[serial]['temp_factor'] = temp_factor

    def renumber(self):
        for model_serial, model in self._data.items():
            for index in range(len(model)):
                model[index]["serial"] = index + 1

    def load(self, file_path):
        if os.path.isfile(file_path) != True:
            logger.critical("file not found: {}".format(file_path))
            raise

        model_serial = 1
        chain_serial = 0
        self._data.setdefault(model_serial, [])

        with open(file_path, "r") as fin:
            while True:
                line = fin.readline()
                if len(line) == 0:
                    break
                line = line.rstrip("\n")

                # if (len(line) != 80):
                #    continue
                if self.debug == True:
                    print(line)

                record_name = line[0:6]
                if record_name == "SSBOND":
                    serNum = line[7:10]
                    chainID1 = str(line[15])
                    seqNum1 = int(line[17:21])
                    icode1 = str(line[21])
                    chainID2 = str(line[29])
                    seqNum2 = int(line[31:35])
                    icode2 = str(line[35])

                    ssbond = ({"chain_id": chainID1, "seq_num": seqNum1}, {"chain_id": chainID2, "seq_num": seqNum2})
                    self._ssbonds.append(ssbond)

                elif (record_name == "ATOM  ") or (record_name == "HETATM"):

                    if len(line) < 80:
                        line = line + (" " * (80 - len(line)))

                    serial = int(line[6:11])
                    name4 = line[12:16]
                    name = name4.strip()
                    alt_loc = line[16]
                    res_name = line[17:20]
                    chain_id = line[21]
                    if (chain_id == " ") and (res_name != "WAT"):
                        chain_id = chr(ord("A") + (chain_serial % 26))
                    res_seq = line[22:26]
                    i_code = line[26]
                    coord_x = line[30:38]
                    coord_y = line[38:46]
                    coord_z = line[46:54]
                    occupancy = line[54:60].strip()
                    temp_factor = line[60:66].strip()
                    element = line[76:78].strip()
                    charge = line[78:80].strip()

                    item = {}
                    item["serial"] = serial
                    item["record_name"] = record_name
                    item["name"] = name
                    item["alt_loc"] = alt_loc

                    # res_name
                    if res_name in ["HID", "HIE", "HIP"]:
                        # rename AMBER residue name dialect
                        res_name = "HIS"
                    item["res_name"] = res_name

                    item["chain_id"] = chain_id
                    item["res_seq"] = int(res_seq)
                    item["i_code"] = i_code
                    item["coord"] = [float(coord_x), float(coord_y), float(coord_z)]

                    if len(occupancy) != 0:
                        item["occupancy"] = float(occupancy)
                    else:
                        item["occupancy"] = 1.0
                    if len(temp_factor) != 0:
                        item["temp_factor"] = float(temp_factor)
                    else:
                        item["temp_factor"] = 0.0

                    if len(element) != 0:
                        # TODO: 原子変換テーブル作成
                        element = element.strip()
                        if element == "D":
                            element = "H"
                        item["element"] = element
                    else:
                        # see https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
                        # TODO: テーブルを持って変換するように変更
                        name4 = name4.upper()
                        name2 = name4[0:2]
                        name2s = name2.strip()

                        if (len(name4.strip()) == 4) and (name2[0] == "H"):
                            element = "H"
                        elif len(name2s) == 2:
                            element = name2[0] + name2[1].lower()
                        else:
                            element = name2[1]

                        # if name2 == "CL":
                        #     element = "Cl"
                        # elif name2 == "NA":
                        #     element = "Na"
                        # elif name2 == "MG":
                        #     element = "Mg"
                        # elif name2 == "FE":
                        #     element = "Fe"
                        # else:
                        #     element = name[0]
                        item["element"] = element

                    if len(charge) != 0:
                        charge_last_char = charge[-1]
                        if (charge_last_char == "+") or (charge_last_char == "-"):
                            charge = charge_last_char + charge[0:-1]
                        item["charge"] = charge
                    else:
                        item["charge"] = "  "

                    self._data[model_serial].append(item)
                    continue
                elif record_name == "MODEL ":
                    serial = int(line[10:14])
                    model_serial = serial
                    self._data.setdefault(model_serial, [])
                    chain_serial = 0
                    continue
                elif record_name == "TER   ":
                    if len(line) < 27:
                        line = line + (" " * (27 - len(line)))
                    serial = int(line[6:11]) if line[6:11].isdigit() else 0
                    resname = line[17:20]
                    chain_id = line[21]
                    if (chain_id == " ") and (res_name != "WAT"):
                        chain_id = chr(ord("A") + (chain_serial % 26))
                        chain_serial += 1
                    res_seq = line[22:26]
                    i_code = line[26]
                    item = {}
                    item["serial"] = serial
                    item["record_name"] = record_name
                    item["res_name"] = res_name
                    item["chain_id"] = chain_id
                    item["res_seq"] = res_seq
                    item["i_code"] = i_code

                    self._data[model_serial].append(item)
                    continue

    def get_atomgroup(self, select_model=None, select_altloc="A"):
        """
        return AtomGroup object
        """
        root = AtomGroup()

        for model_serial, model_items in self._data.items():
            if (select_model == None) or (int(select_model) == int(model_serial)):

                model_name = "model_%d" % (model_serial)
                model = AtomGroup()
                model.name = model_name

                for index in range(len(model_items)):
                    item = model_items[index]
                    record_name = item["record_name"]
                    serial = item["serial"]

                    if (record_name == "ATOM  ") or (record_name == "HETATM"):
                        name = item["name"]
                        alt_loc = item["alt_loc"]
                        res_name = item["res_name"]
                        chain_id = item["chain_id"]
                        res_seq = int(item["res_seq"])
                        i_code = item["i_code"]
                        coord = item["coord"]
                        occupancy = item.get("occupancy", 1.0)
                        temp_factor = item.get("temp_factor", 0.0)
                        element = item.get("element", "X")
                        charge = item.get("charge", 0.0)
                        if charge == "  ":
                            charge = 0.0

                        if chain_id == " ":
                            chain_id = "_"

                        if model.has_group(chain_id) == False:
                            chain = AtomGroup()
                            chain.name = chain_id
                            model.set_group(chain_id, chain)

                        res_key = "%d" % (res_seq)
                        if model[chain_id].has_group(res_key) == False:
                            residue = AtomGroup()
                            residue.name = res_name
                            model[chain_id].set_group(res_key, residue)

                        # create atom object -------------------------------
                        atom = Atom()
                        atom.symbol = element
                        atom.xyz = Position(coord)
                        atom.name = name
                        atom.charge = charge
                        atom_key = "%d_%s" % (serial, name)

                        # set the atom object ------------------------------
                        if (alt_loc == " ") or (alt_loc == select_altloc):
                            model[chain_id][res_key].set_atom(atom_key, atom)
                        else:
                            logger.debug(
                                'skip alt_loc="{alt_loc}" atom: {atom_str}'.format(alt_loc=alt_loc, atom_str=str(atom))
                            )

                for ssbond in self._ssbonds:
                    chain_id1 = ssbond[0]["chain_id"]
                    seq_num1 = ssbond[0]["seq_num"]
                    chain_id2 = ssbond[1]["chain_id"]
                    seq_num2 = ssbond[1]["seq_num"]

                    path1 = "/{chain_id}/{res_key}/SG".format(chain_id=chain_id1, res_key=seq_num1)
                    path2 = "/{chain_id}/{res_key}/SG".format(chain_id=chain_id2, res_key=seq_num2)
                    SG1 = model[chain_id1][seq_num1]["SG"]
                    SG2 = model[chain_id2][seq_num2]["SG"]
                    model.add_bond(SG1, SG2, 1)

                root.set_group(model_name, model)

        return root

    def set_by_atomgroup(self, atomgroup, set_b_factor=None):
        assert isinstance(atomgroup, AtomGroup)
        atomgroup = self.get_modpdb_atomgroup(atomgroup)

        re_model_serial = re.compile("^model_(\d+)")
        re_res_seq = re.compile("^(\d+)")
        re_atom_serial = re.compile("^(\d+)")
        self._data = {}
        item = {}
        model_serial = 1
        for model_key, model in atomgroup.groups():
            match_obj = re_model_serial.match(model_key)
            if match_obj != None:
                model_serial = int(match_obj.group(1))
            self._data.setdefault(model_serial, [])

            serial = 1
            for chain_id, chain in model.groups():
                if chain_id != "_":
                    item["chain_id"] = chain_id
                else:
                    item["chain_id"] = " "

                for res_key, residue in chain.groups():
                    res_seq = 0
                    res_seq_match_obj = re_res_seq.match(res_key)
                    if res_seq_match_obj != None:
                        res_seq = int(res_seq_match_obj.group(1))
                    item["res_seq"] = res_seq

                    # resname
                    item["res_name"] = residue.name

                    has_OXT = False
                    for key, atom in residue.atoms():
                        item["record_name"] = "ATOM  "

                        # serial number
                        # serial_match_obj = re_atom_serial.match(key)
                        # if (serial_match_obj != None):
                        #    serial = int(match_obj.group(1))
                        item["serial"] = serial
                        serial += 1

                        item["alt_loc"] = " "
                        item["i_code"] = " "

                        # name
                        name = atom.name
                        if name.strip() == "OXT":
                            has_OXT = True
                        item["name"] = name
                        item["coord"] = [atom.xyz.x, atom.xyz.y, atom.xyz.z]
                        item["element"] = atom.symbol
                        item["charge"] = atom.charge

                        if set_b_factor == "charge":
                            item["temp_factor"] = atom.charge

                        self._data[model_serial].append(copy.copy(item))

                    if has_OXT == True:
                        item["record_name"] = "TER   "
                        item["serial"] = serial
                        serial += 1
                        self._data[model_serial].append(copy.copy(item))

                # TER
                if len(self._data[model_serial]) > 0 and (self._data[model_serial][-1]["record_name"] != "TER   "):
                    item["record_name"] = "TER   "
                    item["serial"] = serial
                    serial += 1
                    self._data[model_serial].append(copy.copy(item))

        self._sort_by_serial()

    def _sort_by_serial(self):
        for model_serial, model in self._data.items():
            model.sort(key=lambda x: int(x["serial"]))

    def __str__(self):
        occupancy = 1.0
        output = ""
        for model_serial, model in self._data.items():
            output += "MODEL     %4d\n" % (model_serial)
            for index in range(len(model)):
                item = model[index]
                record_name = item["record_name"]
                serial = int(item["serial"])
                if (record_name == "ATOM  ") or (record_name == "HETATM"):
                    name = item["name"]
                    element = item["element"]

                    # see. https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
                    if len(name) < 4:
                        if len(element) == 1:
                            name = " {:<3s}".format(name)
                        else:
                            name = "{:<4s}".format(name)

                    alt_loc = item["alt_loc"]
                    res_name = item["res_name"]
                    chain_id = item["chain_id"]
                    res_seq = int(item["res_seq"])
                    i_code = item["i_code"]
                    coord = item["coord"]
                    occupancy = item.setdefault("occupancy", 1.0)
                    temp_factor = item.setdefault("temp_factor", 1.0)
                    element = item.setdefault("element", "  ")
                    charge = int(item.setdefault("charge", 0))
                    if charge == 0:
                        charge = "  "

                    # line = "ATOM  %5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n" % (
                    #     serial,
                    #     name,
                    #     alt_loc,
                    #     res_name,
                    #     chain_id,
                    #     res_seq,
                    #     i_code,
                    #     coord[0],
                    #     coord[1],
                    #     coord[2],
                    #     occupancy,
                    #     temp_factor,
                    #     element.upper(),
                    #     charge,
                    # )
                    field = {
                        "serial": serial,
                        "name": name,
                        "alt_loc": alt_loc,
                        "res_name": res_name,
                        "chain_id": chain_id,
                        "res_seq": res_seq,
                        "i_code": i_code,
                        "x": coord[0],
                        "y": coord[1],
                        "z": coord[2],
                        "occupancy": occupancy,
                        "temp_factor": temp_factor,
                        "element": element.upper(),
                        "charge": charge,
                    }
                    line = "ATOM  {serial:>5d} {name:>4s}{alt_loc:1s}{res_name:>3s} {chain_id:1s}{res_seq:>4d}{i_code:1s}   {x:>8.3f}{y:>8.3f}{z:>8.3f}{occupancy:6.2f}{temp_factor:6.2f}          {element:>2s}{charge:2s}\n".format(
                        **field
                    )
                    output += line
                elif record_name == "TER   ":
                    # line = "TER   %5d      %3s %c%4d%c\n" % (serial, res_name, chain_id, res_seq, i_code)
                    field = {
                        "serial": serial,
                        "res_name": res_name,
                        "chain_id": chain_id,
                        "res_seq": res_seq,
                        "i_code": i_code,
                    }
                    line = "TER   {serial:>5d}      {res_name:3s} {chain_id:1s}{res_seq:>4d}{i_code:1s}\n".format(
                        **field
                    )
                    output += line
        return output

    def get_modpdb_atomgroup(self, ag_protein):
        """Change the name of a residue or atom.
            Called by 'set_by_atomgroup()'

        Args:
            ag_protein (AtomGroup): AtomGroup object of protein

        Returns:
            AtomGroup: renamed object
        """
        assert isinstance(ag_protein, AtomGroup)
        mode = self._mode

        retval = AtomGroup(ag_protein)
        for model_key, model in retval.groups():
            for chain_key, chain in model.groups():
                for res_key, res in chain.groups():
                    res = self._modpdb_res(res, mode)
                    for atom_key, atom in res.atoms():
                        atom = self._modpdb_resatom(res, atom, mode)
                        atom = self._modpdb_atom(atom, mode)
        return retval

    def _modpdb_atom(self, atom, mode=None):
        new_name = atom.name

        atomname = atom.name.strip().upper()
        symbol = atom.symbol
        if mode == "AMBER":
            for item in self._modpdb_amber_atm_tbl:
                if (item["name"] == atomname) and (item["symbol"] == symbol):
                    new_name = item["rename"]

            # if atomname in self._modpdb_amber_atm_tbl:
            #     new_name = self._modpdb_amber_atm_tbl[atomname]
            # if len(new_name) != 4:
            #     if (len(new_name) < 4) and (atomname[0] == symbol[0]):
            #         new_name = " {}".format(new_name)

            # 原子名対策
            # new_name_lstrip = new_name.lstrip()
            # if len(new_name_lstrip) >= 2:
            #     new_name_2 = new_name_lstrip[0:2]
            #     if new_name_2.upper() == symbol.upper():
            #         new_name = " " * (len(new_name) - len(new_name_lstrip)) + symbol + new_name_lstrip[2:]
        else:
            # "FORMAL"
            for item in self._modpdb_formal_atm_tbl:
                if (item["name"] == atomname) and (item["symbol"] == symbol):
                    new_name = item["rename"]

            # if atomname in self._modpdb_formal_atm_tbl:
            #     new_name = self._modpdb_formal_atm_tbl[atomname]
            # if len(new_name) != 4:
            #     if (0 < len(new_name)) and (len(new_name) < 4) and (atomname[0] == symbol[0]):
            #         new_name = " {}".format(new_name)
            #     else:
            #         new_name = symbol

        atom.name = new_name
        return atom

    def _modpdb_res(self, res, mode=None):
        # rename HIS name in the AMBER mode
        if mode == "AMBER":
            res = self._rename_to_amber_dialect(res)

        # rename res.name by using residue name table
        resname = res.name.upper()
        resname = resname.strip()
        if mode == "AMBER":
            if resname in self._modpdb_amber_res_tbl:
                res.name = self._modpdb_amber_res_tbl[resname]
        else:
            if resname in self._modpdb_formal_res_tbl:
                res.name = self._modpdb_formal_res_tbl[resname]
        return res

    def _modpdb_resatom(self, res, atom, mode=None):
        resname = res.name.upper()
        resname = resname.strip().lstrip()
        atom_name = atom.name.strip().lstrip().upper()

        if mode == "AMBER":
            if resname in self._modpdb_amber_resatom_table:
                if atom_name in self._modpdb_amber_resatom_table[resname]:
                    atom.name = self._modpdb_amber_resatom_table[resname][atom_name]
                    print(":{}@{} -> :{}@{}".format(resname, atom_name, resname, atom.name))
        else:
            if resname in self._modpdb_formal_resatom_table:
                if atom_name in self._modpdb_formal_resatom_table[resname]:
                    atom.name = self._modpdb_formal_resatom_table[resname][atom_name]

        return atom

    def _rename_to_amber_dialect(self, res):
        """translate HIS to HID, HIE or HIP"""
        assert isinstance(res, AtomGroup)
        if res.name == "HIS":
            # check kinds of "HIS"
            has_delta_H = False
            has_epsilon_H = False
            if res.has_atomname("HD1") and res.has_atomname("HD2"):
                has_delta_H = True
            if res.has_atomname("HE1") and res.has_atomname("HE2"):
                has_epsilon_H = True

            if has_delta_H and has_epsilon_H:
                res.name = "HIP"
            elif has_delta_H:
                res.name = "HID"
            elif has_epsilon_H:
                res.name = "HIE"
            logger.debug("found HIS: rename HIS to {}".format(res.name))

        return res


def main():
    # initialize

    # parse args
    parser = optparse.OptionParser(usage="%prog [options] PDB_FILE", version="%prog 1.0")
    parser.add_option("-o", "--output", dest="output_path", help="PDB output file", metavar="FILE")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_false", default=False, help="print message")
    (opts, args) = parser.parse_args()

    if len(args) == 0:
        parser.print_help()
        sys.exit(1)

    # setting
    file_path = args[0]
    verbose = opts.verbose

    #
    pdb_obj = Pdb(file_path)
    print(pdb_obj)

    # end


if __name__ == "__main__":
    main()
