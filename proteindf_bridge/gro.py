#!/usr/bin/env python

import logging
logger = logging.getLogger(__name__)

from .position import Position
from .atom import Atom
from .atomgroup import AtomGroup
from .periodictable import PeriodicTable

class SimpleGro(object):
    """
    Simple .gro file

    sample:
    MD of 2 waters, t= 0.0
        6
        1WATER  OW1    1   0.126   1.624   1.679  0.1227 -0.0580  0.0434
        1WATER  HW2    2   0.190   1.661   1.747  0.8085  0.3191 -0.7791
        1WATER  HW3    3   0.177   1.568   1.613 -0.9045 -2.6469  1.3180
        2WATER  OW1    4   1.275   0.053   0.622  0.2519  0.3140 -0.1734
        2WATER  HW2    5   1.337   0.002   0.680 -1.0641 -1.1349  0.0257
        2WATER  HW3    6   1.326   0.120   0.568  1.9427 -0.8216 -0.0244
    1.82060   1.82060   1.82060
    """
    def __init__(self):
        self._title = ""
        self._num_of_atoms = 0
        self._atoms = []

    def load(self, gro_filepath):
        with open(gro_filepath, "r") as f:
            # line 1
            line = f.readline()
            line = line.rstrip()
            self._title = line

            # line 2
            line = f.readline()
            line = line.rstrip()
            self._num_of_atoms = int(line)

            # atom lines
            print("\n")
            for i in range(self._num_of_atoms):
                line = f.readline()
                line = line.rstrip()

                residue_number = int(line[0:5])
                residue_name = line[5:10].strip()
                atom_name = line[10:15].strip()
                atom_number = int(line[15:20])
                # position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
                position_x = float(line[20:28]) 
                position_y = float(line[28:36]) 
                position_z = float(line[36:44]) 
                # velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with 4 decimal places)
                velocity_x = float(line[44:52])
                velocity_y = float(line[52:60])
                velocity_z = float(line[60:68])

                atom_data = (residue_number, residue_name, atom_name, atom_number,
                             position_x, position_y, position_z,
                            velocity_x, velocity_y, velocity_z)
                #print(atom_data)
                self._atoms.append(atom_data)

            # box vectors
            line = f.readline()
            line = line.strip()
            box_vectors = line.split()
            box_vectors = [float(x) for x in box_vectors]
            for i in range(len(box_vectors), 6):
                box_vectors.append(0.0)
            self._box_vectors = box_vectors


    def get_atomgroup(self):
        ag = AtomGroup()

        pt = PeriodicTable()

        current_res_id = -1
        current_ag = None
        for atom_data in self._atoms:
            (res_id, res_name, name, id, x, y, z, vx, vy, vz) = atom_data
            if current_res_id != res_id:
                if current_ag != None:
                    ag.set_group(current_res_id, current_ag)
                current_res_id = res_id
                current_ag = AtomGroup()
                current_ag.name = res_name

            atom = Atom()
            name = name.strip()
            atom.name = name

            symbol = ""
            if len(name) == 1:
                symbol = name
            else:
                name2 = name[0:2]
                if name2 in pt:
                    symbol = name2
                else:
                    symbol = name[0]
            atom.symbol = symbol
            atom.position = Position(x * 10.0, y * 10.0, z * 10.0) # nm -> angstrom
            current_ag.set_atom(id, atom)

        if current_ag != None:
            ag.set_group(current_res_id, current_ag)

        return ag


    def __str__(self):
        answer = ""
        answer += "{}\n".format(self._title)
        answer += "{}\n".format(self._num_of_atoms)
        for i in range(self._num_of_atoms):
            answer += "{:>5}{:<5}{:>5}{:>5}{:8.3f}{:8.3f}{:8.3f}{:8.4f}{:8.4f}{:8.4f}\n".format(
                self._atoms[i][0], self._atoms[i][1], self._atoms[i][2], self._atoms[i][3],
                self._atoms[i][4], self._atoms[i][5], self._atoms[i][6],
                self._atoms[i][4], self._atoms[i][5], self._atoms[i][6])
        answer += " ".join(map(str, self._box_vectors))

        return answer

