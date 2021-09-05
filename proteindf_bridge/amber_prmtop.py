#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re

from .periodictable import PeriodicTable
from .position import Position
from .atom import Atom
from .atomgroup import AtomGroup

import logging
logger = logging.getLogger(__name__)


class AmberPrmtop(object):
    def __init__(self, prmtop_path, inpcrd_path):
        self._atom_names = None
        self._charges = None
        self._atomic_numbers = None
        self._xyz = None

        self._load_prmtop(prmtop_path)
        self._load_inpcrd(inpcrd_path)
        self._check_data()

    @property
    def atom_names(self):
        return self._atom_names

    @property
    def charges(self):
        return self._charges

    @property
    def atomic_numbers(self):
        return self._atomic_numbers

    @property
    def xyz(self):
        return self._xyz

    def get_atomgroup(self):
        atomgroup = AtomGroup()
        num_of_atoms = len(self.xyz)
        for i in range(num_of_atoms):
            atom = Atom()
            atom.symbol = PeriodicTable.get_symbol(self._atomic_numbers[i])
            atom.xyz = self._xyz[i]
            atom.charge = self._charges[i]
            atom.name = self._atom_names[i]
            atomgroup.set_atom(i, atom)

        return atomgroup

    def _load_prmtop(self, prmtop_path):
        with open(prmtop_path, "r") as f:
            for line in f:
                line = line.rstrip()
                # print(line)
                if line == '%FLAG ATOM_NAME':
                    line = self._read_atom_name(f)
                if line == '%FLAG CHARGE':
                    line = self._read_charges(f)
                if line == '%FLAG ATOMIC_NUMBER':
                    line = self._read_atomic_number(f)

    def _check_data(self):
        print("# atomic_numbers: {}".format(len(self._atomic_numbers)))
        print("# atom_names: {}".format(len(self._atom_names)))
        print("# charges: {}".format(len(self._charges)))
        print("# xyz: {}".format(len(self._xyz)))

    def _read_atom_name(self, f):
        logger.debug('read atomname')
        atom_names = []
        for line in f:
            line = line.rstrip()

            if line[0:7] == '%FORMAT':
                continue
            if line[0:5] == '%FLAG':
                break

            while len(line) > 0:
                name = line[0:4]
                name.rstrip()
                line = line[4:]
                atom_names.append(name)

        self._atom_names = atom_names
        print('line: ', line)
        return line

    def _read_charges(self, f):
        print('read charge')
        logger.debug('read charge')
        charges = []
        for line in f:
            line = line.rstrip()

            if line[0:7] == '%FORMAT':
                continue
            if line[0:5] == '%FLAG':
                break

            values = line.split()
            for v in values:
                # print(v, float(v))
                charges.append(float(v) / 18.2223)

        self._charges = charges
        print('line: ', line)
        return line

    def _read_atomic_number(self, f):
        logger.debug('read atomic number')
        atomic_numbers = []
        for line in f:
            line = line.rstrip()

            if line[0:7] == '%FORMAT':
                continue
            if line[0:5] == '%FLAG':
                break

            values = line.split()
            for v in values:
                atomic_numbers.append(int(v))

        self._atomic_numbers = atomic_numbers
        print('line: ', line)
        return line

    def _load_inpcrd(self, inpcrd_path):
        title = ""
        num_of_atoms = 0
        xyz = []
        with open(inpcrd_path, "r") as f:
            for line in f:
                line = line.rstrip()

                if len(title) == 0:
                    title = line
                    continue
                if num_of_atoms == 0:
                    values = line.split()
                    num_of_atoms = int(values[0])
                    continue

                while len(line) > 0:
                    x = line[0:12]
                    y = line[12:24]
                    z = line[24:36]
                    xyz.append(Position(float(x), float(y), float(z)))
                    line = line[36:]

        self._xyz = xyz
