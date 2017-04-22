#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re

class AmberPrmtop(object):
    def __init__(self, prmtop_path):
        self._charges = None
        self._load(prmtop_path)

    @property
    def charges(self):
        return self._charges
        
    def _load(self, prmtop_path):
        with open(prmtop_path, "r") as f:
            for line in f:
                line = line.rstrip()
                #print(line)
                if line == '%FLAG CHARGE':
                    self._read_charges(f)

    def _read_charges(self, f):
        print('>>>> read charge')
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
