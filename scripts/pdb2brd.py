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

import sys
import argparse
import pprint

import proteindf_bridge as bridge


def main():
    # initialize

    # parse args
    parser = argparse.ArgumentParser(
        description='translate from PDB to bridge file')
    parser.add_argument('PDB_FILE',
                        nargs=1,
                        help='pdb file')
    parser.add_argument('BRD_FILE',
                        nargs=1,
                        help='brd file')
    parser.add_argument('-m', '--model',
                        nargs=1,
                        default=[1],
                        type=int,
                        help='select model')
    parser.add_argument('-l', '--alt_loc',
                        nargs=1,
                        default='A',
                        help='select alc_loc')
    parser.add_argument('-d', '--debug',
                        action='store_true',
                        default=False)
    parser.add_argument('-v', '--verbose',
                        action="store_true",
                        default=False)
    args = parser.parse_args()

    # setting
    pdb_file_path = args.PDB_FILE[0]
    output_path = args.BRD_FILE[0]
    select_model = args.model[0]
    select_altloc = args.alt_loc[0]
    debug = args.debug
    verbose = args.verbose

    # load PDB file
    if (verbose == True):
        print("reading: %s\n" % (pdb_file_path))
    pdb_obj = bridge.Pdb()
    pdb_obj.debug = debug
    pdb_obj.load(pdb_file_path)
    # print(pdb_obj)
    atom_group = pdb_obj.get_atomgroup(select_model=select_model,
                                       select_altloc=select_altloc)

    # output DfData
    data = atom_group.get_raw_data()
    # pprint.pprint(data)

    # output file
    if (verbose == True):
        print("writing: %s\n" % (output_path))
    bridge.save_msgpack(data, output_path)

    # end


if __name__ == '__main__':
    main()
