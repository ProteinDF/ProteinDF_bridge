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
try:
    import msgpack
except:
    import msgpack_pure as msgpack

import proteindf_bridge as bridge

def main():
    # parse args
    parser = argparse.ArgumentParser(description='transform XYZ file to bridge file')
    parser.add_argument('XYZ_PATH',
                        nargs=1,
                        help='xyz file path')
    parser.add_argument('BRD_PATH',
                        nargs=1,
                        help='bridge file path')
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default = False)
    args = parser.parse_args()

    # setting
    xyz_file_path = args.XYZ_PATH[0]
    brd_file_path = args.BRD_PATH[0]
    verbose = args.verbose

    # reading
    if (verbose == True):
        print("reading: %s\n" % (xyz_file_path))
    xyz = bridge.Xyz()
    xyz.load(xyz_file_path)
    atomgroup = xyz.get_atom_group()

    brd_file = open(brd_file_path, "wb");
    mpac = msgpack.packb(atomgroup.get_raw_data())
    brd_file.write(mpac)
    brd_file.close()

if __name__ == '__main__':
    main()
