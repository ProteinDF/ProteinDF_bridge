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

import proteindf_bridge as bridge


def main():
    # parse args
    parser = argparse.ArgumentParser(
        description='transform bridge file to XYZ file')
    parser.add_argument('FILE',
                        nargs=1,
                        help='bridge file')
    parser.add_argument('-o', '--output',
                        nargs='?',
                        help='XYZ output file')
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False)
    args = parser.parse_args()

    # setting
    mpac_file_path = args.FILE[0]
    output = args.output
    verbose = args.verbose

    # reading
    if (verbose == True):
        print("reading: %s\n" % (mpac_file_path))
    mpac_data = bridge.load_msgpack(mpac_file_path)

    # prepare atomgroup
    atom_group = bridge.AtomGroup(mpac_data)
    # print(atom_group)

    # prepare XYZ object
    xyz_obj = bridge.Xyz()
    xyz_obj.set_by_atomgroup(atom_group)

    # output XYZ
    if output:
        fout = open(output, "w")
        contents = str(xyz_obj)
        fout.write(contents)
        fout.close()
    else:
        print(xyz_obj)

    # end


if __name__ == '__main__':
    main()
