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
    parser = argparse.ArgumentParser(description='bridge file selecter')
    parser.add_argument('FILE',
                        nargs=1,
                        help='bridge file')
    parser.add_argument('-o', '--output',
                        nargs=1,
                        type=str,
                        default=[''],
                        help='output brd file')
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default = False)
    parser.add_argument('-q', '--query',
                        nargs=1,
                        type=str,
                        default=["*"],
                        help='select by using path string')
    args = parser.parse_args()

    # setting
    mpac_file_path = args.FILE[0]
    output_path = args.output[0]
    verbose = args.verbose

    query = args.query[0]

    # reading
    if verbose:
        print("reading: {}".format(mpac_file_path))
    mpac_file = open(mpac_file_path, "rb")
    mpac_data =msgpack.unpackb(mpac_file.read())
    mpac_file.close()

    # prepare atomgroup
    atomgroup = bridge.AtomGroup(mpac_data)
    #print(atom_group)

    # selecter
    if verbose:
        print('query=\"{}\"'.format(query))
    path_selecter = bridge.Select_Path_wildcard(query)
    selected = atomgroup.select(path_selecter)

    # output
    if len(output_path) > 0:
        if (verbose == True):
            print("writing: %s\n" % (output_path))
        with open(output_path, "wb") as output_file:
            output_data = selected.get_raw_data()
            output_mpac = msgpack.packb(output_data)
            output_file.write(output_mpac)
    else:
        print(selected)
    

if __name__ == '__main__':
    main()
