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

import pdfbridge

def main():
    # parse args
    parser = argparse.ArgumentParser(description='setup bonds')
    parser.add_argument('FILE',
                        nargs=1,
                        help='bridge file')
    parser.add_argument('-o', '--output',
                        nargs='?',
                        help='output brd file')
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default = False)
    args = parser.parse_args()

    # setting
    mpac_file_path = args.FILE[0]
    output = args.output
    if output == '':
        output = 'output.brd'
    verbose = args.verbose

    # reading
    if verbose:
        print("reading: %s\n" % (mpac_file_path))
    mpac_file = open(mpac_file_path, "rb")
    mpac_data =msgpack.unpackb(mpac_file.read())
    mpac_file.close()
        
    # prepare atomgroup
    atom_group = pdfbridge.AtomGroup(mpac_data)
    #print(atom_group)

    bonds = pdfbridge.Bond()
    bonds.setup(atom_group)
    
    # output
    mpac = msgpack.packb(atom_group.get_dict_data())
    if (verbose == True):
        print("writing: %s\n" % (output))
    savefile = open(output, "wb");
    savefile.write(mpac)
    savefile.close()


if __name__ == '__main__':
    main()
