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

import bridge

def load_brd(path, verbose=False):
    if verbose:
        print("reading: %s\n" % (path))
    mpac_file = open(path, "rb")
    mpac_data =msgpack.unpackb(mpac_file.read())
    mpac_file.close()
    atomgroup = bridge.AtomGroup(mpac_data)
    return atomgroup

def main():
    # parse args
    parser = argparse.ArgumentParser(description='superpose')
    parser.add_argument('FILE1',
                        nargs=1,
                        help='bridge file')
    parser.add_argument('FILE2',
                        nargs=1,
                        help='bridge file')
    parser.add_argument('-q', '--quaternion',
                        action="store_true",
                        default = False,
                        help='use quaternion')
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default = False)
    args = parser.parse_args()

    # setting
    mpac_file_path1 = args.FILE1[0]
    mpac_file_path2 = args.FILE2[0]
    use_quaternion = args.quaternion
    verbose = args.verbose

    # reading
    atomgroup1 = load_brd(mpac_file_path1, verbose)
    atomgroup2 = load_brd(mpac_file_path2, verbose)

    if use_quaternion:
        sp = bridge.Superposer_quaternion(atomgroup1, atomgroup2)
    else:
        sp = bridge.Superposer(atomgroup1, atomgroup2)
    rmsd = sp.rmsd
    print('rmsd: {}'.format(rmsd))

    #print('>>>> ag1')
    #print(atomgroup1)
    #print('>>>> ag2')
    #print(atomgroup2)
    ag = sp.superimpose(atomgroup1)
    #print('>>>> after')
    print(ag)
    
if __name__ == '__main__':
    main()
