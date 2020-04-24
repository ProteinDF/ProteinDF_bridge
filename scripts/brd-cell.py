#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import math

try:
    import msgpack
except:
    import msgpack_pure as msgpack

import proteindf_bridge as bridge


def main():
    # parse args
    parser = argparse.ArgumentParser(description='calc cell size')
    parser.add_argument('FILE',
                        nargs=1,
                        help='bridge file')
    parser.add_argument('num_of_molecules',
                        nargs=1,
                        type=int,
                        default=100,
                        help='the number of molecules')
    parser.add_argument('density',
                        nargs=1,
                        type=float,
                        default=1.0,
                        help='density (g cm^-3')
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False)
    args = parser.parse_args()

    # setting
    mpac_file_path = args.FILE[0]
    num_of_mols = args.num_of_molecules[0]
    density = args.density[0]
    verbose = args.verbose

    # reading
    if (verbose == True):
        print("reading: {}}".format(mpac_file_path))
        print("molecules: {}".format(num_of_molecules))
        print("density: {}".format(density))
    mpac_file = open(mpac_file_path, "rb")
    mpac_data = msgpack.unpackb(mpac_file.read())
    mpac_file.close()

    # prepare atomgroup
    atom_group = bridge.AtomGroup(mpac_data)
    weights = atom_group.weight * num_of_mols / bridge.AVOGADRO_CONST

    volume = weights / (density / 1.0e-6)
    cell_size = math.pow(volume, 1.0/3.0) / (1.0e-10)

    # output
    if (verbose == True):
        print("volume: {}".format(volume))
    print("cell size: {} Angstrom".format(cell_size))


if __name__ == '__main__':
    main()
