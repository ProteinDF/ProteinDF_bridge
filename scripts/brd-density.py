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
    parser = argparse.ArgumentParser(description='calc density')
    parser.add_argument('FILE',
                        nargs=1,
                        help='bridge file')
    parser.add_argument('num_of_molecules',
                        nargs=1,
                        type=int,
                        default=100,
                        help='the number of molecules')
    parser.add_argument('--box',
                        nargs=3,
                        type=float,
                        help='box size')
    parser.add_argument("--use-nm",
                        action="store_true",
                        default=False,
                        help='output length by nm unit (default: angstrom)')
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False)
    args = parser.parse_args()

    # setting
    mpac_file_path = args.FILE[0]
    num_of_mols = args.num_of_molecules[0]
    (box_x, box_y, box_z) = args.box
    is_use_nm = args.use_nm
    verbose = args.verbose

    if (is_use_nm):
        # nm -> angstrom
        box_x *= 10.0
        box_y *= 10.0
        box_z *= 10.0
    # reading
    if (verbose == True):
        print("reading: {}".format(mpac_file_path))
        print("molecules: {}".format(num_of_mols))
        print("box [A]: {} {} {}".format(box_x, box_y, box_z))
        print("use nm: {}".format(is_use_nm))

    mpac_file = open(mpac_file_path, "rb")
    mpac_data = msgpack.unpackb(mpac_file.read())
    mpac_file.close()

    # prepare atomgroup
    atom_group = bridge.AtomGroup(mpac_data)
    weights = atom_group.weight * num_of_mols / bridge.AVOGADRO_CONST
    volume = box_x * box_y * box_z * 1.0E-24
    if (verbose == True):
        print("volume [cm^-3]: {}".format(volume))

    density = weights / volume

    print("density [g / cm^-3]: {}".format(density))


if __name__ == '__main__':
    main()
