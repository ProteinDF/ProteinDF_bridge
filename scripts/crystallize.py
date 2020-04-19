#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
try:
    import msgpack
except:
    import msgpack_pure as msgpack

import proteindf_bridge as bridge


def main():
    # parse args
    parser = argparse.ArgumentParser(description='crystallize molecules')
    parser.add_argument('INPUT_BRD_PATH',
                        nargs=1,
                        help='input bridge file path')
    parser.add_argument('OUTPUT_BRD_PATH',
                        nargs=1,
                        help='output bridge file path')

    parser.add_argument("--num_x",
                        nargs=1,
                        type=int,
                        default=[1],
                        help="number of molecules for x")
    parser.add_argument("--num_y",
                        nargs=1,
                        type=int,
                        default=[1],
                        help="number of molecules for y")
    parser.add_argument("--num_z",
                        nargs=1,
                        type=int,
                        default=[1],
                        help="number of molecules for z")

    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False)
    args = parser.parse_args()

    print(args)
    input_path = args.INPUT_BRD_PATH[0]
    output_path = args.OUTPUT_BRD_PATH[0]
    num_x = args.num_x[0]
    num_y = args.num_y[0]
    num_z = args.num_z[0]
    verbose = args.verbose
    if verbose:
        print("INPUT: {}".format(input_path))
        print("#mols in cell: {} x {} x {}".format(num_x, num_y, num_z))

    # load input file
    ag = bridge.AtomGroup()
    with open(input_path, "rb") as f:
        mpac_data = msgpack.unpackb(f.read())
        ag = bridge.AtomGroup(mpac_data)
    # print(ag)
    num_of_res = ag.get_number_of_groups()
    (pos_min, pos_max) = ag.box()
    cell_size = pos_max - pos_min
    if verbose:
        print("residues: {}".format(num_of_res))
        print("min: {}".format(pos_min))
        print("max: {}".format(pos_max))
        print("cell size: {}".format(cell_size))

    # make cell
    cell = bridge.AtomGroup()
    shift = bridge.Position()
    count = 1
    for count_x in range(num_x):
        shift.x = count_x * cell_size.x
        for count_y in range(num_y):
            shift.y = count_y * cell_size.y
            for count_z in range(num_z):
                shift.z = count_z * cell_size.z

                new_ag = bridge.AtomGroup(ag)
                new_ag.shift_by(shift)

                for resid, res in new_ag.groups():
                    cell.set_group(count, res)
                    count += 1
    # print(cell)

    # output
    if verbose:
        print("OUTPUT: {}".format(output_path))
    with open(output_path, "wb") as f:
        data = cell.get_raw_data()
        mpac = msgpack.packb(data)
        f.write(mpac)


if __name__ == '__main__':
    main()
