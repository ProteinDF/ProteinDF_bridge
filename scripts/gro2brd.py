#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pprint

import proteindf_bridge as bridge


def main():
    # parse args
    parser = argparse.ArgumentParser(
        description='translate from gro to bridge file')
    parser.add_argument('GRO_FILE',
                        nargs=1,
                        help='gro file')
    parser.add_argument('BRD_FILE',
                        nargs=1,
                        help='brd file')
    parser.add_argument('-v', '--verbose',
                        action="store_true",
                        default=False)
    args = parser.parse_args()

    # setting
    gro_file_path = args.GRO_FILE[0]
    output_path = args.BRD_FILE[0]
    verbose = args.verbose

    # load gro file
    if (verbose == True):
        print("reading: %s\n" % (gro_file_path))
    gro_obj = bridge.SimpleGro()
    gro_obj.load(gro_file_path)
    # print(gro_obj)
    atom_group = gro_obj.get_atomgroup()

    # output DfData
    data = atom_group.get_raw_data()
    # pprint.pprint(data)

    # output file
    if (verbose == True):
        print("writing: %s\n" % (output_path))
    bridge.save_msgpack(data, output_path)


if __name__ == '__main__':
    main()
