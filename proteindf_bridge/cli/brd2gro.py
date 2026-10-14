#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import proteindf_bridge.rs as rs


def main():
    # parse args
    parser = argparse.ArgumentParser(
        description='output gro format from bridge file')
    parser.add_argument('FILE',
                        nargs=1,
                        help='bridge file')
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False)
    args = parser.parse_args()

    # setting
    mpac_file_path = args.FILE[0]
    verbose = args.verbose

    # reading
    if (verbose == True):
        print("reading: %s\n" % (mpac_file_path))
    atom_group = rs.load_atomgroup(mpac_file_path)

    gro = rs.SimpleGro()
    gro.set_by_atomgroup(atom_group)

    # output
    print(gro)


if __name__ == '__main__':
    main()
