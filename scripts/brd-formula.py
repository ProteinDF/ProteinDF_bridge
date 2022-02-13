#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse

import proteindf_bridge as bridge


def main():

    # parse args
    parser = argparse.ArgumentParser(description='print molecular formula')
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
    mpac_data = bridge.load_msgpack(mpac_file_path)

    # prepare atomgroup
    atom_group = bridge.AtomGroup(mpac_data)

    # output
    formula = atom_group.formula()

    print("{}: MW={}".format(formula, atom_group.weight))


if __name__ == '__main__':
    main()
