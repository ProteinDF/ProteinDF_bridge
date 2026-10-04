#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import proteindf_bridge as bridge


def divide_lowest_atomgroup(atomgroup):
    group_list = []
    for subgrp_key, subgrp in atomgroup.groups():
        from_subgrp = divide_lowest_atomgroup(subgrp)
        group_list.extend(from_subgrp)

    subgroup = []
    for atom_key, atom in atomgroup.atoms():
        subgroup.append(atom)
    if len(subgroup) > 0:
        group_list.append(subgroup)

    return group_list


def main():
    # parse args
    parser = argparse.ArgumentParser(description='bridge file divider')
    parser.add_argument('brd_path',
                        nargs=1,
                        help='bridge file')
    parser.add_argument('-o', '--output',
                        nargs=1,
                        type=str,
                        default=[''],
                        help='output brd file')
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False)
    args = parser.parse_args()

    # setting
    brd_path = args.brd_path[0]
    output_path = args.output[0]
    verbose = args.verbose

    #
    if verbose:
        print("load: {}".format(brd_path))
    atomgroup = bridge.load_atomgroup(brd_path)

    #
    group_list = divide_lowest_atomgroup(atomgroup)

    # output
    if output_path:
        blocks = bridge.AtomGroup()
        block_index = 1
        for grp in group_list:
            block = bridge.AtomGroup()
            atom_index = 1
            for atom in grp:
                block.set_atom(atom_index, atom)
                atom_index += 1
            blocks.set_group(block_index, block)
            block_index += 1

        if verbose:
            print("save: {}".format(output_path))
        bridge.save_atomgroup(blocks, output_path)

    else:
        for grp in group_list:
            print("----")
            for atom in grp:
                print(atom)


if __name__ == '__main__':
    main()
