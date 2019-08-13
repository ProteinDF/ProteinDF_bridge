#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import proteindf_bridge as bridge

def divide_lowest_atomgroup(atomgroup):
    group_list = []
    for subgrp_key, subgrp in atomgroup.groups():
        from_subgrp = divide_lowest_atomgroup(subgrp)
        group_list.extend(from_subgrp)

    subgroup = bridge.AtomGroup()
    for atom_key, atom in atomgroup.atoms():
        subgroup.add(atom_key, atom)
    if len(subgroup) > 0:
        subgroup.name = atomgroup.name
        group_list.append(subgroup)

    return group_list

def divide_mainchain(atomgroup):
    main_chain_atoms = ['N', 'H', 'CA', 'HA', 'C', 'O']

    group_list = []
    for subgrp_key, subgrp in atomgroup.groups():
        from_subgrp = divide_mainchain(subgrp)
        group_list.extend(from_subgrp)

    main_chain = bridge.AtomGroup()
    side_chain = bridge.AtomGroup()
    for atom_key, atom in atomgroup.atoms():
        if atom.name in main_chain_atoms:
            main_chain.set_atom(atom_key, atom)
        else:
            side_chain.set_atom(atom_key, atom)

    if main_chain.get_number_of_atoms() > 0:
        main_chain.name = atomgroup.name + "_mainchain"
        group_list.append(main_chain)
    if side_chain.get_number_of_atoms() > 0:
        side_chain.name = atomgroup.name + "_sidechain"
        group_list.append(side_chain)

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
                        default = False)
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
    group_list = divide_mainchain(atomgroup)

    # output
    if output_path:
        blocks = bridge.AtomGroup()
        block_index = 1
        for atomgroup in group_list:
            blocks.set_group(block_index, atomgroup)
            block_index += 1

        if verbose:
            print("save: {}".format(output_path))
        bridge.save_atomgroup(blocks, output_path)

    else:
        for atomgroup in group_list:
            print("---- {} ----".format(atomgroup.name))
            for atom_key, atom in atomgroup.atoms():
                print("{atom}: {atom_key}".format(atom_key=atom_key, atom=str(atom)))


if __name__ == '__main__':
    main()
