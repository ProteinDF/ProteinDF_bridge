#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import pprint

import proteindf_bridge as bridge

import logging
import logging.config


def get_rest_of_frame_molecule(frame_molecule, selected_molecule):
    # calc the rest
    selector = bridge.Select_AtomGroup(selected_molecule)
    selected = frame_molecule.select(selector)
    rest_molecule = frame_molecule ^ selected

    return rest_molecule


def assign_rest_molecule(rest_molecule, output_atom_group,
                         model_id="model_1", chain_id="Z", res_name="UNK"):
    chain = bridge.AtomGroup()
    res = bridge.AtomGroup()
    res.name = res_name
    atom_id = 1
    for atom in rest_molecule.get_atom_list():
        res.set_atom(atom_id, atom)
        atom_id += 1
    chain.set_group(1, res)

    output_atom_group[model_id].set_group(chain_id, chain)


def main():
    parser = argparse.ArgumentParser(
        description='restructure brd file by reference file')
    parser.add_argument('target_brd_path',
                        nargs=1,
                        help='target brd file')
    parser.add_argument('ref_brd_path',
                        nargs=1,
                        help='reference brd file')
    parser.add_argument('-o', '--output_path',
                        nargs=1,
                        default=["output.brd"])
    parser.add_argument('-r', '--range',
                        nargs=1,
                        default=[1.0E-5])
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False)
    args = parser.parse_args()
    # print(args)

    target_brd_path = args.target_brd_path[0]
    ref_brd_path = args.ref_brd_path[0]
    output_path = args.output_path[0]
    range = float(args.range[0])
    verbose = args.verbose

    if verbose:
        print("target: {}".format(target_brd_path))
        print("reference: {}".format(ref_brd_path))

    # load
    target_ag = bridge.load_atomgroup(target_brd_path)
    ref_ag = bridge.load_atomgroup(ref_brd_path)

    # matching
    #target_selector = bridge.Select_AtomGroup(target_ag)
    #restructured = ref_ag.select(target_selector)

    # calc the rest
    #rest_of_target = get_rest_of_frame_molecule(target_ag, restructured)
    #assign_rest_molecule(rest_of_target, restructured)

    restructured = target_ag.restructure(ref_ag, range)

    if output_path:
        if verbose:
            print("output brd file: {}".format(output_path))
        bridge.save_atomgroup(restructured, output_path)


if __name__ == '__main__':
    #import cProfile
    #pr = cProfile.Profile()
    # pr.enable()
    main()
    # pr.disable()
    # pr.dump_stats('program.profile')
