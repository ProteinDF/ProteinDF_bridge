#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse

import proteindf_bridge as bridge


def main():
    # parse args
    parser = argparse.ArgumentParser(
        description='renumber resid in bridge file')
    parser.add_argument('FILE',
                        nargs=1,
                        help='bridge file')
    parser.add_argument('increment',
                        nargs=1,
                        type=int,
                        help='increment number')
    parser.add_argument('-o', '--output',
                        nargs=1,
                        type=str,
                        default=[''],
                        help='output brd file')
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False)
    parser.add_argument('-q', '--query',
                        nargs=1,
                        type=str,
                        default=["*"],
                        help='select by using path string')
    args = parser.parse_args()

    # setting
    mpac_file_path = args.FILE[0]
    increment = args.increment[0]
    output_path = args.output[0]
    verbose = args.verbose

    query = args.query[0]

    # reading
    if verbose:
        print("reading: {}".format(mpac_file_path))
    models = bridge.load_atomgroup(mpac_file_path)
    # print(models)

    # renumber
    if verbose:
        print("increment: {}".format(increment))
    new_models = bridge.AtomGroup()
    for model_key, model in models.groups():
        new_model = bridge.AtomGroup()
        new_model.name = model.name
        for chain_key, chain in model.groups():
            new_chain = bridge.AtomGroup()
            new_chain.name = chain.name
            for res_id, res in chain.groups():
                new_chain.set_group(int(res_id) + increment, res)

            new_model.set_group(chain_key, new_chain)
        new_models.set_group(model_key, new_model)

    # output
    if len(output_path) > 0:
        if (verbose == True):
            print("writing: %s\n" % (output_path))
        bridge.save_atomgroup(new_models, output_path)
    else:
        print(new_models)


if __name__ == '__main__':
    main()
