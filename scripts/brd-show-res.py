#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import proteindf_bridge as bridge


def show_protein(model):
    for chain_key, chain in model.groups():
        for res_key, res in chain.groups():
            print("{}, {}, {}".format(chain_key, res_key, res.name))


def main():
    # parse args
    parser = argparse.ArgumentParser(description='print residues in the bridge file')
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
    atomgroup = bridge.AtomGroup(mpac_data)

    if bridge.Format.is_models(atomgroup):
        for model_key, model in atomgroup.groups():
            print(">>>> model: {}".format(model_key))
            show_protein(model)
    elif bridge.Format.is_protein(atomgroup):
        show_protein(atomgroup)
    else:
        print("unknown bridge file format.")


if __name__ == '__main__':
    main()
