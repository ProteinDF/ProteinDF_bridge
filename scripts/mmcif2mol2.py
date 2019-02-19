#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import pprint
import proteindf_bridge as bridge

def main():
    parser = argparse.ArgumentParser(description='parse mmCIF file to mol2 file')
    parser.add_argument('mmCIF_FILE',
                        nargs=1,
                        help="mmCIF file")
    #parser.add_argument('-r', '--read',
    #                    nargs=1,
    #                    help="read message pack file")
    parser.add_argument("-o", "--output",
                        nargs=1,
                        help="mol2 output directory")
    parser.add_argument('-v', '--verbose',
                        action="store_true",
                        default=False)
    args = parser.parse_args()

    mmcif_file_path = args.mmCIF_FILE[0]
    #read_mpac_file_path = ""
    #if args.read:
    #    read_mpac_file_path = args.write[0]
    output_path = "."
    if args.output:
        output_path = args.output[0]
    verbose = args.verbose

    # load mmCIF file
    if verbose:
        print("reading: {}".format(mmcif_file_path))
        print("output dir: {}".format(output_path))
    mmcif = bridge.SimpleMmcif(mmcif_file_path)

    #print("=" * 80)
    #print(mmcif)
    data_ids = mmcif.get_molecule_names()
    for data_id in data_ids:
        print("{}...".format(data_id))
        ag = mmcif.get_atomgroup(data_id)
        #print(ag)

        mol2_path = os.path.join(output_path, "{}.mol2".format(data_id))
        if verbose:
            print("output: {}".format(mol2_path))
        mol2 = bridge.SimpleMol2(ag)
        mol2.save(mol2_path)


if __name__ == '__main__':
    main()
