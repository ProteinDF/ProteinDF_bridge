#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import proteindf_bridge as bridge

def main():
    parser = argparse.ArgumentParser(description='parse mmCIF file to mol2 file')
    parser.add_argument('mmCIF_FILE',
                        nargs=1,
                        help="mmCIF file")
    #parser.add_argument('-r', '--read',
    #                    nargs=1,
    #                    help="read message pack file")
    parser.add_argument('-v', '--verbose',
                        action="store_true",
                        default=False)
    args = parser.parse_args()

    mmcif_file_path = args.mmCIF_FILE[0]
    #read_mpac_file_path = ""
    #if args.read:
    #    read_mpac_file_path = args.write[0]
    verbose = args.verbose

    # load mmCIF file
    if verbose:
        print("reading: {}".format(mmcif_file_path))
    mmcif = bridge.SimpleMmcif(mmcif_file_path)

    print("=" * 80)
    print(mmcif)

if __name__ == '__main__':
    main()
