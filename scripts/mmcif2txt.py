#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import proteindf_bridge as bridge

def main():
    parser = argparse.ArgumentParser(description='parse mmCIF file')
    parser.add_argument('mmCIF_FILE',
                        nargs=1,
                        help="mmCIF file")
    parser.add_argument('-w', '--write',
                        nargs=1,
                        help="save message pack file")
    parser.add_argument('-v', '--verbose',
                        action="store_true",
                        default=False)
    args = parser.parse_args()

    mmcif_file_path = args.mmCIF_FILE[0]
    write_mpac_file_path = ""
    if args.write:
        write_mpac_file_path = args.write[0]
    verbose = args.verbose

    # load mmCIF file
    if verbose:
        print("reading: {}".format(mmcif_file_path))
    mmcif = bridge.SimpleMmcif(mmcif_file_path)

    if len(write_mpac_file_path) > 0:
        if verbose:
            print("write mpac file: {}".format(write_mpac_file_path))
        mmcif.save_msgpack(write_mpac_file_path)

    # print(mmcif)

if __name__ == '__main__':
    main()
