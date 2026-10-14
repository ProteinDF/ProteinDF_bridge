#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import msgpack

import proteindf_bridge.rs as rs


def main():
    parser = argparse.ArgumentParser(description='parse mmCIF file to mol2 file')
    parser.add_argument('mmCIF_FILE',
                        nargs=1,
                        help="mmCIF file")
    parser.add_argument("-o", "--output",
                        nargs=1,
                        help="mol2 output directory")
    parser.add_argument("-w", "--write",
                        nargs=1,
                        help="write scratch file")
    parser.add_argument("-r", "--read",
                        nargs=1,
                        help="read scratch file")
    parser.add_argument('-v', '--verbose',
                        action="store_true",
                        default=False)
    args = parser.parse_args()

    mmcif_file_path = args.mmCIF_FILE[0]
    output_path = ""
    if args.output:
        output_path = args.output[0]

    write_scratch_path = ""
    if args.write:
        write_scratch_path = args.write[0]
    read_scratch_path = ""
    if args.read:
        read_scratch_path = args.read[0]

    verbose = args.verbose

    mmcif = None
    if len(read_scratch_path) > 0:
        # read scratch file
        with open(read_scratch_path, "rb") as f:
            scratch_data = msgpack.unpackb(f.read(), raw=False)
        mmcif = rs.SimpleMmcif.from_data(scratch_data)
    else:
        # load mmCIF file
        if verbose:
            print("reading: {}".format(mmcif_file_path))
            print("output dir: {}".format(output_path))
        mmcif = rs.SimpleMmcif(mmcif_file_path)

        if len(write_scratch_path) > 0:
            if verbose:
                print("write scratch file: {}".format(write_scratch_path))
            with open(write_scratch_path, "wb") as f:
                f.write(msgpack.packb(mmcif.data))

    #print(mmcif)

    if len(output_path) > 0:
        data_ids = mmcif.get_molecule_names()
        for data_id in data_ids:
            print("{}...".format(data_id))
            ag = mmcif.get_atomgroup(data_id)
            #print(ag)

            mol2_path = os.path.join(output_path, "{}.mol2".format(data_id))
            if verbose:
                print("make mol2: {}".format(mol2_path))
            mol2 = rs.SimpleMol2(ag)
            if verbose:
                print("save: {}".format(mol2_path))
            mol2.save(mol2_path)


if __name__ == '__main__':
    main()
