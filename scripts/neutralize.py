#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import logging
import logging.config
try:
    import msgpack
except:
    import msgpack_pure as msgpack

import pdfbridge

def main():
    # parse args
    parser = argparse.ArgumentParser(description='neutralize protein')
    parser.add_argument('INPUT_FILE',
                        nargs=1,
                        help='bridge file')
    parser.add_argument('OUTPUT_FILE',
                        nargs=1,
                        help='output file (bridge format)')
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default = False)
    parser.add_argument("-d", "--debug",
                        action="store_true",
                        default = False)
    args = parser.parse_args()
        

    # setting
    input_path = args.INPUT_FILE[0]
    output_path = args.OUTPUT_FILE[0]
    is_verbose = args.verbose

    # logging
    logging.basicConfig()
    logger = logging.getLogger(__name__)
    handler = logging.StreamHandler()
    if args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    logger.addHandler(handler)
        
    # reading
    protein = None
    with open(input_path, "rb") as mpac_file:
        mpac_data = msgpack.unpackb(mpac_file.read())
        protein = pdfbridge.AtomGroup(mpac_data)

    # neutralize
    neutralizer = pdfbridge.Neutralize(protein)
    mod_protein = neutralizer.neutralized

    # output
    with open(output_path, "wb") as mpac_file:
        raw_data = mod_protein.get_raw_data()
        mpac_file.write(msgpack.packb(raw_data))
    
if __name__ == '__main__':
    main()
