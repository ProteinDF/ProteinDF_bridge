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

def get_ions(atomgroup, logger):
    ions = []
    for key, subgrp in atomgroup.groups():
        new_ions = get_ions(subgrp, logger)
        ions.extend(new_ions)
    for key, atom in atomgroup.atoms():
        if atom.symbol not in ('Na', 'Cl'):
            logger.debug("pass>'{}':'{}'@{}".format(key, atom.symbol, atomgroup.name))
        else:
            logger.debug("FOUND>'{}':'{}'@{}".format(key, atom.symbol, atomgroup.name))
            atomgroup.erase_atom(key)
            ions.append(atom)
    return ions

def get_max_resid(model):
    assert isinstance(model, pdfbridge.AtomGroup)

    max_resid = 0
    for chain_name, chain in model.groups():
        for resid, res in chain.groups():
            max_resid = max(int(resid), max_resid)
    return max_resid
    

def reorder_ions_for_amber(protein, logger):
    assert isinstance(protein, pdfbridge.AtomGroup)
    
    protein = pdfbridge.AtomGroup(protein)

    # atomgroup for ions
    for model_name, model in protein.groups():
        max_resid = get_max_resid(model)
        current_resid = max_resid +1

        for chain_name, chain in model.groups():
            ions = get_ions(chain, logger)

            for ion in ions:
                res = pdfbridge.AtomGroup()
                res.name = ion.name
                res.set_atom(ion.symbol, ion)
                chain.set_group(current_resid, res)
                current_resid += 1

    return protein

def main():
    # parse args
    parser = argparse.ArgumentParser(description='reorder protein')
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

    # reorder
    mod_protein = reorder_ions_for_amber(protein, logger)

    # output
    with open(output_path, "wb") as mpac_file:
        raw_data = mod_protein.get_raw_data()
        mpac_file.write(msgpack.packb(raw_data))


if __name__ == '__main__':
    main()
