#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import proteindf_bridge as bridge


def atomgroup2pdf(atomgroup):
    '''convert atomgroup data to ProteinDF coordinate format
    '''
    # print(atomgroup)

    atoms = []
    for atom in atomgroup.get_atom_list():
        item = {}
        item['symbol'] = atom.symbol
        p = atom.xyz / 0.529177249
        item['xyz'] = [p.x, p.y, p.z]
        item['charge'] = atom.charge
        item['label'] = ''
        atoms.append(item)
    # print(atoms)

    data = {}
    data.setdefault('coordinates', {})
    data['coordinates']['atoms'] = atoms

    return data


def main():
    # initialize
    parser = argparse.ArgumentParser(description='parse Amber prmtop')
    parser.add_argument('prmtop_file',
                        nargs=1,
                        help='Amber prmtop file')
    parser.add_argument('inpcrd_file',
                        nargs=1,
                        help='Amber inpcrd file')
    args = parser.parse_args()

    verbose = True
    prmtop_path = args.prmtop_file[0]
    inpcrd_path = args.inpcrd_file[0]
    output_path = "output.brd"

    ap = bridge.AmberPrmtop(prmtop_path, inpcrd_path)
    # charges = ap.charges
    # for c in charges:
    #     print(c)
    atomgroup = ap.get_atomgroup()
    # print(atomgroup)

    raw_data = None
    # raw_data = atomgroup.get_raw_data()
    raw_data = atomgroup2pdf(atomgroup)

    # output file
    if (verbose == True):
        print("writing: %s\n" % (output_path))
    bridge.save_msgpack(raw_data, output_path)


if __name__ == '__main__':
    main()
