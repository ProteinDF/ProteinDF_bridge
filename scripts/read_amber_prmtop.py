#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import proteindf_bridge as bridge


def main():
    # initialize
    parser = argparse.ArgumentParser(description='parse Amber prmtop')
    parser.add_argument('FILE',
                        nargs=1,
                        help='Amber prmtop file')
    args = parser.parse_args()

    path = args.FILE[0]

    ap = bridge.AmberPrmtop(path)
    charges = ap.charges
    for c in charges:
        print(c)


if __name__ == '__main__':
    main()
