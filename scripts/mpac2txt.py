#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2014 The ProteinDF development team.
# see also AUTHORS and README if provided.
#
# This file is a part of the ProteinDF software package.
#
# The ProteinDF is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The ProteinDF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

import sys
import argparse
import pprint

import proteindf_bridge as bridge


def main():
    # initialize
    parser = argparse.ArgumentParser(
        description='display file formatted by MsgPack using YAML.')
    parser.add_argument('FILE',
                        nargs=1,
                        help='Amber PDB file')
    args = parser.parse_args()

    file_path = args.FILE[0]

    data = bridge.load_msgpack(file_path)
    pprint.pprint(data)


if __name__ == '__main__':
    main()
