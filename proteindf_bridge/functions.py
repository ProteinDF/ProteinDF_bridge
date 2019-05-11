#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
logger = logging.getLogger(__name__)

try:
    import msgpack
except:
    import msgpack_pure as msgpack

from .atomgroup import AtomGroup

def load_atomgroup(brd_file_path):
    mpac_file = open(brd_file_path, "rb")
    mpac_data = msgpack.unpackb(mpac_file.read())
    mpac_file.close()
    atomgroup = AtomGroup(mpac_data)

    return atomgroup

def save_atomgroup(atomgroup, file_path):
    data = atomgroup.get_raw_data()
    mpac = msgpack.packb(data)

    with open(file_path, "wb") as f:
        f.write(mpac)

