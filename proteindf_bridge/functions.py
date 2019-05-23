#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
logger = logging.getLogger(__name__)

try:
    import msgpack
except:
    import msgpack_pure as msgpack

from .atomgroup import AtomGroup


def load_msgpack(mpac_path):
    """load message pack file
    """
    mpac_data = None
    with open(mpac_path, "rb") as f:
        mpac_data = msgpack.unpackb(f.read())

    return mpac_data

def save_msgpack(mpac_path, data):
    mpac_data = msgpack.packb(data)
    with open(mpac_path, "wb") as f:
        f.write(mpac_data)


def load_atomgroup(brd_path):
    """load bridge file (msgpack format)
    """
    mpac_data = load_msgpack(brd_path)
    atomgroup = AtomGroup(mpac_data)

    return atomgroup


def save_atomgroup(atomgroup, file_path):
    """save bridge file (msgpack format)
    """
    data = atomgroup.get_raw_data()
    save_msgpack(file_path, data)
