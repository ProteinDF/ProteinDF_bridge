#!/usr/bin/env python
# -*- coding: utf-8 -*-

from .atomgroup import AtomGroup
from .utils import Utils

try:
    import msgpack
except:
    import msgpack_pure as msgpack

import logging
logger = logging.getLogger(__name__)


def load_msgpack(mpac_path):
    """load message pack file
    """
    assert(isinstance(mpac_path, str))

    mpac_data = None
    with open(mpac_path, "rb") as f:
        mpac_data = msgpack.unpackb(f.read(), strict_map_key=False)

        if isinstance(mpac_data, list):
            mpac_data = Utils.to_unicode_list(mpac_data)
        elif isinstance(mpac_data, dict):
            mpac_data = Utils.to_unicode_dict(mpac_data)

    return mpac_data


def save_msgpack(data, mpac_path):
    assert(isinstance(mpac_path, str))

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
    save_msgpack(data, file_path)
