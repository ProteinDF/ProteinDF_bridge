#!/usr/bin/env python
# -*- coding: utf-8 -*-

import inspect
import os
import yaml
try:
    import msgpack
except:
    import msgpack_pure as msgpack

from .atomgroup import AtomGroup
from .str_processing import StrUtils


import logging
logger = logging.getLogger(__name__)


def locate():
    """
    return tuple of (file, function, line number)
    cf.) https://qiita.com/ymko/items/b46d32b98f013f06d805
    """
    frame = inspect.currentframe().f_back
    return os.path.basename(frame.f_code.co_filename), frame.f_code.co_name, frame.f_lineno


def load_yaml(yaml_path):
    assert(isinstance(yaml_path, str))

    data = None
    with open(yaml_path) as f:
        contents = f.read()
        # contents = contents.decode('utf8')  # for Japanese
        data = parse_yaml(contents)

    return data


def parse_yaml(yaml_data):
    assert(isinstance(yaml_data, str))
    data = list(yaml.load_all(yaml_data, Loader=yaml.SafeLoader))

    return data


def get_yaml(data):
    yaml_str = ""
    if data is not None:
        yaml_str = yaml.dump(data,
                             encoding='utf8',
                             allow_unicode=True,
                             default_flow_style=False,
                             line_break='\n')
        yaml_str = StrUtils.to_unicode(yaml_str)

    return yaml_str


def save_yaml(data, yaml_path):
    assert(isinstance(yaml_path, str))

    yaml_str = get_yaml(data)
    with open(yaml_path, "wb") as f:
        f.write(yaml_str)


def load_msgpack(mpac_path):
    """load message pack file and return python data
    """
    assert(isinstance(mpac_path, str))

    mpac_data = None
    with open(mpac_path, "rb") as f:
        mpac_data = msgpack.unpackb(f.read(), strict_map_key=False)

        if isinstance(mpac_data, list):
            mpac_data = StrUtils.to_unicode_list(mpac_data)
        elif isinstance(mpac_data, dict):
            mpac_data = StrUtils.to_unicode_dict(mpac_data)

    return mpac_data


def save_msgpack(data, mpac_path):
    """save python data to msgpack file
    """
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
