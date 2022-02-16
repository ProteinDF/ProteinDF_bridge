#!/usr/bin/env python
# -*- coding: utf-8 -*-

from .functions import locate
from .atomgroup import AtomGroup
from .atom import Atom

import logging

logger = logging.getLogger(__name__)


class Format(object):
    @classmethod
    def is_residue(cls, res):
        assert isinstance(res, AtomGroup)
        answer = True

        if res.get_number_of_groups() > 0:
            logger.debug("not allowed groups in residue: name={}".format(res.name))
            answer = False

        if res.get_number_of_atoms() == 0:
            logger.debug("no atoms found in residue: name={}".format(res.name))

        return answer

    @classmethod
    def is_chain(cls, chain):
        assert isinstance(chain, AtomGroup)
        answer = True

        if chain.get_number_of_groups() > 0:
            for res_key, res in chain.groups():
                answer &= cls.is_residue(res)
        else:
            logger.debug("no groups found in chain: name={}".format(chain.name))

        if chain.get_number_of_atoms() != 0:
            logger.debug("not allowed any atoms in chain: name={}".format(chain.name))
            answer = False

        return answer

    @classmethod
    def is_protein(cls, model):
        assert isinstance(model, AtomGroup)
        answer = True

        if model.get_number_of_groups() > 0:
            for chain_key, chain in model.groups():
                answer &= cls.is_chain(chain)
        else:
            locate = locate()
            logger.debug(
                "no groups found in model: name={} at {}/{}/{}".format(
                    model.name, locate[0], locate[1], locate[2]
                )
            )

        if model.get_number_of_atoms() != 0:
            locate = locate()
            logger.critical(
                "not allowed existing any atoms in model: name={} at {}/{}/{}".format(
                    model.name, locate[0], locate[1], locate[2]
                )
            )
            answer = False

        return answer

    @classmethod
    def is_models(cls, models):
        assert isinstance(models, AtomGroup)
        answer = True

        if models.get_number_of_groups() > 0:
            for model_id, model in models.groups():
                answer &= cls.is_protein(model)
        else:
            # logger.warning("no groups found in models: name={}".format(models.name))
            pass

        if models.get_number_of_atoms() != 0:
            # logger.error("not allowed any atoms in models: name={}".format(models.name))
            answer = False
        return answer
