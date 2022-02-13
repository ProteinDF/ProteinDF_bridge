#!/usr/bin/env python
# -*- coding: utf-8 -*-

from .select import Select_Path
from .atomgroup import AtomGroup
from .atom import Atom
from .biopdb import Pdb

import logging
logger = logging.getLogger(__name__)


class SSBond(object):
    """ find disulfide bonds in protein models

    >>> tmp_pdb = Pdb('./data/1hls.pdb')
    >>> models = tmp_pdb.get_atomgroup()
    >>> model = models.get_group('model_1')
    >>> ssb = SSBond(model)
    >>> print(ssb.get_bonds())
    [('/model_1/A/6/', '/model_1/A/11/'), ('/model_1/A/7/', '/model_1/B/7/'), ('/model_1/A/20/', '/model_1/B/19/')]
    """
    # /model_1/A/6/ /model_1/A/11/ 2.0203361106508986
    # /model_1/A/7/ /model_1/B/7/ 2.0182586553759654
    # /model_1/A/20/ /model_1/B/19/ 2.0179301276307857

    _ss_bond_max_length = 2.1 * 1.1

    def __init__(self, model):
        self._model = AtomGroup(model)
        self._ssbonds = []
        self._isChecked = False

    def get_bonds(self):
        if self._isChecked == False:
            self._check()
        return self._ssbonds

    def _check(self):
        SGs = []
        for chain_id, chain in self._model.groups():
            for reskey, res in chain.groups():
                if res.name == 'CYS' or res.name == 'CYX':
                    logger.debug('found CYS. path={}'.format(res.path))
                    if res.has_atom('SG'):
                        SG = Atom(res.get_atom('SG'))
                        SGs.append((res.path, SG))
                    else:
                        logger.warn('not found SG atom in {}'.format(res.path))
        self._check_SGs(SGs)
        self._isChecked = True

    def _check_SGs(self, SGs):
        assert(isinstance(SGs, list))

        num_of_SGs = len(SGs)
        for i in range(num_of_SGs):
            (path1, SG1) = SGs[i]
            for j in range(i + 1, num_of_SGs):
                (path2, SG2) = SGs[j]
                distance = SG1.xyz.distance_from(SG2.xyz)
                if distance < self._ss_bond_max_length:
                    logger.debug("found SS bond: {} <-> {}; {}".format(path1, path2, distance))
                    self._ssbonds.append((path1, path2))


if __name__ == "__main__":
    from .biopdb import Pdb

    import doctest
    doctest.testmod()
