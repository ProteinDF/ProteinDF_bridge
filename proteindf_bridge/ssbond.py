#!/usr/bin/env python
# -*- coding: utf-8 -*-

from .select import Select_Path
from .atomgroup import AtomGroup
from .atom import Atom
import logging
logger = logging.getLogger(__name__)


class SSBond(object):
    """ check disulfide bond in protein models

    >>> tmp_pdb = Pdb('./data/1hls.pdb')
    >>> models = tmp_pdb.get_atomgroup()
    >>> ssbond = SSBond()
    >>> ssbond.check(models)

    """
    _ss_bond_max_length = 2.1 * 1.1

    def __init__(self):
        self._ssbonds = []

    def check(self, protein_models):
        assert(isinstance(protein_models, AtomGroup))

        for model_id, model in protein_models.groups():
            cysteines = []
            for chain_id, chain in model.groups():
                for reskey, res in chain.groups():
                    if res.name == 'CYS' or res.name == 'CYX':
                        logger.debug('found CYS. path={}'.format(res.path))
                        if res.has_atom('SG'):
                            SG = Atom(res.get_atom('SG'))
                            cysteines.append((res.path, SG))
                        else:
                            logger.warn(
                                'not found SG atom in {}'.format(res.path))
            self._check_SSbonds(model, cysteines)

    def _check_SSbonds(self, model, cysteines):
        assert(isinstance(model, AtomGroup))
        assert(isinstance(cysteines, list))

        num_of_cysteins = len(cysteines)
        for i in range(num_of_cysteins):
            (path1, SG1) = cysteines[i]
            for j in range(i + 1, num_of_cysteins):
                (path2, SG2) = cysteines[j]
                distance = SG1.xyz.distance_from(SG2.xyz)
                if distance < self._ss_bond_max_length:
                    print(path1, path2, distance)
                    logger.info('found disulfide bond: {path1} - {path2}'.format(path1=path1,
                                                                                 path2=path2))
                    self._ssbonds.append((path1, path2))


if __name__ == "__main__":
    from .biopdb import Pdb

    import doctest
    doctest.testmod()
