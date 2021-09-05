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

from .atomgroup import AtomGroup
from .format import Format

import logging
logger = logging.getLogger(__name__)


class Utils(object):
    @classmethod
    def remove_WAT(cls, atomgroup):
        """ remove water(WAT or HOH) residues
        """
        assert(isinstance(atomgroup, AtomGroup))

        answer = AtomGroup(atomgroup)
        wat_keys = ["HOH", "WAT"]
        remove_groups = []
        for key, grp in answer.groups():
            grp_name = grp.name
            logger.debug("check group name: {}".format(grp_name))
            if grp_name in wat_keys:
                logger.debug("remove name: {}".format(grp_name))
                remove_groups.append(key)
                continue

            answer.set_group(key, cls.remove_WAT(grp))

        for key in remove_groups:
            answer.remove_group(key)

        return answer

    @classmethod
    def get_sequential_residue_id(cls, protein, req_chain_id, req_res_id):
        """Returns the amino acid residue number of the entire sequence.

        Args:
            protein (AtomGroup): AtomGroup representing a protein
            req_chain_id (str): chain id
            req_res_id (str): residue id

        Returns:
            int: Return serial number of the amino acid residue. If it does not apply, 0 is returned.
        """
        assert(Format.is_protein(protein))
        req_chain_id = str(req_chain_id)
        req_res_id = str(req_res_id)

        counter = 0
        for chain_id, chain in protein.groups():
            if req_chain_id == chain_id:
                for res_id, res in chain.groups():
                    counter += 1
                    if req_res_id == res_id:
                        return counter
            else:
                num_of_groups = chain.get_num_of_groups()
                counter += num_of_groups

        return 0


if __name__ == "__main__":
    import doctest
    doctest.testmod()
