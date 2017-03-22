#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 The ProteinDF development team.
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

import pdfbridge
import logging
logger = logging.getLogger(__name__)

class IonPair(object):
    def __init__(self, model):
        self._model = pdfbridge.AtomGroup(model)

    def get_ion_pairs(self):
        ion_pairs = []
        (anion_list, cation_list) = self._get_ion_list()
        for anion_path, anion_pos_array in anion_list.items():
            for (anion_pos, anion_type) in anion_pos_array:
                
                for cation_path, cation_pos_array in cation_list.items():
                    for (cation_pos,cation_type) in cation_pos_array:

                        d = anion_pos.distance_from(cation_pos)
                        if d < 4.0:
                            logger.info("found ion pair: {}".format(anion_path))
                            logger.info("              : {}".format(cation_path))
                            ion_pairs.append((anion_path, cation_path, anion_type, cation_type))

        return ion_pairs
                    
    def _get_ion_list(self):
        anion_list = {}
        cation_list = {}

        for chain_key, chain in self._model.groups():
            for res_key, res in chain.groups():
                name = res.name

                if name == 'GLU':
                    anion_list[res.path] = [(self._get_center_GLU(res), 'GLU')]
                elif name == 'ASP':
                    anion_list[res.path] = [(self._get_center_ASP(res), 'ASP')]
                elif name == 'LYS':
                    cation_list[res.path] = [(self._get_center_LYS(res), 'LYS')]
                elif name == 'ARG':
                    cation_list[res.path] = [(self._get_center_ARG(res, 0), 'ARG'),
                                             (self._get_center_ARG(res, 1), 'ARG1'),
                                             (self._get_center_ARG(res, 2), 'ARG2')]

                if pdfbridge.AminoAcid.is_aminoacid(res):
                    if res.has_atom('H3'):
                        cation_list.setdefault(res.path, [])
                        cation_list[res.path].append((self._get_center_Nterm(res), 'NTM'))
                    if res.has_atom('OXT'):
                        anion_list.setdefault(res.path, [])
                        anion_list[res.path].append((self._get_center_Cterm(res), 'CTM'))

        return (anion_list, cation_list)
                        
        
    def _get_center_Nterm(self, res):
        """
        N末端のイオン対判定用座標を返す
        """
        return res['N'].xyz

        
    def _get_center_Cterm(self, res):
        """
        C末端のイオン対判定用座標を返す
        """
        ag = pdfbridge.AtomGroup()
        ag.set_atom('C', res['C'])
        ag.set_atom('O1', res['O'])
        ag.set_atom('O2', res['OXT'])
        return ag.center()
        

    def _get_center_GLU(self, res):
        """
        GLUのイオン対判定用座標を返す
        """
        ag = pdfbridge.AtomGroup()
        ag.set_atom('C', res['CD'])
        ag.set_atom('O1', res['OE1'])
        ag.set_atom('O2', res['OE2'])
        return ag.center()

        
    def _get_center_ASP(self, res):
        """
        ASPのイオン対判定用座標を返す
        """
        ag = pdfbridge.AtomGroup()
        ag.set_atom('C', res['CG'])
        ag.set_atom('O1', res['OD1'])
        ag.set_atom('O2', res['OD2'])
        return ag.center()
        

    def _get_center_LYS(self, res):
        """
        LYSのイオン対判定用座標を返す
        """
        return res['NZ'].xyz


    def _get_center_ARG(self, res, case=0):
        """
        case: 0; 中央
        case: 1; NH1側
        case: 2; NH2側
        """
        case = int(case)

        answer = pdfbridge.Position()
        ag = pdfbridge.AtomGroup()
        if case == 0:
            ag.set_atom('NH1', res['NH1'])
            ag.set_atom('NH2', res['NH2'])
            ag.set_atom('CZ', res['CZ'])
            answer = ag.center()
        elif case == 1:
            answer = res['NH1'].xyz
        elif case == 2:
            answer = res['NH2'].xyz
        else:
            logger.warning("unknown mode={}".format(mode))
            
        return answer


if __name__ == "__main__":
    main()


