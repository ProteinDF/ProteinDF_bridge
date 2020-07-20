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

from .superposer import Superposer
from .matrix import Matrix
from .atomgroup import AtomGroup
from .atom import Atom
from .functions import load_msgpack
from .position import Position
from .error import BrInputError
# from .xyz import Xyz

import os
import math
import re

import logging
logger = logging.getLogger(__name__)


class Modeling:
    _ACE_ALA_NME_path_base = os.path.join(
        os.environ.get('PDF_HOME', '.'),
        'data',
        "ACE_ALA_NME_{}.brd")
    _ACE_ALA_NME_comformers = ["trans1", "trans2", "cis1", "cis2"]

    def __init__(self):
        self._ACE_ALA_NME = {}
        for comformer in self._ACE_ALA_NME_comformers:
            brd_path = self._ACE_ALA_NME_path_base.format(comformer)
            # print(comformer, brd_path)
            atomgroup = AtomGroup(load_msgpack(brd_path))
            assert(atomgroup.get_number_of_all_atoms() > 0)
            self._ACE_ALA_NME[comformer] = atomgroup

    def _get_ACE_ALA_NME(self, comformer):
        assert(comformer in self._ACE_ALA_NME_comformers)
        return self._ACE_ALA_NME[comformer]

    # -----------------------------------------------------------------
    def get_ACE_simple(self, next_aa):
        """
        隣のC-alphaの位置をメチル基にする。
        """
        answer = AtomGroup()

        CAs = next_aa.pickup_atoms('CA')
        if len(CAs) > 0:
            answer.set_atom('CA', CAs[0])
        else:
            raise BrInputError(next_aa,
                               'cannot found "CA" atom on building ACE.')

        Cs = next_aa.pickup_atoms('C')
        if len(Cs) > 0:
            answer.set_atom('C', Cs[0])
        else:
            raise BrInputError(next_aa,
                               'cannot found "C" atom on building ACE.')

        Os = next_aa.pickup_atoms('O')
        if len(Os) > 0:
            answer.set_atom('O', Os[0])
        else:
            raise BrInputError(next_aa,
                               'cannot found "O" atom on building ACE.')

        answer |= self.add_methyl(answer['CA'], answer['C'])
        answer.path = '/ACE'
        return answer

    def get_NME_simple(self, next_aa):
        """
        隣のC-alphaの位置をメチル基にする。
        """
        answer = AtomGroup()

        CAs = next_aa.pickup_atoms('CA')
        if len(CAs) > 0:
            answer.set_atom('CA', CAs[0])
        else:
            raise BrInputError(next_aa,
                               'cannot found "CA" atom on building NME.')

        Ns = next_aa.pickup_atoms('N')
        if len(Ns) > 0:
            answer.set_atom('N', Ns[0])
        else:
            raise BrInputError(next_aa,
                               'cannot found "N" atom on building NME.')

        Hs = next_aa.pickup_atoms('H')
        if len(Hs) > 0:
            answer.set_atom('H', Hs[0])
        else:
            # for proline
            CDs = next_aa.pickup_atoms('CD')
            if len(CDs) > 0:
                dummy_H = Atom(CDs[0])
                dummy_H.symbol = 'H'
                answer.set_atom('H', dummy_H)
            else:
                raise BrInputError(next_aa,
                                   'cannot found "H" or "CD" atom(for proline) on building NME.')

        answer |= self.add_methyl(answer['CA'], answer['N'])
        answer.path = '/NME'
        return answer

    # -----------------------------------------------------------------
    def get_ACE(self, res, next_aa=None):
        """
        template (ACE-ALA-NME) format:
        HH3[1-3]-CH3-C -  N-CA(HA)-C-    N-CH3-HH3[1-3]
                     ||   | |      ||    |
                     O    H CB     O     H
        """
        AAN = None
        rmsd_min = 1000.0
        for comformer in self._ACE_ALA_NME_comformers:
            ref_AAN = self._get_ACE_ALA_NME(comformer)
            (matched, rmsd) = self._match_ACE(ref_AAN, res, next_aa)
            # print(comformer, rmsd)
            if rmsd < rmsd_min:
                rmsd_min = rmsd
                AAN = matched

        if rmsd_min > 1.0:
            logger.warn("RMSD value is too large: {}".format(rmsd))

        answer = AtomGroup(AAN['1'])
        answer.path = '/ACE'

        return answer

    def _match_ACE(self, AAN, res, next_aa):
        '''AAN (ACE-ALA-NME)
        '''
        assert(isinstance(AAN, AtomGroup))
        assert(isinstance(res, AtomGroup))
        (AAN_part, res_part) = self._match_residues(AAN['2'], res)

        # for ACE
        if next_aa is not None:
            if next_aa.has_atom('N'):
                AAN_part.set_atom('N2', AAN['3']['N'])
                res_part.set_atom('N2', next_aa['N'])
            if next_aa.has_atom('H'):
                AAN_part.set_atom('NH2', AAN['3']['H'])
                res_part.set_atom('NH2', next_aa['H'])
            if next_aa.has_atom('CA'):
                AAN_part.set_atom('CH3', AAN['3']['CH3'])
                res_part.set_atom('CH3', next_aa['CA'])

        sp = Superposer(AAN_part, res_part)
        rmsd = sp.rmsd

        matched_AAN = sp.superimpose(AAN)

        return (matched_AAN, rmsd)

    def get_NME(self, res, next_aa=None):
        """
        template (ACE-ALA-NME) format:
        HH3[1-3]-CH3-C -  N-CA(HA)-C-    N-CH3-HH3[1-3]
                     ||   | |      ||    |
                     O    H CB     O     H
        """
        AAN = None
        rmsd_min = 1000.0
        for comformer in self._ACE_ALA_NME_comformers:
            ref_AAN = self._get_ACE_ALA_NME(comformer)
            (matched, rmsd) = self._match_NME(ref_AAN, res, next_aa)
            # print(comformer, rmsd)
            if rmsd < rmsd_min:
                rmsd_min = rmsd
                AAN = matched

        if rmsd_min > 1.0:
            logger.warn("RMSD value is too large: {}".format(rmsd))

        answer = AtomGroup(AAN['3'])
        answer.path = '/NME'

        return answer

    def _match_NME(self, AAN, res, next_aa):
        '''AAN (ACE-ALA-NME)
        '''
        assert(isinstance(AAN, AtomGroup))
        assert(isinstance(res, AtomGroup))
        (AAN_part, res_part) = self._match_residues(AAN['2'], res)

        # for NME
        if next_aa is not None:
            if next_aa.has_atom('C'):
                AAN_part.set_atom('C2', AAN['1']['C'])
                res_part.set_atom('C2', next_aa['C'])
            if next_aa.has_atom('O'):
                AAN_part.set_atom('O2', AAN['1']['O'])
                res_part.set_atom('O2', next_aa['O'])
            if next_aa.has_atom('CA'):
                AAN_part.set_atom('CH3', AAN['1']['CH3'])
                res_part.set_atom('CH3', next_aa['CA'])

        sp = Superposer(AAN_part, res_part)
        rmsd = sp.rmsd

        matched_AAN = sp.superimpose(AAN)

        return (matched_AAN, rmsd)

    def _match_residues(self, res1, res2, max_number_of_atoms=-1):
        """
        2つのアミノ酸残基のN, H, CA, HA, C, Oの原子を突き合わせる。
        アミノ酸残基がプロリンだった場合は、CDの炭素をHに命名する。
        GLYはHA1, HA2とあるので突き合せない。
        """
        atom_names = ['CA', 'O', 'C', 'N', 'CB', 'HA']
        if max_number_of_atoms == -1:
            max_number_of_atoms = len(atom_names)
        ans_res1 = AtomGroup()
        ans_res2 = AtomGroup()

        for atom_name in atom_names:
            pickup_atoms1 = res1.pickup_atoms(atom_name)
            if len(pickup_atoms1) > 0:
                pickup_atoms2 = res2.pickup_atoms(atom_name)
                if len(pickup_atoms2) > 0:
                    ans_res1.set_atom(atom_name, pickup_atoms1[0])
                    ans_res2.set_atom(atom_name, pickup_atoms2[0])

            if ans_res1.get_number_of_atoms() >= max_number_of_atoms:
                break

        # match amino-'H'
        if ans_res1.get_number_of_atoms() < max_number_of_atoms:
            res1_H = None
            res2_H = None
            if res1.has_atom('H'):
                res1_H = res1['H']
            elif res1.has_atom('CD'):
                # for proline
                res1_H = res1['CD']
            if res2.has_atom('H'):
                res2_H = res2['H']
            elif res2.has_atom('CD'):
                res2_H = res2['CD']
            if ((res1_H is not None) and (res2_H is not None)):
                ans_res1.set_atom('H', res1_H)
                ans_res2.set_atom('H', res2_H)

        return (ans_res1, ans_res2)

    # -----------------------------------------------------------------
    def add_methyl(self, C1, C2):
        """
        -CH3の水素を付加

        C1に水素を付加
        """
        assert(isinstance(C1, Atom))
        assert(isinstance(C2, Atom))

        ethane = AtomGroup()
        ethane.set_atom('C1', Atom(symbol='C', name='C1',
                                   position=Position(0.00000, 0.00000, 0.00000)))
        ethane.set_atom('H11', Atom(symbol='H', name='H11',
                                    position=Position(-0.85617, -0.58901, -0.35051)))
        ethane.set_atom('H12', Atom(symbol='H', name='H12',
                                    position=Position(-0.08202, 1.03597, -0.35051)))
        ethane.set_atom('H13', Atom(symbol='H', name='H13',
                                    position=Position(0.93818, -0.44696, -0.35051)))
        ethane.set_atom('C2', Atom(symbol='C', name='C2',
                                   position=Position(0.00000, 0.00000, 1.47685)))
        ethane.set_atom('H21', Atom(symbol='H', name='H21',
                                    position=Position(-0.93818, 0.44696, 1.82736)))
        ethane.set_atom('H22', Atom(symbol='H', name='H22',
                                    position=Position(0.85617, 0.58901, 1.82736)))
        ethane.set_atom('H23', Atom(symbol='H', name='H23',
                                    position=Position(0.08202, -1.03597, 1.82736)))

        inC21 = C2.xyz - C1.xyz
        refC21 = ethane['C2'].xyz - ethane['C1'].xyz

        shift = C1.xyz - ethane['C1'].xyz
        rot = self.arbitary_rotate_matrix(inC21, refC21)

        ethane.rotate(rot)
        ethane.shift_by(shift)
        assert(C1.xyz == ethane['C1'].xyz)

        answer = AtomGroup()
        answer.set_atom('H11', ethane['H11'])
        answer.set_atom('H12', ethane['H12'])
        answer.set_atom('H13', ethane['H13'])

        return answer

    # -----------------------------------------------------------------
    def get_NH3(self, angle=0.5 * math.pi, length=1.0):
        pi23 = math.pi * 2.0 / 3.0  # (pi * 2/3)
        sin23 = math.sin(pi23)
        cos23 = math.cos(pi23)
        # pi43 = math.pi * 4.0 / 3.0  # (pi * 4/3)
        # sin43 = math.sin(pi43)
        # cos43 = math.cos(pi43)
        sin_input = math.sin(angle)
        cos_input = math.cos(angle)

        # z軸まわりに120度回転
        # z1_rot = Matrix(3, 3)
        # z1_rot.set(0, 0,  cos23)
        # z1_rot.set(0, 1, -sin23)
        # z1_rot.set(1, 0,  sin23)
        # z1_rot.set(1, 1,  cos23)
        # z1_rot.set(2, 2,  1.0)
        # z軸まわりに240度回転
        # z2_rot = Matrix(3, 3)
        # z2_rot.set(0, 0,  cos43)
        # z2_rot.set(0, 1, -sin43)
        # z2_rot.set(1, 0,  sin43)
        # z2_rot.set(1, 1,  cos43)
        # z2_rot.set(2, 2,  1.0)
        # y軸まわりに回転
        # y_rot = Matrix(3, 3)
        # y_rot.set(0, 0,  cos_input)
        # y_rot.set(0, 2, -sin_input)
        # y_rot.set(2, 0,  sin_input)
        # y_rot.set(2, 2,  cos_input)
        # y_rot.set(1, 1,  1.0)

        # pos_H1 = Position(1.0, 0.0, 0.0)
        # pos_H1.rotate(y_rot)
        # pos_H1 *= length
        # pos_H2 = Position(1.0, 0.0, 0.0)
        # pos_H2.rotate(y_rot)
        # pos_H2.rotate(z1_rot)
        # pos_H2 *= length
        # pos_H3 = Position(1.0, 0.0, 0.0)
        # pos_H3.rotate(y_rot)
        # pos_H3.rotate(z2_rot)
        # pos_H3 *= length

        # X-Z平面上、Y軸に対してangle度開く
        xz_rot = Matrix(3, 3)
        xz_rot.set(0, 0, cos_input)
        xz_rot.set(0, 2, -sin_input)
        xz_rot.set(2, 0, sin_input)
        xz_rot.set(2, 2, cos_input)
        xz_rot.set(1, 1, 1.0)

        # X-Y平面上、Z軸に対して120度開く
        xy_rot = Matrix(3, 3)
        xy_rot.set(0, 0, cos23)
        xy_rot.set(0, 1, -sin23)
        xy_rot.set(1, 0, sin23)
        xy_rot.set(1, 1, cos23)
        xy_rot.set(2, 2, 1.0)

        pos_H1 = Position(0.0, 0.0, 1.0)
        pos_H1.rotate(xz_rot)

        pos_H2 = Position(0.0, 0.0, 1.0)
        pos_H2.rotate(xz_rot)
        pos_H2.rotate(xy_rot)

        pos_H3 = Position(0.0, 0.0, 1.0)
        pos_H3.rotate(xz_rot)
        pos_H3.rotate(xy_rot)
        pos_H3.rotate(xy_rot)

        pos_H1 *= length
        pos_H2 *= length
        pos_H3 *= length

        NH3 = AtomGroup()
        N = Atom(symbol='N',
                 position=Position(0.0, 0.0, 0.0))
        H1 = Atom(symbol='H',
                  position=pos_H1)
        H2 = Atom(symbol='H',
                  position=pos_H2)
        H3 = Atom(symbol='H',
                  position=pos_H3)
        # X1 = Atom(symbol = 'X',
        #                position = Position(1.0, 0.0, 0.0))
        # X2 = Atom(symbol = 'X',
        #                position = Position(0.0, 1.0, 0.0))
        # X3 = Atom(symbol = 'X',
        #                position = Position(0.0, 0.0, 1.0))

        NH3.set_atom('N', N)
        NH3.set_atom('H1', H1)
        NH3.set_atom('H2', H2)
        NH3.set_atom('H3', H3)
        # NH3.set_atom('X1', X1)
        # NH3.set_atom('X2', X2)
        # NH3.set_atom('X3', X3)

        return NH3

    # -----------------------------------------------------------------
    def select_residues(self, chain, from_resid, to_resid):
        '''
        連続したアミノ酸残基を返す
        '''
        answer = AtomGroup()
        for resid, res in chain.groups():
            resid = int(resid)
            if from_resid <= resid <= to_resid:
                answer |= res

        return answer

    # -----------------------------------------------------------------
    def arbitary_rotate_matrix(self, in_a, in_b):
        """
        ベクトルaをbへ合わせる回転行列(3x3)を返す
        """
        assert(isinstance(in_a, Position))
        assert(isinstance(in_b, Position))

        a = Position(in_a)
        b = Position(in_b)
        a.norm()
        b.norm()

        cos_theta = a.dot(b)
        sin_theta = math.sqrt(1 - cos_theta * cos_theta)

        n = a.cross(b)
        n.norm()

        nx = n.x
        ny = n.y
        nz = n.z

        rot = Matrix(3, 3)
        rot.set(0, 0, nx * nx * (1.0 - cos_theta) + cos_theta)
        rot.set(0, 1, nx * ny * (1.0 - cos_theta) + nz * sin_theta)
        rot.set(0, 2, nx * nz * (1.0 - cos_theta) - ny * sin_theta)
        rot.set(1, 0, nx * ny * (1.0 - cos_theta) - nz * sin_theta)
        rot.set(1, 1, ny * ny * (1.0 - cos_theta) + cos_theta)
        rot.set(1, 2, nx * nz * (1.0 - cos_theta) + nx * sin_theta)
        rot.set(2, 0, nx * nz * (1.0 - cos_theta) + ny * sin_theta)
        rot.set(2, 1, ny * nz * (1.0 - cos_theta) - nx * sin_theta)
        rot.set(2, 2, nz * nz * (1.0 - cos_theta) + cos_theta)

        return rot
    # -----------------------------------------------------------------

    def get_last_index(self, res):
        answer = 0
        re_obj = re.compile('([0-9]+)')
        for key, atom in res.atoms():
            m = re_obj.search(key)
            if m is not None:
                num = m.group(0)
                num = int(num)
                answer = max(num, answer)
        return answer

    # -----------------------------------------------------------------
    def neutralize_Nterm(self, res):
        answer = None
        if res.name == "PRO":
            answer = self._neutralize_Nterm_PRO(res)
        else:
            answer = self._neutralize_Nterm(res)

        return answer

    def _neutralize_Nterm(self, res):
        """
        N末端側を中性化するためにCl-(AtomGroup)を返す

        H1, N2, HXT(or H3)が指定されている必要があります。
        """
        ag = AtomGroup()
        ag.set_atom('N', res['N'])
        ag.set_atom('H1', res['H1'])
        ag.set_atom('H2', res['H2'])
        if res.has_atom('HXT'):
            ag.set_atom('H3', res['HXT'])
        elif res.has_atom('H3'):
            ag.set_atom('H3', res['H3'])
        pos = self._get_neutralize_pos_NH3_type(ag)

        answer = AtomGroup()
        Cl = Atom(symbol='Cl',
                  name='Cl',
                  position=pos)
        answer.set_atom('Cl', Cl)
        return answer

    def _neutralize_Nterm_PRO(self, res):
        """in case of 'PRO', neutralize N-term
        """
        ag = AtomGroup()
        ag.set_atom('N', res['N'])
        ag.set_atom('H2', res['H2'])
        if res.has_atom('HXT'):
            ag.set_atom('H1', res['HXT'])
        elif res.has_atom('H3'):
            ag.set_atom('H1', res['H3'])
        pos = self._get_neutralize_pos_NH2_type(ag)

        answer = AtomGroup()
        Cl = Atom(symbol='Cl',
                  name='Cl',
                  position=pos)
        answer.set_atom('Cl', Cl)
        return answer

    def neutralize_Cterm(self, res):
        """
        C末端側を中性化するためにNa+(AtomGroup)を返す
        """
        ag = AtomGroup()
        ag.set_atom('C', res['C'])
        ag.set_atom('O1', res['O'])
        ag.set_atom('O2', res['OXT'])
        pos = self._get_neutralize_pos_COO_type(ag)

        answer = AtomGroup()
        Na = Atom(symbol='Na',
                  name='Na',
                  position=pos)
        answer.set_atom('Na', Na)
        return answer

    # -----------------------------------------------------------------
    def neutralize_GLU(self, res):
        ag = AtomGroup()
        ag.set_atom('C', res['CD'])
        ag.set_atom('O1', res['OE1'])
        ag.set_atom('O2', res['OE2'])
        pos = self._get_neutralize_pos_COO_type(ag)

        answer = AtomGroup()
        Na = Atom(symbol='Na',
                  name='Na',
                  position=pos)
        key = self.get_last_index(res)
        answer.set_atom('{}_Na'.format(key + 1), Na)
        return answer

    def neutralize_ASP(self, res):
        ag = AtomGroup()
        ag.set_atom('C', res['CG'])
        ag.set_atom('O1', res['OD1'])
        ag.set_atom('O2', res['OD2'])
        pos = self._get_neutralize_pos_COO_type(ag)

        answer = AtomGroup()
        Na = Atom(symbol='Na',
                  name='Na',
                  position=pos)
        key = self.get_last_index(res)
        answer.set_atom('{}_Na'.format(key + 1), Na)
        return answer

    def neutralize_LYS(self, res):
        ag = AtomGroup()
        ag.set_atom('N', res['NZ'])
        ag.set_atom('H1', res['HZ1'])
        ag.set_atom('H2', res['HZ2'])
        ag.set_atom('H3', res['HZ3'])
        pos = self._get_neutralize_pos_NH3_type(ag)

        answer = AtomGroup()
        Cl = Atom(symbol='Cl',
                  name='Cl',
                  position=pos)
        key = self.get_last_index(res)
        answer.set_atom('{}_Cl'.format(key + 1), Cl)
        return answer

    def neutralize_ARG(self, res, case=0):
        """
        case: 0; 中央
        case: 1; NH1側
        case: 2; NH2側
        """
        case = int(case)
        pos = Position()
        if case == 0:
            length = 3.0
            NH1 = res['NH1']
            NH2 = res['NH2']
            CZ = res['CZ']
            M = Position(0.5 * (NH1.xyz.x + NH2.xyz.x),
                         0.5 * (NH1.xyz.y + NH2.xyz.y),
                         0.5 * (NH1.xyz.z + NH2.xyz.z))
            vCM = M - CZ.xyz
            vCM.norm()
            pos = CZ.xyz + length * vCM
        elif case == 1:
            length = 2.0
            HH11 = res['HH11']
            HH12 = res['HH12']
            N = res['NH1']
            M = Position(0.5 * (HH11.xyz.x + HH12.xyz.x),
                         0.5 * (HH11.xyz.y + HH12.xyz.y),
                         0.5 * (HH11.xyz.z + HH12.xyz.z))
            vNM = M - N.xyz
            vNM.norm()
            pos = N.xyz + length * vNM
        elif case == 2:
            length = 2.0
            HH21 = res['HH21']
            HH22 = res['HH22']
            N = res['NH2']
            M = Position(0.5 * (HH21.xyz.x + HH22.xyz.x),
                         0.5 * (HH21.xyz.y + HH22.xyz.y),
                         0.5 * (HH21.xyz.z + HH22.xyz.z))
            vNM = M - N.xyz
            vNM.norm()
            pos = N.xyz + length * vNM
        else:
            pass

        answer = AtomGroup()
        Cl = Atom(symbol='Cl',
                  name='Cl',
                  position=pos)
        key = self.get_last_index(res)
        answer.set_atom('{}_Cl'.format(key + 1), Cl)
        return answer

    # ------------------------------------------------------------------
    def neutralize_FAD(self, ag):
        print("neutralize_FAD")
        print(ag)
        answer = AtomGroup()

        POO1 = AtomGroup()
        POO1.set_atom('P', ag['P'])

        # amber format: OP1, pdb: O1P
        if ag.has_atom('O1P'):
            POO1.set_atom('O1', ag['O1P'])
        elif ag.has_atom('OP1'):
            POO1.set_atom('O1', ag['OP1'])
        else:
            raise

        # amber format: OP2, pdb: O2P
        if ag.has_atom('O2P'):
            POO1.set_atom('O2', ag['O2P'])
        elif ag.has_atom('OP2'):
            POO1.set_atom('O2', ag['OP2'])
        else:
            raise

        Na1 = Atom(symbol='Na',
                   name='Na',
                   position=self._get_neutralize_pos_POO_type(POO1))

        POO2 = AtomGroup()
        POO2.set_atom('P', ag['PA'])
        POO2.set_atom('O1', ag['O1A'])  # amber format: OA1, pdb: O1A
        POO2.set_atom('O2', ag['O2A'])  # amber format: OA2, pdb: O2A
        Na2 = Atom(symbol='Na',
                   name='Na',
                   position=self._get_neutralize_pos_POO_type(POO2))

        key = self.get_last_index(ag)
        answer.set_atom('{}_Na1'.format(key + 1), Na1)
        answer.set_atom('{}_Na2'.format(key + 1), Na2)
        return answer

    # ------------------------------------------------------------------
    def _get_neutralize_pos_NH3_type(self, ag):
        length = 3.187
        H1 = ag['H1']
        H2 = ag['H2']
        H3 = ag['H3']
        N = ag['N']

        # 重心を計算
        M = Position((H1.xyz.x + H2.xyz.x + H3.xyz.x) / 3.0,
                     (H1.xyz.y + H2.xyz.y + H3.xyz.y) / 3.0,
                     (H1.xyz.z + H2.xyz.z + H3.xyz.z) / 3.0)
        vNM = M - N.xyz
        vNM.norm()

        return N.xyz + length * vNM

    def _get_neutralize_pos_NH2_type(self, ag):
        length = 3.187
        H1 = ag['H1']
        H2 = ag['H2']
        N = ag['N']

        vNH1 = H1.xyz - N.xyz
        vNH2 = H2.xyz - N.xyz
        vM = 0.5 * (vNH1 + vNH2)

        vM.norm()

        answer = N.xyz + length * vM
        return answer

    def _get_neutralize_pos_COO_type(self, ag):
        length = 2.521
        O1 = ag['O1']
        O2 = ag['O2']
        C = ag['C']

        # 中点を計算
        M = Position(0.5 * (O1.xyz.x + O2.xyz.x),
                     0.5 * (O1.xyz.y + O2.xyz.y),
                     0.5 * (O1.xyz.z + O2.xyz.z))
        vCM = M - C.xyz
        vCM.norm()

        return C.xyz + length * vCM

    # -----------------------------------------------------------------
    def _get_neutralize_pos_POO_type(self, ag):
        length = 2.748
        O1 = ag['O1']
        O2 = ag['O2']
        P = ag['P']

        M = Position(0.5 * (O1.xyz.x + O2.xyz.x),
                     0.5 * (O1.xyz.y + O2.xyz.y),
                     0.5 * (O1.xyz.z + O2.xyz.z))

        vPM = M - P.xyz
        vPM.norm()

        return P.xyz + length * vPM


if __name__ == "__main__":
    import doctest
    doctest.testmod()
