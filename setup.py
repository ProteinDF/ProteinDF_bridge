#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2019 The ProteinDF development team.
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

import sys
import os
from setuptools import setup
from imp import reload

# sys.path.append('./pdfbridge')
# sys.path.append('./pdftests')
exec(open("proteindf_bridge/_version.py").read())

setup(name='proteindf_bridge',
      version=__version__,
      description='bridge scripts the ProteinDF package and other data/package',
      author='Toshiyuki HIRANO',
      author_email='hiracchi@gmail.com',
      url='http://proteindf.github.io/',
      license='GPLv3',
      packages=['proteindf_bridge'],
      scripts=[
          'scripts/doctest_runner.py',

          'scripts/mpac2yml.py',
          'scripts/mpac2txt.py',
          'scripts/yml2mpac.py',

          'scripts/brd-restructure.py',
          'scripts/brd-select.py',
          'scripts/brd-divide.py',
          'scripts/brd-divide-mainchain.py',
          'scripts/brd-select-path.py',
          'scripts/brd-renumber-resid.py',

          'scripts/brd2txt.py',

          'scripts/pdb2brd.py',
          'scripts/brd2pdb.py',

          'scripts/mmcif2txt.py',
          'scripts/mmcif2mol2.py',

          'scripts/gro2brd.py',
          'scripts/brd2gro.py',

          'scripts/xyz2brd.py',
          'scripts/brd2xyz.py',

          'scripts/module_inspect.py',
          'scripts/db2txt.py',

          'scripts/brd-setup-bond.py',
          'scripts/brd-show-bonds.py',
          'scripts/superposer.py',

          'scripts/neutralize.py',
          'scripts/reorder.py',
          'scripts/crystallize.py',

          'scripts/read_amber_prmtop.py'

          # 'scripts/remove_wat.py',
          # 'scripts/relax_protein.py',
          # 'scripts/relax_protein.sh'
      ],

      install_requires=[
          'configparser',
          'msgpack-python',
          'pyyaml',
          'numpy'
      ],

      data_files=[('data', ['data/ACE_ALA_NME.brd'])],
      # test_suite='tests'
      )
