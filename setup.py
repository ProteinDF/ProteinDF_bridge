#!/usr/bin/env python
# -*- coding: utf-8 -*-


from setuptools import setup

setup()

# -------
# import sys
# import os
# from setuptools import setup
# from imp import reload

# # sys.path.append('./pdfbridge')
# # sys.path.append('./pdftests')
# exec(open("proteindf_bridge/_version.py").read())

# setup(
#     name="proteindf_bridge",
#     version=__version__,
#     description="bridge scripts the ProteinDF package and other data/package",
#     author="Toshiyuki HIRANO",
#     author_email="hiracchi@gmail.com",
#     url="http://proteindf.github.io/",
#     license="GPLv3",
#     packages=["proteindf_bridge"],
#     scripts=[
#         "scripts/doctest_runner.py",
#         "scripts/mpac2yml.py",
#         "scripts/mpac2txt.py",
#         "scripts/yml2mpac.py",
#         "scripts/brd-box.py",
#         "scripts/brd-density.py",
#         "scripts/brd-formula.py",
#         "scripts/brd-restructure.py",
#         "scripts/brd-select.py",
#         "scripts/brd-divide.py",
#         "scripts/brd-divide-mainchain.py",
#         "scripts/brd-select-path.py",
#         "scripts/brd-show-res.py",
#         "scripts/brd-renumber-resid.py",
#         "scripts/brd2txt.py",
#         "scripts/pdb2brd.py",
#         "scripts/brd2pdb.py",
#         "scripts/mmcif2txt.py",
#         "scripts/mmcif2mol2.py",
#         "scripts/gro2brd.py",
#         "scripts/brd2gro.py",
#         "scripts/xyz2brd.py",
#         "scripts/brd2xyz.py",
#         "scripts/module_inspect.py",
#         "scripts/db2txt.py",
#         "scripts/brd-setup-bond.py",
#         "scripts/brd-show-bonds.py",
#         "scripts/superposer.py",
#         "scripts/neutralize.py",
#         "scripts/reorder.py",
#         "scripts/crystallize.py",
#         "scripts/read_amber_prmtop.py"
#         # 'scripts/remove_wat.py',
#         # 'scripts/relax_protein.py',
#         # 'scripts/relax_protein.sh'
#     ],
#     install_requires=["configparser", "msgpack", "pyyaml", "numpy"],
#     package_data={
#         "proteindf_bridge": [
#             "data/ACE_ALA_NME.brd",
#             "data/ACE_ALA_NME_trans1.brd",
#             "data/ACE_ALA_NME_trans2.brd",
#             "data/ACE_ALA_NME_cis1.brd",
#             "data/ACE_ALA_NME_cis2.brd",
#         ]
#     },  # test_suite='tests'
# )
