#!/bin/bash

MODELS="\
  trans1
  trans2
  cis1
  cis2
"

for model in ${MODELS}; do
    ${PDF_HOME}/bin/yml2mpac.py ACE_ALA_NME_${model}.yml ACE_ALA_NME_${model}.brd
done
