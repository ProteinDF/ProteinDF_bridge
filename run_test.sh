#!/bin/bash

# ======================================================================
# doctest
# ======================================================================
DOCTESTS=" \
  aminoacid \
  common \
  ionpair \
  mail \
  matrix \
  modeling \
  neutralize \
  periodictable \
  ssbond \
  superposer \
  utils \
  vector \
  xyz \
  "

# ======================================================================
# unit test
# ======================================================================
UNITTESTS=" \
  position \
  atom \
  atomgroup \
  select \
  biopdb \
  gro \
  "

rm -rf test_*

for i in ${UNITTESTS}; do
    echo ">>>> unittest: ${i}"
    python -m unittest -v tests.test_${i}
    echo
done
