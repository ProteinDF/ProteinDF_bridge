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

for i in ${DOCTESTS}; do
    echo ">>>> doctest: ${i}"
    python -m proteindf_bridge.${i}
    echo "done."
done

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
