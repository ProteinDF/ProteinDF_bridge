#!/bin/bash

# ======================================================================
# doctest
# ======================================================================
DOCTESTS=" \
  aminoacid \
  atom \
  atomgroup \
  common \
  ionpair \
  mail \
  matrix \
  modeling \
  neutralize \
  periodictable \
  position \
  select \
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
  "

rm -rf test_*

for i in ${UNITTESTS}; do
    echo ">>>> unittest: ${i}"
    python -m unittest -v tests.test_${i}
    echo "done."
    echo
done
