#!/bin/bash

TESTS=" \
  position \
  atom \
  atomgroup \
  select \
  "

rm -rf test_*

for i in ${TESTS}; do
    echo ">>>> test: ${i}"
    python -m unittest tests.test_${i}
    echo "done."
    echo
done

