#!/bin/bash

TESTS=" \
  position \
  atom \
  atomgroup \
  select \
  biopdb \
  "

rm -rf test_*

for i in ${TESTS}; do
    echo ">>>> test: ${i}"
    python -m unittest -v tests.test_${i}
    echo "done."
    echo
done

