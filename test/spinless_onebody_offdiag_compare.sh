#!/bin/sh -e

mkdir -p spinless_onebody_offdiag_compare/
cd spinless_onebody_offdiag_compare

python3 "$1/test/compare_spinless_hubbard.py" ../../src/HPhi
