#!/bin/sh -e

# Regression test for the PairLift terms of SpinGC in FullDiag.
#
# makeHam() used to loop over NPairLift/2 pairs of terms, so that the last
# term was dropped when the number of PairLift terms was odd, while the
# matrix-vector product of Lanczos/CG/TPQ applied all terms.
#
# Heisenberg chain (L = 4, J = 1, h = 0.3) with three PairLift terms
#   H += J_ij (S^+_i S^+_j + S^-_i S^-_j).
# The spectrum is compared with an independent exact diagonalization, and with
# the same terms written as InterAll.

testname="fulldiag_spingc_pairlift"
tolerance="0.00000001"

mkdir -p "${testname}"
cd "${testname}"

fail() {
  echo "FAILED (${testname}): $1" >&2
  exit 1
}

# compare <label> <file of eigenvalues> <file of eigenvalues>
compare() {
  n1=$(wc -l < "$2" | tr -d ' ')
  n2=$(wc -l < "$3" | tr -d ' ')
  [ "x${n1}" = "x${n2}" ] || fail "$1: ${n1} eigenvalues, expected ${n2}"
  paste "$2" "$3" > paste.dat
  maxdiff=$(awk 'BEGIN { m = 0 } { d = $1 - $2; if (d < 0) d = -d; if (d > m) m = d } END { printf "%.12f", m }' paste.dat)
  ok=$(awk -v d="${maxdiff}" -v tol="${tolerance}" 'BEGIN { print (d < tol) ? 1 : 0 }')
  [ "${ok}" -eq 1 ] || {
    cat paste.dat
    fail "$1: max diff ${maxdiff}"
  }
  echo "$1: ${n1} eigenvalues agree (max diff ${maxdiff})"
}

cat > stan.in <<EOF
L = 4
model = "SpinGC"
method = "FullDiag"
lattice = "chain"
J = 1.0
h = 0.3
2S = 1
EOF

../../src/HPhi -sdry stan.in > gen.log 2>&1 || { cat gen.log; fail "standard input generation failed"; }
grep -v -e "OneBodyG" -e "TwoBodyG" namelist.def > namelist_base.def

cat > reference.dat <<EOF
-2.2736031063
-2.0655575509
-1.0000000000
-0.7393232042
-0.6781690191
-0.6605027964
-0.1602936511
-0.0834265261
0.0000000000
0.0834265261
0.4129123674
0.6781690191
0.7393232042
1.1401722757
2.0655575509
2.5413149107
EOF

# (1) PairLift file with an odd number of terms
cat > pairlift.def <<EOF
======================
NPairLift 3
======================
========PairLift======
======================
    0     1  0.7
    1     2  0.4
    2     3  0.9
EOF
cp namelist_base.def namelist_pairlift.def
echo "PairLift pairlift.def" >> namelist_pairlift.def

rm -rf output
../../src/HPhi -e namelist_pairlift.def > run_pairlift.log 2>&1 || { cat run_pairlift.log; fail "HPhi failed with pairlift.def"; }
awk 'NF >= 2 { print $2 }' output/Eigenvalue.dat | sort -n > eigen_pairlift.dat

# (2) The same terms as InterAll:
#     S^+_i S^+_j = c^+_{i up} c_{i down} c^+_{j up} c_{j down} and its conjugate
cat > interall.def <<EOF
======================
NInterAll 6
======================
========InterAll======
======================
    0 0 0 1 1 0 1 1  0.7 0.0
    1 1 1 0 0 1 0 0  0.7 0.0
    1 0 1 1 2 0 2 1  0.4 0.0
    2 1 2 0 1 1 1 0  0.4 0.0
    2 0 2 1 3 0 3 1  0.9 0.0
    3 1 3 0 2 1 2 0  0.9 0.0
EOF
cp namelist_base.def namelist_interall.def
echo "InterAll interall.def" >> namelist_interall.def

rm -rf output
../../src/HPhi -e namelist_interall.def > run_interall.log 2>&1 || { cat run_interall.log; fail "HPhi failed with interall.def"; }
awk 'NF >= 2 { print $2 }' output/Eigenvalue.dat | sort -n > eigen_interall.dat

compare "PairLift vs exact diagonalization" eigen_pairlift.dat reference.dat
compare "PairLift vs InterAll" eigen_pairlift.dat eigen_interall.dat

echo "PairLift terms of SpinGC in FullDiag are correct."
