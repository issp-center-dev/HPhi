#!/bin/sh -e

mkdir -p fulldiag_spingc_gamma/
cd fulldiag_spingc_gamma

cat > stan.in <<EOF
model = "SpinGC"
method = "FullDiag"
lattice = "chain"
L = 8
J = 1.0
Gamma = 0.5
EOF

${MPIRUNFC} ../../src/HPhi -s stan.in

# Check value
#
# This test pins the SpinGC transverse-field (Gamma != 0) term to its
# pre-phase-2 historical behavior. It is the only behavior-touched hunk of
# the phase 2 distributed-Hamiltonian-generation change (makeHam.c): the
# old code computed the row index and the out-param update of `off` in the
# same statement (Ham[off + 1][j] += tmp_trans * child_SpinGC_CisAit(...,
# &off)), whose relative evaluation order the C standard leaves unspecified
# whenever a side effect and a value computation on the same object are
# unsequenced. Phase 2 splits this into two statements (call, then use of
# the now-fully-updated `off`) to make the order explicit.
#
# The reference below was produced by building HPhi at commit d0cde764
# (immediately before the phase-2 implementation commits) in a separate
# git worktree and running this exact stan.in with it. The current build
# reproduced ALL 256 eigenvalues in output/Eigenvalue.dat, and the full
# output/zvo_phys.dat, bit-for-bit against that reference (see
# .superpowers/sdd/p2-final-fix-report.md for the full comparison
# evidence), so only the first 10 (lowest) eigenvalues are embedded here
# to keep the test compact.
cat > reference_eigenvalue.dat <<EOF
 0 -3.6510934089
 1 -3.6284190638
 2 -3.1284190638
 3 -2.9587385089
 4 -2.9587385089
 5 -2.8019377358
 6 -2.6996281483
 7 -2.6451483739
 8 -2.6451483739
 9 -2.6284190638
EOF
head -10 output/Eigenvalue.dat > eigenvalue.dat
paste eigenvalue.dat reference_eigenvalue.dat > paste.dat
diff=`awk 'BEGIN{max=0}{d=$2-$4; if(d<0)d=-d; if(d>max)max=d}END{print max}' paste.dat`
test "`echo "$diff < 0.000001" | bc`" = "1"

exit $?
