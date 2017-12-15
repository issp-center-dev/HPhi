#!/bin/sh -e

mkdir -p lobcg_spin_kagome/
cd lobcg_spin_kagome

cat > stan.in <<EOF
a0w = 1
a0l = 1
a1w = -1
a1l = 2
model = "Spin"
method = "CG"
lattice = "kagome"
J0 = 1.0
J1 = 0.5
J2 = 0.7
J'=0.2
2Sz = 1
exct = 2
EOF

${MPIRUN} ../../src/HPhi -s stan.in

# Check value

cat > reference.dat <<EOF
 0
  -3.6523534273467710
  0.0000000000000000
  0.5000000000000000

 1
  -3.3180023889678378
  0.0000000000000000
  0.5000000000000000
EOF
paste output/zvo_energy.dat reference.dat > paste.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste.dat`
test "${diff}" = "0.000000"

exit $?
