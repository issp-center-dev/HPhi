#!/bin/sh -e

mkdir -p lobcg_genspin_ladder/
cd lobcg_genspin_ladder

cat > stan.in <<EOF
L = 3
W = 2
model = "Spin"
method = "CG"
lattice = "ladder"
J0 = 1.0
J1 = 1.0
D=0.2
2Sz = 0
2S = 3
exct = 2
EOF

${MPIRUN} ../../src/HPhi -s stan.in

# Check value

cat > reference.dat <<EOF
 0
  -16.6896771153847965
  0.0000000000000000
  0.0000000000000000

 1
  -15.8422975895661722
  0.0000000000000000
  0.0000000000000000
EOF
paste output/zvo_energy.dat reference.dat > paste.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste.dat`
test "${diff}" = "0.000000"

exit $?
