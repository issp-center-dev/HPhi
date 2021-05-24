#!/bin/sh -e

mkdir -p lobcg_genspingc_ladder/
cd lobcg_genspingc_ladder

cat > stan.in <<EOF
L = 2
W = 2
model = "SpinGC"
method = "CG"
lattice = "ladder"
J0 = 1.0
J1 = 1.0
2S = 3
h=-0.01
exct = 3
EOF

${MPIRUN} ../../src/HPhi -s stan.in

# Check value

cat > reference.dat <<EOF
 0
  -18.3623620055432717
  0.0000000000000000
  -0.0000000000000000

 1
  -16.8265118261059037
  0.0000000000000000
  0.9999999999999948

 2
  -16.8165118261057955
  0.0000000000000000
  -0.0000000000000047
EOF
paste output/zvo_energy.dat reference.dat > paste.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste.dat`
test "${diff}" = "0.000000"

exit $?
