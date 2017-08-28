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
exct = 3
EOF

${MPIRUN} ../../src/HPhi -s stan.in

# Check value

cat > reference.dat <<EOF
 0
  -18.3623620055430514
  0.0000000000000000
  -0.0000000000000086

 1
  -16.8165118261059305
  0.0000000000000000
  0.7065246882937294

 2
  -16.8165118261057280
    0.0000000000000000
  -0.6256238968226678
EOF
paste output/zvo_energy.dat reference.dat > paste.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)**2)} END{printf "%8.6f", diff}' paste.dat`
test "${diff}" = "0.000000"

exit $?
