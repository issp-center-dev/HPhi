#!/bin/sh -e

mkdir -p lobcg_kondogc_chain/
cd lobcg_kondogc_chain

cat > stan.in <<EOF
L = 4
model = "KondoGC"
method = "CG"
lattice = "chain"
t = 1.0
J = 4.0
mu = 0.6
h = -0.2
exct = 3
EOF

${MPIRUN} ../../src/HPhi -s stan.in

# Check value

cat > reference.dat <<EOF
 0
  -15.0776213781764330
  0.1072445462409950
  -0.0000000000000000

 1
  -13.6371064976052097
  1.0413195150782721
  -0.4999999999999999

 2
  -13.4371064976051855
  1.0413195150782717
  0.4999999999999998
EOF
paste output/zvo_energy.dat reference.dat > paste.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste.dat`
test "${diff}" = "0.000000"

exit $?
