#!/bin/sh -e

mkdir -p lanczos_kondogc_chain/
cd lanczos_kondogc_chain

cat > stan.in <<EOF
L = 4
model = "KondoGC"
method = "Lanczos"
lattice = "chain"
t = 1.0
J = 4.0
EOF

${MPIRUN} ../../src/HPhi -s stan.in

# Check value

cat > reference.dat <<EOF
  -12.6776213781764451
    0.1072445462409946
    0.0000000000000000
EOF
paste output/zvo_energy.dat reference.dat > paste.dat
diff=`awk 'BEGIN{diff=0.0} $1=="Energy"{diff+=sqrt(($2-$3)**2)} END{printf "%8.6f", diff}' paste.dat`
test "${diff}" = "0.000000"

exit $?
