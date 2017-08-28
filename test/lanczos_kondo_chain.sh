#!/bin/sh

mkdir -p lanczos_kondo_chain/
cd lanczos_kondo_chain

cat > stan.in <<EOF
L = 4
model = "Kondo"
method = "Lanczos"
lattice = "chain"
t = 1.0
J = 4.0
nelec = 4
2Sz = 0
EOF

${MPIRUN} ../../src/HPhi -s stan.in
if test $? -ne 0 ; then
    exit 1
fi

echo "\nCheck value\n"

cat > reference.dat <<EOF
  -12.6776213781764451
    0.1072445462409946
    0.0000000000000000
EOF

diff=`paste output/zvo_energy.dat reference.dat | awk '
BEGIN{diff=0.0} {diff+=sqrt(($2-$3)**2)} END{printf "%8.6f", diff}'`

echo ${diff}
test "${diff}" = "0.000000"

exit $?
