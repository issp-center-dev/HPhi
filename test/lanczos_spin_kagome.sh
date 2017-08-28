#!/bin/sh

mkdir -p lanczos_spin_kagome/
cd lanczos_spin_kagome

cat > stan.in <<EOF
a0w = 1
a0l = 1
a1w = -1
a1l = 2
model = "Spin"
method = "Lanczos"
lattice = "kagome"
J = 1.0
2Sz = 1
EOF

${MPIRUN} ../../src/HPhi -s stan.in
if test $? -ne 0 ; then
    exit 1
fi

echo "\nCheck value\n"

cat > reference.dat <<EOF
  -3.9690017499285153
   0.0000000000000000
   0.5000000000000000
EOF

diff=`paste output/zvo_energy.dat reference.dat | awk '
BEGIN{diff=0.0} {diff+=sqrt(($2-$3)**2)} END{printf "%8.6f", diff}'`

echo ${diff}
test "${diff}" = "0.000000"

exit $?
