#!/bin/sh

mkdir -p lanczos_hubbardgc_tri/
cd lanczos_hubbardgc_tri

cat > stan.in <<EOF
a0w = 3
a0l = 0
a1w = -1
a1l = 2
model = "HubbardGC"
method = "Lanczos"
lattice = "Triangular"
t = 1.0
U = 4.0
EOF

${MPIRUN} ../../src/HPhi -s stan.in
if test $? -ne 0 ; then
    exit 1
fi

# Check value

cat > reference.dat <<EOF
  -12.3140744040600509
    0.3472086205791791
   -0.3753503113209436
EOF

diff=`paste output/zvo_energy.dat reference.dat | awk '
BEGIN{diff=0.0} {diff+=sqrt(($2-$3)**2)} END{printf "%8.6f", diff}'`

test "${diff}" = "0.000000"

exit $?
