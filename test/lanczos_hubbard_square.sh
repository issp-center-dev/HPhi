#!/bin/sh -e

mkdir -p lanczos_hubbard_square/
cd lanczos_hubbard_square

cat > stan.in <<EOF
W = 4
L = 2
model = "FermionHubbard"
method = "Lanczos"
lattice = "Tetragonal"
t = 1.0
U = 4.0
nelec = 8
2Sz = 0
EOF

${MPIRUN} ../../src/HPhi -s stan.in

# Check value

cat > reference.dat <<EOF
  -10.2529529552637637
    1.0444442739280932
    0.0000000000000000
EOF
paste output/zvo_energy.dat reference.dat > paste.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)**2)} END{printf "%8.6f", diff}' paste.dat`
test "${diff}" = "0.000000"

exit $?
