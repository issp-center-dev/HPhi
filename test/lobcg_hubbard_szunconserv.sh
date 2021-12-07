#!/bin/sh -e

mkdir -p lobcg_hubbard_szconserv/
cd lobcg_hubbard_szconserv

cat > stan.in <<EOF
W = 4
L = 2
model = "FermionHubbard"
method = "CG"
lattice = "Tetragonal"
t = 1.0
U = 4.0
nelec = 8
exct = 3
h = 2.0
Gamma = 1.0
EOF

${MPIRUN} ../../src/HPhi -s stan.in

# Check value

cat > reference.dat <<EOF
 0
    -12.0689666528263331
    1.1830984216584166
    0.8944271910345338

 1
    -12.0188079575967599
    0.8548755838399427
    1.7888543820161609

 2
    -11.7934541620331608
    0.9348617586588275
    1.7888543739897627
EOF
paste output/zvo_energy.dat reference.dat > paste.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste.dat`
test "${diff}" = "0.000000"

exit $?
