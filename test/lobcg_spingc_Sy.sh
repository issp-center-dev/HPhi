#!/bin/sh -e

mkdir -p lobcg_spingc_Sy/
cd lobcg_spingc_Sy

cat > stan.in <<EOF
model = SpinGC
method = CG
lattice = Honeycomb
W = 2
L = 2
J0y = -3
J0xy = 3
J1xz = 1
J2yz = 2
h = 0.2
Gamma =  0.5
exct = 1
LanczosEps = 16
EOF

${MPIRUN} ../../src/HPhi -sdry stan.in
cat > trans.def <<EOF
===
N 2
===
===
===
0 0 0 1  0.5  0.5 
0 1 0 0  0.5 -0.5 
EOF

${MPIRUN} ../../src/HPhi -e namelist.def

# Check value

cat > reference.dat <<EOF
  0
  -5.5492800776068707 
  0.0000000000000000 
  -0.2767359441047553
EOF
paste output/zvo_energy.dat reference.dat > paste.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste.dat`
test "${diff}" = "0.000000"

exit $?
