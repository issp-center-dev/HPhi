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
  -5.5492800776055233 
  0.0000000000000000 
  -0.2767335797116458 
EOF
paste output/zvo_energy.dat reference.dat > paste.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.5f", diff}' paste.dat`
test "${diff}" = "0.00000"
echo "${diff}"

exit $?