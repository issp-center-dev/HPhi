#!/bin/sh -e

mkdir -p lobcg_hubbardgc_tri/
cd lobcg_hubbardgc_tri

cat > stan.in <<EOF
a0w = 3
a0l = 0
a1w = -1
a1l = 2
model = "HubbardGC"
method = "CG"
lattice = "Triangular"
t = 1.0
U = 4.0
h = -3.0
exct = 2
EOF

${MPIRUN} ../../src/HPhi -s stan.in

# Check value

cat > reference.dat <<EOF
   0
 -15.3140744040600296
   0.3472086205791794
  -1.0000000000000004

   1
 -14.2336534858459416
   0.4631929500975011
  -1.4999999999999998
EOF
paste output/zvo_energy.dat reference.dat > paste.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste.dat`
test "${diff}" = "0.000000"

exit $?
