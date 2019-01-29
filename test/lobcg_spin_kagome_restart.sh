#!/bin/sh -e

mkdir -p lobcg_spin_kagome_restart/
cd lobcg_spin_kagome_restart

cat > stan.in <<EOF
a0w = 1
a0l = 1
a1w = -1
a1l = 2
model = "Spin"
method = "CG"
lattice = "kagome"
J0 = 1.0
J1 = 0.5
J2 = 0.7
J'=0.2
2Sz = 1
exct = 2
Lanczos_max = 10
outputmode = "None"
Restart = "Restart_out"
EOF

${MPIRUN} ../../src/HPhi -s stan.in

# Check value for Restart_in
cp ./stan.in stan2.in
sed -i -e "s/Restart = \"Restart_out\"/Restart = \"Restart_in\"/g" stan2.in
sed -i -e "s/Lanczos_max = 10/Lanczos_max = 1000/g" stan2.in 

${MPIRUN} ../../src/HPhi -s stan2.in

cat > reference.dat <<EOF
 0
  -3.6523534273467710
  0.0000000000000000
  0.5000000000000000

 1
  -3.3180023889678378
  0.0000000000000000
  0.5000000000000000
EOF
paste output/zvo_energy.dat reference.dat > paste.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste.dat`
rm output/zvo_energy.dat

# Check value for Restart_inout
cp ./stan.in stan3.in
sed -i -e "s/Restart = \"Restart_out\"/Restart = \"Restart\"/g" stan3.in
sed -i -e "s/Lanczos_max = 10/Lanczos_max = 1000/g" stan3.in 

${MPIRUN} ../../src/HPhi -s stan3.in

cat > reference.dat <<EOF
 0
  -3.6523534273467710
  0.0000000000000000
  0.5000000000000000

 1
  -3.3180023889678378
  0.0000000000000000
  0.5000000000000000
EOF

paste output/zvo_energy.dat reference.dat > paste2.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste2.dat`

test "${diff}" = "0.000000"

exit $?
