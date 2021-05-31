#!/bin/sh -e

mkdir -p lobcg_genspin_ladder_restart/
cd lobcg_genspin_ladder_restart

cat > stan.in <<EOF
L = 3
W = 2
model = "Spin"
method = "CG"
lattice = "ladder"
J0 = 1.0
J1 = 1.0
D=0.2
2Sz = 0
2S = 3
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

# Check value

cat > reference.dat <<EOF
 0
  -16.6896771153847894
  0.0000000000000000
  0.0000000000000000

 1
  -15.8422975895659750
  0.0000000000000000
  0.0000000000000000
EOF
paste output/zvo_energy.dat reference.dat > paste.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste.dat`

# Check value for Restart_inout
rm output/zvo_energy.dat
cp ./stan.in stan3.in
sed -i -e "s/Restart = \"Restart_out\"/Restart = \"Restart\"/g" stan3.in
sed -i -e "s/Lanczos_max = 10/Lanczos_max = 1000/g" stan3.in 

${MPIRUN} ../../src/HPhi -s stan3.in

cat > reference.dat <<EOF
 0
  -16.6896771153847894
  0.0000000000000000
  0.0000000000000000

 1
  -15.8422975895659750
  0.0000000000000000
  0.0000000000000000
EOF

paste output/zvo_energy.dat reference.dat > paste2.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste2.dat`


test "${diff}" = "0.000000"

exit $?
