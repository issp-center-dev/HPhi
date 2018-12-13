#!/bin/sh -e

mkdir -p lobcg_kondo_chain_restart/
cd lobcg_kondo_chain_restart

cat > stan.in <<EOF
L = 4
model = "Kondo"
method = "CG"
lattice = "chain"
t = 1.0
J = 4.0
nelec = 4
2Sz = 0
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
  -12.6776213781764220 
  0.1072445462409922 
  0.0000000000000000 

 1
  -9.8347989641820277 
  0.4289751043491778 
  0.0000000000000000 
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
  -12.6776213781764220 
  0.1072445462409922 
  0.0000000000000000 

 1
  -9.8347989641820277 
  0.4289751043491778 
  0.0000000000000000 
EOF

paste output/zvo_energy.dat reference.dat > paste2.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste2.dat`

test "${diff}" = "0.000000"

exit $?
