#!/bin/sh -e

mkdir -p lobcg_spinless/
cd lobcg_spinless
python "$1/test/testSpinlessCalc.py" -p "$2/src/HPhi" -mpi "${MPIRUN}" -m "SpinlessFermion" -s 8 -t "LOBCG"

# Check value: flct
cat > reference.dat <<EOF
   0
   -4.6769911637521853
   0.0000000000000000
   0.0000000000000000 

   1
   -3.2141371350547225
   0.0000000000000000
   0.0000000000000000 

   2
   -3.0976937063656376
   0.0000000000000000
   0.0000000000000000 

   3
   -2.9358904049883559
   0.0000000000000000
   0.0000000000000000 

   4
   -2.9358904049883359
   0.0000000000000000
   0.0000000000000000 
EOF
paste output/zvo_energy.dat reference.dat > paste.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste.dat`
test "${diff}" = "0.000000"
exit $?
