#!/bin/sh -e

mkdir -p lobcg_spinless_GC/
cd lobcg_spinless_GC
python "$1/test/testSpinlessCalc.py" -p "$2/src/HPhi" -mpi "${MPIRUN}" -m "SpinlessFermionGC" -s 8 -t "LOBCG"

# Check value: flct
cat > reference.dat <<EOF
   0
   -4.6769911637521870 
   0.0000000000000000
   0.0000000000000000 

   1
   -4.2804775927510308
   0.0000000000000000
   0.0000000000000000 

   2
   -4.280477592751029
   0.0000000000000000
   0.0000000000000000 

   3
   -3.6635737201988126
   0.0000000000000000
   0.0000000000000000 

   4
   -3.280477592750769
   0.0000000000000000
   0.0000000000000000 
EOF
paste output/zvo_energy.dat reference.dat > paste.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste.dat`
test "${diff}" = "0.000000"
exit $?
