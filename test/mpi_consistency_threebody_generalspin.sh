#!/bin/sh -e
# MPI consistency for the general-spin (SpinGC, 2S=2) three-body Green function
# when the 5th/6th operators (c5 a6) form a same-site TRANSVERSE pair
# (sigma5 != sigma6) on an INTER-PROCESS site.
#
# Prerelease finding H-2: mltplyGeneralSpinGC_mini called the inter-process
# routine child_GC_CisAit_GeneralSpin_MPIdouble under plain M_MLTPLY, so the
# Hermitian-conjugate branch was added too (the S=1/2 path switches to
# M_MLTPLY2 to suppress it). The three-body Green of a transverse inter-process
# operator then came out as the Hermitian-pair sum instead of the one-directional
# value -- wrong, with no error. A three-body Green expectation <GS|A|GS> is
# invariant under the global phase of the (non-degenerate) ground state, so a
# direct serial-vs-MPI comparison of zvo_ThreeBody is well defined.
#
# Requires np that is a power of 3 (spin-1 MPI site split); use exact:3 so the
# single inter-process site is site 5 of the L=6 chain.

if [ -z "${MPIRUN}" ]; then echo "MPIRUN not set. Skipping."; exit 0; fi
MPI_NP=$(printf "%s\n" "${MPIRUN}" | awk '{for(i=1;i<=NF;i++){if($i=="-np"||$i=="-n"){print $(i+1); exit}}}')
if ! printf "%s\n" "${MPI_NP}" | grep -Eq "^[0-9]+$" || [ "${MPI_NP}" -le 1 ]; then
  echo "Error: MPIRUN must include -np/-n with an integer > 1 (got '${MPIRUN}')."; exit 1
fi

# $1 = CMAKE_SOURCE_DIR (add_hphi_mpi_test_with_srcdir)
SRCDIR="$1"

mkdir -p mpi_consistency_threebody_generalspin
cd mpi_consistency_threebody_generalspin
cp "${SRCDIR}/test/SpinOneThreeBody.py" .

# Generate the SpinGC 2S=2 L=6 chain definition files (mag/interall/green1/
# green3/namelist/stan), then the standard-mode defs.
python3 SpinOneThreeBody.py generate > log_generate.txt 2>&1
../../src/HPhi -sdry stan.in > log_sdry.txt 2>&1

# Override green3.def: 5th/6th operators are a TRANSVERSE pair on site 5
# (the inter-process site at np=3). The first four operators are transverse
# pairs on local sites 0/1 so the three-body Green is non-trivial.
cat > green3.def <<EOF
=====
NumGreen 6
=====
=====
=====
0 0 0 1 1 0 1 1 5 0 5 1
0 1 0 0 1 1 1 0 5 0 5 1
1 0 1 1 2 0 2 1 5 1 5 2
0 0 0 1 1 1 1 2 5 0 5 1
1 1 1 2 0 0 0 1 5 1 5 2
0 1 0 2 1 0 1 1 5 0 5 1
EOF

run_side() {  # $1 = tag, $2 = MPI prefix (empty for serial)
  tag="$1"; pfx="$2"
  rm -rf output
  ${pfx} ../../src/HPhi -e open_namelist.def > "log_${tag}.txt" 2>&1
  cp output/zvo_ThreeBody_eigen0.dat "tb_${tag}.dat"
}

run_side serial ""
# Confirm site 5 is in the inter-process region for this np.
run_side mpi "${MPIRUN}"

# Compare the three-body Green (columns 13 = Re, 14 = Im) row by row.
d=$(paste tb_serial.dat tb_mpi.dat | awk '
  { ds=$13-$27; di=$14-$28; if(ds<0)ds=-ds; if(di<0)di=-di;
    if(ds>mr)mr=ds; if(di>mi)mi=di }
  END{ m=(mr>mi)?mr:mi; printf "%.3e", m+0 }')
echo "[threebody-genspin] max |serial - mpi| (Re/Im) = ${d}"
awk -v d="${d}" 'BEGIN{ exit (d < 1e-6) ? 0 : 1 }' || {
  echo "[threebody-genspin] MISMATCH (transverse inter-process 3-body Green differs serial vs MPI)"; exit 1; }
echo "General-spin transverse inter-process three-body Green: serial == MPI within tolerance."
