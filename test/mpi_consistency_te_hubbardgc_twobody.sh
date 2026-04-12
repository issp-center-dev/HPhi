#!/bin/sh -e
# Regression test for H-1 (PR #216 review):
# HubbardGC TimeEvolution + TETwoBody MPI consistency.
# Compares serial vs MPI Flct output for step-dependent InterAll terms.

TOLERANCE="0.000001"
HPHI=../src/HPhi
TE_STEPS=20

run_hphi() {
    label="$1"
    logfile="$2"
    shift 2
    if ! "$@" > "${logfile}" 2>&1; then
        echo "${label} failed"
        cat "${logfile}"
        exit 1
    fi
}

if [ -z "${MPIRUN}" ]; then echo "Error: MPIRUN is not set."; exit 1; fi
MPI_NP=$(printf "%s\n" "${MPIRUN}" | awk '{for(i=1;i<=NF;i++){if($i=="-np"||$i=="-n"){print $(i+1); exit}}}')
if ! printf "%s\n" "${MPI_NP}" | grep -Eq "^[0-9]+$"; then echo "Error: bad MPIRUN"; exit 1; fi
if [ "${MPI_NP}" -le 1 ]; then echo "Error: need >1 rank"; exit 1; fi

TOPDIR=mpi_consistency_te_hubbardgc_twobody
rm -rf ${TOPDIR}
mkdir -p ${TOPDIR}

write_tetwobody_def() {
cat > tetwobody.def <<EOF
========================
NTimeSteps    ${TE_STEPS}
========================
========================
========================
EOF
printf '%s\n' '0.00  0' >> tetwobody.def
i=1
while [ "${i}" -lt "${TE_STEPS}" ]; do
    printf '0.%02d  2\n' "${i}" >> tetwobody.def
    printf '%s\n' '0  0  1  0  1  1  0  1  0.5  0.0' >> tetwobody.def
    printf '%s\n' '0  1  1  1  1  0  0  0  0.5  0.0' >> tetwobody.def
    i=$((i + 1))
done
}

# ---- Serial run ----
mkdir -p ${TOPDIR}/serial
cd ${TOPDIR}/serial

cat > stan_gs.in <<'EOF'
L = 2
model = "HubbardGC"
method = "Lanczos"
lattice = "chain"
t = 1.0
U = 4.0
initial_iv = 1
EigenvecIO = "out"
EOF

echo "  [serial] Input generation..."
run_hphi "  [serial] Input generation" hphi_sdry.log ../../${HPHI} -sdry stan_gs.in

echo "  [serial] Ground state..."
run_hphi "  [serial] Ground state" hphi_gs.log ../../${HPHI} -e namelist.def

# Create TETwoBody: step 0 empty, step 1-19 with inter-PE InterAll
write_tetwobody_def

# Patch for TE mode
sed -e 's/^CalcType.*/CalcType   4/' \
    -e 's/^InputEigenVec.*/InputEigenVec   1/' \
    -e 's/^OutputEigenVec.*/OutputEigenVec   0/' \
    calcmod.def > _t && mv _t calcmod.def
sed -e "s/^Lanczos_max.*/Lanczos_max    ${TE_STEPS}/" \
    modpara.def > _t && mv _t modpara.def
echo "   TETwoBody  tetwobody.def" >> namelist.def

echo "  [serial] Time evolution..."
run_hphi "  [serial] Time evolution" hphi_te.log ../../${HPHI} -e namelist.def
cd ../..

# ---- MPI run ----
mkdir -p ${TOPDIR}/mpi
cd ${TOPDIR}/mpi

cat > stan_gs.in <<'EOF'
L = 2
model = "HubbardGC"
method = "Lanczos"
lattice = "chain"
t = 1.0
U = 4.0
initial_iv = 1
EigenvecIO = "out"
EOF

echo "  [mpi] Input generation..."
run_hphi "  [mpi] Input generation" hphi_sdry.log ../../${HPHI} -sdry stan_gs.in

echo "  [mpi] Ground state..."
run_hphi "  [mpi] Ground state" hphi_gs.log ${MPIRUN} ../../${HPHI} -e namelist.def

write_tetwobody_def

sed -e 's/^CalcType.*/CalcType   4/' \
    -e 's/^InputEigenVec.*/InputEigenVec   1/' \
    -e 's/^OutputEigenVec.*/OutputEigenVec   0/' \
    calcmod.def > _t && mv _t calcmod.def
sed -e "s/^Lanczos_max.*/Lanczos_max    ${TE_STEPS}/" \
    modpara.def > _t && mv _t modpara.def
echo "   TETwoBody  tetwobody.def" >> namelist.def

echo "  [mpi] Time evolution..."
run_hphi "  [mpi] Time evolution" hphi_te.log ${MPIRUN} ../../${HPHI} -e namelist.def
cd ../..

# ---- Compare ----
paste ${TOPDIR}/serial/output/Flct.dat ${TOPDIR}/mpi/output/Flct.dat > ${TOPDIR}/compare.dat
diff=$(awk '
BEGIN { maxdiff = 0.0 }
NR > 1 {
    d = sqrt(($4 - $12) * ($4 - $12));  if (d > maxdiff) maxdiff = d
    d = sqrt(($5 - $13) * ($5 - $13));  if (d > maxdiff) maxdiff = d
}
END { printf "%.10f", maxdiff }
' ${TOPDIR}/compare.dat)

echo "Max Flct difference: ${diff}"

result=$(awk -v diff=${diff} -v tol=${TOLERANCE} 'BEGIN { print (diff < tol) ? "PASS" : "FAIL" }')

if [ "${result}" = "PASS" ]; then
    echo "MPI consistency test PASSED for HubbardGC TE TETwoBody"
    exit 0
else
    echo "MPI consistency test FAILED for HubbardGC TE TETwoBody"
    echo "Expected difference < ${TOLERANCE}, got ${diff}"
    echo "--- Serial ---"
    cat ${TOPDIR}/serial/output/Flct.dat
    echo "--- MPI ---"
    cat ${TOPDIR}/mpi/output/Flct.dat
    exit 1
fi
