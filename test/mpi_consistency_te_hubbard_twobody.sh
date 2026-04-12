#!/bin/sh -e
# Regression test for §9.2 (PR #216 review):
# canonical Hubbard TimeEvolution + expert TEOneBody MPI consistency.
# Verifies that step-dependent Transfer terms are correctly processed
# when batching falls back to per-term MPI path.

TOLERANCE="0.000001"
HPHI=../src/HPhi
TE_STEPS=20

if [ -z "${MPIRUN}" ]; then echo "Error: MPIRUN is not set."; exit 1; fi
MPI_NP=$(printf "%s\n" "${MPIRUN}" | awk '{for(i=1;i<=NF;i++){if($i=="-np"||$i=="-n"){print $(i+1); exit}}}')
if ! printf "%s\n" "${MPI_NP}" | grep -Eq "^[0-9]+$"; then echo "Error: bad MPIRUN"; exit 1; fi
if [ "${MPI_NP}" -le 1 ]; then echo "Error: need >1 rank"; exit 1; fi

TOPDIR=mpi_consistency_te_hubbard_twobody
rm -rf ${TOPDIR}
mkdir -p ${TOPDIR}

write_teonebody_def() {
cat > teonebody.def <<EOF
========================
NTimeSteps    ${TE_STEPS}
========================
========================
========================
EOF
printf '%s\n' '0.00  0' >> teonebody.def
i=1
while [ "${i}" -lt "${TE_STEPS}" ]; do
    printf '0.%02d  2\n' "${i}" >> teonebody.def
    printf '%s\n' '0  0  8  0  0.3  0.0' >> teonebody.def
    printf '%s\n' '8  0  0  0  0.3  0.0' >> teonebody.def
    i=$((i + 1))
done
}

# ---- Serial run ----
mkdir -p ${TOPDIR}/serial
cd ${TOPDIR}/serial

cat > stan_gs.in <<'EOF'
a0w = 3
a0l = 0
a1w = 0
a1l = 3
model = "Hubbard"
method = "Lanczos"
lattice = "square"
t = 1.0
U = 4.0
nelec = 8
2Sz = 0
EigenvecIO = "out"
EOF

echo "  [serial] Ground state..."
../../${HPHI} -s stan_gs.in > /dev/null 2>&1

# TEOneBody: step 0 empty, step 1-19 add inter-PE transfer (site 0 <-> site 8)
# With 4 MPI ranks on 9-site square, site 8 is inter-process.
write_teonebody_def

sed -e 's/^CalcType.*/CalcType   4/' \
    -e 's/^InputEigenVec.*/InputEigenVec   1/' \
    calcmod.def > _t && mv _t calcmod.def
sed -e "s/^Lanczos_max.*/Lanczos_max    ${TE_STEPS}/" \
    modpara.def > _t && mv _t modpara.def
echo "   TEOneBody  teonebody.def" >> namelist.def

echo "  [serial] Time evolution..."
../../${HPHI} -e namelist.def > /dev/null 2>&1
cd ../..

# ---- MPI run ----
mkdir -p ${TOPDIR}/mpi
cd ${TOPDIR}/mpi

cat > stan_gs.in <<'EOF'
a0w = 3
a0l = 0
a1w = 0
a1l = 3
model = "Hubbard"
method = "Lanczos"
lattice = "square"
t = 1.0
U = 4.0
nelec = 8
2Sz = 0
EigenvecIO = "out"
EOF

echo "  [mpi] Ground state..."
${MPIRUN} ../../${HPHI} -s stan_gs.in > /dev/null 2>&1

write_teonebody_def

sed -e 's/^CalcType.*/CalcType   4/' \
    -e 's/^InputEigenVec.*/InputEigenVec   1/' \
    calcmod.def > _t && mv _t calcmod.def
sed -e "s/^Lanczos_max.*/Lanczos_max    ${TE_STEPS}/" \
    modpara.def > _t && mv _t modpara.def
echo "   TEOneBody  teonebody.def" >> namelist.def

echo "  [mpi] Time evolution..."
${MPIRUN} ../../${HPHI} -e namelist.def > /dev/null 2>&1
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
    echo "MPI consistency test PASSED for Hubbard TE TEOneBody"
    exit 0
else
    echo "MPI consistency test FAILED for Hubbard TE TEOneBody"
    echo "Expected difference < ${TOLERANCE}, got ${diff}"
    echo "--- Serial ---"
    cat ${TOPDIR}/serial/output/Flct.dat
    echo "--- MPI ---"
    cat ${TOPDIR}/mpi/output/Flct.dat
    exit 1
fi
