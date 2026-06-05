#!/bin/sh -e
# Regression test for the calculate_jb_Spin_Old uninitialized-jb bug.
# CalcHS=0 routes Spin / SpinlessFermion basis construction through
# calculate_jb_Spin_Old(), which declared `jb` without initializing it before
# `list_jb[ib] = jb` (the sibling calculate_jb_Spin_Hacker sets jb=0). The
# garbage offset corrupted the sz-list, giving heap corruption / a crash in
# Lanczos. This asserts CalcHS=0 reproduces the default (CalcHS=1) ground-state
# energy for an 8-site SpinlessFermion chain.
# $1 = CMAKE_SOURCE_DIR (for test/testSpinlessCalc.py).

SRCDIR="$1"
HPHI=../../src/HPhi

mkdir -p lanczos_spinless_calchs0
cd lanczos_spinless_calchs0

# Generate an 8-site, 3-particle SpinlessFermion input (default CalcHS=1).
python3 "${SRCDIR}/test/testSpinlessCalc.py" -p "${HPHI}" -m SpinlessFermion -s 8 -n 3 > gen.log 2>&1

rm -rf output
${HPHI} -e namelist.def > log_calchs1.txt 2>&1
e1=$(awk 'NR==1{print $2}' output/zvo_energy.dat)

# Same input, basis built via the CalcHS=0 (calculate_jb_Spin_Old) path.
echo "CalcHS         0" >> modpara.def
rm -rf output
${HPHI} -e namelist.def > log_calchs0.txt 2>&1
e0=$(awk 'NR==1{print $2}' output/zvo_energy.dat)

echo "E(CalcHS=1)=${e1}  E(CalcHS=0)=${e0}"
awk -v a="${e1}" -v b="${e0}" 'BEGIN{
    if (b == "") { print "CalcHS=0 produced no ground-state energy"; exit 1 }
    d = a - b; if (d < 0) d = -d
    printf "max |E(CalcHS=0) - E(CalcHS=1)| = %.3e\n", d
    exit (d < 1e-10) ? 0 : 1
}' || { echo "MISMATCH"; exit 1; }
echo "SpinlessFermion CalcHS=0 == CalcHS=1: OK"
