#!/bin/sh -e

# Regression test for CalcHS=2 on HubbardNConserved. Uses ncond (no 2Sz
# specified) to route through the HubbardNConserved branch of
# sz_hacker_for_large_systems. Guards against reintroducing the
# printf-in-main-loop spam and the snoob-overshoot off-by-one that
# used to write beyond list_1_.

mkdir -p lanczos_hubbard_square_ncond_calchs2/
cd lanczos_hubbard_square_ncond_calchs2

cat > stan.in <<EOF
W = 4
L = 2
model = "FermionHubbard"
method = "Lanczos"
lattice = "Tetragonal"
t = 1.0
U = 4.0
ncond = 8
outputmode = "all"
EOF

${MPIRUN} ../../src/HPhi -s stan.in
echo "CalcHS         2" >> modpara.def

rm -rf output
mkdir -p output
${MPIRUN} ../../src/HPhi -e namelist.def > run.log 2>&1

cat > reference.dat <<EOF
  -10.2529529552636234
    1.0444442739281228
    0.0000000000000000
EOF
paste output/zvo_energy.dat reference.dat > paste1.dat
diff=`awk 'BEGIN{diff=0.0} {diff+=sqrt(($2-$3)*($2-$3))} END{printf "%8.6f", diff}' paste1.dat`

# Guard against the raw "i=<num> i_cnt=<num>" debug printf that used to
# spam stdout 1 line per state in the HubbardNConserved branch.
printf_count=`grep -c "^i=[0-9]" run.log || true`
if [ "${printf_count}" != "0" ]; then
  echo "ERROR: raw debug printf leaked into stdout (${printf_count} lines)" >&2
  exit 1
fi

test "${diff}" = "0.000000"
exit $?
