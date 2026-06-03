#!/bin/sh
#
# Regression test for the grand-canonical t-J model (CalcModel = tJGC = 10)
# in expert mode.
#
# Guards two bugs that were fixed together:
#
#   BUG 1 (serial / default CalcHS): the default CalcHS=1 (read_hacker=1) routed
#         tJGC through the snoob-based omp_sz_hacker() path, which can only
#         enumerate a fixed particle-number sector and produced imax=1 -> the run
#         aborted with "Error: in sz". This test deliberately does NOT set CalcHS
#         (i.e. it uses the default) so a regression re-triggers the abort.
#
#   BUG 2 (MPI / invalid ranks): a process whose inter-process configuration
#         contains a doublon "11" must have local dimension 0. The old guard
#         "if (Nup < 0 ...)" in check.c was dead code because Nup/Ndown/Ne are
#         unsigned, so invalid ranks kept full 3^Nlocal dimension and the total
#         dimension / energy were wrong (e.g. 108 instead of 81 for Nsite=4,np=4).
#
# Because the CI matrix runs the whole suite at np = 1, 4, 16 (all valid 4^m for
# tJGC), this single test exercises:
#   np = 1  -> BUG 1 (serial, default CalcHS must not abort)
#   np = 4  -> BUG 2 (one inter-process site; the "11" rank must be zero-dim)
#   np = 16 -> BUG 2 (two inter-process sites; odd/even split)
#
# Physical setup: open chain, nearest-neighbour hopping t = 1, no interactions
# (J = 0). Expected Hilbert-space dimension is 3^Nsite. Reference ground-state
# energies are the converged Lanczos values; MPI runs must reproduce them.

testname="lanczos_tjgc_dimension"
hphi="../../../src/HPhi"   # cwd will be <build>/test/${testname}/N<Nsite>
tol="0.00000001"

mkdir -p ${testname}
cd ${testname}

fail() {
  echo "FAILED (${testname}): $1"
  exit 1
}

# Generate expert-mode input files for a tJGC open chain of <Nsite> sites.
# Note: NO "CalcHS" line is written, so the default (read_hacker=1) is used.
gen_input() {
  n=$1
  cat > namelist.def <<EOF
         ModPara  modpara.def
         LocSpin  locspn.def
           Trans  trans.def
    CoulombIntra  coulombintra.def
        OneBodyG  greenone.def
        TwoBodyG  greentwo.def
         CalcMod  calcmod.def
EOF

  cat > calcmod.def <<EOF
CalcType   0
CalcModel   10
ReStart   0
CalcSpec   0
CalcEigenVec   0
InitialVecType   0
InputEigenVec   0
OutputEigenVec   0
InputHam   0
OutputHam   0
OutputExVec   0
EOF

  # Grand-canonical: do NOT define 2Sz or Ncond. Fixed initial vector for
  # reproducibility. No CalcHS (default path is exercised on purpose).
  cat > modpara.def <<EOF
--------------------
Model_Parameters   0
--------------------
HPhi_Cal_Parameters
--------------------
CDataFileHead  zvo
CParaFileHead  zqp
--------------------
Nsite          ${n}
Lanczos_max    2000
initial_iv     1
exct           1
LanczosEps     14
LanczosTarget  2
LargeValue     12.0
NumAve         5
ExpecInterval  20
EOF

  # All sites itinerant.
  {
    echo "================================"
    echo "NlocalSpin     0"
    echo "================================"
    echo "========i_1LocSpn_0IteElc ======"
    echo "================================"
    i=0
    while [ ${i} -lt ${n} ]; do
      printf "    %d      0\n" ${i}
      i=$((i + 1))
    done
  } > locspn.def

  # Nearest-neighbour hopping t = 1 on an open chain, both spins, both directions.
  tmptrans=$(mktemp)
  nt=0
  b=0
  while [ ${b} -lt $((n - 1)) ]; do
    j=$((b + 1))
    for s in 0 1; do
      printf "    %d     %d     %d     %d    1.000000000000000    0.000000000000000\n" ${b} ${s} ${j} ${s} >> ${tmptrans}
      printf "    %d     %d     %d     %d    1.000000000000000    0.000000000000000\n" ${j} ${s} ${b} ${s} >> ${tmptrans}
      nt=$((nt + 2))
    done
    b=$((b + 1))
  done
  {
    echo "========================"
    echo "NTransfer      ${nt}"
    echo "========================"
    echo "========i_j_s_tijs======"
    echo "========================"
    cat ${tmptrans}
  } > trans.def
  rm -f ${tmptrans}

  printf "=============================================\nNCoulombIntra          0\n=============================================\n================== CoulombIntra ================\n=============================================\n" > coulombintra.def
  printf "===========\nNCisAjs          0\n===========\n===========\n===========\n" > greenone.def
  printf "===========\nNCisAjsCktAlt          0\n===========\n===========\n===========\n" > greentwo.def
}

# Nsite, expected dimension (3^Nsite), reference ground-state energy.
check_case() {
  n=$1
  expdim=$2
  refene=$3

  rm -rf N${n}
  mkdir -p N${n}
  cd N${n}
  gen_input ${n}

  rm -rf output
  mkdir -p output
  ${MPIRUN} ${hphi} -e namelist.def > run.log 2>&1
  rc=$?

  if grep -q "Error: in sz" run.log; then
    cat run.log
    fail "Nsite=${n}: 'Error: in sz' (sz() enumeration inconsistent with idim_max)"
  fi
  if [ ${rc} -ne 0 ]; then
    cat run.log
    fail "Nsite=${n}: HPhi exited with status ${rc}"
  fi

  totdim=$(grep "Total dimension :" run.log | tail -1 | awk '{print $NF}')
  if [ "x${totdim}" != "x${expdim}" ]; then
    cat run.log
    fail "Nsite=${n}: total dimension ${totdim} != expected ${expdim}"
  fi

  if [ ! -f output/zvo_energy.dat ]; then
    cat run.log
    fail "Nsite=${n}: output/zvo_energy.dat not produced"
  fi
  ene=$(awk '/^Energy/{print $2; exit}' output/zvo_energy.dat)
  ok=$(awk -v e="${ene}" -v r="${refene}" -v t="${tol}" \
       'BEGIN{d=e-r; if(d<0)d=-d; print (d<t)?"yes":"no"}')
  if [ "${ok}" != "yes" ]; then
    fail "Nsite=${n}: ground energy ${ene} differs from reference ${refene} (tol ${tol})"
  fi

  echo "  OK  Nsite=${n}: dim=${totdim} (=3^${n}), E=${ene}"
  cd ..
}

echo "tJGC dimension/energy test (MPIRUN='${MPIRUN}')"
check_case 3 27 -1.4142135623730954
check_case 4 81 -2.2360679774997894

echo "PASSED (${testname})"
exit 0
