#!/bin/sh
#
# Regression test for the J (spin-exchange) term of the grand-canonical t-J
# model (CalcModel = tJGC = 10) in expert mode.
#
# The lanczos_tjgc_dimension test only uses hopping (J = 0); it cannot catch a
# sign / convention / no-double-occupancy bug in the exchange term, which is the
# defining interaction of t-J and goes through the InterAll two-body kernel.
#
# The exchange is written via InterAll as J * (S_i.S_j - n_i n_j / 4):
#
#   J/2 * S^+_i S^-_j  : c^dag_{i up}  c_{i dn}  c^dag_{j dn}  c_{j up}
#   J/2 * S^-_i S^+_j  : c^dag_{i dn}  c_{i up}  c^dag_{j up}  c_{j dn}
#  -J/2 * n_{i up} n_{j dn}
#  -J/2 * n_{i dn} n_{j up}
#
# (S^z S^z and -n n/4 combine to the two diagonal -J/2 n n terms.)
#
# Part 1 - analytic lock (serial):
#   On a single bond at half filling (1 up + 1 down) every hop is Pauli/doublon
#   blocked, so only the exchange acts and the ground state is the singlet with
#   energy exactly -J, independent of t. With J = 2 the grand-canonical ground
#   state is this singlet (-2 < -t = -1), so E_gs = -2 locks the exchange sign
#   and magnitude. Run serially: a 2-site system cannot be split over np = 16.
#
# Part 2 - MPI consistency with J (via ${MPIRUN}):
#   A 4-site chain with t = 1, J = 1 must give the same total dimension (3^4=81)
#   and ground-state energy at np = 1, 4, 16, exercising the InterAll two-body
#   term through the MPI (zero-dim invalid rank) path.

testname="lanczos_tjgc_exchange"
hphi="../../../src/HPhi"
tol="0.00000001"

# Part 1 analytic reference: E = -J (singlet) for J = 2.
J1="2.0"
ref1="-2.0"
# Part 2 reference (converged Lanczos, np-independent) for 4-site chain, J = 1.
J2="1.0"
ref2="-2.6688977799503588"

mkdir -p ${testname}
cd ${testname}

fail() { echo "FAILED (${testname}): $1"; exit 1; }

# Generate expert-mode tJGC input for an open chain of <Nsite> sites with
# nearest-neighbour hopping t = 1 and exchange J = <J> on every bond.
gen_input() {
  n=$1
  jj=$2
  jq=$(awk -v j="${jj}" 'BEGIN{printf "%.15f", j/4.0}')
  jh=$(awk -v j="${jj}" 'BEGIN{printf "%.15f", j/2.0}')

  cat > namelist.def <<EOF
         ModPara  modpara.def
         LocSpin  locspn.def
           Trans  trans.def
        InterAll  interall.def
        OneBodyG  greenone.def
        TwoBodyG  greentwo.def
         CalcMod  calcmod.def
EOF
  printf "CalcType   0\nCalcModel   10\nReStart   0\nCalcSpec   0\nCalcEigenVec   0\nInputEigenVec   0\nOutputEigenVec   0\n" > calcmod.def
  {
    printf -- "--------------------\nModel_Parameters   0\n--------------------\nHPhi_Cal_Parameters\n--------------------\n"
    printf "CDataFileHead  zvo\nCParaFileHead  zqp\n--------------------\n"
    printf "Nsite          %d\nLanczos_max    2000\ninitial_iv     1\nexct           1\nLanczosEps     14\nLanczosTarget  2\nLargeValue     30.0\nNumAve         5\nExpecInterval  20\n" "${n}"
  } > modpara.def
  {
    printf "================================\nNlocalSpin     0\n================================\n========i_1LocSpn ======\n================================\n"
    i=0; while [ ${i} -lt ${n} ]; do printf "    %d      0\n" ${i}; i=$((i + 1)); done
  } > locspn.def

  tt=$(mktemp); ii=$(mktemp); nt=0; ni=0; b=0
  while [ ${b} -lt $((n - 1)) ]; do
    j=$((b + 1))
    for s in 0 1; do
      printf "    %d     %d     %d     %d    1.000000000000000    0.0\n" ${b} ${s} ${j} ${s} >> ${tt}
      printf "    %d     %d     %d     %d    1.000000000000000    0.0\n" ${j} ${s} ${b} ${s} >> ${tt}
      nt=$((nt + 2))
    done
    printf "    %d 0 %d 1 %d 1 %d 0    %s  0.0\n" ${b} ${b} ${j} ${j} "${jq}" >> ${ii}   # J/4 S+_i S-_j
    printf "    %d 0 %d 1 %d 1 %d 0    %s  0.0\n" ${j} ${j} ${b} ${b} "${jq}" >> ${ii}   # Hermitian pair
    printf "    %d 1 %d 0 %d 0 %d 1    %s  0.0\n" ${b} ${b} ${j} ${j} "${jq}" >> ${ii}   # J/4 S-_i S+_j
    printf "    %d 1 %d 0 %d 0 %d 1    %s  0.0\n" ${j} ${j} ${b} ${b} "${jq}" >> ${ii}   # Hermitian pair
    printf "    %d 0 %d 0 %d 1 %d 1    -%s 0.0\n" ${b} ${b} ${j} ${j} "${jh}" >> ${ii}   # -J/2 n_iup n_jdn
    printf "    %d 1 %d 1 %d 0 %d 0    -%s 0.0\n" ${b} ${b} ${j} ${j} "${jh}" >> ${ii}   # -J/2 n_idn n_jup
    ni=$((ni + 6)); b=$((b + 1))
  done
  { printf "========================\nNTransfer      %d\n========================\n========i_j_s_tijs======\n========================\n" ${nt}; cat ${tt}; } > trans.def
  { printf "======================\nNInterAll      %d\n======================\n========zInterAll=====\n======================\n" ${ni}; cat ${ii}; } > interall.def
  rm -f ${tt} ${ii}
  printf "===========\nNCisAjs          0\n===========\n===========\n===========\n" > greenone.def
  printf "===========\nNCisAjsCktAlt          0\n===========\n===========\n===========\n" > greentwo.def
}

energy_ok() { # $1 = measured, $2 = reference
  awk -v e="$1" -v r="$2" -v t="${tol}" 'BEGIN{d=e-r; if(d<0)d=-d; exit (d<t)?0:1}'
}

# ---- Part 1: analytic singlet E = -J (serial; 2 sites cannot split to np=16) ----
echo "Part 1: 2-site singlet, J=${J1}, expect E=${ref1}"
rm -rf part1; mkdir -p part1; cd part1
gen_input 2 ${J1}
rm -rf output; mkdir -p output
${hphi} -e namelist.def > run.log 2>&1 || { cat run.log; fail "Part 1: HPhi failed"; }
grep -q "Error: in sz" run.log && { cat run.log; fail "Part 1: Error in sz"; }
e1=$(awk '/^Energy/{print $2; exit}' output/zvo_energy.dat)
energy_ok "${e1}" "${ref1}" || fail "Part 1: E=${e1} != ${ref1} (=-J singlet)"
echo "  OK  Part 1: E=${e1} (=-J)"
cd ..

# ---- Part 2: 4-site chain with J, MPI consistency (via MPIRUN) ----
echo "Part 2: 4-site chain, J=${J2}, expect dim=81, E=${ref2} (MPIRUN='${MPIRUN}')"
rm -rf part2; mkdir -p part2; cd part2
gen_input 4 ${J2}
rm -rf output; mkdir -p output
${MPIRUN} ${hphi} -e namelist.def > run.log 2>&1 || { cat run.log; fail "Part 2: HPhi failed"; }
grep -q "Error: in sz" run.log && { cat run.log; fail "Part 2: Error in sz"; }
dim=$(grep "Total dimension :" run.log | tail -1 | awk '{print $NF}')
[ "x${dim}" = "x81" ] || { cat run.log; fail "Part 2: total dimension ${dim} != 81"; }
e2=$(awk '/^Energy/{print $2; exit}' output/zvo_energy.dat)
energy_ok "${e2}" "${ref2}" || fail "Part 2: E=${e2} != reference ${ref2}"
echo "  OK  Part 2: dim=${dim}, E=${e2}"
cd ..

echo "PASSED (${testname})"
exit 0
