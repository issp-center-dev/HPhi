#!/bin/sh -e

# Negative tests for CalcHS validation.
#
# Two kinds of failure modes are pinned here so we can't silently
# regress on either:
#
#   1. CalcHS=2 for any non-Hubbard model. sz_hacker_for_large_systems
#      is only wired up for Hubbard and HubbardNConserved; the sz()
#      entry guard must reject CalcHS=2 with a clear error message for
#      every other model (Kondo, KondoGC, KondoNConserved, Spin, ...),
#      instead of either crashing with an opaque imax != idim_max
#      abort or silently ignoring CalcHS=2 and running CalcHS=0/1.
#
#   2. CalcHS values other than 0/1 on Kondo (e.g. CalcHS=3 or -1).
#      The Kondo branch of sz() reads read_hacker and only handles
#      0/1; historically any other value fell through with a stale
#      icnt and produced the opaque "imax=1, idim_max=..." abort. The
#      Kondo branch must reject unknown values with a parameter error.
#
# This script intentionally runs serially even in MPI CI jobs. Its goal
# is to verify entry-guard validation and error messages, not MPI site
# decomposition. Running these negative cases through mpiexec needlessly
# exercises unrelated standard-mode setup paths; Kondo+ncond=2 at 16
# ranks can hang during `HPhi -s`, blocking the whole CI job before the
# CalcHS validation under test is even reached.

mkdir -p calchs2_unsupported_reject/
cd calchs2_unsupported_reject

run_reject_case () {
  name=$1
  subdir=$2
  stanfile=$3
  calchs_value=$4
  expected_msg=$5

  mkdir -p "${subdir}"
  (
    cd "${subdir}"
    printf "%s" "${stanfile}" > stan.in
    ../../../src/HPhi -sdry stan.in > gen.log 2>&1
    echo "CalcHS         ${calchs_value}" >> modpara.def
    rm -rf output
    mkdir -p output

    set +e
    ../../../src/HPhi -e namelist.def > run.log 2>&1
    rc=$?
    set -e

    if [ "${rc}" = "0" ]; then
      echo "[${name}] ERROR: CalcHS=${calchs_value} run succeeded but was expected to fail" >&2
      exit 1
    fi
    if ! grep -q "${expected_msg}" run.log; then
      echo "[${name}] ERROR: expected error message '${expected_msg}' not found" >&2
      tail -20 run.log >&2
      exit 1
    fi
  )
}

ENTRY_GUARD_MSG_CALCHS2="CalcHS=2 is only supported for Hubbard and HubbardNConserved"
ENTRY_GUARD_MSG_INVALID="is not a valid value for this model"

KONDO_IN='model = "Kondo"
method = "CG"
lattice = "chain"
L = 4
t = 1.0
J = 1.0
nelec = 4
2Sz = 0
exct = 1
'

KONDOGC_IN='model = "KondoGC"
method = "CG"
lattice = "chain"
L = 4
t = 1.0
J = 1.0
exct = 1
'

KONDO_NCOND_IN='model = "Kondo"
method = "CG"
lattice = "chain"
L = 4
t = 1.0
J = 1.0
ncond = 2
exct = 1
'

SPIN_IN='model = "Spin"
method = "CG"
lattice = "chain"
L = 4
J = 1.0
2Sz = 0
exct = 1
'

# (1) CalcHS=2 on non-Hubbard models must hit the "only supported for
#     Hubbard and HubbardNConserved" branch of the sz() entry guard.
run_reject_case "Kondo_calchs2"             kondo_calchs2             "${KONDO_IN}"       2  "${ENTRY_GUARD_MSG_CALCHS2}"
run_reject_case "KondoGC_calchs2"           kondogc_calchs2           "${KONDOGC_IN}"     2  "${ENTRY_GUARD_MSG_CALCHS2}"
run_reject_case "KondoNConserved_calchs2"   kondo_ncond_calchs2       "${KONDO_NCOND_IN}" 2  "${ENTRY_GUARD_MSG_CALCHS2}"
run_reject_case "Spin_calchs2"              spin_calchs2              "${SPIN_IN}"        2  "${ENTRY_GUARD_MSG_CALCHS2}"

# (2) CalcHS values that are invalid for the model (e.g. 3 or -1 on
#     Kondo / KondoGC / KondoNConserved) must hit the generic
#     "is not a valid value" branch of the entry guard. This pins the
#     validation so KondoGC / KondoNConserved cannot silently accept
#     garbage CalcHS values, and so Kondo cannot fall through to an
#     opaque "imax != idim_max" abort.
run_reject_case "Kondo_calchs3"             kondo_calchs3             "${KONDO_IN}"       3   "${ENTRY_GUARD_MSG_INVALID}"
run_reject_case "Kondo_calchs_neg1"         kondo_calchs_neg1         "${KONDO_IN}"       -1  "${ENTRY_GUARD_MSG_INVALID}"
run_reject_case "KondoGC_calchs3"           kondogc_calchs3           "${KONDOGC_IN}"     3   "${ENTRY_GUARD_MSG_INVALID}"
run_reject_case "KondoGC_calchs_neg1"       kondogc_calchs_neg1       "${KONDOGC_IN}"     -1  "${ENTRY_GUARD_MSG_INVALID}"
run_reject_case "KondoNConserved_calchs3"   kondo_ncond_calchs3       "${KONDO_NCOND_IN}" 3   "${ENTRY_GUARD_MSG_INVALID}"

exit 0
