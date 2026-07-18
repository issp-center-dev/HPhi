#!/bin/sh
set -e

solver="${HPHI_FULLDIAG_SOLVER:-3}"
case "${solver}" in
  1|3) ;;
  *) echo "ERROR: HPHI_FULLDIAG_SOLVER must be 1 or 3"; exit 1 ;;
esac

# fulldiag_expecmode_equiv (phase 3a Task 7): ExpecMode 0/1/(2) physics
# equivalence. "ExpecMode changes speed only" (Task 6 guarantee) is checked
# by running the SAME FullDiag problem with ExpecMode 0 (serial) and
# ExpecMode 1 (state-task-parallel), then diffing every zvo_phys_* /
# zvo_phys.dat column (S^2/Sz included -- no more excluding them, see the
# companion strengthening of fulldiag_elpa_hubbard_chain.sh) and every
# Green aggregate output file (per-state "*_eigen%d" files AND the
# aggregated "*_eigen" files) byte-for-byte-numerically (tolerance 1e-8).
#
# Registered as separate named ctest cases pointed at this one script (see
# test/CMakeLists.txt): ELPA at exact:2/exact:3 and ScaLAPACK at exact:2.
# The exact:3 case exercises non-divisible state-panel ownership: N states
# over 3 ranks never divides evenly for the small Hilbert spaces used below.
# This script takes its process count entirely from ${MPIRUN} (set by the
# precheck wrapper); it never launches a second, different -np internally.
#
# Requires Solver 1 (ScaLAPACK) or 3 (ELPA) for ExpecMode to be accepted at
# all (src/readdef.c's cErrExpecMode gate) -- so, like
# fulldiag_elpa_hubbard_chain.sh, this test is registration-only wherever
# neither library is compiled in; see test/CMakeLists.txt's if(USE_ELPA)
# guard and .superpowers/sdd/task-7-report.md for the local (non-ELPA,
# non-ScaLAPACK) verification story.
#
# ---------------------------------------------------------------------
# Task 1 inventory-audit conclusion, reflected in what this script does
# and does NOT attempt to cover (spec Sec.2's audit item):
#   - AnomalousG is REJECTED OUTRIGHT for any nproc>1 FullDiag run by
#     validate_anomalous_scope_common() (src/anomalous_pair.c: "Error:
#     AnomalousG does not support MPI FullDiag."), independent of
#     ExpecMode and unchanged by phase 3a (see
#     docs/superpowers/specs/2026-07-11-elpa-fulldiag-phase3-design.md
#     Sec.2: "AnomalousG 等が既存の検証で分散 FullDiag を拒否している場合、
#     その拒否は ExpecMode に関係なくそのまま有効"). It is therefore
#     impossible to reach with nproc>=2 regardless of ExpecMode, and is
#     intentionally NOT exercised here.
#   - NBodyG and ThreeBodyG/FourBodyG/SixBodyG ARE reachable in
#     distributed FullDiag: they are computed from inside
#     expec_cisajscktaltdc()/expec_nbodyg(), which the Mode 1 local loop
#     calls for every owned state exactly like Mode 0 does (see
#     src/phys_distributed.c's call chain). Design-doc Sec.3's "NBodyG /
#     AnomalousG ファイル: 常にフォールバック" describes NBodyG going
#     through the *same* Mode-1 local-loop path as every other quantity
#     (there being no Mode-2 trace kernel for it in 3a, not that it is
#     unreachable) -- Case 1 below exercises exactly that fallback path
#     plus the GreenOutputNBody aggregate kind. Case 4 below exercises
#     ThreeBody/FourBody/SixBody the same way (canonical "Spin" rejects
#     these -- src/expec_cisajscktaltdc.c:183-185 -- so Case 4 uses
#     SpinGC, matching that constraint).
# ---------------------------------------------------------------------
#
# Phase 3b Task 5 (golden-test hardening + capability-table finalization):
# kTraceCap (src/expec_trace.c) now has Hubbard/HubbardGC/Spin(half)/
# SpinGC(half) TRUE for one-body AND two-body GFs, so Cases 1-4's
# ExpecMode-2 sub-cases now assert the literal "use the trace kernel" INFO
# line for both quantities (assert_kernel_plan(), below) instead of a
# generic fallback-or-kernel match. Case 5 is new: a dedicated canonical
# (non-GC) Hubbard one-body+two-body golden case whose greenone.def/
# greentwo.def deliberately walk every reachable operator-family branch the
# §2c call-inventory audit lists for that (model, quantity) pair -- this is
# the only branch-coverage evidence for the canonical Hubbard basis, which
# is not locally unit-testable (test/unit/expec_trace_map_check.c is
# GC-only). The tJ/tJGC/Kondo/KondoGC rows stay FALSE (known is_gc grouping
# mismatch, recorded in kTraceCap's comment) and are not exercised here.
# ---------------------------------------------------------------------

testname="fulldiag_expecmode_equiv"
# Resolve the HPhi binary path ONCE, as an absolute path, before any `cd`.
# ctest's cwd when this script starts is build/test; prep_case later runs
# HPhi from build/test/${testname}/<case> (two levels deeper) and run_mode
# from build/test/${testname}/<case>/<mode> (three levels deeper) -- a
# single relative path cannot serve both depths, so we resolve it here
# instead of relying on a fixed number of "../" hops.
hphi="$(pwd)/../src/HPhi"
tol="0.00000001"

mkdir -p "${testname}"
cd "${testname}"

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

run_hphi() {
  log="$1"
  shift
  "$@" > "${log}" 2>&1 || {
    echo "Command failed: $*" >&2
    tail -100 "${log}" >&2
    exit 1
  }
}

# Compare every output file that exists under dirA/output against the
# same-named file under dirB/output, column-by-column, numeric columns to
# ${tol}, text columns (e.g. header rows' "<H>" labels) exactly. Both
# directories come from the identical problem run with only ExpecMode
# differing, so the file set and every file's line/column count must be
# identical -- any structural mismatch is itself a failure, not just a
# numeric one.
# TimeKeeper files are excluded: they legitimately record wall-clock
# timings that differ run-to-run and carry no physics content. CalcTimer.dat
# is excluded for the same reason (wall-clock timings, written on every
# MPI-build run; would spuriously fail the numeric diff at tol 1e-8).
#
# This comparison relies on Mode 1's rank-order part-file concatenation
# (GreenOutputMergePartials(), see src/include/green_output.h) equaling
# serial (Mode 0) state order: state ownership across ranks is a contiguous
# block (ceil(N/P) states per rank, rank 0 first), so concatenating parts in
# rank order reproduces ascending state order exactly, byte-for-byte
# comparable against Mode 0's own ascending-state-order output.
compare_output_trees() {
  dirA="$1"
  dirB="$2"
  ( cd "${dirA}/output" && find . -type f ! -name '*TimeKeeper*' ! -name 'CalcTimer.dat' | sort ) > _filesA.lst
  ( cd "${dirB}/output" && find . -type f ! -name '*TimeKeeper*' ! -name 'CalcTimer.dat' | sort ) > _filesB.lst
  [ -s _filesA.lst ] || fail "compare_output_trees: ${dirA}/output produced zero output files -- broken run or over-eager exclude filter"
  [ -s _filesB.lst ] || fail "compare_output_trees: ${dirB}/output produced zero output files -- broken run or over-eager exclude filter"
  grep -Eq '_eigen\.dat$' _filesA.lst || fail "compare_output_trees: ${dirA}/output has no aggregate Green/eigen file (*_eigen.dat, see GreenOutputFileName()'s _Eigen_Aggregate formats in src/green_output.c) -- the aggregate/merge path this test exists to cover did not actually run"
  diff _filesA.lst _filesB.lst > /dev/null || {
    echo "Output file sets differ between ${dirA}/output and ${dirB}/output:" >&2
    diff _filesA.lst _filesB.lst >&2
    fail "ExpecMode output file set mismatch"
  }
  while IFS= read -r f; do
    fa="${dirA}/output/${f}"
    fb="${dirB}/output/${f}"
    na=$(wc -l < "${fa}")
    nb=$(wc -l < "${fb}")
    [ "${na}" = "${nb}" ] || fail "line count mismatch for ${f}: ${na} (${dirA}) vs ${nb} (${dirB})"
    paste "${fa}" "${fb}" | awk -v tol="${tol}" -v fname="${f}" '
      {
        if (NF % 2 != 0) { printf "MISMATCH %s line %d: uneven column count (NF=%d)\n", fname, NR, NF; bad=1; next }
        half = NF / 2
        for (i = 1; i <= half; i++) {
          L = $i; R = $(i + half)
          if (L ~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/ && R ~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/) {
            d = L - R; if (d < 0) d = -d
            if (d > tol) { printf "MISMATCH %s line %d col %d: %s vs %s (diff %.3e)\n", fname, NR, i, L, R, d; bad = 1 }
          } else if (L != R) {
            printf "MISMATCH %s line %d col %d (text): \"%s\" vs \"%s\"\n", fname, NR, i, L, R; bad = 1
          }
        }
      }
      END { exit (bad ? 1 : 0) }
    ' || fail "ExpecMode 0/1 output mismatch in ${f} (see MISMATCH lines above)"
  done < _filesA.lst
}

# Phase 3b Task 5: assert that a mode2 run's log shows the plan that is
# expected now that the capability table (src/expec_trace.c's kTraceCap) has
# Hubbard/HubbardGC/Spin(half)/SpinGC(half) flipped TRUE for both one-body
# and two-body Green functions -- i.e. TraceReportPlan() must show the
# ACTUAL "use the trace kernel" line for both quantities (not merely a
# generic "use the ..." fallback-or-kernel match), plus the fixed
# always-fallback line for energy/fluctuation/S2/NBodyG/AnomalousG, which
# never changes regardless of capability. This replaces phase 3a/early-3b's
# generic grep (which only checked that *some* plan line was printed,
# because the table shipped all-FALSE back then).
assert_kernel_plan() {
  logfile="$1"
  grep -q "ExpecMode 2: one-body Green functions use the trace kernel\." \
    "${logfile}" || fail "ExpecMode 2 did not select the trace kernel for one-body GFs in ${logfile}"
  grep -q "ExpecMode 2: two-body Green functions use the trace kernel\." \
    "${logfile}" || fail "ExpecMode 2 did not select the trace kernel for two-body GFs in ${logfile}"
  grep -q "always use the ExpecMode-1 path" \
    "${logfile}" || fail "ExpecMode 2 always-fallback plan INFO line was not printed in ${logfile}"
}

# Variant for a run that ALSO defines ThreeBodyG/FourBodyG/SixBodyG (case4):
# expec_cisajscktaltdc() computes the two-body GF and the multibody GFs in
# ONE evaluator (expec_cisajscktaltdc.c:115 runs it when ANY of
# NCisAjtCkuAlvDC/NTBody/NFBody/NSBody > 0), so TraceBuildPlan() demotes the
# two-body GF to the ExpecMode-1 fallback whenever a multibody GF is defined
# (plan field demoted_shared_evaluator -- see src/include/expec_trace.h).
# The clavius production run of an earlier Task-5 revision caught exactly
# this: with the two-body kernel selected, mode2 silently LOST
# zvo_ThreeBody/FourBody/SixBody_eigen.dat (the fallback loop skips the
# whole shared evaluator). One-body is unaffected (expec_cisajs() has no
# NTBody/NFBody/NSBody reference), so its kernel line is still required.
assert_kernel_plan_shared_evaluator() {
  logfile="$1"
  grep -q "ExpecMode 2: one-body Green functions use the trace kernel\." \
    "${logfile}" || fail "ExpecMode 2 did not select the trace kernel for one-body GFs in ${logfile}"
  grep -q "ExpecMode 2: two-body Green functions use the ExpecMode-1 fallback (they share their evaluator with three-/four-/six-body Green functions)\." \
    "${logfile}" || fail "ExpecMode 2 did not report the shared-evaluator fallback for two-body GFs in ${logfile}"
  grep -q "always use the ExpecMode-1 path" \
    "${logfile}" || fail "ExpecMode 2 always-fallback plan INFO line was not printed in ${logfile}"
}

# Prepare a case directory: write stan.in, run `HPhi -sdry`, then let the
# caller add extra namelist/def files before calling run_mode(). Always
# turns on the aggregate Green output format (spec: reuse
# green_output_format.sh's method) because the merge/manifest code this
# task is about is only exercised via GreenOutputKindUsesAggregate()==true.
prep_case() {
  casedir="$1"
  standata="$2"
  mkdir -p "${casedir}"
  (
    cd "${casedir}"
    printf '%s\n' "${standata}" > stan.in
    run_hphi log_sdry.txt "${hphi}" -sdry stan.in
    printf 'OutputGreenFormat 1\n' >> calcmod.def
  )
}

# Run one ExpecMode value for a prepared case, into its own mode<N>/
# subdirectory (a full copy of the case's *.def files), via ${MPIRUN}.
# Returns (via the log file left behind) so the caller can grep it for the
# ExpecMode-2 downgrade INFO message when needed.
run_mode() {
  casedir="$1"
  mode="$2"
  moddir="${casedir}/mode${mode}"
  mkdir -p "${moddir}"
  cp "${casedir}"/*.def "${moddir}/" 2>/dev/null || true
  (
    cd "${moddir}"
    printf 'Solver %s\n' "${solver}" >> calcmod.def
    if [ "${solver}" = "3" ]; then
      printf 'NGPU 0\n' >> calcmod.def
    fi
    printf 'ExpecMode %d\n' "${mode}" >> calcmod.def
    run_hphi log_run.txt ${MPIRUN} "${hphi}" -e namelist.def
  )
}

# =========================================================================
# Case 1: Hubbard chain L=4 (one-body + two-body GF, aggregate ON), plus an
# NBodyG definition to exercise the Mode-1 fallback path and the
# GreenOutputNBody aggregate kind. Base Hamiltonian parameters match the
# established Hubbard-chain-L4 FullDiag fixture (fulldiag_hubbard_chain.sh /
# fulldiag_elpa_hubbard_chain.sh: t=1.0, U=4.0, nelec=4, 2Sz=0); the NBodyG
# def itself is the one already proven against this exact model/lattice/L
# combination by test/mpi_nbodyg_hubbard.sh (three operators: two diagonal
# density terms plus one two-site term), reused verbatim here rather than
# hand-rederived.
# =========================================================================
case1="case1_hubbard_nbodyg"
prep_case "${case1}" 'model = "Hubbard"
method = "FullDiag"
lattice = "chain"
L = 4
t = 1.0
U = 4.0
nelec = 4
2Sz = 0
outputmode = "correlation"'
(
  cd "${case1}"
  printf '    NBodyG  nbodyg.def\n' >> namelist.def
  cat > nbodyg.def <<EOF
========================
NNBodyG 3
========================
========NBodyG==========
========================
1 3 0 3 0
1 3 0 0 0
2 3 0 0 0 1 1 1 1
EOF
)
run_mode "${case1}" 0
run_mode "${case1}" 1
compare_output_trees "${case1}/mode0" "${case1}/mode1"

# ExpecMode 2 sanity sub-case. Phase 3b Task 1 replaced the old single
# downgrade-INFO grep with a check of every TraceReportPlan() line
# (src/expec_trace.c); Task 5 flipped kTraceCap's Hubbard row to TRUE for
# both one-body and two-body (golden evidence: this case's ExpecMode-2
# sub-case plus case5_hubbard_onebody_twobody_golden below, and the clavius
# forced-kernel checkpoint np=2/3 2026-07-12 -- see kTraceCap's per-row
# comment), so the plan now actually selects the trace kernel for both
# quantities -- assert_kernel_plan() requires the literal "trace kernel"
# text, not just any fallback-or-kernel INFO line. Mode 2's numeric result
# must still be identical to Mode 1's regardless of which path the plan
# picked (the "ExpecMode changes speed only" guarantee).
run_mode "${case1}" 2
assert_kernel_plan "${case1}/mode2/log_run.txt"
compare_output_trees "${case1}/mode1" "${case1}/mode2"

# =========================================================================
# Case 2: SpinGC, Gamma = 0.5 (transverse field), L = 6 chain.
# =========================================================================
case2="case2_spingc_gamma"
prep_case "${case2}" 'model = "SpinGC"
method = "FullDiag"
lattice = "chain"
L = 6
J = 1.0
Gamma = 0.5
outputmode = "correlation"'
run_mode "${case2}" 0
run_mode "${case2}" 1
compare_output_trees "${case2}/mode0" "${case2}/mode1"

# ExpecMode 2 sub-case: SpinGC (half) is TRUE in kTraceCap as of Task 5 (see
# its per-row evidence comment in src/expec_trace.c), so the plan must
# select the trace kernel here too.
run_mode "${case2}" 2
assert_kernel_plan "${case2}/mode2/log_run.txt"
compare_output_trees "${case2}/mode1" "${case2}/mode2"

# =========================================================================
# Case 3: Spin (canonical) chain, L = 8, 2Sz = 0.
# =========================================================================
case3="case3_spin_chain"
prep_case "${case3}" 'model = "Spin"
method = "FullDiag"
lattice = "chain"
L = 8
J = 1.0
2Sz = 0
outputmode = "correlation"'
run_mode "${case3}" 0
run_mode "${case3}" 1
compare_output_trees "${case3}/mode0" "${case3}/mode1"

# ExpecMode 2 sub-case: canonical Spin (half) is TRUE in kTraceCap as of
# Task 5 (see its per-row evidence comment in src/expec_trace.c) -- the
# canonical Spin basis (GetOffComp/list_1) is not locally unit-testable, so
# this case plus the clavius forced-kernel checkpoint np=2/3 2026-07-12 are
# the verification for that row.
run_mode "${case3}" 2
assert_kernel_plan "${case3}/mode2/log_run.txt"
compare_output_trees "${case3}/mode1" "${case3}/mode2"

# =========================================================================
# Case 4: ThreeBody/FourBody/SixBody aggregate-kind coverage. SpinGC on the
# small Honeycomb cell (W=2, L=2 -> 8 sites) already used by
# test/lanczos_spingc_hcor.sh (ThreeBodyG/FourBodyG defs) and
# test/lobcg_spingc_6body.sh (SixBodyG def); those defs are reused verbatim
# (they are pure index-list Green-function definitions, independent of the
# calculation method), only the method is switched to FullDiag here.
# Canonical "Spin" rejects Three/Four/SixBodyG (src/expec_cisajscktaltdc.c
# validation), which is why this case uses SpinGC rather than folding it
# into Case 3.
# =========================================================================
case4="case4_spingc_honeycomb_manybody"
prep_case "${case4}" 'W = 2
L = 2
model = "SpinGC"
method = "FullDiag"
lattice = "Honeycomb"
J0x = -1.0
J0y =  0.0
J0z =  0.0
J1x =  0.0
J1y = -1.0
J1z =  0.0
J2x =  0.0
J2y =  0.0
J2z = -1.0
2S=1
h=0
outputmode = "correlation"'
(
  cd "${case4}"
  printf '    ThreeBodyG  green3.def\n' >> namelist.def
  printf '    FourBodyG   green4.def\n' >> namelist.def
  printf '    SixBodyG    green6.def\n' >> namelist.def
  cat > green3.def <<EOF
===================
num       2
===================
===================
===================
        0         0           0           1           2           1           2           0           3           0           3           0
        0         0           0           1           2           1           2           0           3           1           3           1
EOF
  cat > green4.def <<EOF
===================
num       2
===================
===================
===================
        0  0  0  1  2  1  2  0  3  0  3  0  4  0  4  0
        0  0  0  1  2  1  2  0  3  1  3  1  4  1  4  1
EOF
  cat > green6.def <<EOF
===================
num       2
===================
===================
===================
   5    0    5    1    2    0    2    1    3    0    3    0    4    0    4    1    1    0    1    1    0    0    0    0
   5    0    5    1    2    0    2    1    3    0    3    0    4    0    4    1    1    0    1    1    0    1    0    1
EOF
)
run_mode "${case4}" 0
run_mode "${case4}" 1
compare_output_trees "${case4}/mode0" "${case4}/mode1"

# ExpecMode 2 sub-case: SpinGC (half, 2S=1 stays iFlgGeneralSpin==FALSE --
# src/readdef.c's `X->LocSpn[i]>LOCSPIN` check only sets general-spin when
# 2S>1) is TRUE in kTraceCap as of Task 5, so the ONE-body GF selects the
# trace kernel -- but the TWO-body GF must fall back here because this case
# defines ThreeBodyG/FourBodyG/SixBodyG, which share their Mode-1 evaluator
# with the two-body GF (the shared-evaluator demotion; see
# assert_kernel_plan_shared_evaluator()'s comment above). This case is
# therefore the golden evidence BOTH that kernel and always-fallback
# quantities coexist correctly in one run AND that the shared-evaluator
# demotion preserves the multibody output files (the compare below fails
# with missing zvo_ThreeBody/FourBody/SixBody_eigen.dat if it regresses).
run_mode "${case4}" 2
assert_kernel_plan_shared_evaluator "${case4}/mode2/log_run.txt"
compare_output_trees "${case4}/mode1" "${case4}/mode2"

# =========================================================================
# Case 5: canonical Hubbard chain L=4, one-body + two-body GF only (no
# NBodyG -- Case 1 already covers the NBodyG/always-fallback interaction).
# This is the phase-3b Task 5 golden case dedicated to branch-coverage of
# the canonical (non-GC) Hubbard trace-kernel dispatch, which is NOT locally
# unit-testable (test/unit/expec_trace_map_check.c is GC-only, per its file
# header -- the canonical basis needs GetOffComp/list_1). greenone.def
# exercises GetOffComp's reachable branches in one case: up-spin hop and
# down-spin hop, each in both the forward and reversed site order, plus the
# boundary site pair (site 0, site L-1=3), plus diagonal density terms for
# both spins. greentwo.def exercises all four canonical two-body element
# families the §2c audit lists for Hubbard (non-diagonal CisAjtCkuAlv,
# density-density-diagonal CisAisCisAis, and both "same-index" mixed
# families CisAisCjtAku / CisAjtCkuAku) plus the Sz-violating zero-result
# row (expec_cisajscktaltdc.c:728-734's 0.0 shortcut, mirrored by
# TraceMapExtractTwoBody()'s map.n==0 sentinel at expec_trace.c). Every site
# index stays strictly < Nsite=4 throughout (the phase-3a lesson: an
# out-of-range site index SIGFPEs).
# =========================================================================
case5="case5_hubbard_onebody_twobody_golden"
prep_case "${case5}" 'model = "Hubbard"
method = "FullDiag"
lattice = "chain"
L = 4
t = 1.0
U = 4.0
nelec = 4
2Sz = 0
outputmode = "correlation"'
(
  cd "${case5}"
  cat > greenone.def <<EOF
========================
NCisAjs 8
========================
========GreenOne========
========================
0 0 1 0
1 0 0 0
0 1 1 1
1 1 0 1
0 0 3 0
3 0 0 0
0 0 0 0
0 1 0 1
EOF
  cat > greentwo.def <<EOF
========================
NCisAjsCktAlt 5
========================
========GreenTwo========
========================
0 0 0 0 1 1 1 1
0 0 0 0 1 0 2 0
0 0 1 0 2 1 2 1
0 0 1 0 2 1 3 1
0 0 1 1 2 0 3 0
EOF
)
run_mode "${case5}" 0
run_mode "${case5}" 1
compare_output_trees "${case5}/mode0" "${case5}/mode1"
run_mode "${case5}" 2
assert_kernel_plan "${case5}/mode2/log_run.txt"
compare_output_trees "${case5}/mode1" "${case5}/mode2"

echo "fulldiag_expecmode_equiv: OK"
