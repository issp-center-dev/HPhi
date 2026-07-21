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

# Phase 3c Task 6: compare ONLY the zvo_phys_* / zvo_phys.dat file(s) between
# two run directories. Used by the InputHam negative case (case 7) and the
# spinless energy-demotion cases (cases 8/9), where the observable comparison
# must run on a DISTRIBUTED solver (ScaLAPACK for case 7; either solver for
# 8/9 depending on the ctest registration -- Solver 1 under
# fulldiag_expecmode_scalapack_np2). zvo_phys carries exactly the
# energy/fluctuation quantities (H,N,Sz,S2,D) the ExpecMode demotion paths
# recompute, so it is the right (and brief-specified) comparison target.
#
# DEGENERACY-AWARE COMPARISON (blocker-2 fix). Both files hold the identical
# energy spectrum in ascending-energy order (same H matrix, deterministic
# eigenVALUES), but a distributed eigensolver may return ANY orthonormal
# rotation within a degenerate eigenspace across two SEPARATE HPhi launches
# (ScaLAPACK's degenerate-subspace basis is not reproducible launch-to-launch,
# unlike a single deterministic ELPA run). Energies then match state-by-state
# but basis-dependent per-state quantities (S2, doublon, ...) need NOT. A
# naive per-line float diff therefore spuriously fails on degeneracy (the
# maintainer reproduced up to ~0.70 S2 differences on case 7 at ScaLAPACK
# np=3). So instead of comparing per line, we:
#   1. Group consecutive rows into DEGENERATE LEVELS: rows whose energy (col 1)
#      agrees with the previous row within etol (1e-6 absolute; energies are
#      O(1..10) here and real level gaps are O(0.1..1) >> etol). The grouping
#      is derived from file A and independently re-derived from file B; a
#      differing level count or boundary is a genuine spectral mismatch -> fail.
#   2. Compare the per-level SUM of each column between A and B. The per-level
#      sum of an observable is Tr(P.A.P) over the degenerate-eigenspace
#      projector P, which is basis-INVARIANT -- equal for any two valid
#      diagonalizations -- while individual per-state values are not. For a
#      non-degenerate level (size 1) this reduces to the exact per-state
#      comparison, unchanged. Tolerance is scaled to tol*max(1,level_size) so a
#      k-term floating sum is not held to a tighter bound than a single value.
# A real regression (e.g. Mode-2 energy kernel wrongly reactivating for a
# spinless model, shifting the <N> column by O(1)) still fails: that per-level
# SUM diverges by O(level_size) >> the scaled tolerance.
compare_phys() {
  dirA="$1"
  dirB="$2"
  found=0
  for fa in "${dirA}"/output/zvo_phys*.dat; do
    [ -e "${fa}" ] || continue
    found=1
    base=$(basename "${fa}")
    fb="${dirB}/output/${base}"
    [ -e "${fb}" ] || fail "compare_phys: ${fb} missing (present in ${dirA})"
    na=$(wc -l < "${fa}")
    nb=$(wc -l < "${fb}")
    [ "${na}" = "${nb}" ] || fail "compare_phys line count mismatch for ${base}: ${na} (${dirA}) vs ${nb} (${dirB})"
    paste "${fa}" "${fb}" | awk -v tol="${tol}" -v etol="1e-6" -v fname="${base}" '
      function abs(x){ return x<0?-x:x }
      {
        # Skip the header / any non-numeric-first-column line ("  <H> <N> ...").
        if ($1 !~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/) next
        if (NF % 2 != 0) { printf "MISMATCH %s line %d: uneven column count (NF=%d)\n", fname, NR, NF; bad=1; next }
        h = NF / 2
        if (half == 0) half = h
        else if (h != half) { printf "MISMATCH %s line %d: column count changed (%d vs %d)\n", fname, NR, h, half; bad=1; next }
        n++
        eA[n] = $1 + 0; eB[n] = $(1 + half) + 0
        for (i = 1; i <= half; i++) { A[n, i] = $i + 0; B[n, i] = $(i + half) + 0 }
      }
      END {
        if (bad) exit 1
        if (n == 0) { printf "MISMATCH %s: no numeric data rows found\n", fname; exit 1 }
        # Group ascending-energy rows into degenerate levels (from file A).
        gA = 1; startA[1] = 1
        for (k = 2; k <= n; k++) { if (abs(eA[k] - eA[k-1]) > etol) { gA++; startA[gA] = k } }
        for (g = 1; g < gA; g++) endA[g] = startA[g+1] - 1
        endA[gA] = n
        # Independently group file B and require identical boundaries.
        gB = 1; startB[1] = 1
        for (k = 2; k <= n; k++) { if (abs(eB[k] - eB[k-1]) > etol) { gB++; startB[gB] = k } }
        for (g = 1; g < gB; g++) endB[g] = startB[g+1] - 1
        endB[gB] = n
        if (gA != gB) { printf "MISMATCH %s: degenerate-level count differs (A=%d B=%d) -- different spectra\n", fname, gA, gB; exit 1 }
        for (g = 1; g <= gA; g++) {
          if (startA[g] != startB[g] || endA[g] != endB[g]) {
            printf "MISMATCH %s: level %d boundary differs (A rows %d..%d, B rows %d..%d) -- different spectra\n", fname, g, startA[g], endA[g], startB[g], endB[g]; exit 1
          }
        }
        # Per-level SUM comparison (basis-invariant Tr(P.A.P)); size-1 levels
        # reduce to the exact per-state check.
        for (g = 1; g <= gA; g++) {
          sz = endA[g] - startA[g] + 1
          lt = tol * (sz > 1 ? sz : 1)
          for (i = 1; i <= half; i++) {
            sa = 0; sb = 0
            for (k = startA[g]; k <= endA[g]; k++) { sa += A[k, i]; sb += B[k, i] }
            d = abs(sa - sb)
            if (d > lt) { printf "MISMATCH %s level %d (rows %d..%d) col %d: sumA=%.10g vs sumB=%.10g (diff %.3e, tol %.3e)\n", fname, g, startA[g], endA[g], i, sa, sb, d, lt; bad = 1 }
          }
        }
        exit (bad ? 1 : 0)
      }
    ' || fail "ExpecMode zvo_phys mismatch in ${base} (see MISMATCH lines above)"
  done
  [ "${found}" = "1" ] || fail "compare_phys: no zvo_phys*.dat found in ${dirA}/output"
}

# Blocker-2: the degeneracy-SAFE subset of compare_output_trees -- assert the
# output FILE SET is identical between two run directories and that the
# aggregate _eigen.dat merge output was produced, WITHOUT diffing per-state /
# per-eigen values (those are basis-dependent and NOT degeneracy-safe on a
# non-deterministic distributed solver). Used by the spinless cases 8/9, whose
# near-diagonal, highly degenerate Hamiltonians make the per-eigen Green-file
# diff of compare_output_trees fragile under ScaLAPACK; their observable
# agreement is asserted separately via the degeneracy-aware compare_phys. File
# existence (unlike eigenvector basis) does not depend on the degenerate
# rotation, so this check stays valid on any solver.
compare_output_presence() {
  dirA="$1"
  dirB="$2"
  ( cd "${dirA}/output" && find . -type f ! -name '*TimeKeeper*' ! -name 'CalcTimer.dat' | sort ) > _pfilesA.lst
  ( cd "${dirB}/output" && find . -type f ! -name '*TimeKeeper*' ! -name 'CalcTimer.dat' | sort ) > _pfilesB.lst
  [ -s _pfilesA.lst ] || fail "compare_output_presence: ${dirA}/output produced zero output files -- broken run or over-eager exclude filter"
  [ -s _pfilesB.lst ] || fail "compare_output_presence: ${dirB}/output produced zero output files -- broken run or over-eager exclude filter"
  grep -Eq '_eigen\.dat$' _pfilesA.lst || fail "compare_output_presence: ${dirA}/output has no aggregate Green/eigen file (*_eigen.dat) -- the aggregate/merge path did not run"
  diff _pfilesA.lst _pfilesB.lst > /dev/null || {
    echo "Output file sets differ between ${dirA}/output and ${dirB}/output:" >&2
    diff _pfilesA.lst _pfilesB.lst >&2
    fail "compare_output_presence: output file set mismatch between ${dirA} and ${dirB}"
  }
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
  assert_energy_kernel "${logfile}"
  grep -q "always use the ExpecMode-1 path" \
    "${logfile}" || fail "ExpecMode 2 always-fallback plan INFO line was not printed in ${logfile}"
}

# Phase 3c Task 6: the energy/fluctuation family now runs through its OWN CSR
# trace kernel (TRACE_Q_ENERGY, src/expec_trace.c TraceReportPlan()), eligible
# for EVERY makeHam-reachable model -- broader than the GF capability table.
# Every case below builds its Hamiltonian through makeHam (no InputHam) with a
# tiny Hilbert space, so the energy slot always ends kernel-active (never
# memory-demoted); this asserts that verbatim INFO line. Kept as a separate
# helper so both assert_kernel_plan() (GF kernel cases) and
# assert_kernel_plan_shared_evaluator() (case 4) reuse it, and so the tJ
# energy-only case (case 6) can assert the SAME energy line without the GF
# kernel lines. The substring is byte-exact from expec_trace.c.
assert_energy_kernel() {
  logfile="$1"
  grep -q "the energy/fluctuation family uses the trace kernel\." \
    "${logfile}" || fail "ExpecMode 2 did not select the trace kernel for the energy/fluctuation family in ${logfile}"
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
  assert_energy_kernel "${logfile}"
  grep -q "always use the ExpecMode-1 path" \
    "${logfile}" || fail "ExpecMode 2 always-fallback plan INFO line was not printed in ${logfile}"
}

# Phase 3c Task 6: taxonomy-orthogonality assertion for a model that IS
# energy-kernel-eligible (makeHam-reachable) but is NOT in the GF capability
# table (kTraceCap: tJ/tJGC/Kondo/KondoGC/general-spin rows are FALSE). The
# energy family selects the trace kernel while BOTH GF quantities take the
# "unsupported model" ExpecMode-1 fallback -- the energy kernel and the GF
# kernels are dispatched independently (TraceBuildPlan() gates the GF slots on
# cap_q[] but the energy slot only on InputHam). Substrings byte-exact from
# expec_trace.c's TraceReportPlan().
assert_energy_kernel_gf_unsupported() {
  logfile="$1"
  assert_energy_kernel "${logfile}"
  grep -q "ExpecMode 2: one-body Green functions use the ExpecMode-1 fallback (unsupported model)\." \
    "${logfile}" || fail "ExpecMode 2 one-body GF did not report the unsupported-model fallback in ${logfile}"
  grep -q "ExpecMode 2: two-body Green functions use the ExpecMode-1 fallback (unsupported model)\." \
    "${logfile}" || fail "ExpecMode 2 two-body GF did not report the unsupported-model fallback in ${logfile}"
}

# Blocker-1 (energy-kernel CORRECTNESS fix): a model that is NOT energy-kernel-
# eligible -- SpinlessFermion / SpinlessFermionGC. Before the fix TraceBuildPlan()
# enabled the energy trace kernel for EVERY non-InputHam model, so canonical
# SpinlessFermion (which trace_model_n_diag() maps to n_diag==0, same as
# canonical Spin) silently activated the kernel and wrote canonical-Spin
# CONSTANT fluctuation fields (num=NsiteMPI, Sz=0.5*Total2SzMPI) -- WRONG for a
# spinless model. TraceModelEnergySupported() now whitelists only the models the
# kernel actually implements, so the energy family here MUST report the
# "unsupported model" ExpecMode-1 fallback (kernel demoted, NOT active), exactly
# like both GF quantities already do. Substrings byte-exact from
# expec_trace.c's TraceReportPlan(). The energy line is the direct regression
# guard: if the blocker regressed, this run would instead print "the
# energy/fluctuation family uses the trace kernel." and the grep below would
# fail. The Mode-0/1/2 zvo_phys equivalence (compare below) is the value guard:
# with the kernel wrongly active, Mode 2's <N> column would read NsiteMPI
# instead of the fallback value, diverging from Mode 0/1.
assert_energy_and_gf_unsupported() {
  logfile="$1"
  grep -q "the energy/fluctuation family uses the ExpecMode-1 fallback (unsupported model)\." \
    "${logfile}" || fail "ExpecMode 2 did not report the unsupported-model energy fallback in ${logfile}"
  grep -q "ExpecMode 2: one-body Green functions use the ExpecMode-1 fallback (unsupported model)\." \
    "${logfile}" || fail "ExpecMode 2 one-body GF did not report the unsupported-model fallback in ${logfile}"
  grep -q "ExpecMode 2: two-body Green functions use the ExpecMode-1 fallback (unsupported model)\." \
    "${logfile}" || fail "ExpecMode 2 two-body GF did not report the unsupported-model fallback in ${logfile}"
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

# =========================================================================
# Case 6 (phase 3c Task 6): tJ chain L=4 -- taxonomy-orthogonality golden
# case. tJ is NOT a GF-kernel model (kTraceCap tJ row is FALSE, is_gc grouping
# mismatch) but IS energy-kernel-eligible (any model reaching makeHam is), so
# ExpecMode 2 selects the energy trace kernel while BOTH GF quantities take the
# "unsupported model" ExpecMode-1 fallback -- the orthogonality this case
# exists to prove.
#
# tJ is UNAVAILABLE in the Standard mode used by prep_case (StdFace has no tJ
# lattice generator -- src/StdFace/src/StdFace_main.c dispatches only
# hubbard/spin/kondo/spinlessfermion), so this case is written as an
# Expert-mode def set (CalcModel 9). It is a genuine tJ Hamiltonian:
# nearest-neighbour hopping t=1 on the periodic 4-site chain (both spins, both
# directions) plus the J=1 tJ exchange J*(S_i.S_j - n_i n_j/4) on every bond,
# written as InterAll with each transverse partner in HPhi's reversed
# conjugate-operator ordering (the form its NonHermite checker requires). The
# whole fixture -- def validity, makeHam reachability, and the aggregate
# zvo_cisajs_eigen.dat / zvo_cisajscktalt_eigen.dat production the
# compare_output_trees() _eigen.dat guard needs -- was validated locally with
# the serial LAPACK build (build_noMPI, ExpecMode 0); only the distributed
# ExpecMode 1/2 dispatch (model-independent, proven by cases 1-5) runs first on
# clavius. Ncond 2 on 4 sites = 2 holes; 2Sz 0 -> zvo_phys_Nup1_Ndown1.dat.
# All site indices stay strictly < Nsite=4 (the phase-3a out-of-range lesson).
# =========================================================================
case6="case6_tj_chain_energy_only"
mkdir -p "${case6}"
(
  cd "${case6}"
  cat > namelist.def <<EOF
         ModPara  modpara.def
         LocSpin  locspn.def
           Trans  trans.def
        InterAll  interall.def
        OneBodyG  greenone.def
        TwoBodyG  greentwo.def
         CalcMod  calcmod.def
EOF
  cat > calcmod.def <<EOF
CalcType        2
CalcModel       9
ReStart         0
CalcSpec        0
CalcEigenVec    0
InitialVecType  0
InputEigenVec   0
OutputEigenVec  0
InputHam        0
OutputHam       0
OutputGreenFormat 1
EOF
  cat > modpara.def <<EOF
--------------------
Model_Parameters   0
--------------------
HPhi_Cal_Parameters
--------------------
CDataFileHead  zvo
CParaFileHead  zqp
--------------------
Nsite             4
Ncond             2
2Sz               0
Lanczos_max       2000
initial_iv        -1
exct              1
LanczosEps        14
LanczosTarget     2
LargeValue        12.0
NumAve            5
ExpecInterval     20
EOF
  cat > locspn.def <<EOF
================================
NlocalSpin     0
================================
========i_1LocSpn_0IteElc ======
================================
    0      0
    1      0
    2      0
    3      0
EOF
  cat > trans.def <<EOF
========================
NTransfer      16
========================
========i_j_s_tijs======
========================
0 0 1 0 -1.0 0.0
1 0 0 0 -1.0 0.0
0 1 1 1 -1.0 0.0
1 1 0 1 -1.0 0.0
1 0 2 0 -1.0 0.0
2 0 1 0 -1.0 0.0
1 1 2 1 -1.0 0.0
2 1 1 1 -1.0 0.0
2 0 3 0 -1.0 0.0
3 0 2 0 -1.0 0.0
2 1 3 1 -1.0 0.0
3 1 2 1 -1.0 0.0
3 0 0 0 -1.0 0.0
0 0 3 0 -1.0 0.0
3 1 0 1 -1.0 0.0
0 1 3 1 -1.0 0.0
EOF
  cat > interall.def <<EOF
======================
NInterAll      16
======================
========zInterAll=====
======================
0 0 0 0 1 1 1 1 -0.5 0.0
0 1 0 1 1 0 1 0 -0.5 0.0
0 0 0 1 1 1 1 0 0.5 0.0
1 0 1 1 0 1 0 0 0.5 0.0
1 0 1 0 2 1 2 1 -0.5 0.0
1 1 1 1 2 0 2 0 -0.5 0.0
1 0 1 1 2 1 2 0 0.5 0.0
2 0 2 1 1 1 1 0 0.5 0.0
2 0 2 0 3 1 3 1 -0.5 0.0
2 1 2 1 3 0 3 0 -0.5 0.0
2 0 2 1 3 1 3 0 0.5 0.0
3 0 3 1 2 1 2 0 0.5 0.0
0 0 0 0 3 1 3 1 -0.5 0.0
0 1 0 1 3 0 3 0 -0.5 0.0
0 0 0 1 3 1 3 0 0.5 0.0
3 0 3 1 0 1 0 0 0.5 0.0
EOF
  cat > greenone.def <<EOF
========================
NCisAjs 8
========================
========GreenOne========
========================
0 0 0 0
0 1 0 1
1 0 1 0
1 1 1 1
0 0 1 0
1 0 0 0
2 0 2 0
2 1 2 1
EOF
  cat > greentwo.def <<EOF
========================
NCisAjsCktAlt 3
========================
========GreenTwo========
========================
0 0 0 0 1 1 1 1
0 0 0 0 2 0 2 0
1 0 1 0 2 1 2 1
EOF
)
run_mode "${case6}" 0
run_mode "${case6}" 1
compare_output_trees "${case6}/mode0" "${case6}/mode1"
run_mode "${case6}" 2
assert_energy_kernel_gf_unsupported "${case6}/mode2/log_run.txt"
compare_output_trees "${case6}/mode1" "${case6}/mode2"

# =========================================================================
# Case 7 (phase 3c Task 6): InputHam negative case -- the energy/fluctuation
# family DEMOTES to the ExpecMode-1 path when the Hamiltonian is read from a
# file (re-running makeHam in the CSR collector would build a matrix DIFFERENT
# from the one diagonalized; src/expec_trace.c TraceBuildPlan()'s InputHam
# static demotion). This is the ONLY static energy demotion.
#
# Solver choice (deliberate, not the script's ${solver}): this case forces
# Solver 1 (ScaLAPACK) at the ${MPIRUN} process count, because:
#   * ExpecMode 2 requires a distributed solver (ScaLAPACK/ELPA) AND nproc>1;
#     at nproc==1 readdef reverts ExpecMode 2->0 (src/readdef.c ~line 593), so
#     TraceReportPlan() never runs and the InputHam plan line never prints --
#     np=1 therefore CANNOT exercise this path.
#   * ELPA (Solver 3) + OutputHam/InputHam + nproc>1 is rejected at startup
#     (src/readdef.c cErrElpaHamIO: the ELPA panel is not the full replicated
#     matrix that inputHam()/outputHam() need). ScaLAPACK keeps the full
#     replicated Ham (src/lapack_diag.c diag_scalapack_cmp takes Ham[][]) and
#     is NOT gated, so it is the only solver on which InputHam + ExpecMode 2
#     coexist. USE_ELPA forces USE_SCALAPACK (top-level CMakeLists.txt), so
#     Solver 1 is compiled into the ELPA build too -- this case runs under both
#     the ELPA and ScaLAPACK ctest registrations.
# The OutputHam->InputHam file roundtrip and the InputHam==makeHam zvo_phys
# equality were validated locally with the serial LAPACK build (build_noMPI,
# Solver 0, np=1); only the ScaLAPACK/ExpecMode-2 dispatch (model-independent)
# runs first on clavius.
# =========================================================================
case7="case7_inputham_negative"
prep_case "${case7}" 'model = "Hubbard"
method = "FullDiag"
lattice = "chain"
L = 4
t = 1.0
U = 4.0
nelec = 4
2Sz = 0
outputmode = "correlation"'

# Step A: produce output/<head>_Ham.dat once (OutputHam=1 makes CalcByFullDiag
# return right after outputHam(), before any diagonalization/observables).
genham7="${case7}/genham"
mkdir -p "${genham7}"
cp "${case7}"/*.def "${genham7}/" 2>/dev/null || true
(
  cd "${genham7}"
  printf 'Solver 1\n' >> calcmod.def
  printf 'OutputHam 1\n' >> calcmod.def
  run_hphi log_run.txt ${MPIRUN} "${hphi}" -e namelist.def
)
[ -f "${genham7}/output/zvo_Ham.dat" ] \
  || fail "case7: OutputHam did not produce ${genham7}/output/zvo_Ham.dat (see cFileNamePhys_FullDiag_Ham)"

# Step B: Mode-0 reference -- normal makeHam FullDiag, ExpecMode 0, Solver 1.
ref7="${case7}/mode0"
mkdir -p "${ref7}"
cp "${case7}"/*.def "${ref7}/" 2>/dev/null || true
(
  cd "${ref7}"
  printf 'Solver 1\n' >> calcmod.def
  printf 'ExpecMode 0\n' >> calcmod.def
  run_hphi log_run.txt ${MPIRUN} "${hphi}" -e namelist.def
)

# Step C: Mode-2 InputHam run -- reads the Ham file, ExpecMode 2, Solver 1.
# The energy family must demote (InputHam); assert the verbatim reason line and
# that zvo_phys still equals the makeHam Mode-0 reference (the loaded matrix IS
# the one makeHam built and OutputHam wrote, so the physics is identical).
test7="${case7}/mode2"
mkdir -p "${test7}/output"
cp "${case7}"/*.def "${test7}/" 2>/dev/null || true
cp "${genham7}/output/zvo_Ham.dat" "${test7}/output/"
(
  cd "${test7}"
  printf 'Solver 1\n' >> calcmod.def
  printf 'InputHam 1\n' >> calcmod.def
  printf 'ExpecMode 2\n' >> calcmod.def
  run_hphi log_run.txt ${MPIRUN} "${hphi}" -e namelist.def
)
grep -q "the energy/fluctuation family uses the ExpecMode-1 fallback (the Hamiltonian was read from InputHam)\." \
  "${test7}/log_run.txt" \
  || fail "case7: ExpecMode 2 did not report the InputHam energy demotion in ${test7}/log_run.txt"
compare_phys "${ref7}" "${test7}"

# =========================================================================
# Case 8 (blocker-1): canonical SpinlessFermion (CalcModel 7) -- the energy
# trace kernel must be DEMOTED (unsupported model), not activated.
#
# SpinlessFermion is NOT energy-kernel-eligible: trace_model_n_diag() maps it
# (via its bare default:) to n_diag==0, the SAME as canonical Spin, so the old
# "enable for every non-InputHam model" logic silently activated the kernel and
# wrote canonical-Spin CONSTANT num/Sz fields -- wrong for a spinless model.
# TraceModelEnergySupported() now excludes it, so ExpecMode 2 falls back to the
# Mode-1 path (asserted below), and Mode 0/1/2 must agree.
#
# Written as an Expert-mode def set (like case 6): StdFace's `model =
# "spinlessfermion"` DOES emit a canonical spinless def set, but it ignores V
# (CoulombInter) and emits OneBodyG/TwoBodyG rows with spin index 1, which
# readdef rejects for a spinless model ("spin index must be 0") -- so a
# hand-authored set is cleaner and fully controlled. It is a genuine spinless
# Hamiltonian: t=1 nearest-neighbour hopping on the periodic 4-site chain
# (spin index 0 only) plus a V=2 nearest-neighbour CoulombInter (the only
# two-body term spinless supports besides hopping; no CoulombIntra/Hund/
# Exchange/PairHop). Ncond 2 on 4 sites -> zvo_phys_Nup2_Ndown0.dat. All site
# indices stay strictly < Nsite=4 (the phase-3a out-of-range lesson). Green
# defs use spin index 0 only (spinless), so aggregate zvo_cisajs_eigen.dat /
# zvo_cisajscktalt_eigen.dat are produced for the compare_output_trees
# _eigen.dat guard. The def validity, makeHam reachability, FullDiag completion
# and zvo_phys production were validated locally with the serial LAPACK build
# (build_noMPI, ExpecMode 0); only the distributed ExpecMode 1/2 dispatch
# (model-independent, proven by cases 1-7) runs first on clavius.
# =========================================================================
case8="case8_spinlessfermion_energy_demote"
mkdir -p "${case8}"
(
  cd "${case8}"
  cat > namelist.def <<EOF
         ModPara  modpara.def
         LocSpin  locspn.def
           Trans  trans.def
    CoulombInter  coulombinter.def
        OneBodyG  greenone.def
        TwoBodyG  greentwo.def
         CalcMod  calcmod.def
EOF
  cat > calcmod.def <<EOF
CalcType        2
CalcModel       7
ReStart         0
CalcSpec        0
CalcEigenVec    0
InitialVecType  0
InputEigenVec   0
OutputEigenVec  0
InputHam        0
OutputHam       0
OutputGreenFormat 1
EOF
  cat > modpara.def <<EOF
--------------------
Model_Parameters   0
--------------------
HPhi_Cal_Parameters
--------------------
CDataFileHead  zvo
CParaFileHead  zqp
--------------------
Nsite             4
Ncond             2
Lanczos_max       2000
initial_iv        -1
exct              1
LanczosEps        14
LanczosTarget     2
LargeValue        12.0
NumAve            5
ExpecInterval     20
EOF
  cat > locspn.def <<EOF
================================
NlocalSpin     0
================================
========i_1LocSpn_0IteElc ======
================================
    0      0
    1      0
    2      0
    3      0
EOF
  cat > trans.def <<EOF
========================
NTransfer      8
========================
========i_j_s_tijs======
========================
0 0 1 0 -1.0 0.0
1 0 0 0 -1.0 0.0
1 0 2 0 -1.0 0.0
2 0 1 0 -1.0 0.0
2 0 3 0 -1.0 0.0
3 0 2 0 -1.0 0.0
3 0 0 0 -1.0 0.0
0 0 3 0 -1.0 0.0
EOF
  cat > coulombinter.def <<EOF
========================
NCoulombInter 4
========================
========CoulombInter====
========================
0 1 2.0
1 2 2.0
2 3 2.0
3 0 2.0
EOF
  cat > greenone.def <<EOF
========================
NCisAjs 4
========================
========GreenOne========
========================
0 0 0 0
1 0 1 0
2 0 2 0
3 0 3 0
EOF
  cat > greentwo.def <<EOF
========================
NCisAjsCktAlt 3
========================
========GreenTwo========
========================
0 0 0 0 1 0 1 0
1 0 1 0 2 0 2 0
2 0 2 0 3 0 3 0
EOF
)
run_mode "${case8}" 0
run_mode "${case8}" 1
# Blocker-2: spinless FullDiag builds a near-diagonal, highly degenerate H, so
# the per-eigen Green-file value diff of compare_output_trees is NOT
# degeneracy-safe under the ScaLAPACK registration (Solver 1, non-deterministic
# degenerate basis across separate launches). Assert the file set / aggregate
# merge output with the degeneracy-safe presence check, and the Mode-0/1/2
# agreement of the degeneracy-INVARIANT observables with the degeneracy-aware
# compare_phys (per-level Tr(P.A.P) sums) instead.
compare_output_presence "${case8}/mode0" "${case8}/mode1"
compare_phys "${case8}/mode0" "${case8}/mode1"
run_mode "${case8}" 2
assert_energy_and_gf_unsupported "${case8}/mode2/log_run.txt"
compare_output_presence "${case8}/mode1" "${case8}/mode2"
compare_phys "${case8}/mode1" "${case8}/mode2"

# =========================================================================
# Case 9 (blocker-1): SpinlessFermionGC (CalcModel 8) -- same energy-kernel
# demotion, grand-canonical basis (2^Nsite). StdFace has no GC alias for
# spinless (src/StdFace/src/StdFace_main.c only maps "spinlessfermion"/
# "spinless" -> canonical CalcModel 7), so this MUST be an Expert-mode def set.
# Same t=1 hopping + V=2 CoulombInter as case 8; no Ncond/2Sz (GC), so the phys
# file is the GC-named zvo_phys.dat. 16 basis states. Same local validation
# story as case 8.
# =========================================================================
case9="case9_spinlessfermiongc_energy_demote"
mkdir -p "${case9}"
(
  cd "${case9}"
  cat > namelist.def <<EOF
         ModPara  modpara.def
         LocSpin  locspn.def
           Trans  trans.def
    CoulombInter  coulombinter.def
        OneBodyG  greenone.def
        TwoBodyG  greentwo.def
         CalcMod  calcmod.def
EOF
  cat > calcmod.def <<EOF
CalcType        2
CalcModel       8
ReStart         0
CalcSpec        0
CalcEigenVec    0
InitialVecType  0
InputEigenVec   0
OutputEigenVec  0
InputHam        0
OutputHam       0
OutputGreenFormat 1
EOF
  cat > modpara.def <<EOF
--------------------
Model_Parameters   0
--------------------
HPhi_Cal_Parameters
--------------------
CDataFileHead  zvo
CParaFileHead  zqp
--------------------
Nsite             4
Lanczos_max       2000
initial_iv        -1
exct              1
LanczosEps        14
LanczosTarget     2
LargeValue        12.0
NumAve            5
ExpecInterval     20
EOF
  cat > locspn.def <<EOF
================================
NlocalSpin     0
================================
========i_1LocSpn_0IteElc ======
================================
    0      0
    1      0
    2      0
    3      0
EOF
  cat > trans.def <<EOF
========================
NTransfer      8
========================
========i_j_s_tijs======
========================
0 0 1 0 -1.0 0.0
1 0 0 0 -1.0 0.0
1 0 2 0 -1.0 0.0
2 0 1 0 -1.0 0.0
2 0 3 0 -1.0 0.0
3 0 2 0 -1.0 0.0
3 0 0 0 -1.0 0.0
0 0 3 0 -1.0 0.0
EOF
  cat > coulombinter.def <<EOF
========================
NCoulombInter 4
========================
========CoulombInter====
========================
0 1 2.0
1 2 2.0
2 3 2.0
3 0 2.0
EOF
  cat > greenone.def <<EOF
========================
NCisAjs 4
========================
========GreenOne========
========================
0 0 0 0
1 0 1 0
2 0 2 0
3 0 3 0
EOF
  cat > greentwo.def <<EOF
========================
NCisAjsCktAlt 3
========================
========GreenTwo========
========================
0 0 0 0 1 0 1 0
1 0 1 0 2 0 2 0
2 0 2 0 3 0 3 0
EOF
)
run_mode "${case9}" 0
run_mode "${case9}" 1
# Blocker-2: same degeneracy-robust comparison as case 8 (spinless GC, 16
# basis states, near-diagonal highly degenerate H -- see case 8's note).
compare_output_presence "${case9}/mode0" "${case9}/mode1"
compare_phys "${case9}/mode0" "${case9}/mode1"
run_mode "${case9}" 2
assert_energy_and_gf_unsupported "${case9}/mode2/log_run.txt"
compare_output_presence "${case9}/mode1" "${case9}/mode2"
compare_phys "${case9}/mode1" "${case9}/mode2"

echo "fulldiag_expecmode_equiv: OK"
