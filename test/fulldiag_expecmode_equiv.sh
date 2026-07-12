#!/bin/sh
set -e

# fulldiag_expecmode_equiv (phase 3a Task 7): ExpecMode 0/1/(2) physics
# equivalence. "ExpecMode changes speed only" (Task 6 guarantee) is checked
# by running the SAME FullDiag problem with ExpecMode 0 (serial) and
# ExpecMode 1 (state-task-parallel), then diffing every zvo_phys_* /
# zvo_phys.dat column (S^2/Sz included -- no more excluding them, see the
# companion strengthening of fulldiag_elpa_hubbard_chain.sh) and every
# Green aggregate output file (per-state "*_eigen%d" files AND the
# aggregated "*_eigen" files) byte-for-byte-numerically (tolerance 1e-8).
#
# Registered as two SEPARATE named ctest cases pointed at this one script
# (see test/CMakeLists.txt): fulldiag_expecmode_equiv_np2 (exact:2) and
# fulldiag_expecmode_equiv_np3 (exact:3, exercising non-divisible
# state-panel ownership: N states over 3 ranks never divides evenly for
# the small Hilbert spaces used below). This script takes its process
# count entirely from ${MPIRUN} (set by the precheck wrapper); it never
# launches a second, different -np internally.
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

testname="fulldiag_expecmode_equiv"
hphi="../../src/HPhi"
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
# timings that differ run-to-run and carry no physics content.
compare_output_trees() {
  dirA="$1"
  dirB="$2"
  ( cd "${dirA}/output" && find . -type f ! -name '*TimeKeeper*' | sort ) > _filesA.lst
  ( cd "${dirB}/output" && find . -type f ! -name '*TimeKeeper*' | sort ) > _filesB.lst
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
    printf 'Solver 3\nNGPU 0\nExpecMode %d\n' "${mode}" >> calcmod.def
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

# ExpecMode 2 sanity sub-case (only exercised once, per spec: 3a always
# downgrades trace-kernel Mode 2 to Mode 1, so this both confirms the INFO
# message fires and that the result is bit-identical to Mode 1's).
run_mode "${case1}" 2
grep -q "ExpecMode 2 kernels are not available in this build; running as ExpecMode 1" \
  "${case1}/mode2/log_run.txt" || fail "ExpecMode 2 downgrade INFO message was not printed"
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
  11    0   11    1    2    0    2    1    3    0    3    0    4    0    4    1    1    0    1    1    0    0    0    0
  11    0   11    1    2    0    2    1    3    0    3    0    4    0    4    1    1    0    1    1    0    1    0    1
EOF
)
run_mode "${case4}" 0
run_mode "${case4}" 1
compare_output_trees "${case4}/mode0" "${case4}/mode1"

echo "fulldiag_expecmode_equiv: OK"
