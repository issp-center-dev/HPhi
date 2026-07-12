#!/bin/sh
# HPhi -- ExpecLocal call-inventory guard.
#
# Fails if any of the FullDiag expec_* observable-evaluation files below
# gain a new raw MPI-flavoured call (a real MPI_* library call, exitMPI, or
# any wrapperMPI.c helper) that is not on the frozen ExpecLocal allow-list
# and is not wrapped in an EXPEC_LOCAL_GUARDED_BEGIN/END marker region.
#
# Background: docs/superpowers/specs/2026-07-11-expec-call-inventory.md
# (phase 3a, Task 1) freezes which MPI-touching calls are reachable from the
# expec_* layer used by FullDiag observable evaluation. This script is the
# regression guard for that inventory. Task 3 inserts the defensive guards
# (`if (iExpecLocal) return ...;`) around the partner_rank != myrank branches
# in nbody_correlation.c / anomalous_pair.c, wraps them in
# EXPEC_LOCAL_GUARDED_BEGIN/END markers, and empties TEMP_UNGUARDED_FILES
# below. Task 6 adds the new MPI-free local-loop layer file to FILES.
#
# Usage: check_expec_local_calls.sh <CMAKE_SOURCE_DIR>

set -eu

srcdir="${1:-}"
if [ -z "$srcdir" ]; then
  echo "usage: check_expec_local_calls.sh <CMAKE_SOURCE_DIR>" >&2
  exit 1
fi

STRIPPER="$srcdir/test/strip_c_comments.py"
if [ ! -f "$STRIPPER" ]; then
  echo "ERROR: stripper not found at $STRIPPER" >&2
  exit 1
fi

# Files reachable from the FullDiag expec_* observable-evaluation layer,
# exactly as enumerated in docs/superpowers/specs/2026-07-11-expec-call-inventory.md
# Step 1's "expec 到達層" file set. Do NOT add the underlying mltply*Core.c
# element-function layer here -- those are shared, general-purpose MPI
# site-decomposition machinery used by Lanczos/TPQ as well, not something
# this expec-focused guard polices; their reachability was audited instead
# in the inventory document itself.
FILES="src/expec_energy_flct.c src/expec_cisajs.c src/expec_cisajscktaltdc.c src/expec_totalspin.c src/nbody_correlation.c src/anomalous_pair.c src/phys_distributed_local.c src/expec_trace.c"
# The MPI-free local-loop layer (Task 6) is scanned here; the orchestration
# layer src/phys_distributed.c holds the legitimate raw MPI collectives and is
# permanently NOT scanned.
# src/expec_trace.c (phase 3b Task 1: the ExpecMode 2 trace-kernel plan
# construction / INFO reporting / dispatch skeleton) is added to this list in
# the same commit that creates it -- it must stay MPI-free (no mpi.h, no raw
# MPI_*/exitMPI, only the allow-listed wrapperMPI helpers below) for the same
# reason phys_distributed_local.c does.

# Task 3 emptied this once ExpecLocal defensive guards were inserted and
# marker-wrapped around the partner_rank != myrank branches in both files
# (see EXPEC_LOCAL_GUARDED_BEGIN/END regions in src/nbody_correlation.c and
# src/anomalous_pair.c).
TEMP_UNGUARDED_FILES=""

# ExpecLocal permitted call matrix (must match the inventory document's
# table verbatim).
ALLOW="SumMPI_dc SumMPI_d SumMPI_li SumMPI_i fopenMPI childfopenMPI stdoutMPI"

# Token-recognition pattern for MPI-flavoured calls: same call vocabulary as
# the pattern used to generate the Step 1 inventory in
# docs/superpowers/specs/2026-07-11-expec-call-inventory.md, plus two
# additions made during the final whole-branch review:
#   - `fgetsMPI` added to the vocabulary: it does a real MPI_Bcast
#     (src/wrapperMPI.c) but its name doesn't match any of the other
#     alternatives, so an unguarded future call in a scanned file would
#     otherwise slip past this guard undetected.
#   - A `(^|[^A-Za-z0-9_])` left-boundary group prepended to the whole
#     alternation, so a call name can never match as a substring of a
#     longer identifier (e.g. `foo_MPI_Barrier(` no longer falsely matches
#     `MPI_Barrier(`). The name-extraction step below strips this boundary
#     character back off before comparing against ALLOW.
PATTERN='(^|[^A-Za-z0-9_])(SumMPI_[a-z]+|MaxMPI_[a-z]+|BcastMPI_[a-z]+|BarrierMPI|NormMPI_dc|VecProdMPI|MPI_[A-Za-z_]+|exitMPI|fopenMPI|childfopenMPI|fgetsMPI)\('

workdir=$(mktemp -d "${TMPDIR:-/tmp}/check_expec_local_calls.XXXXXX")
trap 'rm -rf "$workdir"' EXIT INT TERM

is_in_list() {
  needle="$1"
  shift
  for item in "$@"; do
    if [ "$item" = "$needle" ]; then
      return 0
    fi
  done
  return 1
}

violations=0
skipped_temp=0

for f in $FILES; do
  path="$srcdir/$f"
  if [ ! -f "$path" ]; then
    echo "ERROR: $f is listed in check_expec_local_calls.sh FILES but does not exist" >&2
    violations=$((violations + 1))
    continue
  fi

  # Always run the stripper so a broken/empty-output stripper is caught even
  # for temporarily-unguarded files.
  stripped="$workdir/stripped.c"
  python3 "$STRIPPER" "$path" > "$stripped"

  if is_in_list "$f" $TEMP_UNGUARDED_FILES; then
    skipped_temp=$((skipped_temp + 1))
    echo "SKIP (TEMP_UNGUARDED_FILES): $f"
    continue
  fi

  # Pass 1 (on the ORIGINAL, pre-strip source): collect
  # EXPEC_LOCAL_GUARDED_BEGIN/END line ranges. The marker convention is
  # defined on raw source line numbers, and strip_c_comments.py preserves
  # line numbers (it keeps every newline), so ranges collected here line up
  # exactly with the line numbers of matches found in the stripped pass.
  ranges="$workdir/ranges.txt"
  awk '
    /\/\*[ \t]*EXPEC_LOCAL_GUARDED_BEGIN[ \t]*\*\// { begin = NR; next }
    /\/\*[ \t]*EXPEC_LOCAL_GUARDED_END[ \t]*\*\// {
      if (begin != "") { print begin, NR; begin = "" }
    }
  ' "$path" > "$ranges"

  # Pass 2 (on the comment-stripped source): find MPI-flavoured call tokens.
  matches="$workdir/matches.txt"
  grep -noE "$PATTERN" "$stripped" > "$matches" || true

  if [ ! -s "$matches" ]; then
    continue
  fi

  while IFS=: read -r lineno token; do
    # Strip the trailing "(" first, then the single leading boundary
    # character captured by the PATTERN's "(^|[^A-Za-z0-9_])" group (there is
    # none to strip when the match starts at the true beginning of the line).
    name=$(printf '%s' "$token" | sed -E 's/\($//; s/^[^A-Za-z0-9_]//')

    # Excluded if inside a defensive-guard marker region.
    guarded=0
    if [ -s "$ranges" ]; then
      while read -r b e; do
        if [ "$lineno" -ge "$b" ] && [ "$lineno" -le "$e" ]; then
          guarded=1
          break
        fi
      done < "$ranges"
    fi
    if [ "$guarded" -eq 1 ]; then
      continue
    fi

    if is_in_list "$name" $ALLOW; then
      continue
    fi

    echo "VIOLATION: $f:$lineno: disallowed call '${name}(' not in ExpecLocal allow-list and not EXPEC_LOCAL_GUARDED" >&2
    violations=$((violations + 1))
  done < "$matches"
done

if [ "$violations" -gt 0 ]; then
  echo "FAILED: $violations disallowed MPI-flavoured call(s) found reachable from the expec_* layer." >&2
  echo "See docs/superpowers/specs/2026-07-11-expec-call-inventory.md for the frozen allow-list" >&2
  echo "and the EXPEC_LOCAL_GUARDED_BEGIN/END marker convention." >&2
  exit 1
fi

echo "PASSED: check_expec_local_calls (skipped via TEMP_UNGUARDED_FILES: $skipped_temp file(s))"
exit 0
