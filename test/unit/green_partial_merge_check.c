/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 The University of Tokyo */

/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>. */

/*
 * Direct, in-process unit test of GreenOutputMergePartials() and the
 * partial-output session API (phase 3a, Task 4/green_output.c), driven
 * against the REAL green_output objects (no fakes) with real MPI ranks
 * writing real files under a dedicated scratch output directory.
 *
 * Registered min:2 (see test/CMakeLists.txt): every scenario below is
 * written to generalize to any nprocs >= 2 -- rank 0 plays the
 * "always contributes" role and the highest rank plays the
 * "special/failing" role, so the same binary is meaningful whether run
 * with exactly 2 ranks or more.
 *
 * Six scenarios. The first five use separate GreenOutputKind values; the
 * sixth deliberately writes two kinds in one session to exercise
 * cross-kind publication rollback:
 *   (a) GreenOutputOneBody  -- normal success with an old final and a stale
 *       recovery-name collision: every rank writes one
 *       identifying row; Merge succeeds on every rank, the final file
 *       contains every rank's row IN RANK ORDER, and every part file is
 *       gone afterwards.
 *   (b) GreenOutputTwoBody  -- legitimately empty: every rank opens and
 *       closes without necessarily writing (only rank 0 writes a row);
 *       the non-writing ranks still get attempted=1/bytes=0, which the
 *       merge must treat as valid, not as a failure. Merge succeeds and
 *       the final file contains exactly rank 0's row.
 *   (c) GreenOutputThreeBody -- failure injection: every rank opens,
 *       writes, and closes successfully, THEN (after Close, exactly as
 *       required by the spec -- not by making the API itself fail) the
 *       test unlinks the highest rank's part file out from under the
 *       manifest. Merge must return nonzero on EVERY rank, must not
 *       publish any final file (verified via stat, after removing any
 *       pre-existing final so this really tests "no NEW final"), and
 *       every OTHER rank's part file must still be on disk afterwards
 *       (nothing is deleted on failure).
 *   (d) GreenOutputFourBody -- lifecycle failure: the highest rank opens
 *       and writes but never closes (manifest: opened=1, closed_ok=0);
 *       merge validation must reject the manifest on every rank and
 *       publish nothing.
 *   (e) GreenOutputSixBody  -- size mismatch: all ranks open/write/close
 *       cleanly, then the highest rank's part is truncated from outside
 *       the API; the merge's pass-2 probe must detect that the on-disk
 *       length differs from the manifest's recorded byte count and fail
 *       on every rank, publishing nothing and deleting nothing.
 *   (f) GreenOutputOneBody + GreenOutputTwoBody -- both have pre-existing
 *       finals and valid new parts. The test injects a failure while the
 *       second new final is being published, after the first was already
 *       published. Merge must restore BOTH old finals, leave all parts,
 *       remove temp/backup artifacts, and fail on every rank.
 *
 * Non-vacuousness check (documented per the phase-3a plan's TDD
 * requirement -- this test is not driven by a prior failing test, so its
 * own discriminating power must be checked by hand): during development,
 * scenario (c)'s two key assertions were each temporarily inverted --
 * `rc == 0` in place of `rc != 0`, and the final-file stat check flipped
 * to require existence -- and the test was rebuilt and rerun; both
 * mutations made it fail (see task-7-report.md for the exact commands and
 * output), confirming the assertions are load-bearing rather than
 * vacuously true. The mutations were reverted before committing.
 */
#include <mpi.h>
#include <errno.h>
#include <glob.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "green_output.h"
#include "wrapperMPI.h"
#include "global.h"
#include "DefCommon.h"
#include "struct.h"

/* wrapperMPI.c's InitializeMPI() calls splash() (src/splash.c); this test
   never calls InitializeMPI(), but the wrapperMPI.o translation unit's
   undefined reference to splash() must still be resolved at link time
   (same reasoning as unit/elpa_eigen_check.c's `int myrank = 0;` stub for
   src/matrixscalapack.o). We deliberately do not link src/splash.c here
   (it pulls in version_major.h/version_minor.h/version_patch.h for no
   benefit to this test), so provide a trivial stub instead. */
void splash(void) { }

static int g_ok = 1;
static int g_rank = 0;
static int g_rename_call = 0;
static int g_rename_fail_on_call = 0;

/* src/green_output.c calls this instead of rename() only in this unit-test
   target (GREEN_OUTPUT_TESTING). The injected error is one-shot so rollback
   renames can proceed normally after the targeted publish fails. */
int GreenOutputTestRename(const char *old_path, const char *new_path) {
  g_rename_call++;
  if (g_rename_fail_on_call == g_rename_call) {
    g_rename_fail_on_call = 0;
    errno = EIO;
    return -1;
  }
  return rename(old_path, new_path);
}

#define CHECK(cond, ...) \
  do { \
    if (!(cond)) { \
      fprintf(stderr, "[rank %d] CHECK FAILED at %s:%d: ", g_rank, __FILE__, __LINE__); \
      fprintf(stderr, __VA_ARGS__); \
      fprintf(stderr, "\n"); \
      g_ok = 0; \
    } \
  } while (0)

static int FileExists(const char *path) {
  struct stat st;
  return stat(path, &st) == 0;
}

static int BackupArtifactExists(const char *joined) {
  char pattern[D_FileNameMax + 128];
  glob_t matches;
  int found;
  memset(&matches, 0, sizeof(matches));
  snprintf(pattern, sizeof(pattern), "%s.bak_merge*", joined);
  found = (glob(pattern, 0, NULL, &matches) == 0 && matches.gl_pathc > 0);
  globfree(&matches);
  return found;
}

static void RemoveBackupArtifacts(const char *joined) {
  char pattern[D_FileNameMax + 128];
  glob_t matches;
  size_t i;
  memset(&matches, 0, sizeof(matches));
  snprintf(pattern, sizeof(pattern), "%s.bak_merge*", joined);
  if (glob(pattern, 0, NULL, &matches) == 0) {
    for (i = 0; i < matches.gl_pathc; i++) remove(matches.gl_pathv[i]);
  }
  globfree(&matches);
}

static int FileHasExactContents(const char *path, const char *expected) {
  FILE *fp = fopen(path, "rb");
  char buf[256];
  size_t got;
  size_t expected_len = strlen(expected);
  int same;
  if (fp == NULL) return 0;
  got = fread(buf, 1, sizeof(buf), fp);
  same = (got == expected_len && memcmp(buf, expected, expected_len) == 0 &&
          !ferror(fp) && feof(fp));
  fclose(fp);
  return same;
}

/* Join cParentOutputFolder + relative path, exactly as GreenOutputJoinOutputPath()
   does internally (that helper is static to green_output.c, so this test
   reimplements the same join rule against the public cParentOutputFolder). */
static void JoinOutputPath(const char *rel, char *out, size_t outsz) {
  out[0] = '\0';
  strncat(out, cParentOutputFolder, outsz - 1);
  strncat(out, rel, outsz - 1 - strlen(out));
}

static void PartPath(const char *final_rel, int rank, char *out, size_t outsz) {
  char rel[512];
  snprintf(rel, sizeof(rel), "%s.part%d", final_rel, rank);
  JoinOutputPath(rel, out, outsz);
}

/* Remove a kind's final file and every rank-0..nprocs_hint-1 part file, so
   each scenario starts from a clean slate regardless of leftovers from a
   previous run of this binary. Only rank 0 needs to do this (shared
   filesystem for all ranks of a single-node MPI test), followed by a
   barrier before any rank proceeds. */
static void PreClean(const struct BindStruct *X, GreenOutputKind kind, int nprocs_hint) {
  char final_rel[D_FileNameMax];
  char joined[D_FileNameMax + 64];
  int r;
  if (g_rank != 0) return;
  if (GreenOutputFileName(X, kind, final_rel) != 0) return;
  JoinOutputPath(final_rel, joined, sizeof(joined));
  remove(joined);
  {
    char auxiliary[D_FileNameMax + 96];
    snprintf(auxiliary, sizeof(auxiliary), "%s.tmp_merge", joined);
    remove(auxiliary);
    RemoveBackupArtifacts(joined);
  }
  for (r = 0; r < nprocs_hint; r++) {
    char part_joined[D_FileNameMax + 64];
    PartPath(final_rel, r, part_joined, sizeof(part_joined));
    remove(part_joined);
  }
}

int main(int argc, char **argv) {
  int nprocs = 1;
  struct BindStruct X;
  char final_rel[D_FileNameMax];
  char joined[D_FileNameMax + 64];

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  /* We deliberately do NOT call InitializeMPI() (it would MPI_Init() a
     second time, print the splash banner/parallelization info, and pull in
     HPHI_MPI_NOBATCH parsing this test has no use for) -- but two of the
     globals it sets are load-bearing for the code under test and must be
     replicated here exactly as InitializeMPI() does (src/wrapperMPI.c):
     `myrank` (childfopenMPI()/fopenMPI()'s non-ExpecLocal rank-0 gate reads
     this global, not the local MPI rank) and `stdoutMPI` (every error path
     in childfopenMPI()/GreenOutputMergePartials() fprintf()s to it; left at
     its global.c default of NULL, any error path -- including the one this
     test's own scenario (c) deliberately triggers -- segfaults instead of
     printing a diagnostic). */
  myrank = g_rank;
  stdoutMPI = (g_rank == 0) ? stdout : fopen("/dev/null", "w");

  if (nprocs < 2) {
    fprintf(stderr, "green_partial_merge_check requires at least 2 MPI ranks (got %d).\n", nprocs);
    MPI_Finalize();
    return 1;
  }

  /* Dedicated scratch output directory (spec requirement: do not touch
     any other test's output/ directory). cParentOutputFolder is a plain
     `const char*` global (not `const char* const`), so reassigning it
     here is exactly how the real driver would point childfopenMPI() et
     al. at a different tree -- no green_output.c change needed. */
  cParentOutputFolder = "./green_partial_merge_check_out/";
  if (g_rank == 0) {
    mkdir("green_partial_merge_check_out", 0777); /* ok if it already exists */
  }
  MPI_Barrier(MPI_COMM_WORLD);

  memset(&X, 0, sizeof(X));
  X.Def.iOutputGreenFormat = OUTPUTGREENFORMAT_AGGREGATE;
  X.Def.iCalcType = FullDiag; /* family 3: "_eigen" aggregate file names */
  X.Def.CDataFileHead = "zvo";

  /* ---- Pre-clean every kind used below (spec: clear final/part before
     each scenario so "not published"/"still present" stats are meaningful,
     not accidental leftovers from an earlier run of this binary). ---- */
  PreClean(&X, GreenOutputOneBody, nprocs);
  PreClean(&X, GreenOutputTwoBody, nprocs);
  PreClean(&X, GreenOutputThreeBody, nprocs);
  PreClean(&X, GreenOutputFourBody, nprocs);
  PreClean(&X, GreenOutputSixBody, nprocs);
  MPI_Barrier(MPI_COMM_WORLD);

  /* =====================================================================
   * (a) Normal success with a stale backup-name collision: every rank writes
   *     one identifying row into the OneBody kind; Merge must succeed
   *     everywhere, the final file must contain every rank's row IN RANK
   *     ORDER, and every part file must be gone afterwards.
   * ===================================================================*/
  {
    FILE *fp = NULL;
    int rc_open, rc_close, rc_merge;
    char stale_backup[D_FileNameMax + 128] = {0};

    if (g_rank == 0) {
      FILE *seed;
      CHECK(GreenOutputFileName(&X, GreenOutputOneBody, final_rel) == 0,
            "(a) GreenOutputFileName failed while seeding old output");
      JoinOutputPath(final_rel, joined, sizeof(joined));
      seed = fopen(joined, "wb");
      CHECK(seed != NULL, "(a) could not seed old OneBody final");
      if (seed != NULL) { fputs("old-one\n", seed); fclose(seed); }
      snprintf(stale_backup, sizeof(stale_backup), "%s.bak_merge.%ld.0",
               joined, (long)getpid());
      seed = fopen(stale_backup, "wb");
      CHECK(seed != NULL, "(a) could not seed stale recovery backup");
      if (seed != NULL) { fputs("stale-backup\n", seed); fclose(seed); }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    GreenOutputSetPartialSuffix(g_rank);
    ExpecLocalEnter();
    rc_open = GreenOutputOpenAggregate(&X, GreenOutputOneBody, &fp);
    CHECK(rc_open == 0, "(a) GreenOutputOpenAggregate failed for OneBody");
    if (rc_open == 0) {
      fprintf(fp, "rank%d\n", g_rank);
    }
    rc_close = GreenOutputCloseAggregate(GreenOutputOneBody, fp);
    CHECK(rc_close == 0, "(a) GreenOutputCloseAggregate failed for OneBody");
    ExpecLocalLeave();

    rc_merge = GreenOutputMergePartials(&X);
    CHECK(rc_merge == 0, "(a) GreenOutputMergePartials returned %d, expected 0 (normal success)", rc_merge);
    GreenOutputClearPartialSuffix();

    if (g_rank == 0) {
      CHECK(GreenOutputFileName(&X, GreenOutputOneBody, final_rel) == 0,
            "(a) GreenOutputFileName failed for OneBody");
      JoinOutputPath(final_rel, joined, sizeof(joined));
      {
        FILE *fin = fopen(joined, "r");
        char line[256];
        int r = 0;
        CHECK(fin != NULL, "(a) final OneBody file '%s' was not published", joined);
        if (fin != NULL) {
          while (fgets(line, sizeof(line), fin) != NULL) {
            char expect[64];
            size_t len = strlen(line);
            if (len > 0 && line[len - 1] == '\n') line[len - 1] = '\0';
            snprintf(expect, sizeof(expect), "rank%d", r);
            CHECK(strcmp(line, expect) == 0,
                  "(a) row %d was '%s', expected '%s' (rows must be in rank order)", r, line, expect);
            r++;
          }
          fclose(fin);
          CHECK(r == nprocs, "(a) final OneBody file had %d rows, expected %d (one per rank)", r, nprocs);
        }
      }
      /* Every part file must be gone: success deletes all parts. */
      {
        int r;
        for (r = 0; r < nprocs; r++) {
          char part_joined[D_FileNameMax + 64];
          PartPath(final_rel, r, part_joined, sizeof(part_joined));
          CHECK(!FileExists(part_joined), "(a) part file '%s' should have been deleted after success", part_joined);
        }
      }
      CHECK(FileHasExactContents(stale_backup, "stale-backup\n"),
            "(a) stale recovery backup was overwritten or removed");
      remove(stale_backup);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  /* =====================================================================
   * (b) Legitimately empty: every rank opens and closes the TwoBody kind
   *     (attempted=1 for all), but only rank 0 writes a row -- every
   *     other rank's part file is a real, successfully-closed, 0-byte
   *     file. Merge must treat that as valid (not a failure) and the
   *     final file must contain exactly rank 0's row.
   * ===================================================================*/
  {
    FILE *fp = NULL;
    int rc_open, rc_close, rc_merge;

    GreenOutputSetPartialSuffix(g_rank);
    ExpecLocalEnter();
    rc_open = GreenOutputOpenAggregate(&X, GreenOutputTwoBody, &fp);
    CHECK(rc_open == 0, "(b) GreenOutputOpenAggregate failed for TwoBody");
    if (rc_open == 0 && g_rank == 0) {
      fprintf(fp, "only-rank0-row\n");
    }
    /* Other ranks: open then close with zero bytes written -- the
       "legitimately empty" contribution the spec calls out. */
    rc_close = GreenOutputCloseAggregate(GreenOutputTwoBody, fp);
    CHECK(rc_close == 0, "(b) GreenOutputCloseAggregate failed for TwoBody");
    ExpecLocalLeave();

    rc_merge = GreenOutputMergePartials(&X);
    CHECK(rc_merge == 0, "(b) GreenOutputMergePartials returned %d, expected 0 (legitimately-empty ranks are not a failure)", rc_merge);
    GreenOutputClearPartialSuffix();

    if (g_rank == 0) {
      CHECK(GreenOutputFileName(&X, GreenOutputTwoBody, final_rel) == 0,
            "(b) GreenOutputFileName failed for TwoBody");
      JoinOutputPath(final_rel, joined, sizeof(joined));
      {
        FILE *fin = fopen(joined, "r");
        char line[256];
        int r = 0;
        CHECK(fin != NULL, "(b) final TwoBody file '%s' was not published", joined);
        if (fin != NULL) {
          while (fgets(line, sizeof(line), fin) != NULL) r++;
          fclose(fin);
          CHECK(r == 1, "(b) final TwoBody file had %d rows, expected exactly 1 (only rank 0 wrote)", r);
        }
      }
      {
        int r;
        for (r = 0; r < nprocs; r++) {
          char part_joined[D_FileNameMax + 64];
          PartPath(final_rel, r, part_joined, sizeof(part_joined));
          CHECK(!FileExists(part_joined), "(b) part file '%s' should have been deleted after success", part_joined);
        }
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  /* =====================================================================
   * (c) Failure injection: every rank opens, writes, and closes the
   *     ThreeBody kind successfully (a clean manifest on every rank).
   *     THEN, after Close, the test itself unlinks the HIGHEST rank's
   *     part file (simulating it vanishing from disk after a successful
   *     close -- the exact case GreenOutputMergePartials()'s "pass 2"
   *     re-open verification exists to catch). Merge must return nonzero
   *     on EVERY rank, must not publish a final file, and must leave
   *     every surviving part file in place.
   * ===================================================================*/
  {
    FILE *fp = NULL;
    int rc_open, rc_close, rc_merge;
    int victim = nprocs - 1;

    GreenOutputSetPartialSuffix(g_rank);
    ExpecLocalEnter();
    rc_open = GreenOutputOpenAggregate(&X, GreenOutputThreeBody, &fp);
    CHECK(rc_open == 0, "(c) GreenOutputOpenAggregate failed for ThreeBody");
    if (rc_open == 0) {
      fprintf(fp, "rank%d\n", g_rank);
    }
    rc_close = GreenOutputCloseAggregate(GreenOutputThreeBody, fp);
    CHECK(rc_close == 0, "(c) GreenOutputCloseAggregate failed for ThreeBody");
    ExpecLocalLeave();

    /* The failure is injected by the TEST, from outside the API, strictly
       after Close -- never by making Open/Close themselves fail. */
    if (g_rank == victim) {
      CHECK(GreenOutputFileName(&X, GreenOutputThreeBody, final_rel) == 0,
            "(c) GreenOutputFileName failed for ThreeBody");
      {
        char part_joined[D_FileNameMax + 64];
        PartPath(final_rel, victim, part_joined, sizeof(part_joined));
        CHECK(unlink(part_joined) == 0, "(c) failed to unlink victim part file '%s'", part_joined);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    rc_merge = GreenOutputMergePartials(&X);
    CHECK(rc_merge != 0, "(c) GreenOutputMergePartials returned 0, expected nonzero on every rank after a part file vanished");
    GreenOutputClearPartialSuffix();

    if (g_rank == 0) {
      CHECK(GreenOutputFileName(&X, GreenOutputThreeBody, final_rel) == 0,
            "(c) GreenOutputFileName failed for ThreeBody (post-merge check)");
      JoinOutputPath(final_rel, joined, sizeof(joined));
      CHECK(!FileExists(joined), "(c) final ThreeBody file '%s' must NOT be published after a failed merge", joined);
      {
        int r;
        for (r = 0; r < nprocs; r++) {
          char part_joined[D_FileNameMax + 64];
          PartPath(final_rel, r, part_joined, sizeof(part_joined));
          if (r == victim) {
            CHECK(!FileExists(part_joined), "(c) victim part file '%s' should not reappear", part_joined);
          } else {
            CHECK(FileExists(part_joined), "(c) surviving part file '%s' must be retained after a failed merge", part_joined);
          }
        }
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  /* =====================================================================
   * (d) Lifecycle-failure injection: the HIGHEST rank opens and writes
   *     the FourBody kind but NEVER closes it, so its manifest record has
   *     opened=1 / closed_ok=0 (buffered data may not be on disk even
   *     though the part file exists and is readable). Merge must reject
   *     the manifest (pass 1 lifecycle check) with nonzero rc on EVERY
   *     rank and must not publish a final file.
   * ===================================================================*/
  {
    FILE *fp = NULL;
    int rc_open, rc_merge;
    int victim = nprocs - 1;

    GreenOutputSetPartialSuffix(g_rank);
    ExpecLocalEnter();
    rc_open = GreenOutputOpenAggregate(&X, GreenOutputFourBody, &fp);
    CHECK(rc_open == 0, "(d) GreenOutputOpenAggregate failed for FourBody");
    if (rc_open == 0) {
      fprintf(fp, "rank%d\n", g_rank);
    }
    if (g_rank == victim) {
      /* deliberately do NOT close: closed_ok stays 0 in the manifest */
    } else {
      int rc_close = GreenOutputCloseAggregate(GreenOutputFourBody, fp);
      CHECK(rc_close == 0, "(d) GreenOutputCloseAggregate failed for FourBody");
      fp = NULL;
    }
    ExpecLocalLeave();
    MPI_Barrier(MPI_COMM_WORLD);

    rc_merge = GreenOutputMergePartials(&X);
    CHECK(rc_merge != 0, "(d) GreenOutputMergePartials returned 0, expected nonzero on every rank for an opened-but-never-closed part");
    GreenOutputClearPartialSuffix();
    if (g_rank == victim && fp != NULL) fclose(fp); /* test hygiene only */

    if (g_rank == 0) {
      CHECK(GreenOutputFileName(&X, GreenOutputFourBody, final_rel) == 0,
            "(d) GreenOutputFileName failed for FourBody (post-merge check)");
      JoinOutputPath(final_rel, joined, sizeof(joined));
      CHECK(!FileExists(joined), "(d) final FourBody file '%s' must NOT be published after a lifecycle-failed merge", joined);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  /* =====================================================================
   * (e) Size-mismatch injection: every rank opens, writes, and closes the
   *     SixBody kind successfully, THEN the HIGHEST rank truncates its
   *     own part file from outside the API (rewrites it shorter than the
   *     manifest's recorded byte count). Merge's pass-2 probe must detect
   *     that the on-disk length no longer matches the manifest and fail
   *     on EVERY rank without publishing; all parts stay in place.
   * ===================================================================*/
  {
    FILE *fp = NULL;
    int rc_open, rc_close, rc_merge;
    int victim = nprocs - 1;

    GreenOutputSetPartialSuffix(g_rank);
    ExpecLocalEnter();
    rc_open = GreenOutputOpenAggregate(&X, GreenOutputSixBody, &fp);
    CHECK(rc_open == 0, "(e) GreenOutputOpenAggregate failed for SixBody");
    if (rc_open == 0) {
      fprintf(fp, "rank%d payload payload payload\n", g_rank);
    }
    rc_close = GreenOutputCloseAggregate(GreenOutputSixBody, fp);
    CHECK(rc_close == 0, "(e) GreenOutputCloseAggregate failed for SixBody");
    ExpecLocalLeave();

    if (g_rank == victim) {
      CHECK(GreenOutputFileName(&X, GreenOutputSixBody, final_rel) == 0,
            "(e) GreenOutputFileName failed for SixBody");
      {
        char part_joined[D_FileNameMax + 64];
        FILE *trunc = NULL;
        PartPath(final_rel, victim, part_joined, sizeof(part_joined));
        trunc = fopen(part_joined, "wb");
        CHECK(trunc != NULL, "(e) failed to reopen victim part file '%s' for truncation", part_joined);
        if (trunc != NULL) {
          fputs("x\n", trunc); /* shorter than the manifest's byte count */
          fclose(trunc);
        }
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    rc_merge = GreenOutputMergePartials(&X);
    CHECK(rc_merge != 0, "(e) GreenOutputMergePartials returned 0, expected nonzero on every rank for a size-mismatched part");
    GreenOutputClearPartialSuffix();

    if (g_rank == 0) {
      CHECK(GreenOutputFileName(&X, GreenOutputSixBody, final_rel) == 0,
            "(e) GreenOutputFileName failed for SixBody (post-merge check)");
      JoinOutputPath(final_rel, joined, sizeof(joined));
      CHECK(!FileExists(joined), "(e) final SixBody file '%s' must NOT be published after a size-mismatch merge failure", joined);
      {
        int r;
        for (r = 0; r < nprocs; r++) {
          char part_joined[D_FileNameMax + 64];
          PartPath(final_rel, r, part_joined, sizeof(part_joined));
          CHECK(FileExists(part_joined), "(e) part file '%s' must be retained after a failed merge", part_joined);
        }
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  /* =====================================================================
   * (f) Cross-kind publish failure: seed old OneBody and TwoBody finals,
   *     build valid new parts for both, then fail the fourth rename on rank
   *     0. The first two renames create backups, the third publishes the
   *     new OneBody file, and the fourth attempts to publish TwoBody. The
   *     merge must roll the whole output set back to the old generation.
   * ===================================================================*/
  {
    FILE *fp_one = NULL;
    FILE *fp_two = NULL;
    int rc_merge;
    char one_rel[D_FileNameMax] = {0}, two_rel[D_FileNameMax] = {0};
    char one_joined[D_FileNameMax + 64] = {0};
    char two_joined[D_FileNameMax + 64] = {0};

    PreClean(&X, GreenOutputOneBody, nprocs);
    PreClean(&X, GreenOutputTwoBody, nprocs);
    if (g_rank == 0) {
      FILE *seed;
      CHECK(GreenOutputFileName(&X, GreenOutputOneBody, one_rel) == 0,
            "(f) GreenOutputFileName failed for OneBody");
      CHECK(GreenOutputFileName(&X, GreenOutputTwoBody, two_rel) == 0,
            "(f) GreenOutputFileName failed for TwoBody");
      JoinOutputPath(one_rel, one_joined, sizeof(one_joined));
      JoinOutputPath(two_rel, two_joined, sizeof(two_joined));
      seed = fopen(one_joined, "wb");
      CHECK(seed != NULL, "(f) could not seed old OneBody final");
      if (seed != NULL) { fputs("old-one\n", seed); fclose(seed); }
      seed = fopen(two_joined, "wb");
      CHECK(seed != NULL, "(f) could not seed old TwoBody final");
      if (seed != NULL) { fputs("old-two\n", seed); fclose(seed); }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    GreenOutputSetPartialSuffix(g_rank);
    ExpecLocalEnter();
    CHECK(GreenOutputOpenAggregate(&X, GreenOutputOneBody, &fp_one) == 0,
          "(f) GreenOutputOpenAggregate failed for OneBody");
    if (fp_one != NULL) fprintf(fp_one, "new-one-rank%d\n", g_rank);
    CHECK(GreenOutputCloseAggregate(GreenOutputOneBody, fp_one) == 0,
          "(f) GreenOutputCloseAggregate failed for OneBody");
    CHECK(GreenOutputOpenAggregate(&X, GreenOutputTwoBody, &fp_two) == 0,
          "(f) GreenOutputOpenAggregate failed for TwoBody");
    if (fp_two != NULL) fprintf(fp_two, "new-two-rank%d\n", g_rank);
    CHECK(GreenOutputCloseAggregate(GreenOutputTwoBody, fp_two) == 0,
          "(f) GreenOutputCloseAggregate failed for TwoBody");
    ExpecLocalLeave();

    if (g_rank == 0) {
      g_rename_call = 0;
      g_rename_fail_on_call = 4;
    }
    rc_merge = GreenOutputMergePartials(&X);
    CHECK(rc_merge != 0,
          "(f) GreenOutputMergePartials returned 0 after injected publish failure");
    GreenOutputClearPartialSuffix();

    if (g_rank == 0) {
      int r;
      char auxiliary[D_FileNameMax + 96];
      CHECK(g_rename_call == 6,
            "(f) rename wrapper saw %d calls, expected 6 (2 backup, 2 publish attempts, 2 restores)",
            g_rename_call);
      CHECK(FileHasExactContents(one_joined, "old-one\n"),
            "(f) old OneBody final was not restored exactly");
      CHECK(FileHasExactContents(two_joined, "old-two\n"),
            "(f) old TwoBody final was not restored exactly");
      snprintf(auxiliary, sizeof(auxiliary), "%s.tmp_merge", one_joined);
      CHECK(!FileExists(auxiliary), "(f) OneBody temp remained after rollback");
      CHECK(!BackupArtifactExists(one_joined),
            "(f) OneBody backup remained after rollback");
      snprintf(auxiliary, sizeof(auxiliary), "%s.tmp_merge", two_joined);
      CHECK(!FileExists(auxiliary), "(f) TwoBody temp remained after rollback");
      CHECK(!BackupArtifactExists(two_joined),
            "(f) TwoBody backup remained after rollback");
      for (r = 0; r < nprocs; r++) {
        char part_joined[D_FileNameMax + 64];
        PartPath(one_rel, r, part_joined, sizeof(part_joined));
        CHECK(FileExists(part_joined),
              "(f) OneBody part '%s' must remain after rollback", part_joined);
        PartPath(two_rel, r, part_joined, sizeof(part_joined));
        CHECK(FileExists(part_joined),
              "(f) TwoBody part '%s' must remain after rollback", part_joined);
      }
    }
  }

  {
    int global_ok = 0;
    MPI_Allreduce(&g_ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (g_rank == 0) {
      printf("green_partial_merge_check: %s\n", global_ok ? "OK" : "FAILED");
    }
    MPI_Finalize();
    return global_ok ? 0 : 1;
  }
}
