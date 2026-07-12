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
 * Three scenarios, one per GreenOutputKind so they cannot interfere with
 * each other's manifest bookkeeping within a single process:
 *   (a) GreenOutputOneBody  -- normal success: every rank writes one
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
  MPI_Barrier(MPI_COMM_WORLD);

  /* =====================================================================
   * (a) Normal success: every rank writes one identifying row into the
   *     OneBody kind; Merge must succeed everywhere, the final file must
   *     contain every rank's row IN RANK ORDER, and every part file must
   *     be gone afterwards.
   * ===================================================================*/
  {
    FILE *fp = NULL;
    int rc_open, rc_close, rc_merge;

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
