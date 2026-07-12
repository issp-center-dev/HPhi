/* HPhi  -  Quantum Lattice Model Simulator */
#ifndef HPHI_GREEN_OUTPUT_H
#define HPHI_GREEN_OUTPUT_H

#include <stdio.h>
#include "struct.h"

typedef enum {
  GreenOutputOneBody = 0,
  GreenOutputTwoBody,
  GreenOutputThreeBody,
  GreenOutputFourBody,
  GreenOutputSixBody,
  GreenOutputNBody,
  GreenOutputAnomalous
} GreenOutputKind;

typedef enum {
  GreenOutputTPQDataSS = 0,
  GreenOutputTPQDataNorm,
  GreenOutputTPQDataFlct
} GreenOutputTPQDataKind;

int GreenOutputUsesAggregate(const struct BindStruct *X);
int GreenOutputKindUsesAggregate(const struct BindStruct *X, GreenOutputKind kind);
const char *GreenOutputOpenMode(const struct BindStruct *X);
int GreenOutputFileName(const struct BindStruct *X, GreenOutputKind kind, char *sdt);
int GreenOutputWriteIndexPrefix(FILE *fp, const struct BindStruct *X);
int GreenOutputInitializeAggregateFiles(struct BindStruct *X);

/**
 * @brief Start a partial-output session for rank-local Green aggregate
 * writes (phase 3a, Mode 1 / state-task-parallel FullDiag).
 *
 * Zeroes every per-kind manifest record -- including sticky open_error and
 * first-open ("opened") state -- so a session started by this call never
 * inherits bookkeeping from a previous session in the same process (e.g. a
 * second phys() call). Must be paired with GreenOutputClearPartialSuffix().
 * Not meaningful outside MPI, but safe to call in a non-MPI build (there is
 * simply only one rank/session).
 * @param rank the calling rank's MPI rank (used to build unique part-file
 * names, e.g. "<final>.part<rank>").
 *
 * NOTE: partial-mode opens still go through childfopenMPI()/fopenMPI(), so
 * on ranks other than 0 they only reach the real filesystem while ExpecLocal
 * mode is active (see ExpecLocalEnter()/ExpecLocalLeave() in wrapperMPI.h).
 * The Mode 1 driver activates both together, in the actual order the driver
 * runs them (src/phys_distributed_local.c + src/phys_distributed.c):
 *   ExpecLocalEnter(); GreenOutputSetPartialSuffix(myrank); ... state loop ...
 *   GreenOutputClearPartialSuffix(); ExpecLocalLeave(); ... (rank-0 gather,
 *   collective) ... GreenOutputMergePartials(X);
 * GreenOutputClearPartialSuffix() runs before the collective
 * GreenOutputMergePartials() (not after, as an earlier draft of this comment
 * showed) -- this is manifest-neutral: Clear only stops future opens from
 * targeting the just-ended session's part-file suffix, it does not touch
 * the manifest or delete any part file, so ending the session before Merge
 * runs does not affect what Merge later reads.
 */
void GreenOutputSetPartialSuffix(int rank);

/**
 * @brief End the current partial-output session. After this call,
 * GreenOutputOpenAggregate()/GreenOutputCloseAggregate() behave exactly as
 * they do outside any session (the non-partial, childfopenMPI +
 * GreenOutputOpenMode() path). Does NOT clear the manifest -- the manifest
 * from the just-ended session remains readable (e.g. by
 * GreenOutputMergePartials()) until the next GreenOutputSetPartialSuffix()
 * call zeroes it.
 */
void GreenOutputClearPartialSuffix(void);

/**
 * @brief Open the aggregate Green-output file for the given kind.
 *
 * Outside a partial-output session: identical to the pre-existing
 * childfopenMPI(sdt, GreenOutputOpenMode(X), fp) pattern used at every
 * expec_* aggregate call site (sdt is derived internally via
 * GreenOutputFileName()).
 *
 * Inside a partial-output session (see GreenOutputSetPartialSuffix()): opens
 * this rank's private part file for `kind` instead of the shared aggregate
 * file. The FIRST open of a given kind in the session unlinks any stale part
 * file and creates it fresh ("w"); every subsequent open of the same kind in
 * the same session appends ("a") -- expec_* functions open/close a kind
 * repeatedly across the state loop, once per evaluated eigenstate, and all
 * of those writes must land in one contiguous per-rank file. Records
 * attempted=1 on first attempt for the kind; open_error is sticky (once set,
 * every later open of the same kind in the session fails immediately without
 * touching the filesystem again).
 *
 * @return 0 on success, -1 on failure (including: kind is not aggregate for
 * X, or the derived part/final path would exceed the manifest record's
 * 256-byte path fields).
 */
int GreenOutputOpenAggregate(struct BindStruct *X, GreenOutputKind kind, FILE **fp);

/**
 * @brief Close a FILE* previously returned by GreenOutputOpenAggregate().
 *
 * Outside a partial-output session: plain fclose().
 *
 * Inside a session: records bytes = ftell(fp) taken immediately before
 * fclose() (the running total for this rank/kind, since the file is opened
 * in append mode after the first write) and closed_ok (sticky: once a close
 * fails for a kind, closed_ok stays false for the rest of the session even
 * if a later close of the same kind succeeds), then fclose()s fp.
 *
 * @return 0 on success, -1 on failure.
 */
int GreenOutputCloseAggregate(GreenOutputKind kind, FILE *fp);

/**
 * @brief Rank 0 merges every rank's partial Green-output files into the
 * final aggregate files, gated on every rank's manifest reporting success.
 *
 * Collective over MPI_COMM_WORLD (every rank must call this). Rank 0
 * gathers all ranks' manifests (MPI_Gather of the fixed-shape in-memory
 * manifest record, see the .c file for why this is not a portable wire
 * format), determines a single pass/fail verdict, and MPI_Bcasts it to
 * every rank so the return value is identical everywhere.
 *
 * A run fails if ANY (rank, kind) record has attempted && open_error, OR if
 * rank 0 cannot re-open a part file that a manifest claims succeeded (this
 * catches a part file that vanished or was corrupted after a successful
 * close -- manifest contents alone are never trusted as proof of a
 * concatenatable file). A record with attempted && !open_error && bytes==0
 * (e.g. a zero-owner rank that legitimately never wrote a row) is a
 * perfectly valid empty contribution.
 *
 * On success: for every kind attempted by at least one rank, rank 0
 * concatenates every attempting rank's part file in rank order into a
 * private "<final>.tmp_merge" file (every fread/fwrite/fclose checked), and
 * only after ALL kinds concatenated cleanly rename()s each temp onto its
 * final name (this publish step is the Mode-1 replacement for
 * GreenOutputInitializeAggregateFiles(), which Mode 1 must NOT call
 * directly), then deletes all part files. On failure: temps are removed,
 * nothing is published under a final name, no part file is deleted
 * (surviving parts are left in place for diagnosis), and every rank
 * returns non-zero. A truncated file can therefore never appear under a
 * final aggregate name, even on a mid-merge disk-full/write error.
 *
 * In a build without MPI this is a no-op that returns 0 (no cross-rank
 * partial files can exist to merge).
 */
int GreenOutputMergePartials(struct BindStruct *X);

int GreenOutputUsesTPQDataAggregate(const struct BindStruct *X);
int GreenOutputTPQDataFileName(const struct BindStruct *X, GreenOutputTPQDataKind kind, char *sdt);
int GreenOutputInitializeTPQDataAggregateFiles(struct BindStruct *X);
void GreenOutputWriteTPQSSRow(FILE *fp, const struct BindStruct *X, int step, double inv_temp);
void GreenOutputWriteTPQNormRow(FILE *fp, const struct BindStruct *X, int step, double inv_temp,
                                double norm, double first_norm);
void GreenOutputWriteTPQFlctRow(FILE *fp, const struct BindStruct *X, int step, double inv_temp);

#endif
