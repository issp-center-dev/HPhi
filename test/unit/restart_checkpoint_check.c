/* Raw common-module test: real file I/O and collectives, no solver lifecycle. */
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "Common.h"
#include "restart_checkpoint.h"
#include "restart_manifest.h"
#include "symmetry_basis.h"
#include "symmetry_checkpoint.h"
#include "wrapperMPI.h"
#include "restart_check.h"
#ifdef MPI
#include <mpi.h>
#endif

/* Sector identities are exercised by integration tests; this target only
 * constructs a raw runtime, and must fail if the sector identity helper is called. */
int SymmetryCheckpointIdentityWords(const struct BindStruct *X, uint64_t *words)
{ (void)X; (void)words; abort(); }

/* Old sector imports belong to the integration suite, not this raw harness. */
int ReadSymmetryCheckpoint(const struct BindStruct *X, const char *name,
                           double complex *vector, struct SymmetryCheckpointInfo *info)
{ (void)X; (void)name; (void)vector; (void)info; abort(); }

int SymmetryCheckpointLegacyFileHash(FILE *fp, uint64_t rank, uint64_t offset,
                                     uint64_t local_dim, uint64_t *hash,
                                     double *norm2)
{
  (void)fp; (void)rank; (void)offset; (void)local_dim; (void)hash; (void)norm2;
  abort();
}

/* The standalone common-module target does not link the application timer. */
void StartTimer(int id) { (void)id; }
void StopTimer(int id) { (void)id; }

static int all_ok(int ok) { return SumMPI_i(!ok) == 0; }
static char root_state[D_FileNameMax];

static int replace_word(int index, uint64_t value, uint64_t *old)
{
  int bad = 0;
  if (myrank == 0) {
    FILE *fp = fopen(root_state, "r+b");
    bad = fp == NULL;
    if (!bad) bad = fseek(fp, 8L * index, SEEK_SET) != 0 || RestartGetWord(fp, old) != 0 ||
                    fseek(fp, 8L * index, SEEK_SET) != 0 || RestartPutWord(fp, value) != 0;
    if (fp != NULL) bad |= fclose(fp) != 0;
  }
  return all_ok(!bad);
}

static int rejected_read(struct BindStruct *X, const struct RestartManifestRow *row)
{
  double complex value[3] = {7 + 7*I, 7 + 7*I, 7 + 7*I};
  struct RestartVectorInfo info;
  int rc = RestartReadVector(X, row, 0, RESTART_READ_CONTINUATION, value, &info);
  return all_ok(rc != 0 && value[0] == 7 + 7*I && value[1] == 7 + 7*I && value[2] == 7 + 7*I);
}

int main(int argc, char **argv)
{
#ifdef MPI
  CHECK(MPI_Init(&argc, &argv) == MPI_SUCCESS);
  CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &myrank) == MPI_SUCCESS);
  CHECK(MPI_Comm_size(MPI_COMM_WORLD, &nproc) == MPI_SUCCESS);
#else
  (void)argc; (void)argv; myrank = 0; nproc = 1;
#endif
  stdoutMPI = stdout;
  CHECK(ExpecLocalActive() == 0);
  struct BindStruct X = {0};
  X.Def.CDataFileHead = "module";
  X.Def.iCalcType = CG; X.Def.k_exct = 1; X.Def.iReStart = RESTART_OUT;
  X.Def.iCalcModelInput = SpinGC;
  X.Check.idim_max = nproc > 1 && myrank == nproc - 1 ? 0 : 2;
  X.Check.idim_maxMPI = SumMPI_li(X.Check.idim_max);
  CHECK(all_ok(RestartLocalOffset(&X) == (uint64_t)myrank * 2));
  char folder[D_FileNameMax] = {0};
  if (myrank == 0) snprintf(folder, sizeof(folder), "restart_checkpoint_check.%ld/", (long)getpid());
#ifdef MPI
  CHECK(MPI_Bcast(folder, sizeof(folder), MPI_CHAR, 0, MPI_COMM_WORLD) == MPI_SUCCESS);
#endif
  CHECK(all_ok(myrank != 0 || mkdir(folder, 0700) == 0));
  cParentOutputFolder = folder;
  snprintf(root_state, sizeof(root_state), "%smodule_restart_set0_comp0_gen0_rank_0.dat", folder);
  double complex state[3] = {0}, direction[3] = {0}, loaded[3] = {0};
  if (myrank == 0) state[1] = 1;
  double complex *const components[] = {state, direction};
  double complex *const *const vectors[] = {components};
  uint64_t sets[] = {0};
  struct RestartUnit unit = {1, sets, 2};
  struct RestartVectorInfo info = {.method = CG, .set = 0, .step = 4, .nset = 1, .scalar_a = -2};
  CHECK(RestartWriteUnit(&X, &unit, vectors, &info, 0) == 0);
  struct RestartManifest m = {0}; int exists = 0;
  CHECK(RestartManifestLoad(&X, &m, &exists) == 0 && exists);
  const struct RestartManifestRow *row = RestartManifestFind(&m, 0);
  CHECK(row != NULL && row->gen == 0);
  /* A rank-local missing row must agree before the reader dereferences it. */
  CHECK(rejected_read(&X, myrank == 0 ? NULL : row));
  enum RestartSource source; int legacy_step;
  CHECK(RestartClassifySet(&X, &m, 0, NULL, &source, &legacy_step) == 0 && source == RESTART_SOURCE_V3);
  struct RestartVectorInfo result;
  CHECK(RestartReadVector(&X, row, 0, RESTART_READ_CONTINUATION, loaded, &result) == 0);
  CHECK(all_ok(memcmp(state, loaded, (X.Check.idim_max + 1) * sizeof(*state)) == 0));
  CHECK(result.step == 4 && result.scalar_a == -2 && result.nset == 1);
  CHECK(RestartReadVector(&X, row, 1, RESTART_READ_CONTINUATION, loaded, &result) == 0);
  /* Import is a state-only contract even through the shared public reader. */
  CHECK(RestartReadVector(&X, row, 0, RESTART_READ_IMPORT, loaded, &result) == 0);
  CHECK(all_ok(memcmp(state, loaded, (X.Check.idim_max + 1) * sizeof(*state)) == 0));
  CHECK(RestartReadVector(&X, row, 1, RESTART_READ_IMPORT, loaded, &result) != 0);

  /* Header mutations are rejected during classification before recovery. */
  const int words[] = {RW_VERSION, RW_PHASE, RW_COMPONENT, RW_STATE, RW_STEP, RW_FLAGS,
                       RW_RANK, RW_OFFSET, RW_LOCAL_DIM, RW_SCALAR_A};
  const uint64_t values[] = {99, 1, 1, 1, INT_MAX, 2, 1, 1, UINT64_MAX,
                             UINT64_C(0x7ff8000000000000)};
  for (size_t i = 0; i < sizeof(words)/sizeof(words[0]); ++i) {
    uint64_t old = 0, unused = 0;
    CHECK(replace_word(words[i], values[i], &old));
    CHECK(RestartClassifySet(&X, &m, 0, NULL, &source, &legacy_step) != 0);
    CHECK(rejected_read(&X, row));
    CHECK(replace_word(words[i], old, &unused));
  }
  /* Flip the occupied amplitude's sign: norm stays exactly 1, digest changes. */
  uint64_t old = 0, unused = 0;
  CHECK(replace_word(RESTART_V3_WORDS, RestartDoubleWord(-1), &old));
  CHECK(rejected_read(&X, row));
  CHECK(replace_word(RESTART_V3_WORDS, old, &unused));
  /* A trailing byte is rejected even on the empty owner's file. */
  int bad = 0;
  if (myrank == nproc - 1) {
    char path[D_FileNameMax];
    snprintf(path, sizeof(path), "%smodule_restart_set0_comp0_gen0_rank_%d.dat", folder, myrank);
    FILE *fp = fopen(path, "ab");
    bad = fp == NULL;
    if (fp != NULL) { bad |= fputc(42, fp) == EOF; bad |= fclose(fp) != 0; }
  }
  CHECK(all_ok(!bad));
  CHECK(RestartClassifySet(&X, &m, 0, NULL, &source, &legacy_step) != 0);
  if (myrank == nproc - 1) {
    char path[D_FileNameMax];
    snprintf(path, sizeof(path), "%smodule_restart_set0_comp0_gen0_rank_%d.dat", folder, myrank);
    bad = truncate(path, RESTART_V3_HEADER_BYTES + 16 * X.Check.idim_max) != 0;
  }
  CHECK(all_ok(!bad));
  /* A rank-local bad direction must collectively fail before publication;
   * the previously published manifest and files remain usable. */
  if (myrank == 0) direction[1] = NAN;
  info.step = 5;
  CHECK(RestartWriteUnit(&X, &unit, vectors, &info, 0) != 0);
  direction[1] = 0;
  CHECK(RestartReadVector(&X, row, 0, RESTART_READ_CONTINUATION, loaded, &result) == 0);
  CHECK(result.step == 4);
  CHECK(RestartWriteUnit(&X, &unit, vectors, &info, 0) == 0);
  RestartManifestFree(&m);
  CHECK(RestartManifestLoad(&X, &m, &exists) == 0 && exists);
  row = RestartManifestFind(&m, 0);
  CHECK(row != NULL && row->gen == 1 && row->step == 5);
  CHECK(RestartClassifySet(&X, &m, 0, NULL, &source, &legacy_step) == 0);
  RestartManifestFree(&m);
  /* Retain the tiny artifacts for diagnostics. */
  if (myrank == 0) puts("restart_checkpoint_check: ok (expected rejection diagnostics above)");
#ifdef MPI
  CHECK(MPI_Finalize() == MPI_SUCCESS);
#endif
  return 0;
}
