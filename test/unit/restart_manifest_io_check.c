#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "Common.h"
#include "restart_manifest.h"
#include "restart_check.h"
#include "wrapperMPI.h"

#ifdef MPI
#include <mpi.h>
#endif

static int mpi_ready = 0;

static int fatal_check(int condition, const char *expression,
                       const char *file, int line)
{
  if (condition) return 1;
  fprintf(stderr, "FAIL %s:%d: %s\n", file, line, expression);
  fflush(stderr);
#ifdef MPI
  if (mpi_ready) MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  return 0;
}

/* A rank-local unit assertion must stop every MPI peer rather than return
 * from one rank while the others enter the next collective. */
#undef CHECK
#define CHECK(condition) do { \
  if (!fatal_check((condition), #condition, __FILE__, __LINE__)) return 1; \
} while (0)

static int all_ok(int condition)
{
  return SumMPI_i(condition ? 0 : 1) == 0;
}

static void set_failure(const char *stage, int rank)
{
  char text[32];
  snprintf(text, sizeof(text), "%d", rank);
  setenv("HPHI_TEST_RESTART_MANIFEST_FAIL", stage, 1);
  setenv("HPHI_TEST_RESTART_MANIFEST_FAIL_RANK", text, 1);
}

static void clear_failure(void)
{
  unsetenv("HPHI_TEST_RESTART_MANIFEST_FAIL");
  unsetenv("HPHI_TEST_RESTART_MANIFEST_FAIL_RANK");
}

static int root_write(const char *path, const char *text)
{
  int failed = 0;
  if (myrank == 0) {
    FILE *fp = fopen(path, "w");
    failed = fp == NULL || fputs(text, fp) == EOF;
    if (fp != NULL && fclose(fp) != 0) failed = 1;
  }
  return all_ok(!failed);
}

int main(int argc, char **argv)
{
  struct BindStruct X;
  struct RestartManifest m;
  char directory[D_FileNameMax] = {0};
  char folder[D_FileNameMax], path[D_FileNameMax], partial[D_FileNameMax];
  int exists, failed = 0;
#ifdef MPI
  CHECK(MPI_Init(&argc, &argv) == MPI_SUCCESS);
  mpi_ready = 1;
  CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &myrank) == MPI_SUCCESS);
  CHECK(MPI_Comm_size(MPI_COMM_WORLD, &nproc) == MPI_SUCCESS);
  if (argc == 2 && strcmp(argv[1], "--rank-local-assertion") == 0) {
    if (myrank == nproc - 1)
      fprintf(stderr, "intentional rank-local assertion on rank %d\n", myrank);
    CHECK(myrank != nproc - 1);
  }
#else
  (void)argc;
  (void)argv;
  myrank = 0;
  nproc = 1;
#endif
  stdoutMPI = stdout;
  memset(&X, 0, sizeof(X));
  X.Def.CDataFileHead = "manifest_io";
  if (myrank == 0) {
    snprintf(directory, sizeof(directory), "restart_manifest_io_check.%ld",
             (long)getpid());
  }
#ifdef MPI
  CHECK(MPI_Bcast(directory, (int)sizeof(directory), MPI_CHAR, 0,
                  MPI_COMM_WORLD) == MPI_SUCCESS);
#endif
  if (myrank == 0 && mkdir(directory, 0700) != 0) failed = 1;
  CHECK(all_ok(!failed));
  CHECK(all_ok(snprintf(folder, sizeof(folder), "%s/", directory) > 0));
  cParentOutputFolder = folder;
  CHECK(all_ok(snprintf(path, sizeof(path), "%s/manifest_io_restart_manifest.dat",
                        directory) > 0));
  CHECK(all_ok(snprintf(partial, sizeof(partial), "%s.part", path) > 0));

  CHECK(RestartManifestLoad(&X, &m, &exists) == 0);
  CHECK(all_ok(exists == 0));
  RestartManifestFree(&m);

  CHECK(root_write(path,
      "hphi_restart_manifest 1\n"
      "method cTPQ\n"
      "basis sector\n"
      "set 3 ranks 4 gen 1 step 9 ncomp 2\n"));
  BarrierMPI();
  CHECK(RestartManifestLoad(&X, &m, &exists) == 0);
  CHECK(all_ok(exists == 1 && m.method == cTPQ && m.basis == 1 && m.nrows == 1));
  const struct RestartManifestRow *found = RestartManifestFind(&m, 3);
  CHECK(all_ok(found != NULL));
  CHECK(all_ok(found->ranks == 4));
  CHECK(all_ok(found->ncomp == 2));
  RestartManifestFree(&m);

  memset(&m, 0, sizeof(m));
  m.method = CG;
  m.basis = 0;
  struct RestartManifestRow row = {7, (uint64_t)nproc, 0, 12, 1};
  CHECK(all_ok(RestartManifestUpsert(&m, &row) == 0));
  CHECK(RestartManifestPublish(&X, &m) == 0);
  RestartManifestFree(&m);
  BarrierMPI();
  if (myrank == 0 && access(partial, F_OK) == 0) failed = 1;
  CHECK(all_ok(!failed));
  CHECK(RestartManifestLoad(&X, &m, &exists) == 0);
  CHECK(all_ok(exists == 1 && m.method == CG && m.basis == 0 && m.nrows == 1));
  found = RestartManifestFind(&m, 7);
  CHECK(all_ok(found != NULL));
  CHECK(all_ok(found->ranks == (uint64_t)nproc));
  RestartManifestFree(&m);

  /* Non-ENOENT input open failure and receiving-rank allocation failure are
   * collective: every rank returns, then a clean load still succeeds. */
  set_failure("load_open", 0);
  CHECK(RestartManifestLoad(&X, &m, &exists) != 0);
  clear_failure();
  CHECK(RestartManifestLoad(&X, &m, &exists) == 0);
  CHECK(all_ok(exists && RestartManifestFind(&m, 7) != NULL));
  RestartManifestFree(&m);
#ifdef MPI
  if (nproc > 1) {
    set_failure("receive_alloc", nproc - 1);
    CHECK(RestartManifestLoad(&X, &m, &exists) != 0);
    clear_failure();
    CHECK(RestartManifestLoad(&X, &m, &exists) == 0);
    CHECK(all_ok(exists && RestartManifestFind(&m, 7) != NULL));
    RestartManifestFree(&m);
  }
#endif

  /* Each root-side publication failure preserves the old manifest and
   * removes its temporary file. */
  const char *publish_failures[] = {"flush", "close", "publish"};
  for (size_t i = 0; i < sizeof(publish_failures)/sizeof(publish_failures[0]); ++i) {
    memset(&m, 0, sizeof(m));
    m.method = cTPQ; m.basis = 1;
    struct RestartManifestRow newer = {9, (uint64_t)nproc, 1, 20, 2};
    CHECK(all_ok(RestartManifestUpsert(&m, &newer) == 0));
    set_failure(publish_failures[i], 0);
    CHECK(RestartManifestPublish(&X, &m) != 0);
    clear_failure();
    RestartManifestFree(&m);
    CHECK(all_ok(access(partial, F_OK) != 0));
    CHECK(RestartManifestLoad(&X, &m, &exists) == 0);
    CHECK(all_ok(exists && m.method == CG && m.basis == 0 &&
                 RestartManifestFind(&m, 7) != NULL &&
                 RestartManifestFind(&m, 9) == NULL));
    RestartManifestFree(&m);
  }

  memset(&m, 0, sizeof(m));
  m.method = CG;
  m.basis = 9;
  CHECK(RestartManifestPublish(&X, &m) != 0);
  RestartManifestFree(&m);
  CHECK(RestartManifestLoad(&X, &m, &exists) == 0);
  CHECK(all_ok(exists == 1 && m.method == CG && m.basis == 0 && m.nrows == 1));
  RestartManifestFree(&m);

  CHECK(root_write(path,
      "hphi_restart_manifest 1\nmethod CG\nbasis raw\ninvalid\n"));
  BarrierMPI();
  CHECK(RestartManifestLoad(&X, &m, &exists) != 0);
  CHECK(all_ok(exists == 0 && m.rows == NULL && m.nrows == 0));
  RestartManifestFree(&m);

  if (myrank == 0) {
    if (remove(path) != 0 || rmdir(directory) != 0) failed = 1;
  }
  CHECK(all_ok(!failed));
  if (myrank == 0) puts("restart_manifest_io_check: ok");
#ifdef MPI
  int mpi_result = MPI_Finalize();
  mpi_ready = 0;
  CHECK(mpi_result == MPI_SUCCESS);
#endif
  return 0;
}
