#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "Common.h"
#include "restart_manifest.h"
#include "wrapperMPI.h"

#ifdef MPI
#include <mpi.h>
#endif

enum { MANIFEST_OK = 0, MANIFEST_ABSENT = 1, MANIFEST_ERROR = 2 };

#ifdef HPHI_TEST_RESTART_IO
static int manifest_test_fail(const char *stage)
{
  const char *wanted = getenv("HPHI_TEST_RESTART_MANIFEST_FAIL");
  const char *rank = getenv("HPHI_TEST_RESTART_MANIFEST_FAIL_RANK");
  return wanted != NULL && strcmp(wanted, stage) == 0 &&
         (rank == NULL || atoi(rank) == myrank);
}
#endif

static void manifest_init(struct RestartManifest *m)
{
  memset(m, 0, sizeof(*m));
  m->method = -1;
  m->basis = -1;
}

static int manifest_path(
    const struct BindStruct *X, char path[D_FileNameMax], int partial)
{
  char name[D_FileNameMax];
  int length;
  if (X == NULL || X->Def.CDataFileHead == NULL ||
      cParentOutputFolder == NULL || cFileNameRestartManifest == NULL)
    return -1;
  length = snprintf(name, sizeof(name), cFileNameRestartManifest,
                    X->Def.CDataFileHead);
  if (length < 0 || (size_t)length >= sizeof(name)) return -1;
  length = snprintf(path, D_FileNameMax, partial ? "%s%s.part" : "%s%s",
                    cParentOutputFolder, name);
  return length < 0 || length >= D_FileNameMax ? -1 : 0;
}

static void manifest_error(const char *stage)
{
  FILE *out = stdoutMPI != NULL ? stdoutMPI : stderr;
  fprintf(out, "Error: restart manifest %s failed.\n", stage);
  fflush(out);
}

static int root_load(
    const struct BindStruct *X, struct RestartManifest *m, int *status)
{
  char path[D_FileNameMax];
  FILE *fp;
  if (manifest_path(X, path, 0) != 0) {
    *status = MANIFEST_ERROR;
    return -1;
  }
  errno = 0;
#ifdef HPHI_TEST_RESTART_IO
  if (manifest_test_fail("load_open")) {
    errno = EACCES;
    fp = NULL;
  } else
#endif
  fp = fopen(path, "r");
  if (fp == NULL) {
    *status = errno == ENOENT ? MANIFEST_ABSENT : MANIFEST_ERROR;
    return *status == MANIFEST_ABSENT ? 0 : -1;
  }
  *status = RestartManifestParse(fp, m) == 0 ? MANIFEST_OK : MANIFEST_ERROR;
  if (fclose(fp) != 0) *status = MANIFEST_ERROR;
  return *status == MANIFEST_OK ? 0 : -1;
}

#ifdef MPI
static int broadcast_words(uint64_t *words, uint64_t count)
{
  if (count > INT_MAX) return -1;
  return MPI_Bcast(words, (int)count, MPI_UINT64_T, 0, MPI_COMM_WORLD) ==
                 MPI_SUCCESS
             ? 0
             : -1;
}

static int pack_manifest(
    const struct RestartManifest *m, uint64_t **wire, uint64_t *count)
{
  size_t nwords;
  if (m->nrows > ((size_t)INT_MAX - 3) / 5) return -1;
  nwords = 3 + 5 * m->nrows;
  *wire = malloc(nwords * sizeof(**wire));
  if (*wire == NULL) return -1;
  (*wire)[0] = (uint64_t)m->method;
  (*wire)[1] = (uint64_t)m->basis;
  (*wire)[2] = (uint64_t)m->nrows;
  for (size_t i = 0; i < m->nrows; ++i) {
    (*wire)[3 + 5 * i] = m->rows[i].set;
    (*wire)[4 + 5 * i] = m->rows[i].ranks;
    (*wire)[5 + 5 * i] = m->rows[i].gen;
    (*wire)[6 + 5 * i] = m->rows[i].step;
    (*wire)[7 + 5 * i] = m->rows[i].ncomp;
  }
  *count = (uint64_t)nwords;
  return 0;
}

static int unpack_manifest(
    const uint64_t *wire, uint64_t count, struct RestartManifest *m)
{
  struct RestartManifest received;
  int invalid;
  manifest_init(&received);
  invalid = count < 3 || (count - 3) % 5 != 0 ||
            wire[2] != (count - 3) / 5 || wire[0] > INT_MAX ||
            wire[1] > INT_MAX;
  if (!invalid) {
    received.method = (int)wire[0];
    received.basis = (int)wire[1];
    for (uint64_t i = 0; i < wire[2]; ++i) {
      struct RestartManifestRow row = {
        wire[3 + 5 * i], wire[4 + 5 * i], wire[5 + 5 * i],
        wire[6 + 5 * i], wire[7 + 5 * i]
      };
      if (row.ranks < 1 || row.gen > 1 || row.ncomp < 1 || row.ncomp > 2 ||
          RestartManifestFind(&received, row.set) != NULL ||
          RestartManifestUpsert(&received, &row) != 0) {
        invalid = 1;
        break;
      }
    }
  }
  if (SumMPI_i(invalid) != 0) {
    RestartManifestFree(&received);
    return -1;
  }
  RestartManifestFree(m);
  *m = received;
  return 0;
}
#endif

int RestartManifestLoad(
    const struct BindStruct *X, struct RestartManifest *m, int *exists)
{
  int status = MANIFEST_OK;
#ifdef MPI
  if (SumMPI_i(m == NULL || exists == NULL) != 0) {
    manifest_error("load arguments");
    return -1;
  }
#else
  if (m == NULL || exists == NULL) return -1;
#endif
  manifest_init(m);
  *exists = 0;
#ifdef MPI
  if (myrank == 0) root_load(X, m, &status);
  status = BcastMPI_i(0, status);
#else
  root_load(X, m, &status);
#endif
  if (status == MANIFEST_ABSENT) return 0;
  if (status != MANIFEST_OK) {
    RestartManifestFree(m);
    manifest_error("load");
    return -1;
  }

#ifdef MPI
  uint64_t *wire = NULL, count = 0;
  int invalid = 0;
  if (myrank == 0 && pack_manifest(m, &wire, &count) != 0) invalid = 1;
  invalid = BcastMPI_i(0, invalid);
  if (invalid) {
    free(wire);
    RestartManifestFree(m);
    manifest_error("load allocation");
    return -1;
  }
  invalid = MPI_Bcast(&count, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD) != MPI_SUCCESS;
  if (SumMPI_i(invalid) != 0) {
    free(wire);
    RestartManifestFree(m);
    manifest_error("load length");
    return -1;
  }
  invalid = count < 3 || count > INT_MAX ||
            count > SIZE_MAX / sizeof(*wire);
  if (SumMPI_i(invalid) != 0) {
    free(wire);
    RestartManifestFree(m);
    manifest_error("load length");
    return -1;
  }
  if (myrank != 0) {
#ifdef HPHI_TEST_RESTART_IO
    if (!manifest_test_fail("receive_alloc"))
#endif
      wire = malloc((size_t)count * sizeof(*wire));
  }
  if (SumMPI_i(wire == NULL) != 0) {
    free(wire);
    RestartManifestFree(m);
    manifest_error("load allocation");
    return -1;
  }
  invalid = broadcast_words(wire, count) != 0;
  if (SumMPI_i(invalid) != 0) {
    free(wire);
    RestartManifestFree(m);
    manifest_error("load broadcast");
    return -1;
  }
  if (unpack_manifest(wire, count, m) != 0) {
    free(wire);
    RestartManifestFree(m);
    manifest_error("load broadcast");
    return -1;
  }
  free(wire);
#endif
  *exists = 1;
  return 0;
}

int RestartManifestPublish(
    const struct BindStruct *X, const struct RestartManifest *m)
{
  int failed = 0;
#ifdef MPI
  if (myrank == 0) {
#endif
    char path[D_FileNameMax], partial[D_FileNameMax];
    FILE *fp = NULL;
    const char *stage = "path";
    if (manifest_path(X, path, 0) != 0 ||
        manifest_path(X, partial, 1) != 0) {
      failed = 1;
    } else {
      stage = "open";
      fp = fopen(partial, "w");
      if (fp == NULL) {
        failed = 1;
      } else {
        stage = "format";
        if (RestartManifestFormat(m, fp) != 0) failed = 1;
        if (!failed) {
          stage = "flush";
#ifdef HPHI_TEST_RESTART_IO
          if (manifest_test_fail("flush")) failed = 1; else
#endif
          if (fflush(fp) != 0) failed = 1;
        }
        int close_failed = fclose(fp) != 0;
#ifdef HPHI_TEST_RESTART_IO
        close_failed |= manifest_test_fail("close");
#endif
        if (close_failed && !failed) {
          stage = "close";
          failed = 1;
        }
        fp = NULL;
        if (!failed) {
          stage = "publish";
#ifdef HPHI_TEST_RESTART_IO
          if (manifest_test_fail("publish")) failed = 1; else
#endif
          if (rename(partial, path) != 0) failed = 1;
        }
        if (failed) remove(partial);
      }
    }
    if (failed) manifest_error(stage);
#ifdef MPI
  }
  failed = BcastMPI_i(0, failed);
#endif
  return failed ? -1 : 0;
}
