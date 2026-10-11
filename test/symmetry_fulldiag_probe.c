/* Test-only sector FullDiag matrix/state contracts; never linked into HPhi. */
#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
struct BindStruct;
#include "struct.h"
#include "global.h"
#include "wrapperMPI.h"
#include "makeHamSym.h"
#include "symmetry_basis.h"
#include "symmetry_checkpoint.h"
#include "symmetry_hamiltonian_io.h"
#include "symmetry_fulldiag_probe.h"

/* Column split of the ELPA panel, NC = ceil(D/P); used for every storage so
 * that each lower-triangle element is dumped by exactly one rank. */
static void dump_columns(unsigned long dim, unsigned long *first, unsigned long *last)
{
  unsigned long chunk = (dim + (unsigned long)nproc - 1) / (unsigned long)nproc;
  *first = (unsigned long)myrank * chunk + 1;
  *last = *first + chunk - 1;
  if (*last > dim) *last = dim;
}

static int owns(unsigned long col)
{
  return !iHamPanelActive || ((long)col >= HamColBegin && (long)col <= HamColEnd);
}

static double complex *element(unsigned long row, unsigned long col)
{
  if (!iHamPanelActive) return &Ham[row][col];
  return &Ham_local[(size_t)((long)col - HamColBegin) * (size_t)HamPanelLd + (row - 1)];
}

/* Overwrite with an exactly representable Hermitian matrix: diagonal i/8,
 * neighbours .25+.125i below the diagonal and the conjugate above it. */
static void write_known_matrix(unsigned long dim)
{
  for (unsigned long col = 1; col <= dim; ++col) {
    if (!owns(col)) continue;
    for (unsigned long row = 1; row <= dim; ++row) *element(row, col) = 0;
    *element(col, col) = (double)col / 8.0;
    if (col < dim) *element(col + 1, col) = .25 + .125*I;
    if (col > 1) *element(col - 1, col) = .25 - .125*I;
  }
}

static int mutate(const struct BindStruct *X, const char *action)
{
  unsigned long dim = X->Sym->dim;
  if (strcmp(action, "fulldiag-matrix") == 0) return 0;
  if (strcmp(action, "fulldiag-matrix-known") == 0) { write_known_matrix(dim); return 0; }
  if (strcmp(action, "fulldiag-matrix-nan") == 0) {
    if (owns(1)) *element(1, 1) = NAN;
    return 0;
  }
  if (strcmp(action, "fulldiag-matrix-nonhermitian") == 0) {
    if (dim < 2) return -1;
    if (owns(2)) *element(1, 2) += .01;
    return 0;
  }
  if (strcmp(action, "fulldiag-matrix-diagonal") == 0) {
    if (owns(1)) *element(1, 1) += .01*I;
    return 0;
  }
  if (strcmp(action, "fulldiag-matrix-diagonal-small") == 0) {
    if (owns(1)) *element(1, 1) += 1e-14*I;
    return 0;
  }
  if (strcmp(action, "fulldiag-matrix-rank-mismatch") == 0) {
    /* Each rank alone is finite and Hermitian; only the rank agreement fails. */
    if (iHamPanelActive || nproc < 2) return -1;
    if (myrank == 1) *element(1, 1) += .125;
    return 0;
  }
  return -1;
}

static int write_contract(const struct SymmetryHamiltonianInfo *info)
{
  const struct SymmetryHamiltonianIdentity *id = &info->identity;
  FILE *fp;
  int error = 0;
  if (myrank != 0) return 0;
  fp = fopen("matrix-contract.dat", "w");
  if (fp == NULL) return -1;
  if (fprintf(fp, "model=%" PRIu64 "\nnsite=%" PRIu64 "\nnup=%" PRIu64 "\nndown=%" PRIu64
              "\nne=%" PRIu64 "\nlocal_site_mask=%" PRIu64 "\nfixed_flags=%" PRIu64
              "\nphase=%" PRIu64 "\nfull_dim=%" PRIu64 "\nsector_dim=%" PRIu64
              "\ngroup_digest=%" PRIu64 "\nsector_count=%" PRIu64 "\nsector_xor=%" PRIu64
              "\nsector_sum=%" PRIu64 "\nbasis_order=%" PRIu64 "\nmatrix_digest=%" PRIu64 "\n",
              id->model, id->nsite, id->nup, id->ndown, id->ne, id->local_site_mask,
              id->fixed_flags, id->phase, id->full_dim, id->sector_dim, id->group_digest,
              id->sector_count, id->sector_xor, id->sector_sum, id->basis_order,
              info->matrix_digest) < 0) error = 1;
  if (fclose(fp) != 0) error = 1;
  return error ? -1 : 0;
}

/* Plain fopen on every rank: childfopenMPI would give non-root ranks /dev/null. */
static int write_rank_matrix(unsigned long dim)
{
  char name[64];
  unsigned long first, last;
  FILE *fp;
  int error = 0;
  snprintf(name, sizeof(name), "matrix-contract.rank%d.dat", myrank);
  fp = fopen(name, "w");
  if (fp == NULL) return -1;
  dump_columns(dim, &first, &last);
  for (unsigned long col = first; !error && col <= last; ++col) {
    if (!owns(col)) { error = 1; break; }
    for (unsigned long row = col; row <= dim; ++row) {
      double complex value = *element(row, col);
      if (fprintf(fp, "%lu %lu %.17g %.17g\n", row, col, creal(value), cimag(value)) < 0) {
        error = 1;
        break;
      }
    }
  }
  if (fclose(fp) != 0) error = 1;
  return error ? -1 : 0;
}

static int matrix_action(struct BindStruct *X, const char *action)
{
  struct SymmetryHamiltonianInfo info;
  if (SumMPI_i(makeHamSym(X) != 0) != 0) {
    fprintf(stdoutMPI, "Error: FullDiag probe could not build the sector Hamiltonian.\n");
    return -1;
  }
  if (SumMPI_i(mutate(X, action) != 0) != 0) {
    fprintf(stdoutMPI, "Error: unknown or inapplicable FullDiag probe action %s.\n", action);
    return -1;
  }
  if (PrepareSymmetryHamiltonian(X, 0, &info) != 0) return -1;
  if (SumMPI_i(write_contract(&info) != 0) != 0 ||
      SumMPI_i(write_rank_matrix(X->Sym->dim) != 0) != 0) {
    fprintf(stdoutMPI, "Error: FullDiag probe could not write the matrix contract.\n");
    return -1;
  }
  return 1;
}

/* Rank-local view for checkpoint I/O; every pointer is borrowed. */
static struct BindStruct local_view(const struct BindStruct *X)
{
  struct BindStruct local = *X;
  local.Check.idim_max = X->Sym->local_dim;
  local.Check.idim_maxMPI = X->Sym->dim;
  return local;
}

static int state_name(char *name, size_t size, const struct BindStruct *X, int contract)
{
  int length = contract
      ? snprintf(name, size, "contract_rank_%d.dat", myrank)
      : snprintf(name, size, "%s_eigenvec_0_rank_%d.dat", X->Def.CDataFileHead, myrank);
  return SumMPI_i(length < 0 || (size_t)length >= size) != 0 ? -1 : 0;
}

/* Write global basis vector 1 with FullDiag provenance, read it back into a
 * separate array, then also publish it under the regular state-0 name. */
static int checkpoint_action(struct BindStruct *X, const char *action)
{
  struct SymmetryHamiltonianInfo info;
  struct SymmetryCheckpointInfo loaded;
  struct BindStruct local = local_view(X);
  double complex *vector = NULL, *readback = NULL;
  char name[D_FileNameMax];
  int status = -1;
  if (SumMPI_i(makeHamSym(X) != 0) != 0 ||
      SumMPI_i(mutate(X, strcmp(action, "fulldiag-checkpoint-known") == 0
                         ? "fulldiag-matrix-known" : "fulldiag-matrix") != 0) != 0 ||
      PrepareSymmetryHamiltonian(X, 0, &info) != 0) return -1;
  vector = calloc(X->Sym->local_dim + 1, sizeof(*vector));
  readback = calloc(X->Sym->local_dim + 1, sizeof(*readback));
  if (SumMPI_i(vector == NULL || readback == NULL) != 0) goto done;
  if (X->Sym->local_offset == 0 && X->Sym->local_dim > 0) vector[1] = 1;
  struct SymmetryCheckpointInfo source = {
    .source_method = FullDiag, .state_index = 0, .step = 0, .time = 0,
    .hamiltonian_digest = info.matrix_digest,
    .hamiltonian_digest_kind = SYMMETRY_HAMILTONIAN_MATRIX_V1
  };
  if (state_name(name, sizeof(name), X, 1) != 0 ||
      WriteSymmetryCheckpoint(&local, name, vector, &source) != 0) goto done;
  memset(&loaded, 0xff, sizeof(loaded));
  if (ReadSymmetryCheckpoint(&local, name, readback, &loaded) != 0) goto done;
  int mismatch = memcmp(vector, readback, (X->Sym->local_dim + 1) * sizeof(*vector)) != 0 ||
                 loaded.source_method != FullDiag || loaded.state_index != 0 ||
                 loaded.step != 0 || loaded.time != 0 ||
                 loaded.hamiltonian_digest != info.matrix_digest ||
                 loaded.hamiltonian_digest_kind != SYMMETRY_HAMILTONIAN_MATRIX_V1;
  if (SumMPI_i(mismatch) != 0) {
    fprintf(stdoutMPI, "Error: FullDiag probe checkpoint round-trip changed the state or provenance.\n");
    goto done;
  }
  fprintf(stdoutMPI, "FullDiagProbe checkpoint round-trip ok\n");
  if (state_name(name, sizeof(name), X, 0) != 0 ||
      WriteSymmetryCheckpoint(&local, name, vector, &source) != 0) goto done;
  if (SumMPI_i(write_contract(&info) != 0) != 0) goto done;
  status = 1;
done:
  free(vector);
  free(readback);
  return status;
}

/* Read the regular state-0 files; a rejected file must not touch the vector. */
static int checkpoint_read_action(struct BindStruct *X)
{
  struct SymmetryCheckpointInfo loaded;
  struct BindStruct local = local_view(X);
  size_t count = X->Sym->local_dim + 1;
  double complex *vector = malloc(count * sizeof(*vector));
  double complex *sentinel = malloc(count * sizeof(*sentinel));
  char name[D_FileNameMax];
  int status = -1;
  if (SumMPI_i(vector == NULL || sentinel == NULL) != 0) goto done;
  for (size_t i = 0; i < count; ++i) vector[i] = sentinel[i] = 7.0 + 7.0*I;
  if (state_name(name, sizeof(name), X, 0) != 0) goto done;
  if (ReadSymmetryCheckpoint(&local, name, vector, &loaded) != 0) {
    if (SumMPI_i(memcmp(vector, sentinel, count * sizeof(*vector)) != 0) == 0)
      fprintf(stdoutMPI, "FullDiagProbe sentinel intact after rejected checkpoint.\n");
    goto done;
  }
  fprintf(stdoutMPI, "FullDiagProbe checkpoint read ok state=%" PRIu64 " kind=%" PRIu64 "\n",
          loaded.state_index, loaded.hamiltonian_digest_kind);
  status = 1;
done:
  free(vector);
  free(sentinel);
  return status;
}

int SymmetryFullDiagTestFail(const char *stage)
{
  static int fired;
  const char *wanted = getenv("HPHI_TEST_FD_FAIL");
  const char *rank_text = getenv("HPHI_TEST_FD_FAIL_RANK");
  if (fired || wanted == NULL || rank_text == NULL ||
      strcmp(wanted, stage) != 0 || atoi(rank_text) != myrank) return 0;
  fired = 1;
  return 1;
}

int SymmetryFullDiagProbe(struct BindStruct *X, const char *action)
{
  int invalid = X == NULL || action == NULL || X->Sym == NULL ||
                X->Def.iCalcType != FullDiag || !X->Def.iFlgSymmetryBasis;
  if (!invalid)
    invalid = X->Sym->enabled != TRUE ||
              X->Sym->basis_layout != SYMMETRY_BASIS_REPLICATED || X->Sym->dim == 0;
  if (SumMPI_i(invalid) != 0) {
    fprintf(stdoutMPI, "Error: FullDiag probe requires a replicated TransSym FullDiag run.\n");
    return -1;
  }
  if (strncmp(action, "fulldiag-matrix", 15) == 0) return matrix_action(X, action);
  if (strcmp(action, "fulldiag-checkpoint") == 0 ||
      strcmp(action, "fulldiag-checkpoint-known") == 0) return checkpoint_action(X, action);
  if (strcmp(action, "fulldiag-checkpoint-read") == 0) return checkpoint_read_action(X);
  fprintf(stdoutMPI, "Error: unknown FullDiag probe action %s.\n", action);
  return -1;
}
