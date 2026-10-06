/* Test-only observation of the production sector matvec before any solver. */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "struct.h"
#include "global.h"
#include "wrapperMPI.h"
#include "mltply.h"
#include "symmetry_basis.h"
#include "symmetry_observables.h"
#include "symmetry_spingc_probe.h"

int SpinGCProbeBeforeSolver(struct BindStruct *X)
{
  const char *action = getenv("HPHI_TEST_SPINGC_ACTION");
  struct SymmetryBasisRuntime *sym = X->Sym;
  double complex *in = NULL, *out = NULL;
  FILE *matrix = NULL, *info = NULL;
  char filename[128];
  int error = 0;
  if (action == NULL) return 0;
  if (strcmp(action, "moments") == 0) {
    double complex *vector = NULL;
    FILE *input = NULL;
    unsigned long int saved_capacity = 0UL;
    int corrupt_storage = 0;
    const char *corrupt_rank = getenv("HPHI_TEST_SPINGC_INVALID_STORAGE_RANK");
    error = X == NULL || sym == NULL || X->Def.iCalcModel != SpinGC ||
            !X->Def.iFlgSymmetryBasis;
    if (SumMPI_i(error)) goto moments_failure;
    vector = calloc(sym->local_dim + 1UL, sizeof(*vector));
    error = vector == NULL;
    if (SumMPI_i(error)) goto moments_cleanup;
    snprintf(filename, sizeof(filename), "probe-input.rank%d.dat", myrank);
    input = fopen(filename, "r");
    error = input == NULL;
    if (SumMPI_i(error)) goto moments_cleanup;
    for (unsigned long int j = 1; j <= sym->local_dim; ++j) {
      double real, imag;
      if (fscanf(input, "%lf %lf", &real, &imag) != 2) {
        error = 1;
        break;
      }
      vector[j] = real + imag*I;
    }
    if (!error) {
      int next;
      do next = fgetc(input); while (next != EOF && isspace(next));
      if (next != EOF || ferror(input)) error = 1;
    }
    if (fclose(input) != 0) error = 1;
    input = NULL;
    if (SumMPI_i(error)) goto moments_cleanup;
    if (corrupt_rank != NULL && atoi(corrupt_rank) == myrank) {
      if (sym->basis_layout == SYMMETRY_BASIS_DISTRIBUTED) {
        saved_capacity = sym->local_capacity;
        sym->local_capacity = 0UL;
      } else {
        saved_capacity = sym->capacity;
        sym->capacity = 0UL;
      }
      corrupt_storage = 1;
    }
    error = EvaluateSymmetrySpinGCMoments(X, vector) != 0;
    if (corrupt_storage) {
      if (sym->basis_layout == SYMMETRY_BASIS_DISTRIBUTED)
        sym->local_capacity = saved_capacity;
      else
        sym->capacity = saved_capacity;
    }
    if (SumMPI_i(error)) goto moments_cleanup;
    if (myrank == 0) {
      if (fprintf(stdoutMPI, "SpinGCProbe Sz %.17g\nSpinGCProbe Sz2 %.17g\n",
                  X->Phys.Sz, X->Phys.Sz2) < 0) error = 1;
      if (fflush(stdoutMPI) != 0) error = 1;
    }
moments_cleanup:
    if (input != NULL && fclose(input) != 0) error = 1;
    free(vector);
moments_failure:
    if (SumMPI_i(error)) {
      fprintf(stdoutMPI, "Error: SpinGC moments probe failed.\n");
      return -1;
    }
    return 1;
  }
  if (strcmp(action, "matvec") != 0) return 0;
  error = X->Def.iCalcModel != SpinGC || !X->Def.iFlgSymmetryBasis ||
          X->Def.iCalcType == FullDiag || sym == NULL;
  if (!error) error = sym->dim > 256 || sym->dim == 0;
  if (SumMPI_i(error)) return -1;
  in = calloc(sym->local_dim + 1, sizeof(*in));
  out = calloc(sym->local_dim + 1, sizeof(*out));
  error = in == NULL || out == NULL;
  if (SumMPI_i(error)) { error = 1; goto cleanup; }
  snprintf(filename, sizeof(filename), "spingc_probe_rank_%d.dat", myrank);
  matrix = fopen(filename, "w");
  snprintf(filename, sizeof(filename), "spingc_probe_rank_%d.info", myrank);
  info = fopen(filename, "w");
  error = matrix == NULL || info == NULL;
  if (SumMPI_i(error)) { error = 1; goto cleanup; }
  if (fprintf(info, "dim=%lu\noffset=%lu\nlocal_dim=%lu\n"
              "raw_basis_elements=%llu\nraw_diagonal_elements=%llu\n"
              "initial_vector_elements=%llu\nprobe_vector_elements=%lu\n",
              sym->dim, sym->local_offset, sym->local_dim,
              sym->allocation_raw_basis_list_elements,
              sym->allocation_raw_diagonal_elements,
              sym->allocation_initial_vector_elements,
              2*(sym->local_dim+1)) < 0) error = 1;
  if (fclose(info)) error = 1;
  info = NULL;
  if (SumMPI_i(error)) { error = 1; goto cleanup; }
  for (unsigned long col = 1; col <= sym->dim; ++col) {
    memset(in, 0, (sym->local_dim+1)*sizeof(*in));
    memset(out, 0, (sym->local_dim+1)*sizeof(*out));
    if (col > sym->local_offset && col <= sym->local_offset+sym->local_dim)
      in[col-sym->local_offset] = 1;
    error = mltply(X, out, in) != 0;
    if (SumMPI_i(error)) { error = 1; goto cleanup; }
    for (unsigned long row = 1; row <= sym->local_dim; ++row) {
      const struct SymmetryBasisVector *entry = SymmetryBasisLocalEntry(sym, row);
      if (entry == NULL) { error = 1; break; }
      if (fprintf(matrix, "%lu %lu %lu %.17g %.17g\n",
                  sym->local_offset+row, col, entry->rep_state,
                  creal(out[row]), cimag(out[row])) < 0) error = 1;
    }
    if (fflush(matrix)) error = 1;
    if (SumMPI_i(error)) { error = 1; goto cleanup; }
  }
cleanup:
  if (info != NULL && fclose(info)) error = 1;
  if (matrix != NULL && fclose(matrix)) error = 1;
  free(in);
  free(out);
  if (SumMPI_i(error)) return -1;
  return 1;
}

int SpinGCProbeFinalVector(const struct BindStruct *X,
                           const double complex *vec)
{
  const char *action = getenv("HPHI_TEST_SPINGC_ACTION");
  const struct SymmetryBasisRuntime *sym;
  FILE *output = NULL;
  char filename[128];
  int error;
  if (action == NULL || strcmp(action, "lanczos") != 0) return 0;
  error = X == NULL || vec == NULL;
  sym = error ? NULL : X->Sym;
  if (!error)
    error = X->Def.iCalcModel != SpinGC || !X->Def.iFlgSymmetryBasis ||
            SymmetryBasisOwnedStorageReady(sym, sym->local_dim) != TRUE;
  if (SumMPI_i(error)) return -1;
  snprintf(filename, sizeof(filename), "spingc-lanczos.rank%d.dat", myrank);
  output = fopen(filename, "w");
  error = output == NULL;
  if (SumMPI_i(error)) goto cleanup;
  for (unsigned long int j = 1; j <= sym->local_dim; ++j) {
    const struct SymmetryBasisVector *entry = SymmetryBasisLocalEntry(sym, j);
    double real = creal(vec[j]);
    double imag = cimag(vec[j]);
    if (entry == NULL || !isfinite(real) || !isfinite(imag) ||
        fprintf(output, "%lu %lu %.17g %.17g\n",
                sym->local_offset+j, entry->rep_state, real, imag) < 0) {
      error = 1;
      break;
    }
  }
  if (output != NULL && fflush(output) != 0) error = 1;
cleanup:
  if (output != NULL && fclose(output) != 0) error = 1;
  if (SumMPI_i(error)) return -1;
  return 0;
}
