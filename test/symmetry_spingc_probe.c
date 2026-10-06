/* Test-only observation of the production sector matvec before any solver. */
#include <stdlib.h>
#include <string.h>
#include "struct.h"
#include "global.h"
#include "wrapperMPI.h"
#include "mltply.h"
#include "symmetry_basis.h"
#include "symmetry_spingc_probe.h"

int SpinGCProbeBeforeSolver(struct BindStruct *X)
{
  const char *action = getenv("HPHI_TEST_SPINGC_ACTION");
  struct SymmetryBasisRuntime *sym = X->Sym;
  double complex *in = NULL, *out = NULL;
  FILE *matrix = NULL, *info = NULL;
  char filename[128];
  int error = 0;
  if (action == NULL || strcmp(action, "matvec") != 0) return 0;
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
