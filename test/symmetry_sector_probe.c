/* Test-only observation of production sector operations; never linked into HPhi. */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
struct BindStruct;
#include "struct.h"
#include "global.h"
#include "wrapperMPI.h"
#include "mltply.h"
#include "symmetry_basis.h"
#include "symmetry_observables.h"
#include "expec_totalspin.h"
#include "symmetry_sector_probe.h"

static const char *probe_action(const struct BindStruct *X, int *legacy)
{
  const char *action = getenv("HPHI_TEST_SYMMETRY_ACTION");
  *legacy = 0;
  if (action == NULL && X != NULL && X->Def.iCalcModel == SpinGC) {
    action = getenv("HPHI_TEST_SPINGC_ACTION");
    *legacy = action != NULL;
  }
  return action;
}

static int ready(const struct BindStruct *X, const double complex *vector)
{
  return X != NULL && vector != NULL && X->Def.iFlgSymmetryBasis && X->Sym != NULL &&
         SymmetryBasisOwnedStorageReady(X->Sym, X->Sym->local_dim) == TRUE;
}

static int read_vector(const struct BindStruct *X, double complex *vector)
{
  char filename[128];
  FILE *fp;
  int error = 0;
  snprintf(filename, sizeof(filename), "probe-input.rank%d.dat", myrank);
  fp = fopen(filename, "r");
  if (fp == NULL) return -1;
  for (unsigned long j = 1; j <= X->Sym->local_dim; ++j) {
    double re, im;
    if (fscanf(fp, "%lf %lf", &re, &im) != 2 || !isfinite(re) || !isfinite(im)) {
      error = 1;
      break;
    }
    vector[j] = re + I*im;
  }
  if (!error) {
    int next;
    do next = fgetc(fp); while (next != EOF && isspace(next));
    if (next != EOF || ferror(fp)) error = 1;
  }
  if (fclose(fp)) error = 1;
  return error ? -1 : 0;
}

static int write_info(const struct SymmetryBasisRuntime *sym, const char *stem,
                      unsigned long elements, double prenorm)
{
  char filename[256];
  FILE *fp;
  int error = 0;
  snprintf(filename, sizeof(filename), "%s.info", stem);
  fp = fopen(filename, "w");
  if (fp == NULL) return -1;
  if (fprintf(fp, "dim=%lu\noffset=%lu\nlocal_dim=%lu\n"
              "raw_basis_elements=%llu\nraw_diagonal_elements=%llu\n"
              "initial_vector_elements=%llu\nprobe_vector_elements=%lu\nprenorm=%.17g\n",
              sym->dim, sym->local_offset, sym->local_dim,
              sym->allocation_raw_basis_list_elements,
              sym->allocation_raw_diagonal_elements,
              sym->allocation_initial_vector_elements, elements, prenorm) < 0) error = 1;
  if (fclose(fp)) error = 1;
  return error ? -1 : 0;
}

static int write_vector(const struct BindStruct *X, const double complex *vector,
                        const char *stem)
{
  char filename[256];
  FILE *fp;
  int error = 0;
  snprintf(filename, sizeof(filename), "%s.dat", stem);
  fp = fopen(filename, "w");
  if (fp == NULL) return -1;
  for (unsigned long j = 1; j <= X->Sym->local_dim; ++j) {
    const struct SymmetryBasisVector *entry = SymmetryBasisLocalEntry(X->Sym, j);
    if (entry == NULL || !isfinite(creal(vector[j])) || !isfinite(cimag(vector[j])) ||
        fprintf(fp, "%lu %lu %.17g %.17g\n", X->Sym->local_offset+j,
                entry->rep_state, creal(vector[j]), cimag(vector[j])) < 0) {
      error = 1;
      break;
    }
  }
  if (fclose(fp)) error = 1;
  return error ? -1 : 0;
}

int SymmetryProbeWriteVector(const struct BindStruct *X, const double complex *vector,
                             const char *stage, int sample, int step, double prenorm)
{
  const char *capture = getenv("HPHI_TEST_SYMMETRY_CAPTURE");
  char stem[192];
  int error;
  double root_norm;
  if (capture == NULL || strcmp(capture, "1") != 0) return 0;
  error = !ready(X, vector) || stage == NULL || sample < 0 || step < 0 ||
          !isfinite(prenorm) || prenorm <= 0;
  if (!error) error = strcmp(stage, "initial") && strcmp(stage, "final");
  if (SumMPI_i(error)) return -1;
  root_norm = MaxMPI_d(prenorm);
  if (SumMPI_i(prenorm != root_norm)) return -1;
  snprintf(stem, sizeof(stem), "sector_%s_sample%d_step%d.rank%d", stage, sample, step, myrank);
  error = write_vector(X, vector, stem) != 0;
  error |= write_info(X->Sym, stem, 0, prenorm) != 0;
  return SumMPI_i(error) ? -1 : 0;
}

int SymmetryProbeBeforeSolver(struct BindStruct *X)
{
  int legacy, error = 0;
  const char *action = probe_action(X, &legacy);
  const char *poison = getenv("HPHI_TEST_SPINGC_POISON_FIXED");
  struct SymmetryBasisRuntime *sym = X == NULL ? NULL : X->Sym;
  double complex *in = NULL, *out = NULL;
  FILE *matrix = NULL;
  char stem[128], filename[160];
  if (poison != NULL && strcmp(poison, "1") == 0 && X != NULL &&
      X->Def.iCalcModel == SpinGC && X->Def.iFlgSymmetryBasis) {
    X->Def.Nup = 1; X->Def.Ndown = 2; X->Def.Ne = 3;
    fprintf(stdoutMPI, "SpinGCProbe poisoned fixed fields after basis setup.\n");
  }
  if (action == NULL || strcmp(action, "lanczos") == 0) return 0;
  if (strcmp(action, "matvec") && strcmp(action, "apply") && strcmp(action, "moments"))
    return legacy ? 0 : -1;
  error = X == NULL || sym == NULL || !X->Def.iFlgSymmetryBasis ||
          X->Def.iCalcType == FullDiag;
  if (!error) error = sym->dim == 0 || (strcmp(action, "matvec") == 0 && sym->dim > 256);
  if (SumMPI_i(error)) return -1;
  in = calloc(sym->local_dim+1, sizeof(*in));
  out = calloc(sym->local_dim+1, sizeof(*out));
  error = in == NULL || out == NULL;
  if (SumMPI_i(error)) { error = 1; goto cleanup; }
  if (strcmp(action, "matvec")) {
    error = read_vector(X, in) != 0;
    if (SumMPI_i(error)) { error = 1; goto cleanup; }
  }
  if (strcmp(action, "moments") == 0) {
    const char *corrupt = legacy ? getenv("HPHI_TEST_SPINGC_INVALID_STORAGE_RANK") : NULL;
    unsigned long saved = 0;
    unsigned long *capacity = sym->basis_layout == SYMMETRY_BASIS_DISTRIBUTED ?
                              &sym->local_capacity : &sym->capacity;
    if (corrupt != NULL && atoi(corrupt) == myrank) { saved = *capacity; *capacity = 0; }
    error = legacy ? EvaluateSymmetrySpinGCMoments(X, in) != 0 : expec_totalSz(X, in) != 0;
    if (corrupt != NULL && atoi(corrupt) == myrank) *capacity = saved;
    if (SumMPI_i(error)) { error = 1; goto cleanup; }
    if (myrank == 0) {
      fprintf(stdoutMPI, "%s Sz %.17g\n%s Sz2 %.17g\n",
              legacy ? "SpinGCProbe" : "SymmetryProbe", X->Phys.Sz,
              legacy ? "SpinGCProbe" : "SymmetryProbe", X->Phys.Sz2);
      if (!legacy)
        fprintf(stdoutMPI, "SymmetryProbe num %.17g\nSymmetryProbe num2 %.17g\n"
                "SymmetryProbe doublon %.17g\n", X->Phys.num, X->Phys.num2, X->Phys.doublon);
      if (fflush(stdoutMPI)) error = 1;
    }
    if (legacy) goto cleanup;
    memcpy(out, in, (sym->local_dim+1)*sizeof(*in));
  }
  snprintf(stem, sizeof(stem), "%s_probe_rank_%d", legacy ? "spingc" : "sector", myrank);
  error = write_info(sym, stem, 2*(sym->local_dim+1), 1.0) != 0;
  if (SumMPI_i(error)) { error = 1; goto cleanup; }
  if (strcmp(action, "matvec")) {
    if (strcmp(action, "apply") == 0) error = mltply(X, out, in) != 0;
    if (SumMPI_i(error)) { error = 1; goto cleanup; }
    error = write_vector(X, out, stem) != 0;
    goto cleanup;
  }
  snprintf(filename, sizeof(filename), "%s.dat", stem);
  matrix = fopen(filename, "w");
  error = matrix == NULL;
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
      if (entry == NULL || !isfinite(creal(out[row])) || !isfinite(cimag(out[row])) ||
          fprintf(matrix, "%lu %lu %lu %.17g %.17g\n", sym->local_offset+row, col,
                  entry->rep_state, creal(out[row]), cimag(out[row])) < 0) { error = 1; break; }
    }
    if (fflush(matrix)) error = 1;
    if (SumMPI_i(error)) { error = 1; goto cleanup; }
  }
cleanup:
  if (matrix != NULL && fclose(matrix)) error = 1;
  free(in); free(out);
  if (SumMPI_i(error)) {
    if (legacy && strcmp(action, "moments") == 0)
      fprintf(stdoutMPI, "Error: SpinGC moments probe failed.\n");
    return -1;
  }
  return 1;
}

int SymmetryProbeFinalVector(const struct BindStruct *X, const double complex *vector)
{
  int legacy, error;
  const char *action = probe_action(X, &legacy);
  char stem[128];
  if (action == NULL || strcmp(action, "lanczos")) return 0;
  error = !ready(X, vector);
  if (SumMPI_i(error)) return -1;
  snprintf(stem, sizeof(stem), legacy ? "spingc-lanczos.rank%d" : "sector_probe_rank_%d", myrank);
  error = write_vector(X, vector, stem) != 0;
  if (!legacy) error |= write_info(X->Sym, stem, 0, 1.0) != 0;
  if (SumMPI_i(error)) return -1;
  return SymmetryProbeWriteVector(X, vector, "final", 0, 0, 1.0);
}
