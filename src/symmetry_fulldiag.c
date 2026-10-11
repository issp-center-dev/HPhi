/* Sector FullDiag state post-processing, one eigenstate at a time.
 *
 * Every rank receives its owned block of the state, so the existing sector
 * correlation and checkpoint code runs unchanged on a rank-local view of the
 * solver context. The view borrows every pointer; only its scalar
 * dimensions and the state number differ from the solver context. */
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include "Common.h"
#include "CalcTime.h"
#include "expec_cisajs.h"
#include "expec_cisajscktaltdc.h"
#include "green_output.h"
#include "nbody_correlation.h"
#include "symmetry_basis.h"
#include "symmetry_checkpoint.h"
#include "symmetry_distribution.h"
#include "symmetry_fulldiag.h"
#include "symmetry_hamiltonian_io.h"
#include "symmetry_sector.h"
#include "wrapperMPI.h"
#ifdef _SCALAPACK
#include "matrixscalapack.h"
#endif

struct StateBuffers {
  double complex *owned;    /* 1-based rank-local block of the state */
  double complex *root;     /* rank 0: D gathered components, others: one */
  int *counts, *displacements;
};

static int same_identity(const struct SymmetryHamiltonianIdentity *a,
                         const struct SymmetryHamiltonianIdentity *b)
{
  return a->model == b->model && a->nsite == b->nsite && a->nup == b->nup &&
         a->ndown == b->ndown && a->ne == b->ne &&
         a->local_site_mask == b->local_site_mask && a->fixed_flags == b->fixed_flags &&
         a->phase == b->phase && a->full_dim == b->full_dim &&
         a->sector_dim == b->sector_dim && a->group_digest == b->group_digest &&
         a->sector_count == b->sector_count && a->sector_xor == b->sector_xor &&
         a->sector_sum == b->sector_sum && a->basis_order == b->basis_order;
}

static int validate_context(const struct BindStruct *X,
                            const struct SymmetryHamiltonianInfo *info)
{
  const struct SymmetryBasisRuntime *sym;
  unsigned long offset, count;
  if (X == NULL || info == NULL || !X->Def.iFlgSymmetryBasis ||
      X->Def.iCalcType != FullDiag) return -1;
  sym = X->Sym;
  if (sym == NULL || sym->enabled != TRUE ||
      sym->basis_layout != SYMMETRY_BASIS_REPLICATED || sym->dim == 0 ||
      sym->dim > INT_MAX || X->Check.idim_max != sym->dim) return -1;
  if (SymmetryBlockRange(sym->dim, myrank, nproc, &offset, &count) != 0 ||
      offset != sym->local_offset || count != sym->local_dim) return -1;
  if (sym->local_dim >= SIZE_MAX / sizeof(double complex) - 1) return -1;
  return info->identity.sector_dim == sym->dim ? 0 : -1;
}

/* LAPACK keeps L_vec[state][component] (both 0-based) on every rank; the
 * distributed solvers gather the state on rank 0 and scatter the blocks. */
static int copy_state(const struct BindStruct *X, unsigned long state,
                      struct StateBuffers *buffers)
{
  const struct SymmetryBasisRuntime *sym = X->Sym;
  if (sym->dim == 1) {
    /* The normalized eigenvector of a 1x1 Hermitian matrix is 1 up to a
     * phase; LAPACK zheev returns exactly 1. ScaLAPACK pzheev has no N=1
     * quick return and was observed to return 1+1i, so the backend vector
     * is not used for a one-state sector. */
    if (sym->local_dim == 1) buffers->owned[1] = 1.0;
    return 0;
  }
#ifdef _SCALAPACK
  if (use_scalapack) {
    if (myrank == 0)
      for (unsigned long i = 0; i < sym->dim; ++i) buffers->root[i] = NAN;
    /* The existing API always returns 0. A missed transfer leaves NaN on
     * rank 0, which the finite/norm check of the scattered state rejects. */
    (void)GetEigenVectorBlock((long)state, (long)sym->dim, Z_vec, descZ_vec,
                              buffers->root);
    return MPI_Scatterv(buffers->root, buffers->counts, buffers->displacements,
                        MPI_DOUBLE_COMPLEX, buffers->owned + 1, (int)sym->local_dim,
                        MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD) == MPI_SUCCESS ? 0 : -1;
  }
#endif
  if (L_vec == NULL || L_vec[state] == NULL) return -1;
  memcpy(buffers->owned + 1, L_vec[state] + sym->local_offset,
         sym->local_dim * sizeof(*buffers->owned));
  return 0;
}

/* The checkpoint contract: finite components and |norm^2 - 1| <= 1e-8.
 * The solver result is not renormalized. */
static int check_state(const struct BindStruct *X, unsigned long state,
                       const double complex *owned)
{
  double norm = 0.0;
  int nonfinite = 0;
  for (unsigned long i = 1; i <= X->Sym->local_dim; ++i) {
    double re = creal(owned[i]), im = cimag(owned[i]);
    if (!isfinite(re) || !isfinite(im)) { nonfinite = 1; break; }
    norm += re*re + im*im;
  }
  if (SumMPI_i(nonfinite) != 0) {
    fprintf(stdoutMPI, "Error: sector FullDiag eigenstate %lu has a non-finite component.\n", state);
    return -1;
  }
  norm = SumMPI_d(norm);
  if (!isfinite(norm) || fabs(norm - 1.0) > 1e-8) {
    fprintf(stdoutMPI, "Error: sector FullDiag eigenstate %lu is not normalized (norm^2=%.17g).\n",
            state, norm);
    return -1;
  }
  return 0;
}

static int allocate_buffers(const struct BindStruct *X, struct StateBuffers *buffers)
{
  int invalid;
  buffers->owned = calloc(X->Sym->local_dim + 1, sizeof(*buffers->owned));
  invalid = buffers->owned == NULL;
#ifdef _SCALAPACK
  if (use_scalapack) {
    buffers->root = calloc(myrank == 0 ? X->Sym->dim : 1, sizeof(*buffers->root));
    buffers->counts = malloc((size_t)nproc * sizeof(*buffers->counts));
    buffers->displacements = malloc((size_t)nproc * sizeof(*buffers->displacements));
    invalid |= buffers->root == NULL || buffers->counts == NULL ||
               buffers->displacements == NULL;
    for (int rank = 0; !invalid && rank < nproc; ++rank) {
      unsigned long offset, count;
      if (SymmetryBlockRange(X->Sym->dim, rank, nproc, &offset, &count) != 0) {
        invalid = 1;
        break;
      }
      buffers->counts[rank] = (int)count;
      buffers->displacements[rank] = (int)offset;
    }
  }
#endif
  return SumMPI_i(invalid) != 0 ? -1 : 0;
}

static void free_buffers(struct StateBuffers *buffers)
{
  free(buffers->owned);
  free(buffers->root);
  free(buffers->counts);
  free(buffers->displacements);
}

static int write_state(const struct BindStruct *local, unsigned long state,
                       const double complex *owned,
                       const struct SymmetryHamiltonianInfo *info)
{
  char name[D_FileNameMax];
  struct SymmetryCheckpointInfo checkpoint = {
    .source_method = FullDiag, .state_index = state, .step = 0, .time = 0,
    .hamiltonian_digest = info->matrix_digest,
    .hamiltonian_digest_kind = SYMMETRY_HAMILTONIAN_MATRIX_V1
  };
  int length = snprintf(name, sizeof(name), "%s_eigenvec_%lu_rank_%d.dat",
                        local->Def.CDataFileHead, state, myrank);
  if (SumMPI_i(length < 0 ||
               (size_t)length + strlen(cParentOutputFolder) >= sizeof(name)) != 0) {
    fprintf(stdoutMPI, "Error: sector FullDiag eigenstate file name is too long.\n");
    return -1;
  }
  return WriteSymmetryCheckpoint(local, name, owned, &checkpoint);
}

int ProcessSymmetryFullDiagStates(const struct BindStruct *X,
    const struct SymmetryHamiltonianInfo *info)
{
  struct SymmetryHamiltonianIdentity current;
  struct StateBuffers buffers = { NULL, NULL, NULL, NULL };
  struct BindStruct local;
  unsigned int outputs;
  int status = -1;
  if (SumMPI_i(validate_context(X, info) != 0) != 0) {
    fprintf(stdoutMPI, "Error: sector FullDiag state output requires the replicated FullDiag context of the prepared matrix.\n");
    return -1;
  }
  if (BuildSymmetryHamiltonianIdentity(X, &current) != 0) return -1;
  if (SumMPI_i(!same_identity(&current, &info->identity)) != 0) {
    fprintf(stdoutMPI, "Error: sector FullDiag basis identity changed after matrix preparation.\n");
    return -1;
  }
  outputs = SymmetryFullDiagOutputMask(&X->Def);
  local = *X;
  local.Check.idim_max = X->Sym->local_dim;
  local.Check.idim_maxMPI = X->Sym->dim;
  StartTimer(5700);
  if (allocate_buffers(X, &buffers) != 0) {
    fprintf(stdoutMPI, "Error: sector FullDiag state buffer allocation failed.\n");
    goto cleanup;
  }
  /* Aggregate Green files are truncated once; states then append in order. */
  if ((outputs & SYM_FD_CORRELATIONS) && GreenOutputInitializeAggregateFiles(&local) != 0)
    goto cleanup;
  for (unsigned long state = 0; state < X->Sym->dim; ++state) {
    int failed;
    local.Phys.eigen_num = (int)state;
    StartTimer(5710);
    failed = copy_state(X, state, &buffers) != 0;
    StopTimer(5710);
    if (SumMPI_i(failed) != 0) {
      fprintf(stdoutMPI, "Error: sector FullDiag eigenstate %lu could not be distributed.\n", state);
      goto cleanup;
    }
    if (check_state(X, state, buffers.owned) != 0) goto cleanup;
    if (outputs & SYM_FD_CORRELATIONS) {
      /* Each wrapper returns the same result on every rank, so the
       * short-circuit keeps the collective sequence aligned. */
      StartTimer(5720);
      failed = expec_cisajs(&local, buffers.owned) != 0 ||
               expec_cisajscktaltdc(&local, buffers.owned) != 0 ||
               expec_nbodyg(&local, buffers.owned) != 0;
      StopTimer(5720);
      if (failed) goto cleanup;
    }
    /* A state is saved only after its correlations succeeded. */
    if (outputs & SYM_FD_EIGENVECTORS) {
      StartTimer(5730);
      failed = write_state(&local, state, buffers.owned, info) != 0;
      StopTimer(5730);
      if (failed) goto cleanup;
    }
  }
  status = 0;
cleanup:
  StopTimer(5700);
  free_buffers(&buffers);
  return status;
}
