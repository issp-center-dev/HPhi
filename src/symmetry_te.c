#include <limits.h>
#include <math.h>
#include "symmetry_te.h"
#include "symmetry_basis.h"
#include "symmetry_diagonal.h"
#include "symmetry_basis_io.h"
#include "symmetry_matvec_plan.h"
#include "symmetry_sector.h"
#include "HPhiTrans.h"
#include "wrapperMPI.h"

int SymmetryTEIsDynamic(const struct DefineList *def)
{
  return def->iFlgSymmetryBasis && def->iCalcType == TimeEvolution &&
         (def->NLaser || def->NTETransferMax || def->NTEInterAllMax);
}

double SymmetryTETime(const struct DefineList *def, unsigned int step)
{
  return def->NLaser ? def->Param.Tinit + step * def->Param.TimeSlice : def->TETime[step];
}

static void free_arrays(struct SymmetryTEHamiltonian *view)
{
  free(view->transfer); free(view->interall);
  free(view->transfer_data); free(view->interall_data);
  free(view->transfer_values); free(view->interall_values);
  view->transfer = view->interall = NULL;
  view->transfer_data = view->interall_data = NULL;
  view->transfer_values = view->interall_values = NULL;
}

void FreeSymmetryTEHamiltonian(struct SymmetryTEHamiltonian *view)
{
  free_arrays(view);
  free(view->digests);
  view->digests = NULL;
}

int SelectSymmetryTEHamiltonian(struct BindStruct *X, struct SymmetryTEHamiltonian *view,
                                unsigned int step, double time)
{
  const struct DefineList *base = &view->base;
  X->Def = *base;
  free_arrays(view);
  unsigned int nt = base->NLaser ? 0 : base->NTETransfer[step];
  unsigned int nd = base->NLaser ? 0 : base->NTETransferDiagonal[step];
  unsigned int ni = base->NLaser ? 0 : base->NTEInterAllOffDiagonal[step];
  unsigned int nid = base->NLaser ? 0 : base->NTEInterAllDiagonal[step];
  unsigned int nc = base->NLaser ? 0 : base->NTEChemi[step];
  uint64_t transfers = (uint64_t)base->EDNTransfer + nt + nd + nc;
  uint64_t interactions = (uint64_t)base->NInterAll_OffDiagonal + ni + nid;
  int failed = transfers > UINT_MAX || interactions > UINT_MAX ||
               transfers > SIZE_MAX / (4*sizeof(int)) || interactions > SIZE_MAX / (8*sizeof(int));
  if (SumMPI_i(failed) != 0) return -1;
  size_t tn = transfers ? (size_t)transfers : 1;
  size_t in = interactions ? (size_t)interactions : 1;
  view->transfer = calloc(tn, sizeof(*view->transfer));
  view->transfer_data = calloc(tn * 4, sizeof(int));
  view->transfer_values = calloc(tn, sizeof(*view->transfer_values));
  view->interall = calloc(in, sizeof(*view->interall));
  view->interall_data = calloc(in * 8, sizeof(int));
  view->interall_values = calloc(in, sizeof(*view->interall_values));
  failed = !view->transfer || !view->transfer_data || !view->transfer_values ||
           !view->interall || !view->interall_data || !view->interall_values;
  if (SumMPI_i(failed) != 0) return -1;
  for (size_t i = 0; i < tn; ++i) view->transfer[i] = view->transfer_data + 4*i;
  for (size_t i = 0; i < in; ++i) view->interall[i] = view->interall_data + 8*i;
  unsigned int pos = 0;
#define TRANSFER(row, val) do { memcpy(view->transfer[pos], (row), 4*sizeof(int)); \
  view->transfer_values[pos++] = (val); } while (0)
  for (unsigned int i = 0; i < base->EDNTransfer; ++i)
    TRANSFER(base->EDGeneralTransfer[i], base->EDParaGeneralTransfer[i]);
  for (unsigned int i = 0; i < nt; ++i)
    TRANSFER(base->TETransfer[step][i], base->ParaTETransfer[step][i]);
  for (unsigned int i = 0; i < nd; ++i) {
    int a = base->TETransferDiagonal[step][i][0], s = base->TETransferDiagonal[step][i][1];
    int row[] = {a,s,a,s};
    TRANSFER(row, base->ParaTETransferDiagonal[step][i]);
  }
  for (unsigned int i = 0; i < nc; ++i) {
    int a = base->TEChemi[step][i], s = base->SpinTEChemi[step][i];
    int row[] = {a,s,a,s};
    TRANSFER(row, base->ParaTEChemi[step][i]);
  }
#undef TRANSFER
  pos = 0;
#define INTERALL(row, val) do { memcpy(view->interall[pos], (row), 8*sizeof(int)); \
  view->interall_values[pos++] = (val); } while (0)
  for (unsigned int i = 0; i < base->NInterAll_OffDiagonal; ++i)
    INTERALL(base->InterAll_OffDiagonal[i], base->ParaInterAll_OffDiagonal[i]);
  for (unsigned int i = 0; i < ni; ++i)
    INTERALL(base->TEInterAllOffDiagonal[step][i], base->ParaTEInterAllOffDiagonal[step][i]);
  for (unsigned int i = 0; i < nid; ++i) {
    const int *q = base->TEInterAllDiagonal[step][i];
    int row[] = {q[0],q[1],q[0],q[1],q[2],q[3],q[2],q[3]};
    INTERALL(row, base->ParaTEInterAllDiagonal[step][i]);
  }
#undef INTERALL
  X->Def.EDNTransfer = (unsigned int)transfers;
  X->Def.EDGeneralTransfer = view->transfer;
  X->Def.EDParaGeneralTransfer = view->transfer_values;
  X->Def.NInterAll_OffDiagonal = (unsigned int)interactions;
  X->Def.InterAll_OffDiagonal = view->interall;
  X->Def.ParaInterAll_OffDiagonal = view->interall_values;
  X->Def.istep = step;
  if (base->NLaser) {
    /* Use the parsed static coefficients, in exactly the ED index ordering.
     * The raw Peierls driver instead indexes the original transfer array. */
    X->Def.ParaGeneralTransfer = view->transfer_values;
    if (TransferWithPeierls(X, time) != 0) return -1;
  }
  return 0;
}

int ValidateSymmetryTESchedule(struct BindStruct *X, struct SymmetryTEHamiltonian *view)
{
  const struct DefineList *def = &view->base;
  int failed = (def->NTETransferMax && def->NTEInterAllMax) ||
               (def->NLaser && (def->NTETransferMax || def->NTEInterAllMax));
  if (def->NLaser) {
    failed |= def->NLaser != 9 || !def->ParaLaser;
    if (!failed) {
      for (int i = 0; i < 9; ++i) failed |= !isfinite(def->ParaLaser[i]);
      double mode = def->ParaLaser[0], lx = def->ParaLaser[5], ly = def->ParaLaser[6];
      failed |= mode < 0 || mode > 5 || floor(mode) != mode ||
                lx < 1 || lx >= INT_MAX || floor(lx) != lx ||
                ly < 1 || ly >= INT_MAX || floor(ly) != ly;
      if (mode == 0 || mode == 4 || mode == 5) failed |= def->ParaLaser[3] <= 0;
      if (mode == 5) failed |= def->ParaLaser[4] <= 0;
    }
  }
  if (SumMPI_i(failed) != 0) {
    fprintf(stdoutMPI, "Error: invalid or mixed symmetry TE driving definitions.\n");
    return -1;
  }
  view->digests = calloc(def->Lanczos_max, sizeof(*view->digests));
  if (SumMPI_i(view->digests == NULL) != 0) return -1;
  for (unsigned int step = 0; step < def->Lanczos_max; ++step) {
    double time = SymmetryTETime(def, step);
    failed = !isfinite(time) ||
             (step && (!isfinite(time - SymmetryTETime(def, step-1)) ||
                       time < SymmetryTETime(def, step-1)));
    if (SumMPI_i(failed) != 0 || SelectSymmetryTEHamiltonian(X, view, step, time) != 0)
      return -1;
    failed = ValidateSymmetryHamiltonian(X) != 0 ||
             ComputeSymmetryHamiltonianDigest(&X->Def, &view->digests[step]) != 0;
    X->Def = *def;
    if (SumMPI_i(failed) != 0) {
      fprintf(stdoutMPI, "Error: symmetry TE preflight failed at step %u (time %.17g).\n", step, time);
      return -1;
    }
  }
  free_arrays(view);
  fprintf(stdoutMPI, "Symmetry TE: all %u Hamiltonians passed preflight.\n", def->Lanczos_max);
  return 0;
}

int RebuildSymmetryTEPlan(struct BindStruct *X)
{
  struct SymmetryBasisRuntime *sym = X->Sym;
  int replicated = sym->basis_layout == SYMMETRY_BASIS_REPLICATED;
  unsigned long count = replicated ? sym->dim : sym->local_dim;
  struct SymmetryBasisVector *basis = replicated ? sym->basis : sym->local_basis;
  double *diagonal = calloc(count ? count : 1, sizeof(*diagonal));
  int failed = diagonal == NULL || basis == NULL;
  if (SumMPI_i(failed) != 0) { free(diagonal); return -1; }
  for (unsigned long i = 1; i <= count; ++i)
    if (EvaluateSymmetryStateDiagonal(&X->Def, basis[i].rep_state, &diagonal[i-1]) != 0)
      failed = 1;
  if (SumMPI_i(failed) != 0) { free(diagonal); return -1; }
  /* Representatives, normalizations, phases and MPI ownership never change.
   * Reuse them directly instead of enumerating the raw Hilbert space again. */
  for (unsigned long i = 1; i <= count; ++i) basis[i].diagonal = diagonal[i-1];
  free(diagonal);
  FreeSymmetryMatvecPlan(sym->matvec_plan);
  sym->matvec_plan = NULL;
  if (replicated) return BuildSymmetryMatvecPlan(X);

  /* Solver plans release their lookup directory and global ghost indices.
   * Recreate the directory from the retained rank-local basis before rebuilding
   * the graph. This also admits changes in sparsity, including empty rows. */
  FreeSymmetryRepresentativeDirectory(sym->representative_directory);
  sym->representative_directory = NULL;
  sym->representative_directory_stats_ready = FALSE;
  sym->representative_directory_heavy_storage_released = FALSE;
  if (BuildSymmetryRepresentativeDirectory(
          sym->local_basis, sym->dim, sym->local_dim, sym->local_capacity,
          sym->local_offset, sym->rank_offsets, myrank, nproc,
          &sym->representative_directory) != 0) return -1;
  return BuildSymmetryDistributedMatvecPlan(X);
}
