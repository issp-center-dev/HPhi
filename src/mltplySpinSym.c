#include "mltplySpinSym.h"
#ifdef MPI
#include <mpi.h>
#endif
#include "CalcTime.h"
#include "DefCommon.h"
#include "global.h"
#include "struct.h"
#include "symmetry_basis.h"
#include "symmetry_matvec_plan.h"

struct LegacyApplyContext {
  const struct BindStruct *X;
  double complex *tmp_v0;
  const double complex *full_v1;
  double complex input_amp;
  double complex prdct;
};

static int apply_legacy_entry(unsigned long int out_index,
                              double complex coefficient,
                              void *context)
{
  unsigned long int local_index;
  double complex contribution;
  struct LegacyApplyContext *apply = (struct LegacyApplyContext *)context;
  if (SymmetryBasisGlobalToLocal(apply->X->Sym, out_index, &local_index) != TRUE) {
    return 0;
  }
  contribution = coefficient * apply->input_amp;
  apply->tmp_v0[local_index] += contribution;
  apply->prdct += conj(apply->full_v1[out_index]) * contribution;
  return 0;
}

static int apply_legacy_scan(struct BindStruct *X,
                             double complex *tmp_v0,
                             const double complex *full_v1,
                             double complex *prdct)
{
  unsigned long int beta;
  struct LegacyApplyContext context;
  context.X = X;
  context.tmp_v0 = tmp_v0;
  context.full_v1 = full_v1;
  context.prdct = 0.0;
  for (beta = 1UL; beta <= X->Sym->dim; beta++) {
    context.input_amp = full_v1[beta];
    if (cabs(context.input_amp) == 0.0) continue;
    if (SymmetryEnumerateColumn(X, beta, apply_legacy_entry, &context) != 0) {
      return -1;
    }
  }
  *prdct = context.prdct;
  return 0;
}

static const double complex *get_full_input_vector(struct BindStruct *X,
                                                   double complex *tmp_v1)
{
#ifdef MPI
  if (nproc > 1) {
    int ierr;
    const double complex *sendbuf = tmp_v1;
    if (X->Sym->mpi_full_v1 == NULL || X->Sym->mpi_recvcounts == NULL ||
        X->Sym->mpi_displs == NULL) {
      return NULL;
    }
    if (X->Sym->local_dim > 0UL) sendbuf = &tmp_v1[1];
    ierr = MPI_Allgatherv(sendbuf, (int)X->Sym->local_dim,
                          MPI_DOUBLE_COMPLEX,
                          &X->Sym->mpi_full_v1[1], X->Sym->mpi_recvcounts,
                          X->Sym->mpi_displs, MPI_DOUBLE_COMPLEX,
                          MPI_COMM_WORLD);
    if (ierr != MPI_SUCCESS) return NULL;
    if (X->Sym->matvec_plan != NULL) {
      X->Sym->matvec_plan->input_allgather_calls++;
    }
    return X->Sym->mpi_full_v1;
  }
#else
  (void)X;
#endif
  return tmp_v1;
}

int mltplySpinSym(struct BindStruct *X,
                  double complex *tmp_v0,
                  double complex *tmp_v1)
{
  double complex prdct = 0.0;
  const double complex *full_v1;

  StartTimer(1501);
  full_v1 = get_full_input_vector(X, tmp_v1);
  StopTimer(1501);
  if (full_v1 == NULL) return -1;

  if (X->Sym->matvec_mode == SYMMETRY_MATVEC_MODE_PLAN &&
      X->Sym->matvec_plan != NULL &&
      X->Sym->matvec_plan->halo.reference_enabled == TRUE) {
    if (ExchangeSymmetryVectorHaloReference(
            &X->Sym->matvec_plan->halo, tmp_v1, full_v1) != 0) {
      fprintf(stdoutMPI,
              "Error: symmetry halo reference exchange did not match "
              "the gathered input vector.\n");
      return -1;
    }
  }

  if (X->Sym->matvec_mode == SYMMETRY_MATVEC_MODE_LEGACY) {
    StartTimer(1502);
    if (apply_legacy_scan(X, tmp_v0, full_v1, &prdct) != 0) {
      StopTimer(1502);
      return -1;
    }
    StopTimer(1502);
  } else if (X->Sym->matvec_mode == SYMMETRY_MATVEC_MODE_PLAN) {
    StartTimer(1503);
    if (ApplySymmetryMatvecPlan(X, tmp_v0, full_v1, &prdct) != 0) {
      StopTimer(1503);
      return -1;
    }
    StopTimer(1503);
  } else {
    fprintf(stdoutMPI, "Error: invalid symmetry matvec mode.\n");
    return -1;
  }

  if (X->Sym->matvec_plan != NULL) {
    X->Sym->matvec_plan->matvec_calls++;
  }
  X->Large.prdct += prdct;
  return 0;
}
