/* Exercise the real FullDiag/lapack_diag callers with an injected solver
 * failure, including a failure reported on only one MPI rank. */
#include "CalcByFullDiag.h"
#include "matrixlapack.h"
#include "wrapperMPI.h"
#include "FileIO.h"
#ifdef MPI
#include <mpi.h>
#endif

static int build_failed, input_failed, solver_failed;
static int solver_calls, output_opens, phys_calls, output_calls;
void StartTimer(int timer) { (void)timer; }
void StopTimer(int timer) { (void)timer; }
int SumMPI_i(int value)
{
#ifdef MPI
  int sum;
  MPI_Allreduce(&value, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return sum;
#else
  return value;
#endif
}
int makeHam(struct BindStruct *x) { (void)x; return build_failed ? -1 : 0; }
int makeHamSym(const struct BindStruct *x) { (void)x; return build_failed ? -1 : 0; }
int inputHam(struct BindStruct *x) { (void)x; return input_failed ? -1 : 0; }
int outputHam(struct BindStruct *x) { (void)x; return 0; }
void phys(struct BindStruct *x, unsigned long n) { (void)x; (void)n; ++phys_calls; }
int output(struct BindStruct *x) { (void)x; ++output_calls; return 0; }
int ZHEEVall(int n, double complex **a, double complex *r, double complex **v)
{
  int i;
  (void)a; (void)v;
  ++solver_calls;
  if (solver_failed) return 0;
  for (i = 0; i < n; ++i) r[i] = i + 0.25;
  return 1;
}
int childfopenMPI(const char *name, const char *mode, FILE **fp)
{
  (void)name; (void)mode;
  ++output_opens;
  *fp = tmpfile();
  return *fp == NULL ? -1 : 0;
}
static void require(int condition, const char *label)
{
  if (SumMPI_i(!condition)) {
    fprintf(stderr, "FullDiag failure contract: %s\n", label);
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    exit(1);
  }
}
int main(int argc, char **argv)
{
  struct EDMainCalStruct x = {0};
  double complex row0[3] = {0}, row1[3] = {0, 1, 0}, row2[3] = {0, 0, 2};
  double complex *matrix[3] = {row0, row1, row2}, eigenvalues[3] = {0};
#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#else
  (void)argc; (void)argv;
#endif
  stdoutMPI = stdout;
  Ham = matrix; L_vec = matrix; v0 = eigenvalues;
  x.Bind.Def.iCalcType = FullDiag;
  x.Bind.Def.iSolver = SOLVER_LAPACK;
  x.Bind.Check.idim_max = 2;
  build_failed = myrank == nproc-1;
  require(CalcByFullDiag(&x) == FALSE && solver_calls == 0 && output_opens == 0,
          "matrix-build failure prevents diagonalization on all ranks");
  build_failed = 0; input_failed = myrank == nproc-1;
  x.Bind.Def.iInputHam = TRUE;
  require(CalcByFullDiag(&x) == FALSE && solver_calls == 0 && output_opens == 0,
          "matrix-input failure prevents diagonalization on all ranks");
  x.Bind.Def.iInputHam = FALSE; input_failed = 0;
  solver_failed = myrank == nproc-1;
  require(CalcByFullDiag(&x) == FALSE && solver_calls == 1 && output_opens == 0 &&
          phys_calls == 0 && output_calls == 0,
          "ZHEEVall returns zero: no output or observables on any rank");
  solver_failed = 0;
  require(CalcByFullDiag(&x) == TRUE && solver_calls == 2 && output_opens == 1 &&
          phys_calls == 1 && output_calls == 1,
          "ZHEEVall returns one: normal output contract");
  if (!myrank) puts("FullDiag failure propagation PASS");
#ifdef MPI
  MPI_Finalize();
#endif
  return 0;
}
