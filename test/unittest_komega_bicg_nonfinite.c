/* Unit test for BiCG non-finite input detection.
 * The pre-fix path propagated NaN through rho/alpha/resnorm and reported it
 * only indirectly after contaminating the update state. */
#include "komega/komega.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>

#ifdef MPI
#include <mpi.h>
#endif

int main(int argc, char **argv) {
  int ndim = 1;
  int nl = 1;
  int nz = 1;
  int itermax = 2;
  int comm = 0;
  int status[3] = {0, 0, 0};
  int diag_stage = 0;
  int diag_iter = 0;
  int diag_seed = 0;
  double threshold = 1.0e-12;
  double complex diag_values[KOMEGA_BICG_DIAG_VALUE_COUNT];
  double complex x[1] = {0.0};
  double complex z[1] = {1.0 + 0.1 * I};
  double complex v12[1] = {0.0};
  double complex v2[1] = {NAN + 0.0 * I};
  double complex v14[1] = {0.0};
  double complex v4[1] = {1.0};
  double complex r_l[1] = {1.0};

#ifdef MPI
  MPI_Init(&argc, &argv);
  comm = (int)MPI_Comm_c2f(MPI_COMM_WORLD);
#else
  (void)argc;
  (void)argv;
#endif

  komega_bicg_init(&ndim, &nl, &nz, x, z, &itermax, &threshold, &comm);
  komega_bicg_update(v12, v2, v14, v4, x, r_l, status);
  komega_bicg_getdiag(&diag_stage, &diag_iter, &diag_seed, diag_values);
  komega_bicg_finalize();

#ifdef MPI
  MPI_Finalize();
#endif

  if (status[0] < 0 && status[1] == 5 &&
      diag_stage == 2 && diag_iter == 1 && diag_seed == 1 &&
      (!isfinite(creal(diag_values[0])) || !isfinite(cimag(diag_values[0]))) &&
      (!isfinite(creal(diag_values[1])) || !isfinite(cimag(diag_values[1])))) {
    printf("UNIT TEST PASSED\n");
    return 0;
  }

  printf("UNIT TEST FAILED: status=(%d,%d,%d) diag=(stage=%d,iter=%d,seed=%d,"
         "rho_local=(%g,%g),rho_global=(%g,%g))\n",
         status[0], status[1], status[2], diag_stage, diag_iter, diag_seed,
         creal(diag_values[0]), cimag(diag_values[0]),
         creal(diag_values[1]), cimag(diag_values[1]));
  return 1;
}
