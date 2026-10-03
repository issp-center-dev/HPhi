#include "Common.h"
#include "matrixscalapack.h"
#include <math.h>

static int fail_stage;
/* Interpose the backend call, while retaining real MPI/BLACS distribution.
 * A failure on the final rank must stop every rank at the same boundary. */
void pzheev_(char *jobz, char *uplo, const long *n, double complex *a,
    const int *ia, const int *ja, int *da, double *w, double complex *z,
    const int *iz, const int *jz, int *dz, double complex *work,
    const long *lw, double *rwork, const long *lrw, int *info)
{
  (void)jobz; (void)uplo; (void)a; (void)ia; (void)ja; (void)da;
  (void)z; (void)iz; (void)jz; (void)dz; (void)lrw;
  *info = 0;
  if (*lw == -1) {
    *work = 256; *rwork = 256;
    if (fail_stage == 1 && myrank == nproc-1) *info = -1;
    if (fail_stage == 3 && myrank == nproc-1) *work = NAN;
  } else {
    if (fail_stage == 2 && myrank == nproc-1) *info = 1;
    for (long i = 0; i < *n; ++i) w[i] = i + 0.25;
  }
}
static void require(int condition, const char *label)
{
  int local = !condition, failed;
  MPI_Allreduce(&local, &failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (failed) {
    if (!myrank) fprintf(stderr, "ScaLAPACK failure contract: %s\n", label);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}
int main(int argc, char **argv)
{
  double complex storage[4][4] = {{0}}, *a[4], values[4], *z;
  int desc[9] = {0};
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  stdoutMPI = stdout;
  for (int i = 0; i < 4; ++i) { a[i] = storage[i]; storage[i][i] = i + 1; }
  z = calloc(16, sizeof(*z));
  require(z != NULL, "test allocation");
  for (fail_stage = 1; fail_stage <= 3; ++fail_stage) {
    for (int i = 0; i < 4; ++i) values[i] = 987;
    require(diag_scalapack_cmp(4, a, values, z, desc) == -1,
            "query, solve and invalid workspace failures propagate");
    require(use_scalapack == 0 && desc[1] == -1 && values[0] == 987,
            "no eigenvalues published and BLACS grid released");
  }
  fail_stage = 0;
  require(diag_scalapack_cmp(4, a, values, z, desc) == 0 && use_scalapack == 1,
          "success after failed calls");
  FreeDistributedEigenvectors(&z, desc, &use_scalapack);
  require(z == NULL && use_scalapack == 0 && desc[1] == -1, "normal cleanup");
  if (!myrank) puts("ScaLAPACK collective failure propagation PASS");
  MPI_Finalize();
  return 0;
}
