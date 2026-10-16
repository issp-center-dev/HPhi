/* Exercise real TPQ initialization/normalization, with an injectable scalar
 * Hamiltonian. Only rank 0 owns vector elements; other ranks still reduce. */
#include <math.h>
#include "Common.h"
#include "FirstMultiply.h"
#include "MakeIniVec.h"
#include "Multiply.h"
#include "wrapperMPI.h"
#ifdef MPI
#include <mpi.h>
#endif

double complex *v0, *v1, *v2;
double LargeValue = 1.0, global_norm, global_1st_norm;
int myrank = 0, nproc = 1, nthreads = 1, step_i = 0;
FILE *stdoutMPI;
const char *cFileNameTimeKeep = "unused", *cTPQStep = "unused", *cTPQStepEnd = "unused";
static int energy_calls, energy_failure, matvec_failure;
static double scalar_h;

void StartTimer(int timer) { (void)timer; }
void StopTimer(int timer) { (void)timer; }
int TimeKeeperWithRandAndStep(struct BindStruct *x, const char *file,
    const char *message, const char *mode, const int step, const int sample)
{
  (void)x; (void)file; (void)message; (void)mode; (void)step; (void)sample;
  return 0;
}
void exitMPI(int code)
{
#ifdef MPI
  MPI_Abort(MPI_COMM_WORLD, code);
#endif
  exit(code);
}
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
double complex SumMPI_dc(double complex value)
{
#ifdef MPI
  double complex sum;
  MPI_Allreduce(&value, &sum, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
  return sum;
#else
  return value;
#endif
}
int mltply(struct BindStruct *x, double complex *out, double complex *in)
{
  unsigned long i;
  if (matvec_failure && myrank == nproc-1) return -1;
  for (i = 1; i <= x->Check.idim_max; ++i) out[i] = scalar_h * in[i];
  return 0;
}
int expec_energy_flct(struct BindStruct *x)
{
  ++energy_calls;
  if (energy_failure) return -1;
  return mltply(x, v0, v1);
}
static void require(int condition, const char *label)
{
  if (SumMPI_i(!condition)) {
    if (!myrank) fprintf(stderr, "TPQ failure contract: %s\n", label);
    exitMPI(1);
  }
}
static double norm(const struct BindStruct *x, const double complex *v)
{
  unsigned long i;
  double value = 0;
  for (i = 1; i <= x->Check.idim_max; ++i) value += creal(conj(v[i])*v[i]);
  return creal(SumMPI_dc(value));
}
int main(int argc, char **argv)
{
  struct BindStruct x = {0};
  double complex a[4] = {123,0,0,0}, b[4] = {456,0,0,0}, c[4] = {0};
#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#else
  (void)argc; (void)argv;
#endif
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  stdoutMPI = stdout;
  /* Test-only oracle input: expose the existing RNG without adding a product
   * vector-output option. Python applies an independently constructed H_q. */
  if (argc == 6 && strcmp(argv[1], "--dump-initial") == 0) {
    unsigned long dim = strtoul(argv[2], NULL, 10), j;
    int sample, samples = atoi(argv[5]);
    x.Check.idim_max = dim / nproc + ((unsigned long)myrank < dim % nproc);
    x.Def.iInitialVecType = atoi(argv[3]);
    x.Def.initial_iv = strtol(argv[4], NULL, 10);
    v0 = calloc(x.Check.idim_max + 1, sizeof(*v0));
    v1 = calloc(x.Check.idim_max + 1, sizeof(*v1));
    require(v0 != NULL && v1 != NULL, "oracle vector allocation");
    for (sample = 0; sample < samples; ++sample) {
      char name[128];
      FILE *fp;
      require(MakeIniVec(sample, &x) == 0, "oracle initialization");
      snprintf(name, sizeof(name), "initial_%d_rank%d.dat", sample, myrank);
      fp = fopen(name, "w");
      require(fp != NULL, "oracle file open");
      fprintf(fp, "%.17g\n", global_1st_norm);
      for (j = 1; j <= x.Check.idim_max; ++j)
        fprintf(fp, "%.17g %.17g\n", creal(v1[j]), cimag(v1[j]));
      require(fclose(fp) == 0, "oracle file close");
    }
    free(v0); free(v1);
#ifdef MPI
    MPI_Finalize();
#endif
    return 0;
  }
  v0 = a; v1 = b; v2 = c;
  x.Def.NsiteMPI = 2;
  x.Def.initial_iv = 7;
  x.Check.idim_max = myrank == 0 ? 3 : 0;
  require(MakeIniVec(0, &x) == 0, "zero-row ranks participate in initialization");
  require(fabs(norm(&x, v0)-1) < 1e-12 && fabs(norm(&x, v1)-1) < 1e-12,
          "initial global norm is one");
  require(a[0] == 123 && b[0] == 456, "one-based sentinels remain untouched");
  require(FirstMultiply(0, NULL) == -1 && energy_calls == 0,
          "invalid initialization propagates before energy evaluation");
  if (myrank == nproc-1) v0 = NULL;
  require(FirstMultiply(0, &x) == -1 && energy_calls == 0,
          "one rank's storage failure reaches every rank");
  v0 = a;
  x.Check.idim_max = 0;
  require(FirstMultiply(0, &x) == -1 && energy_calls == 0,
          "zero global norm rejects instead of dividing by zero");
  x.Check.idim_max = myrank == 0 ? 3 : 0;
  energy_failure = 1;
  require(FirstMultiply(0, &x) == -1, "energy failure propagates");
  energy_failure = 0;
  scalar_h = 2;
  require(FirstMultiply(0, &x) == -1, "annihilated first step rejects");
  scalar_h = 0;
  require(FirstMultiply(0, &x) == 0, "normal first step succeeds");
  require(fabs(norm(&x, v0)-1) < 1e-12, "first step normalizes globally");
  mltply(&x, v0, v1);
  require(Multiply(&x) == 0, "normal later step succeeds");
  scalar_h = 2;
  mltply(&x, v0, v1);
  require(Multiply(&x) == -1, "annihilated later step rejects");
  if (!myrank) v0[1] = NAN;
  require(Multiply(&x) == -1, "nonfinite global norm reaches all ranks");
  require(Multiply(NULL) == -1, "invalid later-step storage rejects");
  require(MakeIniVec(0, &x) == 0, "reset cTPQ input");
  scalar_h = 2;
  x.Def.Param.ExpandCoef = 1;
  mltply(&x, v0, v1);
  require(MultiplyForCanonicalTPQ(&x, 1) == -1, "annihilated cTPQ step rejects");
  require(MultiplyForCanonicalTPQ(&x, NAN) == -1, "nonfinite cTPQ step rejects");
  require(MultiplyForCanonicalTPQ(&x, -1) == -1, "negative cTPQ step rejects");
  require(MultiplyForCanonicalTPQ(NULL, .1) == -1, "invalid cTPQ storage rejects");
  x.Def.Param.ExpandCoef = 0;
  require(MultiplyForCanonicalTPQ(&x, .1) == -1, "invalid cTPQ order rejects");
  x.Def.Param.ExpandCoef = 4;
  if (myrank == nproc-1) v2 = NULL;
  require(MultiplyForCanonicalTPQ(&x, .1) == -1, "one-rank cTPQ workspace failure");
  v2 = c;
  require(MakeIniVec(0, &x) == 0, "reset cTPQ before matvec failure");
  mltply(&x, v0, v1);
  matvec_failure = 1;
  require(MultiplyForCanonicalTPQ(&x, .1) == -1, "one-rank Taylor matvec failure propagates");
  matvec_failure = 0;
  require(MakeIniVec(0, &x) == 0, "reset cTPQ before normal step");
  mltply(&x, v0, v1);
  require(MultiplyForCanonicalTPQ(&x, .1) == 0, "normal cTPQ with empty ranks");
  require(fabs(norm(&x, v0)-1) < 1e-12, "cTPQ global normalization");
  require(fabs(global_norm - (1-.1+.01/2-.001/6+.0001/24)) < 1e-12,
          "cTPQ scalar Taylor polynomial norm");
  if (!myrank) v0[1] = NAN;
  require(MultiplyForCanonicalTPQ(&x, .1) == -1, "nonfinite cTPQ norm reaches all ranks");
  require(MultiplyForTEM(NULL) == -1, "invalid TE storage");
  x.Def.Param.TimeSlice = NAN;
  require(MultiplyForTEM(&x) == -1, "nonfinite TE time step");
  x.Def.Param.TimeSlice = -.1;
  require(MultiplyForTEM(&x) == -1, "negative TE time step");
  x.Def.Param.TimeSlice = .1;
  x.Def.Param.ExpandCoef = 0;
  require(MultiplyForTEM(&x) == -1, "missing TE Taylor order");
  x.Def.Param.ExpandCoef = 4;
  require(MakeIniVec(0, &x) == 0, "reset TE input");
  mltply(&x, v0, v1);
  matvec_failure = 1;
  require(MultiplyForTEM(&x) == -1, "one-rank TE matvec failure");
  matvec_failure = 0;
  require(MakeIniVec(0, &x) == 0, "reset TE before normal step");
  mltply(&x, v0, v1);
  require(MultiplyForTEM(&x) == 0, "normal TE step with empty ranks");
  require(fabs(norm(&x, v0)-1) < 1e-12, "TE normalization");
  require(fabs(global_norm-cabs(1-.2*I-.04/2+.008*I/6+.0016/24)) < 1e-12,
          "TE scalar Taylor polynomial norm");
  if (!myrank) v0[1] = NAN;
  require(MultiplyForTEM(&x) == -1, "nonfinite TE norm reaches all ranks");
  for (unsigned long i = 1; i <= x.Check.idim_max; ++i) v0[i] = v1[i] = 0;
  require(MultiplyForTEM(&x) == -1, "zero TE norm reaches all ranks");
  if (!myrank) puts("TPQ failure propagation and zero-row normalization PASS");
#ifdef MPI
  MPI_Finalize();
#endif
  return 0;
}
