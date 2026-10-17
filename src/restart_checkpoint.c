/* v3 solver checkpoints. File operations are rank-local, decisions collective. */
#include <errno.h>
#include <float.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include "Common.h"
#include "CalcTime.h"
#include "wrapperMPI.h"
#include "restart_checkpoint.h"
#include "restart_manifest.h"
#include "symmetry_checkpoint.h"
#include "symmetry_basis.h"
#include "symmetry_sector.h"
#ifdef HPHI_TEST_SYMMETRY_PROBE
#include "symmetry_fulldiag_probe.h"
#include "restart_probe.h"
#endif
#ifdef MPI
#include <mpi.h>
#endif

static int agree(int invalid, const char *stage)
{
  if (SumMPI_i(invalid != 0) == 0) return 0;
#ifdef HPHI_TEST_SYMMETRY_PROBE
  /* Test-only evidence that every rank crossed the same agreement boundary
   * before the solver's top-level MPI abort can terminate its peers. */
  fprintf(stdout, "restart publication agreed failure stage=%s rank=%d\n", stage, myrank);
  fflush(stdout);
#endif
  if (strcmp(stage, "manifest update") == 0)
    fprintf(stdoutMPI, "Error: restart manifest update failed.\n");
  else
    fprintf(stdoutMPI, "Error: restart checkpoint %s failed.\n", stage);
  return -1;
}
static int basis_kind(const struct BindStruct *X) { return X->Def.iFlgSymmetryBasis ? RESTART_BASIS_SECTOR : RESTART_BASIS_RAW; }
static FILE *open_output(const char *name, const char *mode)
{
  char path[D_FileNameMax];
  int length = snprintf(path, sizeof(path), "%s%s", cParentOutputFolder, name);
  if (length < 0 || length >= (int)sizeof(path)) { errno = ENAMETOOLONG; return NULL; }
  return fopen(path, mode);
}
uint64_t RestartLocalOffset(const struct BindStruct *X)
{
  return X->Def.iFlgSymmetryBasis ? (uint64_t)X->Sym->local_offset : (uint64_t)ExscanMPI_li(X->Check.idim_max);
}
static uint64_t raw_site_digest(const struct DefineList *def)
{
  uint64_t h = RestartHashTag("hphi-raw-sites-v3");
  for (unsigned int i = 0; i < def->NsiteMPI; ++i) {
    h = RestartHashWord(h, (uint64_t)def->SiteToBit[i]);
    h = RestartHashWord(h, (uint64_t)def->LocSpn[i]);
  }
  return h;
}
/* Identity words shared by every file of this run (RANK/OFFSET/LOCAL_DIM/COMPONENT/set words are set per file). */
static int build_identity(const struct BindStruct *X, uint64_t *w)
{
  int invalid = X == NULL || w == NULL || sizeof(double) != 8 ||
                sizeof(double complex) != 16 || DBL_MANT_DIG != 53 || DBL_MAX_EXP != 1024;
  if (agree(invalid, "runtime validation")) return -1;
  memset(w, 0, RESTART_V3_WORDS * sizeof(*w));
  if (X->Def.iFlgSymmetryBasis) {
    if (SymmetryCheckpointIdentityWords(X, w) != 0) return -1;      /* words 0-29 as v1/v2 */
    uint64_t order_xor = 0, order_sum = 0, tag = RestartHashTag("hphi-ckpt-order-v3");
    for (unsigned long i = 1; i <= X->Sym->local_dim; ++i) {
      const struct SymmetryBasisVector *e = SymmetryBasisLocalEntry(X->Sym, i);
      if (e == NULL) { invalid = 1; break; }
      uint64_t h = RestartHashWord(tag, X->Sym->local_offset + i);
      h = RestartHashWord(h, (uint64_t)e->rep_state);
      h = RestartHashWord(h, e->orbit_size);
      h = RestartHashWord(h, e->stabilizer_size);
      order_xor ^= h; order_sum += h;
    }
    if (agree(invalid, "basis order identity")) return -1;
    w[RW_ORDER] = order_xor; w[RW_ORDER_SUM] = order_sum;             /* reduced below */
    w[RW_BASIS_KIND] = RESTART_BASIS_SECTOR;
  } else {
    const struct DefineList *d = &X->Def;
    uint64_t flags = (uint64_t)d->iRestartFixedFlags;                 /* normalized in readdef (Q2) */
    int sz_fixed = (flags & 2) != 0, ne_fixed = (flags & 1) != 0;
    w[RW_MAGIC] = UINT64_C(0x0a31565349485048); w[RW_PHASE] = 0; w[RW_SCALAR] = 128;
    w[RW_MODEL] = (uint64_t)d->iCalcModelInput; w[RW_NSITE] = d->NsiteMPI;
    w[RW_NUP] = sz_fixed && !d->iFlgGeneralSpin ? d->NupMPI : 0;
    w[RW_NDOWN] = sz_fixed && !d->iFlgGeneralSpin ? d->NdownMPI : 0;
    w[RW_NE] = ne_fixed ? d->NeMPI : 0;
    w[RW_DIM] = X->Check.idim_maxMPI; w[RW_LAYOUT] = 2;
    w[RW_FIXED_FLAGS] = flags; w[RW_BASIS_KIND] = RESTART_BASIS_RAW;
    w[RW_TOTAL2SZ] = sz_fixed ? (uint64_t)(int64_t)d->Total2SzMPI : 0;
    w[RW_SITE_DIGEST] = raw_site_digest(d);
    w[RW_HAMILTONIAN] = 0;
  }
  w[RW_VERSION] = 3; w[RW_RANKS] = (uint64_t)nproc;
#ifdef MPI
  if (X->Def.iFlgSymmetryBasis) {
    uint64_t x = w[RW_ORDER], s = w[RW_ORDER_SUM];
    if (MPI_Allreduce(&x, &w[RW_ORDER], 1, MPI_UINT64_T, MPI_BXOR, MPI_COMM_WORLD) != MPI_SUCCESS) exitMPI(-1);
    if (MPI_Allreduce(&s, &w[RW_ORDER_SUM], 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS) exitMPI(-1);
  }
#endif
  return 0;
}

static const char *word_names[RESTART_V3_WORDS] = {"magic", "version", "phase convention", "scalar encoding", "model", "nsite",
  "nup", "ndown", "particle number", "raw dimension", "dimension", "MPI rank count", "MPI rank", "basis layout", "local offset",
  "local dimension", "group/character digest", "sector count", "sector XOR", "sector sum", "basis order", "hamiltonian", "method",
  "state", "step", "time", "data XOR", "data sum", "local site mask", "fixed quantity flags", "basis kind", "component", "total 2Sz",
  "site digest", "basis order sum", "set count", "scalar A", "scalar B", "scalar C", "flags"};
static const int identity_words[] = {RW_MAGIC, RW_VERSION, RW_PHASE, RW_SCALAR, RW_MODEL, RW_NSITE, RW_NUP, RW_NDOWN, RW_NE,
  RW_FULL_DIM, RW_DIM, RW_GROUP, RW_SECTOR_COUNT, RW_SECTOR_XOR, RW_SECTOR_SUM, RW_ORDER, RW_LOCAL_SITE_MASK,
  RW_FIXED_FLAGS, RW_BASIS_KIND, RW_TOTAL2SZ, RW_SITE_DIGEST, RW_ORDER_SUM};
static int check_identity(const uint64_t *expected, const uint64_t *header, const char **bad)
{
  for (size_t i = 0; i < sizeof(identity_words)/sizeof(identity_words[0]); ++i) {
    int w = identity_words[i];
    if (header[w] != expected[w]) { *bad = word_names[w]; return -1; }
  }
  int layout_ok = header[RW_BASIS_KIND] == RESTART_BASIS_RAW ? header[RW_LAYOUT] == 2 : header[RW_LAYOUT] <= 1;
  if (!layout_ok) { *bad = word_names[RW_LAYOUT]; return -1; }
  if (header[RW_FLAGS] & ~UINT64_C(1)) { *bad = word_names[RW_FLAGS]; return -1; }
  return 0;
}
static uint64_t step_lower_bound(int method)
{
  return method == TPQCalc ? 2 : 1;   /* mTPQ bootstrap emits steps 0/1 (N1); cTPQ/CG/TE start at 1 (Q6) */
}

static void allreduce_xor_sum(uint64_t x, uint64_t s, uint64_t *gx, uint64_t *gs)
{
#ifdef MPI
  if (MPI_Allreduce(&x, gx, 1, MPI_UINT64_T, MPI_BXOR, MPI_COMM_WORLD) != MPI_SUCCESS) exitMPI(-1);
  if (MPI_Allreduce(&s, gs, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS) exitMPI(-1);
#else
  *gx = x; *gs = s;
#endif
}
static void bcast_words(uint64_t *words, uint64_t count)
{
#ifdef MPI
  if (count > INT_MAX || MPI_Bcast(words, (int)count, MPI_UINT64_T, 0, MPI_COMM_WORLD) != MPI_SUCCESS) exitMPI(-1);
#else
  (void)words; (void)count;
#endif
}
static void remove_generation(const struct BindStruct *X, const struct RestartUnit *unit,
                              uint64_t gen, int part)
{
  for (int i = 0; i < unit->nset; ++i)
    for (int c = 0; c < unit->ncomp; ++c) {
      char name[D_FileNameMax], path[D_FileNameMax];
      int length = snprintf(name, sizeof(name), cFileNameRestartVector, X->Def.CDataFileHead,
                            (int)unit->sets[i], c, (int)gen, myrank);
      if (length < 0 || length >= (int)sizeof(name)) {
        fprintf(stderr, "Warning: restart generation filename is too long.\n");
        continue;
      }
      int len = snprintf(path, sizeof(path), "%s%s%s", cParentOutputFolder, name, part ? ".part" : "");
      if (len < 0 || len >= (int)sizeof(path) || (remove(path) != 0 && errno != ENOENT))
        fprintf(stderr, "Warning: cannot remove restart generation file %s.\n", name);
    }
}

static void remove_superseded_files(const struct BindStruct *X, const struct RestartManifest *old,
                                    const struct RestartManifest *current)
{
  for (size_t i = 0; i < old->nrows; ++i) {
    const struct RestartManifestRow *a = &old->rows[i], *b = RestartManifestFind(current, a->set);
    for (uint64_t r = (uint64_t)myrank; r < a->ranks; r += (uint64_t)nproc)
      for (uint64_t c = 0; c < a->ncomp; ++c) {
        if (b != NULL && a->gen == b->gen && r < b->ranks && c < b->ncomp) continue;
        char path[D_FileNameMax], name[D_FileNameMax];
        int length = snprintf(name, sizeof(name), cFileNameRestartVector, X->Def.CDataFileHead,
                              (int)a->set, (int)c, (int)a->gen, (int)r);
        if (length < 0 || length >= (int)sizeof(name)) continue;
        length = snprintf(path, sizeof(path), "%s%s", cParentOutputFolder, name);
        if (length < 0 || length >= (int)sizeof(path)) continue;
        if (remove(path) != 0 && errno != ENOENT) fprintf(stderr, "Warning: cannot remove superseded restart file %s.\n", path);
      }
  }
}

static int write_part(struct BindStruct *X, const uint64_t *identity, const struct RestartVectorInfo *info,
                      int comp, uint64_t gen, uint64_t flags, const double complex *vector, char *final_name)
{
  uint64_t h[RESTART_V3_WORDS]; char part[D_FileNameMax];
  memcpy(h, identity, sizeof(h));
  h[RW_RANK] = (uint64_t)myrank; h[RW_OFFSET] = RestartLocalOffset(X); h[RW_LOCAL_DIM] = X->Check.idim_max;
  h[RW_METHOD] = info->method; h[RW_STATE] = info->set; h[RW_STEP] = info->step; h[RW_TIME] = RestartDoubleWord(info->time);
  h[RW_NSET] = info->nset; h[RW_COMPONENT] = (uint64_t)comp; h[RW_FLAGS] = flags;
  h[RW_SCALAR_A] = RestartDoubleWord(info->scalar_a); h[RW_SCALAR_B] = RestartDoubleWord(info->scalar_b);
  h[RW_SCALAR_C] = RestartDoubleWord(info->scalar_c);
  uint64_t x = 0, s = 0; double norm = 0; int invalid = 0;
  for (unsigned long i = 1; i <= X->Check.idim_max; ++i) {
    double re = creal(vector[i]), im = cimag(vector[i]);
    if (!isfinite(re) || !isfinite(im)) invalid = 1;
    norm += re*re + im*im;
  }
  RestartDigestAccumulate(h[RW_OFFSET], vector, X->Check.idim_max, &x, &s);
  norm = SumMPI_d(norm);
  if (comp == RESTART_COMP_STATE && (!isfinite(norm) || fabs(norm - 1.0) > 1e-8)) invalid = 1;
  allreduce_xor_sum(x, s, &h[RW_DATA_XOR], &h[RW_DATA_SUM]);
  if (agree(invalid, "vector validation")) return -1;
  int length = snprintf(final_name, D_FileNameMax, cFileNameRestartVector, X->Def.CDataFileHead,
                        (int)info->set, comp, (int)gen, myrank);
  invalid = length < 0 || length >= D_FileNameMax;
  length = snprintf(part, sizeof(part), "%s%s.part", cParentOutputFolder, final_name);
  invalid |= length < 0 || length >= (int)sizeof(part);
  if (agree(invalid, "path validation")) return -1;
  FILE *fp = fopen(part, "wb");
  invalid = fp == NULL;
#ifdef HPHI_TEST_SYMMETRY_PROBE
  if (!invalid && SymmetryFullDiagTestFail("restart_write")) invalid = 1;
#endif
  for (int i = 0; !invalid && i < RESTART_V3_WORDS; ++i) invalid = RestartPutWord(fp, h[i]) != 0;
  for (unsigned long i = 1; !invalid && i <= X->Check.idim_max; ++i)
    invalid = RestartPutWord(fp, RestartDoubleWord(creal(vector[i]))) != 0 ||
              RestartPutWord(fp, RestartDoubleWord(cimag(vector[i]))) != 0;
  if (fp != NULL) { int flushed = fflush(fp) != 0; int closed = fclose(fp) != 0; invalid |= flushed || closed; }
  if (invalid && fp != NULL) remove(part);
  return agree(invalid, "write");            /* every rank reaches this collective in the same order (Q1) */
}

int RestartWriteUnit(struct BindStruct *X, const struct RestartUnit *unit,
                     double complex *const *const *vectors, const struct RestartVectorInfo *infos, uint64_t flags)
{
  static int fresh_manifest_published = 0;   /* ReStart=1: the first publication replaces the manifest (Q11) */
  struct RestartManifest old = {0}, m = {0}; int exists, invalid = 0;
  uint64_t identity[RESTART_V3_WORDS], gen = 0;
  StartTimer(7100);
  invalid = X == NULL || unit == NULL || vectors == NULL || infos == NULL;
  if (!invalid) invalid = unit->nset <= 0 || unit->ncomp < 1 || unit->ncomp > 2 || unit->sets == NULL ||
                         unit->nset > INT_MAX / unit->ncomp ||
                         (size_t)unit->nset > SIZE_MAX / unit->ncomp / D_FileNameMax;
  if (agree(invalid, "unit validation")) { StopTimer(7100); return -1; }
  for (int i = 0; i < unit->nset; ++i) {
    invalid |= vectors[i] == NULL || unit->sets[i] > INT_MAX || infos[i].set != unit->sets[i];
    if (vectors[i] != NULL) for (int c = 0; c < unit->ncomp; ++c) invalid |= vectors[i][c] == NULL;
  }
  if (agree(invalid, "unit vectors / set indices")) { StopTimer(7100); return -1; }
  char (*names)[D_FileNameMax] = malloc((size_t)unit->nset * unit->ncomp * D_FileNameMax);
  if (RestartManifestLoad(X, &old, &exists) != 0) { free(names); StopTimer(7100); return -1; }
  int fresh = X->Def.iReStart == RESTART_OUT && !fresh_manifest_published;
  m.method = fresh || !exists ? X->Def.iCalcType : old.method;
  m.basis = fresh || !exists ? basis_kind(X) : old.basis;
  /* Keep old intact for collision checks and post-publication cleanup; rows are deep-copied. */
  if (!fresh) for (size_t i = 0; i < old.nrows; ++i)
    if (RestartManifestUpsert(&m, &old.rows[i]) != 0) invalid = 1;
  invalid |= names == NULL;
  if (agree(invalid, "publication allocation")) goto fail;
  if (agree(RestartChooseGeneration(&old, unit->sets, unit->nset, &gen) != 0,
            "no collision-free generation; use a different output directory or CDataFileHead")) goto fail;
  if (agree(build_identity(X, identity) != 0, "identity")) goto fail;
  for (int i = 0; i < unit->nset && !invalid; ++i)
    for (int c = 0; c < unit->ncomp && !invalid; ++c)
      invalid = write_part(X, identity, &infos[i], c, gen, flags, vectors[i][c], names[i * unit->ncomp + c]) != 0;
  /* `invalid` is now identical on every rank (write_part agreed), so breaking here keeps the collective order. */
  if (invalid) { remove_generation(X, unit, gen, 1 /* .part */); goto fail; }
  for (int k = 0; k < unit->nset * unit->ncomp && !invalid; ++k) {
    char part[D_FileNameMax], final[D_FileNameMax];
    snprintf(part, sizeof(part), "%s%s.part", cParentOutputFolder, names[k]);
    snprintf(final, sizeof(final), "%s%s", cParentOutputFolder, names[k]);
#ifdef HPHI_TEST_SYMMETRY_PROBE
    if (SymmetryFullDiagTestFail("restart_rename")) invalid = 1; else
#endif
    invalid = rename(part, final) != 0;
  }
  if (agree(invalid, "rename")) { remove_generation(X, unit, gen, 0); remove_generation(X, unit, gen, 1); goto fail; }
  for (int i = 0; i < unit->nset; ++i) {
    struct RestartManifestRow row = {unit->sets[i], (uint64_t)nproc, gen, infos[i].step, (uint64_t)unit->ncomp};
    if (RestartManifestUpsert(&m, &row) != 0) invalid = 1;
  }
#ifdef HPHI_TEST_SYMMETRY_PROBE
  if (SymmetryFullDiagTestFail("restart_manifest")) invalid = 1;
#endif
  if (agree(invalid, "manifest update")) { remove_generation(X, unit, gen, 0); goto fail; }
  if (RestartManifestPublish(X, &m) != 0) { remove_generation(X, unit, gen, 0); goto fail; }   /* collective on every rank */
  if (fresh) fresh_manifest_published = 1;
  remove_superseded_files(X, &old, &m);   /* published now: delete only old files no longer referenced */
  fprintf(stdoutMPI, "Restart checkpoint published: set %" PRIu64 " step %" PRIu64 " gen %" PRIu64 " flags %" PRIu64 "\n",
          unit->sets[0], infos[0].step, gen, flags);
#ifdef HPHI_TEST_SYMMETRY_PROBE
  if (SumMPI_i(RestartTestAbortAfter(unit->sets[0], infos[0].step)) != 0) {
    fprintf(stdoutMPI, "Test hook: stopping after checkpoint publication.\n");
    RestartManifestFree(&m); RestartManifestFree(&old); free(names); StopTimer(7100); FinalizeMPI(); exit(0);
  }
#endif
  RestartManifestFree(&m); RestartManifestFree(&old); free(names); StopTimer(7100); return 0;
fail:
  RestartManifestFree(&m); RestartManifestFree(&old); free(names); StopTimer(7100); return -1;
}

/* Both continuation manifests and explicitly named state imports use the same
 * header scan and payload checks. A NULL base selects the manifest filenames. */
static int vector_name(char *name, size_t size, const struct BindStruct *X,
                       const struct RestartManifestRow *row, int comp, int rank,
                       const char *base)
{
  int n = base ? snprintf(name, size, "%s_rank_%d.dat", base, rank) :
    snprintf(name, size, cFileNameRestartVector, X->Def.CDataFileHead,
             (int)row->set, comp, (int)row->gen, rank);
  return n < 0 || (size_t)n >= size ? -1 : 0;
}

static int inspect_vector_header(struct BindStruct *X, const struct RestartManifestRow *row, int comp,
                                 enum RestartReadPurpose purpose, const char *base, uint64_t *reference,
                                 uint64_t **saved_offsets, uint64_t **saved_dims)
{
  if (agree(row == NULL, "missing manifest row")) return -1;
  if (purpose == RESTART_READ_IMPORT &&
      agree(comp != RESTART_COMP_STATE, "import requires state component")) return -1;
  const uint64_t P = row->ranks;
  uint64_t expected[RESTART_V3_WORDS];
  uint64_t *offsets = NULL, *dims = NULL; int status = 0;
  if (agree(row->set > INT_MAX || row->gen > 1 || row->ncomp < 1 || row->ncomp > 2 ||
            comp < 0 || comp >= (int)row->ncomp || P > INT_MAX ||
            P > SIZE_MAX / sizeof(uint64_t), "manifest row")) return -1;
  if (agree(X->Def.RestartSourceRanks > 0 && (uint64_t)X->Def.RestartSourceRanks != P, "RestartSourceRanks does not match the saved rank count")) return -1;
  if (agree(build_identity(X, expected) != 0, "identity")) return -1;
  offsets = calloc(P, sizeof(*offsets)); dims = calloc(P, sizeof(*dims));
  if (agree(offsets == NULL || dims == NULL, "allocation")) goto fail;
  if (myrank == 0) {           /* header scan: every old rank, including empty ones */
    uint64_t running = 0;
    for (uint64_t r = 0; r < P && status == 0; ++r) {
      uint64_t h[RESTART_V3_WORDS]; char name[D_FileNameMax];
      if (vector_name(name, sizeof(name), X, row, comp, (int)r, base) != 0) { status = 1; break; }
      FILE *fp = open_output(name, "rb");
      if (fp == NULL) { status = 1; break; }
      for (int i = 0; i < RESTART_V3_WORDS && status == 0; ++i) if (RestartGetWord(fp, &h[i]) != 0) status = 2;
      if (status == 0 && (h[RW_MAGIC] != expected[RW_MAGIC] || h[RW_VERSION] != 3 || h[RW_RANKS] != P || h[RW_RANK] != r)) status = 3;
      if (status == 0 && h[RW_OFFSET] != running) status = 4;
      if (status == 0) {
        long end;
        if (fseek(fp, 0, SEEK_END) != 0 || (end = ftell(fp)) < 0 || h[RW_LOCAL_DIM] > (UINT64_MAX - RESTART_V3_HEADER_BYTES) / 16 ||
            (uint64_t)end != RESTART_V3_HEADER_BYTES + 16 * h[RW_LOCAL_DIM]) status = 5;
      }
      if (status == 0) {
        if (r == 0) memcpy(reference, h, sizeof(h));
        else for (int i = 0; i < RESTART_V3_WORDS; ++i)
          if (i != RW_RANK && i != RW_OFFSET && i != RW_LOCAL_DIM && h[i] != reference[i]) status = 6;
        if (h[RW_LOCAL_DIM] > UINT64_MAX - running) status = 8;
        else { offsets[r] = h[RW_OFFSET]; dims[r] = h[RW_LOCAL_DIM]; running += h[RW_LOCAL_DIM]; }
      }
      fclose(fp);
    }
    if (status == 0 && running != reference[RW_DIM]) status = 7;
  }
  status = BcastMPI_i(0, status);
  if (status != 0) { fprintf(stdoutMPI, "Error: restart checkpoint header scan failed (code %d).\n", status); goto fail; }
  bcast_words(reference, RESTART_V3_WORDS); bcast_words(offsets, P); bcast_words(dims, P);
  const char *bad = NULL;
  int invalid = check_identity(expected, reference, &bad) != 0 || reference[RW_COMPONENT] != (uint64_t)comp;
  if (purpose == RESTART_READ_CONTINUATION) {           /* Q3: method, unit size, step range, int range */
    if (agree(reference[RW_STEP] < step_lower_bound(X->Def.iCalcType) ||
              reference[RW_STEP] > (uint64_t)(INT_MAX - 1), "step range")) goto fail;
    invalid |= reference[RW_STATE] != row->set || reference[RW_STEP] != row->step ||
               reference[RW_METHOD] != (uint64_t)X->Def.iCalcType ||
               (X->Def.iCalcType == CG && reference[RW_NSET] != (uint64_t)X->Def.k_exct);
  }
  if (purpose == RESTART_READ_CONTINUATION) {
    invalid |= X->Def.iCalcType != CG && reference[RW_FLAGS] != 0;
    if (X->Def.iCalcType == cTPQ) {
      double expected_dt = 1.0 / LargeValue, saved_dt = RestartWordDouble(reference[RW_SCALAR_C]);
      invalid |= RestartWordDouble(reference[RW_SCALAR_B]) < 0 || !isfinite(expected_dt) || expected_dt <= 0 ||
                 fabs(saved_dt - expected_dt) > 1e-12 * expected_dt;
    }
  }
  if (!invalid) invalid = !isfinite(RestartWordDouble(reference[RW_TIME])) || !isfinite(RestartWordDouble(reference[RW_SCALAR_A])) ||
                          !isfinite(RestartWordDouble(reference[RW_SCALAR_B])) || !isfinite(RestartWordDouble(reference[RW_SCALAR_C]));
  if (invalid && bad != NULL) fprintf(stderr, "Error: rank %d restart checkpoint %s differs.\n", myrank, bad);
  if (agree(invalid, "identity / unit validation")) goto fail;
  *saved_offsets = offsets; *saved_dims = dims; return 0;
fail:
  free(offsets); free(dims); return -1;
}

static int record_import(const struct BindStruct *X, const struct RestartVectorInfo *info,
                         uint64_t hamiltonian, int matrix_digest)
{
  int invalid = 0;
  FILE *fp = NULL;
  char name[D_FileNameMax];
  if (X->Def.iFlgSymmetryBasis && myrank == 0) {
    int n = snprintf(name, sizeof(name), "%s%ssymmetry_sector.dat",
                     X->Def.iOutputDataHead ? X->Def.CDataFileHead : "",
                     X->Def.iOutputDataHead ? "_" : "");
    invalid = n < 0 || (size_t)n >= sizeof(name);
    if (!invalid) invalid = (fp = open_output(name, "a")) == NULL;
    if (!invalid) {
      fprintf(fp, "source_method=%" PRIu64 "\nsource_state=%" PRIu64 "\nsource_step=%" PRIu64
              "\nsource_time=%.17g\nsource_nset=%" PRIu64 "\nsource_flags=%" PRIu64
              "\nsource_scalar_a=%.17g\nsource_scalar_b=%.17g\nsource_scalar_c=%.17g\n"
              "source_hamiltonian_digest=%016" PRIx64 "\n",
              info->method, info->set, info->step, info->time, info->nset, info->flags,
              info->scalar_a, info->scalar_b, info->scalar_c, hamiltonian);
      if (matrix_digest) fprintf(fp, "source_hamiltonian_kind=sector_matrix_v1\n");
      invalid = ferror(fp) != 0;
      if (fclose(fp) != 0) invalid = 1;
    }
  }
  if (agree(invalid, "import provenance output")) return -1;
  fprintf(stdoutMPI, "Restart checkpoint imported: source_method=%" PRIu64 " source_state=%" PRIu64
          " source_step=%" PRIu64 " source_flags=%" PRIu64 " source_time=%.17g source_nset=%" PRIu64
          " source_scalar_a=%.17g source_scalar_b=%.17g source_scalar_c=%.17g\n",
          info->method, info->set, info->step, info->flags, info->time, info->nset,
          info->scalar_a, info->scalar_b, info->scalar_c);
  return 0;
}

static int read_set_files(struct BindStruct *X, const struct RestartManifestRow *row, int comp,
                          enum RestartReadPurpose purpose, const char *base,
                          double complex *vector, struct RestartVectorInfo *info)
{
  uint64_t reference[RESTART_V3_WORDS], expected[RESTART_V3_WORDS];
  uint64_t *offsets = NULL, *dims = NULL;
  int invalid = 0;
  if (inspect_vector_header(X, row, comp, purpose, base, reference, &offsets, &dims) != 0) return -1;
  const uint64_t P = row->ranks;
  /* expected is used only for the Hamiltonian provenance message below. */
  if (agree(build_identity(X, expected) != 0, "identity")) goto fail;
  /* intersection reads (S3): lo/hi per old file */
  uint64_t a = RestartLocalOffset(X);
  if (agree(a > reference[RW_DIM] || X->Check.idim_max > reference[RW_DIM] - a ||
            X->Check.idim_max >= SIZE_MAX / sizeof(double complex), "local interval / allocation size")) goto fail;
  uint64_t b = a + X->Check.idim_max;
  double complex *tmp = calloc(X->Check.idim_max + 1, sizeof(*tmp));
  unsigned long payload_files_opened = 0;
  invalid = tmp == NULL;
  for (uint64_t j = 0; !invalid && j < P; ++j) {
    uint64_t lo, hi;
    if (!RestartIntersect(a, b, offsets[j], dims[j], &lo, &hi)) continue;
    char name[D_FileNameMax];
    if (vector_name(name, sizeof(name), X, row, comp, (int)j, base) != 0) { invalid = 1; break; }
    FILE *fp = open_output(name, "rb");
    if (fp != NULL) ++payload_files_opened;
    uint64_t byte = RESTART_V3_HEADER_BYTES + 16 * (lo - offsets[j]);
    invalid = fp == NULL || byte > (uint64_t)LONG_MAX || fseek(fp, (long)byte, SEEK_SET) != 0;
    for (uint64_t g = lo; !invalid && g < hi; ++g) {
      uint64_t re, im;
      invalid = RestartGetWord(fp, &re) != 0 || RestartGetWord(fp, &im) != 0;
      if (!invalid) { double c[2] = {RestartWordDouble(re), RestartWordDouble(im)}; memcpy(&tmp[(g - a) + 1], c, sizeof(c)); }
    }
    if (!invalid && hi == offsets[j] + dims[j] && (fgetc(fp) != EOF || ferror(fp))) invalid = 1;
    if (fp != NULL) fclose(fp);
  }
  if (agree(invalid, "payload read")) { free(tmp); goto fail; }
  payload_files_opened = SumMPI_li(payload_files_opened);
  fprintf(stdoutMPI, "Restart read: header files %" PRIu64
          ", payload files opened %lu.\n", P, payload_files_opened);
  uint64_t x = 0, s = 0, gx, gs; double norm = 0;
  for (unsigned long i = 1; i <= X->Check.idim_max; ++i) {
    double re = creal(tmp[i]), im = cimag(tmp[i]);
    if (!isfinite(re) || !isfinite(im)) invalid = 1;
    norm += re*re + im*im;
  }
  RestartDigestAccumulate(a, tmp, X->Check.idim_max, &x, &s);
  norm = SumMPI_d(norm);
  allreduce_xor_sum(x, s, &gx, &gs);
  invalid |= gx != reference[RW_DATA_XOR] || gs != reference[RW_DATA_SUM];
  if (comp == RESTART_COMP_STATE) invalid |= !isfinite(norm) || fabs(norm - 1.0) > 1e-8;
  if (agree(invalid, "payload checksum / norm")) { free(tmp); goto fail; }
  memcpy(vector, tmp, (X->Check.idim_max + 1) * sizeof(*vector)); free(tmp);
  info->method = reference[RW_METHOD]; info->set = reference[RW_STATE]; info->step = reference[RW_STEP];
  info->time = RestartWordDouble(reference[RW_TIME]); info->nset = reference[RW_NSET]; info->flags = reference[RW_FLAGS];
  info->scalar_a = RestartWordDouble(reference[RW_SCALAR_A]); info->scalar_b = RestartWordDouble(reference[RW_SCALAR_B]);
  info->scalar_c = RestartWordDouble(reference[RW_SCALAR_C]);
  fprintf(stdoutMPI, "Restart checkpoint loaded: set=%" PRIu64 " comp=%d step=%" PRIu64 " flags=%" PRIu64 " time=%.17g hamiltonian_changed=%s\n",
          info->set, comp, info->step, info->flags, info->time,
          reference[RW_BASIS_KIND] == RESTART_BASIS_RAW ? "unknown" : (reference[RW_HAMILTONIAN] == expected[RW_HAMILTONIAN] ? "no" : "yes"));
  if (purpose == RESTART_READ_IMPORT &&
      record_import(X, info, reference[RW_HAMILTONIAN], info->method == FullDiag) != 0) goto fail;
  free(offsets); free(dims); return 0;
fail:
  free(offsets); free(dims); return -1;
}

int RestartReadVector(struct BindStruct *X, const struct RestartManifestRow *row, int comp,
                      enum RestartReadPurpose purpose, double complex *vector, struct RestartVectorInfo *info)
{
  int result;
  StartTimer(7110);
  result = read_set_files(X, row, comp, purpose, NULL, vector, info);
  StopTimer(7110);
  return result;
}

/* Retain the old import callers' high-level header diagnostic when root's
 * format probe catches the failure before the legacy reader can report it. */
static int import_probe_failure(const struct BindStruct *X, int invalid, const char *stage)
{
  if (agree(invalid, stage) == 0) return 0;
  if (X->Def.iFlgSymmetryBasis)
    fprintf(stdoutMPI, "Error: symmetry checkpoint header / sector / layout validation failed.\n");
  else if (X->Def.iCalcType == TimeEvolution)
    fprintf(stdoutMPI, "Error: missing, truncated or incompatible TE Inputvector file.\n");
  return -1;
}

static int sector_legacy_runtime_identity(const uint64_t *expected,
                                          const uint64_t *header, int words)
{
  for (int i = RW_MAGIC; i <= RW_DIM; ++i)
    if (header[i] != expected[i]) return -1;
  for (int i = RW_GROUP; i <= RW_SECTOR_SUM; ++i)
    if (header[i] != expected[i]) return -1;
  if (words == 30 &&
      (header[RW_LOCAL_SITE_MASK] != expected[RW_LOCAL_SITE_MASK] ||
       header[RW_FIXED_FLAGS] != expected[RW_FIXED_FLAGS])) return -1;
  return 0;
}

static int sector_legacy_same_unit(const uint64_t *reference,
                                   const uint64_t *header, int words)
{
  for (int i = 0; i < words; ++i) {
    if (i == RW_RANK || i == RW_OFFSET || i == RW_LOCAL_DIM || i == RW_ORDER)
      continue;
    if (header[i] != reference[i]) return -1;
  }
  return 0;
}

/* Read a sector v1/v2 state with a different rank count.  The old DATA words
 * hash each writer-rank partition, so first reproduce that checksum over every
 * old file, then read only the intervals intersecting this run's ownership. */
static int read_sector_cross_rank(struct BindStruct *X, const char *base,
                                  int version, int source_ranks,
                                  double complex *vector,
                                  struct SymmetryCheckpointInfo *info)
{
  const int words = version == RESTART_FORMAT_SECTOR_V2 ? 30 : 28;
  const uint64_t P = (uint64_t)source_ranks;
  uint64_t expected[RESTART_V3_WORDS] = {0};
  uint64_t reference[RESTART_V3_WORDS] = {0};
  uint64_t *offsets = NULL, *dims = NULL;
  int status = 0, invalid = 0;
  if (agree(base == NULL || vector == NULL || info == NULL ||
            (version != RESTART_FORMAT_SECTOR_V1 &&
             version != RESTART_FORMAT_SECTOR_V2) || source_ranks < 1,
            "sector cross-rank arguments")) return -1;
  if (SymmetryCheckpointIdentityWords(X, expected) != 0) return -1;
  offsets = calloc(P, sizeof(*offsets));
  dims = calloc(P, sizeof(*dims));
  if (agree(offsets == NULL || dims == NULL, "sector cross-rank allocation")) goto fail;

  if (myrank == 0) {
    uint64_t running = 0;
    for (uint64_t rank = 0; rank < P && status == 0; ++rank) {
      uint64_t header[RESTART_V3_WORDS] = {0};
      char name[D_FileNameMax];
      FILE *fp = NULL;
      if (vector_name(name, sizeof(name), X, NULL, 0, (int)rank, base) != 0) {
        status = 1;
        break;
      }
      fp = open_output(name, "rb");
      if (fp == NULL) {
        status = 2;
        break;
      }
      for (int i = 0; i < words && status == 0; ++i)
        if (RestartGetWord(fp, &header[i]) != 0) status = 3;
      if (status == 0 &&
          (sector_legacy_runtime_identity(expected, header, words) != 0 ||
           header[RW_RANKS] != P || header[RW_RANK] != rank ||
           header[RW_LAYOUT] > 1 || header[RW_OFFSET] != running)) status = 4;
      if (status == 0) {
        long end = -1;
        uint64_t bytes = (uint64_t)(8 * words);
        if (header[RW_LOCAL_DIM] > (UINT64_MAX - bytes) / 16 ||
            fseek(fp, 0, SEEK_END) != 0 || (end = ftell(fp)) < 0 ||
            (uint64_t)end != bytes + 16 * header[RW_LOCAL_DIM]) status = 5;
      }
      if (status == 0) {
        double time = RestartWordDouble(header[RW_TIME]);
        invalid = !isfinite(time) ||
                  (header[RW_METHOD] != FullDiag && header[RW_METHOD] != CG &&
                   header[RW_METHOD] != TimeEvolution);
        if (header[RW_METHOD] == FullDiag)
          invalid |= header[RW_STEP] != 0 || header[RW_TIME] != 0 ||
                     header[RW_STATE] >= header[RW_DIM];
        if (invalid) status = 6;
      }
      if (status == 0) {
        if (rank == 0) memcpy(reference, header, sizeof(reference));
        else if (sector_legacy_same_unit(reference, header, words) != 0) status = 7;
      }
      if (status == 0) {
        if (header[RW_LOCAL_DIM] > UINT64_MAX - running) status = 8;
        else {
          offsets[rank] = header[RW_OFFSET];
          dims[rank] = header[RW_LOCAL_DIM];
          running += header[RW_LOCAL_DIM];
        }
      }
      if (fclose(fp) != 0 && status == 0) status = 9;
    }
    if (status == 0 && running != reference[RW_DIM]) status = 10;
  }
  status = BcastMPI_i(0, status);
  if (agree(status != 0, "sector cross-rank header scan")) goto fail;
  bcast_words(reference, (uint64_t)words);
  bcast_words(offsets, P);
  bcast_words(dims, P);

  uint64_t local_xor = 0, local_sum = 0, global_xor = 0, global_sum = 0;
  double local_norm = 0;
  for (uint64_t rank = (uint64_t)myrank; !invalid && rank < P; rank += (uint64_t)nproc) {
    char name[D_FileNameMax];
    FILE *fp = NULL;
    uint64_t hash;
    double norm;
    invalid = vector_name(name, sizeof(name), X, NULL, 0, (int)rank, base) != 0;
    if (!invalid) fp = open_output(name, "rb");
    invalid |= fp == NULL || fseek(fp, 8L * words, SEEK_SET) != 0 ||
               SymmetryCheckpointLegacyFileHash(fp, rank, offsets[rank],
                                                dims[rank], &hash, &norm) != 0;
    if (!invalid) {
      local_xor ^= hash;
      local_sum += hash;
      local_norm += norm;
    }
    if (fp != NULL && fclose(fp) != 0) invalid = 1;
  }
  if (agree(invalid, "sector v1/v2 old-partition DATA read")) goto fail;
  allreduce_xor_sum(local_xor, local_sum, &global_xor, &global_sum);
  double norm = SumMPI_d(local_norm);
  invalid = global_xor != reference[RW_DATA_XOR] ||
            global_sum != reference[RW_DATA_SUM] || !isfinite(norm) ||
            fabs(norm - 1.0) > 1e-8;
  if (agree(invalid, "sector v1/v2 payload checksum / norm")) goto fail;

  uint64_t a = RestartLocalOffset(X);
  if (agree(a > reference[RW_DIM] || X->Check.idim_max > reference[RW_DIM] - a ||
            X->Check.idim_max >= SIZE_MAX / sizeof(double complex),
            "sector cross-rank local interval")) goto fail;
  uint64_t b = a + X->Check.idim_max;
  double complex *temporary = calloc(X->Check.idim_max + 1, sizeof(*temporary));
  invalid = temporary == NULL;
  for (uint64_t rank = 0; !invalid && rank < P; ++rank) {
    uint64_t lo, hi;
    if (!RestartIntersect(a, b, offsets[rank], dims[rank], &lo, &hi)) continue;
    char name[D_FileNameMax];
    FILE *fp = NULL;
    uint64_t byte = (uint64_t)(8 * words) + 16 * (lo - offsets[rank]);
    invalid = vector_name(name, sizeof(name), X, NULL, 0, (int)rank, base) != 0;
    if (!invalid) fp = open_output(name, "rb");
    invalid |= fp == NULL || byte > (uint64_t)LONG_MAX ||
               fseek(fp, (long)byte, SEEK_SET) != 0;
    for (uint64_t global = lo; !invalid && global < hi; ++global) {
      uint64_t re_word, im_word;
      invalid = RestartGetWord(fp, &re_word) != 0 ||
                RestartGetWord(fp, &im_word) != 0;
      if (!invalid) {
        double components[2] = {RestartWordDouble(re_word),
                                RestartWordDouble(im_word)};
        if (!isfinite(components[0]) || !isfinite(components[1])) invalid = 1;
        else memcpy(&temporary[(global - a) + 1], components, sizeof(components));
      }
    }
    if (fp != NULL && fclose(fp) != 0) invalid = 1;
  }
  if (agree(invalid, "sector cross-rank payload read")) {
    free(temporary);
    goto fail;
  }
  local_norm = 0;
  for (unsigned long i = 1; i <= X->Check.idim_max; ++i) {
    double re = creal(temporary[i]), im = cimag(temporary[i]);
    local_norm += re*re + im*im;
  }
  norm = SumMPI_d(local_norm);
  if (agree(!isfinite(norm) || fabs(norm - 1.0) > 1e-8,
            "sector cross-rank redistributed norm")) {
    free(temporary);
    goto fail;
  }
  memcpy(vector, temporary,
         (X->Check.idim_max + 1) * sizeof(*temporary));
  free(temporary);
  info->source_method = reference[RW_METHOD];
  info->state_index = reference[RW_STATE];
  info->step = reference[RW_STEP];
  info->time = RestartWordDouble(reference[RW_TIME]);
  info->hamiltonian_digest = reference[RW_HAMILTONIAN];
  info->hamiltonian_digest_kind = reference[RW_METHOD] == FullDiag
      ? SYMMETRY_HAMILTONIAN_MATRIX_V1 : SYMMETRY_HAMILTONIAN_PARSED;
  fprintf(stdoutMPI, "Symmetry checkpoint cross-rank v1/v2 read: payload checksum verified, "
          "basis order guaranteed by sector digest (ORDER word not compared).\n");
  free(offsets);
  free(dims);
  return 0;
fail:
  free(offsets);
  free(dims);
  return -1;
}

int RestartImportVector(struct BindStruct *X, const char *base, double complex *vector,
                        struct RestartVectorInfo *provenance, enum RestartFileFormat *format, int *source_ranks)
{
  int status = 0, kind = RESTART_FORMAT_LEGACY_RAW, ranks = nproc;
  if (agree(base == NULL || vector == NULL || provenance == NULL || format == NULL ||
            source_ranks == NULL, "import arguments")) return -1;
  if (myrank == 0) {
    char name[D_FileNameMax];
    FILE *fp = NULL;
    uint64_t magic, version, saved_ranks;
    status = vector_name(name, sizeof(name), X, NULL, 0, 0, base) != 0;
    if (!status) status = (fp = open_output(name, "rb")) == NULL;
    if (!status) status = RestartGetWord(fp, &magic) != 0 || RestartGetWord(fp, &version) != 0;
    if (!status && magic == UINT64_C(0x0a31565349485048)) {
      if (version < 1 || version > 3) status = 1;
      else {
        kind = (int)version;
        status = fseek(fp, 8 * RW_RANKS, SEEK_SET) != 0 || RestartGetWord(fp, &saved_ranks) != 0;
        if (!status) {
          status = saved_ranks == 0 || saved_ranks > INT_MAX;
          if (status) fprintf(stdoutMPI, "Error: restart checkpoint MPI rank count differs.\n");
          else ranks = (int)saved_ranks;
        }
      }
    }
    if (fp && fclose(fp) != 0) status = 1;
  }
  status = BcastMPI_i(0, status);
  kind = BcastMPI_i(0, kind); ranks = BcastMPI_i(0, ranks);
  if (import_probe_failure(X, status, "import format / version / source rank count")) return -1;
  *format = (enum RestartFileFormat)kind; *source_ranks = ranks;
  if (agree(X->Def.RestartSourceRanks > 0 && X->Def.RestartSourceRanks != ranks,
            "RestartSourceRanks does not match the saved rank count")) return -1;
  memset(provenance, 0, sizeof(*provenance));
  if (kind == RESTART_FORMAT_V3) {
    const struct RestartManifestRow row = {0, (uint64_t)ranks, 0, 0, 1};
    return read_set_files(X, &row, RESTART_COMP_STATE, RESTART_READ_IMPORT, base, vector, provenance);
  }
  if (kind == RESTART_FORMAT_LEGACY_RAW) {
    if (import_probe_failure(X, X->Def.iFlgSymmetryBasis, "legacy raw import requires a raw basis")) return -1;
    return -2;
  }
  if (agree(!X->Def.iFlgSymmetryBasis, "sector import requires a sector basis")) return -1;
  struct SymmetryCheckpointInfo info;
  if (ranks != nproc) {
    if (read_sector_cross_rank(X, base, kind, ranks, vector, &info) != 0) return -1;
  } else {
    char name[D_FileNameMax];
    if (agree(vector_name(name, sizeof(name), X, NULL, 0, myrank, base) != 0,
              "import filename")) return -1;
    if (ReadSymmetryCheckpoint(X, name, vector, &info) != 0) return -1;
  }
  provenance->method = info.source_method; provenance->set = info.state_index;
  provenance->step = info.step; provenance->time = info.time;
  return record_import(X, provenance, info.hamiltonian_digest,
                       info.hamiltonian_digest_kind == SYMMETRY_HAMILTONIAN_MATRIX_V1);
}

static int validate_set_boundary(struct BindStruct *X, const struct RestartManifest *m, uint64_t set)
{
  const struct RestartManifestRow *row = RestartManifestFind(m, set);
  uint64_t first[RESTART_V3_WORDS], h[RESTART_V3_WORDS];
  int bad = row == NULL || m->method != X->Def.iCalcType || m->basis != basis_kind(X);
  if (!bad) bad = row->ncomp != (uint64_t)(X->Def.iCalcType == CG ? 2 : 1);
  if (agree(bad, "recovery boundary manifest")) return -1;
  for (int c = 0; c < (int)row->ncomp; ++c) {
    uint64_t *offsets = NULL, *dims = NULL;
    if (inspect_vector_header(X, row, c, RESTART_READ_CONTINUATION, NULL, h, &offsets, &dims) != 0) return -1;
    free(offsets); free(dims);
    if (c == 0) memcpy(first, h, sizeof(first));
    else for (int w = 0; w < RESTART_V3_WORDS; ++w)
      if (w != RW_COMPONENT && w != RW_DATA_XOR && w != RW_DATA_SUM && first[w] != h[w]) bad = 1;
  }
  if (agree(bad, "component scalars differ / recovery boundary")) return -1;
  /* A saved CG block has one basis identity, including its recorded layout.
   * Current-vs-saved layout compatibility is a separate reader decision. */
  const struct RestartManifestRow *first_row = RestartManifestFind(m, 0);
  if (X->Def.iCalcType == CG && set != 0 && first_row != NULL) {
    uint64_t *offsets = NULL, *dims = NULL;
    if (inspect_vector_header(X, first_row, RESTART_COMP_STATE, RESTART_READ_CONTINUATION, NULL,
                              h, &offsets, &dims) != 0) return -1;
    free(offsets); free(dims);
    for (size_t i = 0; i < sizeof(identity_words) / sizeof(identity_words[0]); ++i)
      if (first[identity_words[i]] != h[identity_words[i]]) bad = 1;
    bad |= first[RW_LAYOUT] != h[RW_LAYOUT];
  }
  return agree(bad, "CG state identities differ / recovery boundary");
}

int RestartClassifySet(struct BindStruct *X, const struct RestartManifest *manifest, uint64_t set,
                       const char *legacy_pattern, enum RestartSource *source, int *legacy_step)
{
  *source = RESTART_SOURCE_NONE; *legacy_step = 0;
  if (agree(set > INT_MAX || (manifest != NULL &&
            (manifest->method != X->Def.iCalcType || manifest->basis != basis_kind(X))),
            "manifest method / basis / set")) return -1;
  if (manifest != NULL && RestartManifestFind(manifest, set) != NULL) {
    if (validate_set_boundary(X, manifest, set) != 0) return -1;
    *source = RESTART_SOURCE_V3; return 0;
  }
  if (legacy_pattern == NULL) return 0;
  if (X->Def.RestartSourceRanks > 0 && X->Def.RestartSourceRanks != nproc) {
    fprintf(stdoutMPI, "Error: legacy restart files with a different rank count arrive in a later release (RestartSourceRanks).\n");
    return -1;   /* PR-3 lifts this */
  }
  char name[D_FileNameMax];
  if (X->Def.iCalcType == TimeEvolution) snprintf(name, sizeof(name), legacy_pattern, myrank);
  else snprintf(name, sizeof(name), legacy_pattern, (int)set, myrank);
  FILE *fp = open_output(name, "rb"); int present = fp != NULL, step = 0, bad = fp == NULL && errno != ENOENT;
  uint64_t dim = 0;
  if (present) {
    long end = -1;
    bad |= fread(&step, sizeof(int), 1, fp) != 1 || fread(&dim, sizeof(dim), 1, fp) != 1;
    bad |= step < (int)step_lower_bound(X->Def.iCalcType) || step > INT_MAX - 1;
    bad |= dim != (uint64_t)X->Check.idim_max || dim > (UINT64_MAX - 28) / 16;
    bad |= fseek(fp, 0, SEEK_END) != 0 || (end = ftell(fp)) < 0;
    if (!bad) bad = (uint64_t)end != RESTART_LEGACY_HEADER_BYTES + 16 * (dim + 1);
    bad |= fclose(fp) != 0;
  }
  if (agree(bad, "legacy recovery header / step / file length")) return -1;
  int count = SumMPI_i(present);
  if (count == 0) return 0;
  if (count != nproc) { fprintf(stdoutMPI, "Error: legacy restart files for set %" PRIu64 " exist on %d of %d ranks.\n", set, count, nproc); return -1; }
  int root_step = BcastMPI_i(0, step);
  if (agree(bad || step != root_step, "legacy header step agreement")) return -1;
  *source = RESTART_SOURCE_LEGACY; *legacy_step = step;
  return 0;
}
