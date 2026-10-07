#include <float.h>
#include <inttypes.h>
#include <math.h>
#include "Common.h"
#include "symmetry_basis.h"
#include "symmetry_sector.h"
#include "symmetry_checkpoint.h"
#include "wrapperMPI.h"
#ifdef MPI
#include <mpi.h>
#endif

/* Version 1: 28 little-endian uint64 words, followed by local_dim pairs of
 * IEEE binary64 real/imaginary components (index zero is not serialized).
 * Phase 1: normalized sum_g conjugate(chi(g)) T_g |smallest representative>,
 * with a positive real coefficient at that representative. */
enum {
  MAGIC, VERSION, PHASE, SCALAR, MODEL, NSITE, NUP, NDOWN, NE,
  FULL_DIM, DIM, RANKS, RANK, LAYOUT, OFFSET, LOCAL_DIM, GROUP,
  SECTOR_COUNT, SECTOR_XOR, SECTOR_SUM, ORDER, HAMILTONIAN,
  METHOD, STATE, STEP, TIME, DATA_XOR, DATA_SUM, HEADER_WORDS
};
static const char *identity_names[HAMILTONIAN] = {
  "magic", "version", "phase convention", "scalar encoding", "model", "nsite",
  "nup", "ndown", "particle number", "raw dimension", "sector dimension",
  "MPI rank count", "MPI rank", "basis layout", "local offset", "local dimension",
  "group/character digest", "sector count", "sector XOR", "sector sum", "basis order"
};
#define CHECKPOINT_MAGIC UINT64_C(0x0a31565349485048)
#define FNV_OFFSET UINT64_C(14695981039346656037)
#define FNV_PRIME UINT64_C(1099511628211)

static uint64_t hash_word(uint64_t hash, uint64_t value)
{
  for (int i = 0; i < 8; ++i) { hash ^= value & 255; hash *= FNV_PRIME; value >>= 8; }
  return hash;
}
static uint64_t double_word(double value)
{
  uint64_t word;
  memcpy(&word, &value, 8);
  return word;
}
static double word_double(uint64_t word)
{
  double value;
  memcpy(&value, &word, 8);
  return value;
}
static int put_word(FILE *fp, uint64_t word)
{
  unsigned char bytes[8];
  for (int i = 0; i < 8; ++i) { bytes[i] = word & 255; word >>= 8; }
  return fwrite(bytes, 1, 8, fp) == 8 ? 0 : -1;
}
static int get_word(FILE *fp, uint64_t *word)
{
  unsigned char bytes[8];
  if (fread(bytes, 1, 8, fp) != 8) return -1;
  *word = 0;
  for (int i = 7; i >= 0; --i) *word = (*word << 8) | bytes[i];
  return 0;
}
static int failure(int invalid, const char *stage)
{
  if (SumMPI_i(invalid != 0) == 0) return 0;
  fprintf(stdoutMPI, "Error: symmetry checkpoint %s failed.\n", stage);
  return -1;
}
static int make_header(const struct BindStruct *X, const char *name,
                       uint64_t *header, char *path)
{
  struct SymmetrySectorDigest sector;
  uint64_t group = 0, hamiltonian = 0, order = FNV_OFFSET;
  int invalid = X == NULL || name == NULL || sizeof(double) != 8 || sizeof(double complex) != 16 ||
                DBL_MANT_DIG != 53 || DBL_MAX_EXP != 1024;
  if (!invalid)
    invalid = !X->Def.iFlgSymmetryBasis || X->Sym == NULL ||
              X->Sym->enabled != TRUE || X->Check.idim_max != X->Sym->local_dim ||
              X->Check.idim_maxMPI != X->Sym->dim;
  if (failure(invalid, "runtime validation")) return -1;
  int length = snprintf(path, D_FileNameMax, "%s%s", cParentOutputFolder, name);
  invalid = length < 0 || length >= D_FileNameMax || name[0] == '\0';
  if (failure(invalid, "path validation")) return -1;
  if (ComputeSymmetrySectorDigest(X->Sym, &sector) != 0) return -1;
  invalid = ComputeSymmetryGroupDigest(&X->Def, &group) != 0 ||
            ComputeSymmetryHamiltonianDigest(&X->Def, &hamiltonian) != 0;
  for (unsigned long i = 1; !invalid && i <= X->Sym->local_dim; ++i) {
    const struct SymmetryBasisVector *entry = SymmetryBasisLocalEntry(X->Sym, i);
    if (entry == NULL) { invalid = 1; break; }
    order = hash_word(order, X->Sym->local_offset + i);
    order = hash_word(order, entry->rep_state);
    order = hash_word(order, entry->orbit_size);
    order = hash_word(order, entry->stabilizer_size);
  }
  if (failure(invalid, "basis identity")) return -1;
  memset(header, 0, HEADER_WORDS * sizeof(*header));
  header[MAGIC] = CHECKPOINT_MAGIC; header[VERSION] = 1;
  header[PHASE] = 1; header[SCALAR] = 128; /* little-endian complex binary64 */
  header[MODEL] = X->Def.iCalcModel; header[NSITE] = X->Def.Nsite;
  header[NUP] = X->Def.iCalcModel == SpinGC ? 0 : X->Def.Nup;
  header[NDOWN] = X->Def.iCalcModel == SpinGC ? 0 : X->Def.Ndown;
  header[NE] = X->Def.iCalcModel == SpinGC ? 0 : X->Def.Ne;
  header[FULL_DIM] = X->Sym->full_dim; header[DIM] = X->Sym->dim;
  header[RANKS] = nproc; header[RANK] = myrank;
  header[LAYOUT] = X->Sym->basis_layout;
  header[OFFSET] = X->Sym->local_offset; header[LOCAL_DIM] = X->Sym->local_dim;
  header[GROUP] = group; header[SECTOR_COUNT] = sector.count;
  header[SECTOR_XOR] = sector.xor_hash; header[SECTOR_SUM] = sector.sum_hash;
  header[ORDER] = order; header[HAMILTONIAN] = hamiltonian;
  return 0;
}
/* Catch mixed checkpoint sets, including different H/time/state provenance.
 * Rank-specific identity fields were separately checked against the runtime. */
static int consistent_header(const uint64_t *header)
{
  uint64_t local[HEADER_WORDS], root[HEADER_WORDS];
  memcpy(local, header, sizeof(local));
  local[RANK] = local[OFFSET] = local[LOCAL_DIM] = local[ORDER] = 0;
  memcpy(root, local, sizeof(root));
#ifdef MPI
  if (MPI_Bcast(root, HEADER_WORDS, MPI_UINT64_T, 0, MPI_COMM_WORLD) != MPI_SUCCESS) return -1;
#endif
  return failure(memcmp(root, local, sizeof(root)) != 0, "cross-rank metadata consistency");
}
static int vector_digest(const struct BindStruct *X, const double complex *vector,
                          uint64_t *xor_hash, uint64_t *sum_hash)
{
  uint64_t hash = hash_word(FNV_OFFSET, (uint64_t)myrank);
  double norm = 0;
  int invalid = vector == NULL;
  hash = hash_word(hash, X->Sym->local_offset);
  hash = hash_word(hash, X->Sym->local_dim);
  for (unsigned long i = 1; !invalid && i <= X->Sym->local_dim; ++i) {
    double re = creal(vector[i]), im = cimag(vector[i]);
    if (!isfinite(re) || !isfinite(im)) { invalid = 1; break; }
    hash = hash_word(hash_word(hash, double_word(re)), double_word(im));
    norm += re*re + im*im;
  }
  if (failure(invalid, "finite vector validation")) return -1;
  norm = SumMPI_d(norm);
  if (failure(!isfinite(norm) || fabs(norm - 1.0) > 1e-8, "global norm validation")) return -1;
  *xor_hash = *sum_hash = hash;
#ifdef MPI
  if (MPI_Allreduce(&hash, xor_hash, 1, MPI_UINT64_T, MPI_BXOR, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Allreduce(&hash, sum_hash, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS) return -1;
#endif
  return 0;
}

int WriteSymmetryCheckpoint(const struct BindStruct *X, const char *name,
                            const double complex *vector,
                            const struct SymmetryCheckpointInfo *info)
{
  uint64_t header[HEADER_WORDS];
  char path[D_FileNameMax], temporary[D_FileNameMax];
  FILE *fp = NULL;
  int invalid, created;
  if (make_header(X, name, header, path) != 0) return -1;
  if (failure(info == NULL || vector == NULL, "output arguments")) return -1;
  if (failure(!isfinite(info->time) || (info->source_method != CG &&
                                      info->source_method != TimeEvolution), "output provenance")) return -1;
  header[METHOD] = info->source_method; header[STATE] = info->state_index;
  header[STEP] = info->step; header[TIME] = double_word(info->time);
  if (vector_digest(X, vector, &header[DATA_XOR], &header[DATA_SUM]) != 0 ||
      consistent_header(header) != 0) return -1;
  int length = snprintf(temporary, sizeof(temporary), "%s.part", path);
  if (failure(length < 0 || (size_t)length >= sizeof(temporary), "temporary path")) return -1;
  fp = fopen(temporary, "wb");
  created = fp != NULL;
  invalid = fp == NULL;
  for (int i = 0; !invalid && i < HEADER_WORDS; ++i) invalid = put_word(fp, header[i]) != 0;
  for (unsigned long i = 1; !invalid && i <= X->Sym->local_dim; ++i)
    invalid = put_word(fp, double_word(creal(vector[i]))) != 0 ||
              put_word(fp, double_word(cimag(vector[i]))) != 0;
  if (fp != NULL && fclose(fp) != 0) invalid = 1;
  if (failure(invalid, "write")) { if (created) remove(temporary); return -1; }
  invalid = rename(temporary, path) != 0;
  if (failure(invalid, "publish")) { remove(temporary); return -1; }
  return 0;
}

int ReadSymmetryCheckpoint(const struct BindStruct *X, const char *name,
                           double complex *vector,
                           struct SymmetryCheckpointInfo *info)
{
  uint64_t expected[HEADER_WORDS], header[HEADER_WORDS] = {0}, xor_hash, sum_hash;
  char path[D_FileNameMax];
  FILE *fp = NULL;
  double complex *temporary = NULL;
  int invalid;
  if (make_header(X, name, expected, path) != 0) return -1;
  if (failure(info == NULL || vector == NULL, "input arguments")) return -1;
  fp = fopen(path, "rb");
  invalid = fp == NULL;
  for (int i = 0; !invalid && i < HEADER_WORDS; ++i) invalid = get_word(fp, &header[i]) != 0;
  /* H is deliberately excluded: this is a quench-capable state import. */
  for (int i = 0; !invalid && i < HAMILTONIAN; ++i) {
    invalid = header[i] != expected[i];
    if (invalid) fprintf(stderr, "Error: rank %d checkpoint %s: %s differs (file=%" PRIu64
                         ", expected=%" PRIu64 ").\n", myrank, name,
                         identity_names[i], header[i], expected[i]);
  }
  if (!invalid) invalid = !isfinite(word_double(header[TIME])) ||
                         (header[METHOD] != CG && header[METHOD] != TimeEvolution);
  if (failure(invalid, "header / sector / layout validation")) goto fail;
  if (consistent_header(header) != 0) goto fail;
  if (X->Sym->local_dim >= SIZE_MAX / sizeof(*temporary)) invalid = 1;
  if (!invalid) temporary = calloc(X->Sym->local_dim + 1, sizeof(*temporary));
  if (failure(invalid || temporary == NULL, "input allocation")) goto fail;
  for (unsigned long i = 1; !invalid && i <= X->Sym->local_dim; ++i) {
    uint64_t re, im;
    if (get_word(fp, &re) || get_word(fp, &im)) { invalid = 1; break; }
    double components[2] = {word_double(re), word_double(im)};
    memcpy(&temporary[i], components, sizeof(components)); /* preserve signed zero */
  }
  if (!invalid && (fgetc(fp) != EOF || ferror(fp))) invalid = 1;
  if (fclose(fp) != 0) invalid = 1;
  fp = NULL;
  if (failure(invalid, "payload length")) goto fail;
  if (vector_digest(X, temporary, &xor_hash, &sum_hash) != 0) goto fail;
  if (failure(xor_hash != header[DATA_XOR] || sum_hash != header[DATA_SUM], "payload checksum")) goto fail;
  memcpy(vector, temporary, (X->Sym->local_dim + 1) * sizeof(*vector));
  info->source_method = header[METHOD]; info->state_index = header[STATE];
  info->step = header[STEP]; info->time = word_double(header[TIME]);
  info->hamiltonian_digest = header[HAMILTONIAN];
  fprintf(stdoutMPI, "Symmetry checkpoint loaded: state=%" PRIu64 " step=%" PRIu64
          " time=%.17g source_hamiltonian=%016" PRIx64 " hamiltonian_changed=%s\n",
          info->state_index, info->step, info->time, info->hamiltonian_digest,
          header[HAMILTONIAN] != expected[HAMILTONIAN] ? "yes" : "no");
  free(temporary);
  return 0;
fail:
  if (fp != NULL) fclose(fp);
  free(temporary);
  return -1;
}
