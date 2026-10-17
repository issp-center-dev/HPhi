/* src/include/restart_checkpoint.h（PR-2 Task 9で作成。以降のtaskはこの宣言を使う） */
#ifndef HPHI_RESTART_CHECKPOINT_H
#define HPHI_RESTART_CHECKPOINT_H
#include <complex.h>
#include <stdint.h>
#include <stdio.h>
struct BindStruct;
struct RestartManifest;
struct RestartManifestRow;

enum { RESTART_V3_WORDS = 40, RESTART_V3_HEADER_BYTES = 320, RESTART_LEGACY_HEADER_BYTES = 12 };
enum RestartBasisKind { RESTART_BASIS_RAW = 0, RESTART_BASIS_SECTOR = 1 };
enum RestartComponent { RESTART_COMP_STATE = 0, RESTART_COMP_DIRECTION = 1 };
enum RestartReadPurpose { RESTART_READ_CONTINUATION = 0, RESTART_READ_IMPORT = 1 };
enum RestartSource { RESTART_SOURCE_NONE = 0, RESTART_SOURCE_V3 = 1, RESTART_SOURCE_LEGACY = 2 };

/* Header word indices (words 0-29 share positions with sector v1/v2). */
enum RestartWord {
  RW_MAGIC, RW_VERSION, RW_PHASE, RW_SCALAR, RW_MODEL, RW_NSITE, RW_NUP, RW_NDOWN, RW_NE,
  RW_FULL_DIM, RW_DIM, RW_RANKS, RW_RANK, RW_LAYOUT, RW_OFFSET, RW_LOCAL_DIM, RW_GROUP,
  RW_SECTOR_COUNT, RW_SECTOR_XOR, RW_SECTOR_SUM, RW_ORDER, RW_HAMILTONIAN,
  RW_METHOD, RW_STATE, RW_STEP, RW_TIME, RW_DATA_XOR, RW_DATA_SUM,
  RW_LOCAL_SITE_MASK, RW_FIXED_FLAGS,
  RW_BASIS_KIND, RW_COMPONENT, RW_TOTAL2SZ, RW_SITE_DIGEST, RW_ORDER_SUM, RW_NSET,
  RW_SCALAR_A, RW_SCALAR_B, RW_SCALAR_C, RW_FLAGS
};

/* Per-set scalar state carried in the header. */
struct RestartVectorInfo {
  uint64_t method;   /* iCalcType */
  uint64_t set;      /* STATE */
  uint64_t step;     /* STEP: next index to execute */
  double time;       /* TIME (TE only) */
  uint64_t nset;     /* NSET */
  uint64_t flags;    /* FLAGS word: bit0 = CG converged terminal checkpoint (Q5) */
  double scalar_a, scalar_b, scalar_c;
};

/* One checkpoint unit (what one manifest update covers). */
struct RestartUnit {
  int nset;                 /* CG: k_exct, TPQ/TE: 1 */
  const uint64_t *sets;     /* set indices */
  int ncomp;                /* 1 or 2 */
};

/* Pure helpers (unit-tested, no MPI). */
void RestartDigestAccumulate(uint64_t offset, const double complex *vector, unsigned long count,
                             uint64_t *xor_hash, uint64_t *sum_hash); /* vector is one-based */
int  RestartIntersect(uint64_t a, uint64_t b, uint64_t offset, uint64_t count,
                      uint64_t *lo, uint64_t *hi);                    /* 1 when hi>lo */

/* Collective. vectors[iset][comp] are one-based, length idim_max+1. */
int RestartWriteUnit(struct BindStruct *X, const struct RestartUnit *unit,
                     double complex *const *const *vectors, const struct RestartVectorInfo *infos, uint64_t flags);
int RestartReadVector(struct BindStruct *X, const struct RestartManifestRow *row, int comp,
                      enum RestartReadPurpose purpose, double complex *vector,
                      struct RestartVectorInfo *info);
/* Classification (§4.6). legacy_pattern: printf pattern with (int set, int rank), TE uses (int rank), or NULL.
   Success verifies every required header/length/identity/step before callers may recover output (V2).
   *legacy_step is a validated next step, agreed across ranks. Payload digest is checked when reading. */
int RestartClassifySet(struct BindStruct *X, const struct RestartManifest *manifest, uint64_t set,
                       const char *legacy_pattern, enum RestartSource *source, int *legacy_step);
/* Import (§4.7): rank 0 inspects "<base>_rank_0.dat" and broadcasts format/ranks.
 * Returns 0 after v3 or sector v1/v2 import, including rank redistribution;
 * returns -2 for the caller's same-rank legacy raw reader, or -1 on collective
 * failure. Source solver fields are provenance; finite scalars and defined
 * FLAGS remain format requirements. */
enum RestartFileFormat { RESTART_FORMAT_LEGACY_RAW = 0, RESTART_FORMAT_SECTOR_V1 = 1, RESTART_FORMAT_SECTOR_V2 = 2, RESTART_FORMAT_V3 = 3 };
int RestartImportVector(struct BindStruct *X, const char *base, double complex *vector,
                        struct RestartVectorInfo *provenance, enum RestartFileFormat *format, int *source_ranks);
uint64_t RestartLocalOffset(const struct BindStruct *X); /* raw: ExscanMPI_li(idim_max); sector: local_offset */
uint64_t RestartHashTag(const char *tag);
uint64_t RestartHashWord(uint64_t hash, uint64_t value);
uint64_t RestartDoubleWord(double value);
double   RestartWordDouble(uint64_t word);
int      RestartPutWord(FILE *fp, uint64_t word);     /* little-endian */
int      RestartGetWord(FILE *fp, uint64_t *word);

#endif
