/* Sector FullDiag Hamiltonian: basis identity, canonical storage and digest. */
#include <float.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include "Common.h"
#include "CalcTime.h"
#include "symmetry_basis.h"
#include "symmetry_hamiltonian_io.h"
#include "symmetry_kondo.h"
#include "symmetry_sector.h"
#include "wrapperMPI.h"
#ifdef HPHI_TEST_SYMMETRY_PROBE
#include "symmetry_fulldiag_probe.h"
#endif
#ifdef MPI
#include <mpi.h>
#endif

#define FNV_OFFSET UINT64_C(14695981039346656037)
/* hphi-sector-matrix-fnv1a64-v1: version word, D, then the lower triangle. */
#define MATRIX_DIGEST_VERSION UINT64_C(1)
#define HERMITIAN_TOLERANCE 1e-12

enum { IDENTITY_WORDS = 15 };

/* Field order of the identity, also the cross-rank comparison order. */
static void pack_identity(const struct SymmetryHamiltonianIdentity *id, uint64_t *words)
{
  const uint64_t fields[IDENTITY_WORDS] = {
    id->model, id->nsite, id->nup, id->ndown, id->ne, id->local_site_mask,
    id->fixed_flags, id->phase, id->full_dim, id->sector_dim, id->group_digest,
    id->sector_count, id->sector_xor, id->sector_sum, id->basis_order
  };
  memcpy(words, fields, sizeof(fields));
}

/* Collective: 0 when every rank holds the root's words. */
static int agree_with_root(const uint64_t *local, int count)
{
  uint64_t root[IDENTITY_WORDS];
  if (count < 1 || count > IDENTITY_WORDS) return -1;
  memcpy(root, local, (size_t)count * sizeof(*root));
#ifdef MPI
  if (MPI_Bcast(root, count, MPI_UINT64_T, 0, MPI_COMM_WORLD) != MPI_SUCCESS) return -1;
#endif
  return SumMPI_i(memcmp(root, local, (size_t)count * sizeof(*root)) != 0) != 0 ? -1 : 0;
}

int BuildSymmetryHamiltonianIdentity(const struct BindStruct *X,
    struct SymmetryHamiltonianIdentity *identity)
{
  struct SymmetryHamiltonianIdentity local;
  struct SymmetrySectorDigest sector;
  uint64_t words[IDENTITY_WORDS];
  int invalid = identity == NULL || X == NULL || X->Sym == NULL ||
                !X->Def.iFlgSymmetryBasis || sizeof(double) != 8 ||
                DBL_MANT_DIG != 53 || DBL_MAX_EXP != 1024;
  if (!invalid)
    invalid = X->Sym->enabled != TRUE ||
              X->Sym->basis_layout != SYMMETRY_BASIS_REPLICATED ||
              X->Sym->basis == NULL || X->Sym->dim == 0 || X->Sym->dim > INT_MAX;
  if (SumMPI_i(invalid) != 0) {
    fprintf(stdoutMPI, "Error: sector Hamiltonian identity requires a replicated, nonempty TransSym basis with binary64 scalars.\n");
    return -1;
  }
  if (ComputeSymmetrySectorDigest(X->Sym, &sector) != 0) {
    fprintf(stdoutMPI, "Error: sector Hamiltonian identity: sector digest failed.\n");
    return -1;
  }
  memset(&local, 0, sizeof(local));
  local.model = (uint64_t)X->Def.iCalcModel;
  local.nsite = X->Def.Nsite;
  if (X->Def.iCalcModel != SpinGC) {
    local.nup = X->Def.Nup;
    local.ndown = X->Def.Ndown;
    local.ne = X->Def.Ne;
  }
  local.phase = 1;
  if (IsSymmetryKondoModel(X->Def.iCalcModel)) {
    struct SymmetryKondoIdentity kondo;
    if (GetSymmetryKondoIdentity(&X->Def, &kondo) != 0) {
      invalid = 1;
    } else {
      local.nup = kondo.nup;
      local.ndown = kondo.ndown;
      local.ne = kondo.ne;
      local.local_site_mask = kondo.local_site_mask;
      local.fixed_flags = kondo.fixed_flags;
      local.phase = kondo.phase;
    }
  }
  local.full_dim = X->Sym->full_dim;
  local.sector_dim = X->Sym->dim;
  if (!invalid && ComputeSymmetryGroupDigest(&X->Def, &local.group_digest) != 0) invalid = 1;
  local.sector_count = sector.count;
  local.sector_xor = sector.xor_hash;
  local.sector_sum = sector.sum_hash;
  /* One global stream, the word encoding of the checkpoint ORDER field. */
  local.basis_order = FNV_OFFSET;
  for (unsigned long i = 1; !invalid && i <= X->Sym->dim; ++i) {
    const struct SymmetryBasisVector *b = SymmetryBasisReplicatedGlobalEntry(X->Sym, i);
    if (b == NULL) { invalid = 1; break; }
    local.basis_order = SymmetryHashInteger(local.basis_order, i, 8);
    local.basis_order = SymmetryHashInteger(local.basis_order, b->rep_state, 8);
    local.basis_order = SymmetryHashInteger(local.basis_order, b->orbit_size, 8);
    local.basis_order = SymmetryHashInteger(local.basis_order, b->stabilizer_size, 8);
  }
  if (SumMPI_i(invalid) != 0) {
    fprintf(stdoutMPI, "Error: sector Hamiltonian identity: group or basis metadata is invalid.\n");
    return -1;
  }
  pack_identity(&local, words);
  if (agree_with_root(words, IDENTITY_WORDS) != 0) {
    fprintf(stdoutMPI, "Error: sector Hamiltonian identity differs across MPI ranks.\n");
    return -1;
  }
  *identity = local;
  return 0;
}

#ifdef MPI
/* Storage adapter: replicated Ham[row][col] or the ELPA owned-column panel,
 * Ham_local[(col-HamColBegin)*HamPanelLd + row-1]. Indices are 1-based.
 * Only the panel path, which needs MPI, uses it. */
static int owns_column(unsigned long col)
{
  return !iHamPanelActive ||
         ((long)col >= HamColBegin && (long)col <= HamColEnd);
}

static double complex read_owned(unsigned long row, unsigned long col)
{
  if (!iHamPanelActive) return Ham[row][col];
  return Ham_local[(size_t)((long)col - HamColBegin) * (size_t)HamPanelLd + (row - 1)];
}

static void set_owned(unsigned long row, unsigned long col, double complex value)
{
  if (!iHamPanelActive) {
    Ham[row][col] = value;
    return;
  }
  Ham_local[(size_t)((long)col - HamColBegin) * (size_t)HamPanelLd + (row - 1)] = value;
}

/* Panel only: buffer[0..dim-1] receives column col from its owner
 * (col-1)/NC, NC = ceil(dim/nproc), as setmem_large() distributes it. */
static int broadcast_column(unsigned long dim, unsigned long col,
                            double complex *buffer)
{
  unsigned long chunk = (dim + (unsigned long)nproc - 1) / (unsigned long)nproc;
  int owner = (int)((col - 1) / chunk);
  if (owner == myrank)
    memcpy(buffer, &Ham_local[(size_t)((long)col - HamColBegin) * (size_t)HamPanelLd],
           dim * sizeof(*buffer));
  return MPI_Bcast(buffer, (int)dim, MPI_DOUBLE_COMPLEX, owner, MPI_COMM_WORLD)
         == MPI_SUCCESS ? 0 : -1;
}
#endif

static int validate_storage(unsigned long dim)
{
  if (!iHamPanelActive) {
    if (Ham == NULL) return -1;
    for (unsigned long row = 0; row <= dim; ++row)
      if (Ham[row] == NULL) return -1;
    return 0;
  }
  long chunk = ((long)dim + nproc - 1) / nproc;
  long first = (long)myrank * chunk + 1, last = ((long)myrank + 1) * chunk;
  if (last > (long)dim) last = (long)dim;
  if (first > (long)dim) { first = 1; last = 0; }
  if (Ham_local == NULL || HamPanelLd != (long)dim ||
      HamColBegin != first || HamColEnd != last) return -1;
  if (last >= first &&
      (size_t)(last - first + 1) > SIZE_MAX / dim / sizeof(*Ham_local)) return -1;
  return 0;
}

static int finite_element(double complex value)
{
  return isfinite(creal(value)) && isfinite(cimag(value));
}

/* |a - conj(b)| <= 1e-12 max(1, |a|, |b|); for a diagonal a == b this
 * bounds the imaginary part with the same scale. */
static int hermitian_pair(double complex lower, double complex upper)
{
  double scale = fmax(1.0, fmax(cabs(lower), cabs(upper)));
  return cabs(lower - conj(upper)) <= HERMITIAN_TOLERANCE * scale;
}

static uint64_t hash_component(uint64_t hash, double value)
{
  uint64_t bits;
  if (value == 0.0) value = 0.0; /* signed zero hashes as +0 */
  memcpy(&bits, &value, sizeof(bits));
  return SymmetryHashInteger(hash, bits, 8);
}

static uint64_t hash_element(uint64_t hash, double complex value)
{
  return hash_component(hash_component(hash, creal(value)), cimag(value));
}

static int report_matrix(int nonfinite, int nonhermitian)
{
  if (nonfinite) {
    fprintf(stdoutMPI, "Error: sector Hamiltonian has a non-finite element; FullDiag stops before diagonalization.\n");
    return -1;
  }
  if (nonhermitian) {
    fprintf(stdoutMPI, "Error: sector Hamiltonian is not Hermitian within 1e-12*max(1,|a|,|b|); FullDiag stops before diagonalization.\n");
    return -1;
  }
  return 0;
}

static void note_element(const char *what, unsigned long row, unsigned long col)
{
  fprintf(stderr, "Error: rank %d sector Hamiltonian element (%lu,%lu) is %s.\n",
          myrank, row, col, what);
}

/* Every rank owns the full matrix: check, canonicalize and hash locally. */
static int prepare_replicated(unsigned long dim, uint64_t *digest)
{
  int nonfinite = 0, nonhermitian = 0;
  for (unsigned long col = 1; col <= dim && !nonfinite && !nonhermitian; ++col) {
    for (unsigned long row = col; row <= dim; ++row) {
      double complex lower = Ham[row][col], upper = Ham[col][row];
      if (!finite_element(lower) || !finite_element(upper)) {
        note_element("not finite", row, col);
        nonfinite = 1;
        break;
      }
      if (!hermitian_pair(lower, upper)) {
        note_element("not Hermitian", row, col);
        nonhermitian = 1;
        break;
      }
    }
  }
  nonfinite = SumMPI_i(nonfinite) != 0;
  nonhermitian = SumMPI_i(nonhermitian) != 0;
  if (report_matrix(nonfinite, nonhermitian) != 0) return -1;
  uint64_t hash = SymmetryHashInteger(FNV_OFFSET, MATRIX_DIGEST_VERSION, 8);
  hash = SymmetryHashInteger(hash, dim, 8);
  for (unsigned long col = 1; col <= dim; ++col) {
    double complex diagonal = creal(Ham[col][col]);
    Ham[col][col] = diagonal;
    hash = hash_element(hash, diagonal);
    for (unsigned long row = col + 1; row <= dim; ++row) {
      Ham[col][row] = conj(Ham[row][col]);
      hash = hash_element(hash, Ham[row][col]);
    }
  }
  *digest = hash;
  return 0;
}

#ifdef MPI
/* ELPA panel: column c is sent once from its owner. The owner of column r
 * compares its H[c,r] with the received H[r,c] before writing the conjugate,
 * and every rank hashes the identical received lower triangle. */
static int prepare_panel(unsigned long dim, uint64_t *digest)
{
  double complex *column = malloc(dim * sizeof(*column));
  uint64_t hash = SymmetryHashInteger(FNV_OFFSET, MATRIX_DIGEST_VERSION, 8);
#ifdef HPHI_TEST_SYMMETRY_PROBE
  if (SymmetryFullDiagTestFail("matrix_alloc")) { free(column); column = NULL; }
#endif
  if (SumMPI_i(column == NULL) != 0) {
    free(column);
    fprintf(stdoutMPI, "Error: sector Hamiltonian column buffer allocation failed.\n");
    return -1;
  }
  hash = SymmetryHashInteger(hash, dim, 8);
  for (unsigned long col = 1; col <= dim; ++col) {
    int flags = 0, any = 0;
    if (SumMPI_i(broadcast_column(dim, col, column) != 0) != 0) {
      free(column);
      fprintf(stdoutMPI, "Error: sector Hamiltonian column broadcast failed.\n");
      return -1;
    }
    for (unsigned long row = col; row <= dim && !flags; ++row) {
      double complex lower = column[row - 1];
      if (!finite_element(lower)) {
        /* Every rank received this element; report it once. */
        flags |= 1;
        if (myrank == 0) note_element("not finite", row, col);
      } else if (owns_column(row)) {
        double complex upper = read_owned(col, row);
        if (!finite_element(upper)) flags |= 1;
        else if (!hermitian_pair(lower, upper)) flags |= 2;
        if (flags) note_element(flags & 1 ? "not finite" : "not Hermitian", row, col);
      }
    }
    if (MPI_Allreduce(&flags, &any, 1, MPI_INT, MPI_BOR, MPI_COMM_WORLD) != MPI_SUCCESS)
      any = 1;
    if (report_matrix(any & 1, any & 2) != 0) {
      free(column);
      return -1;
    }
    column[col - 1] = creal(column[col - 1]);
    if (owns_column(col)) set_owned(col, col, column[col - 1]);
    for (unsigned long row = col; row <= dim; ++row) {
      if (row > col && owns_column(row)) set_owned(col, row, conj(column[row - 1]));
      hash = hash_element(hash, column[row - 1]);
    }
  }
  free(column);
  *digest = hash;
  return 0;
}
#endif

int PrepareSymmetryHamiltonian(const struct BindStruct *X, int from_file,
    struct SymmetryHamiltonianInfo *info)
{
  struct SymmetryHamiltonianInfo local;
  uint64_t digest = 0;
  int status = -1;
  StartTimer(5600);
  memset(&local, 0, sizeof(local));
  if (SumMPI_i(info == NULL || (from_file != 0 && from_file != 1)) != 0) {
    fprintf(stdoutMPI, "Error: invalid sector Hamiltonian preparation arguments.\n");
    goto done;
  }
  if (BuildSymmetryHamiltonianIdentity(X, &local.identity) != 0) goto done;
  if (SumMPI_i(validate_storage(X->Sym->dim) != 0) != 0) {
    fprintf(stdoutMPI, "Error: sector Hamiltonian storage does not match the FullDiag sector dimension.\n");
    goto done;
  }
#ifdef MPI
  if (iHamPanelActive) {
    if (prepare_panel(X->Sym->dim, &digest) != 0) goto done;
  } else
#endif
  if (prepare_replicated(X->Sym->dim, &digest) != 0) goto done;
  if (agree_with_root(&digest, 1) != 0) {
    fprintf(stdoutMPI, "Error: sector Hamiltonian digest differs across MPI ranks.\n");
    goto done;
  }
  local.matrix_digest = digest;
  local.from_file = from_file;
  *info = local;
  status = 0;
done:
  StopTimer(5600);
  return status;
}
