/* Sector FullDiag Hamiltonian: basis identity, canonical storage and digest. */
#include <float.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
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

/* Off diagonal: |a - conj(b)| <= tol max(1, |a|, |b|).
 * Diagonal: |Im(a)| <= tol max(1, |a|), without the factor two from a-conj(a).
 * Inputs have finite components, but their magnitudes/difference may overflow.
 * Scale by a power of two first: norms/difference stay bounded by sqrt(8),
 * without rounding each operand through division before subtracting. */
static int hermitian_pair(double complex lower, double complex upper, int diagonal)
{
  double scale = fmax(1.0, fmax(fmax(fabs(creal(lower)), fabs(cimag(lower))),
                               fmax(fabs(creal(upper)), fabs(cimag(upper)))));
  int exponent = ilogb(scale) + 1;
  double lr = scalbn(creal(lower), -exponent), li = scalbn(cimag(lower), -exponent);
  double ur = scalbn(creal(upper), -exponent), ui = scalbn(cimag(upper), -exponent);
  double bound = HERMITIAN_TOLERANCE *
      fmax(scalbn(1.0, -exponent), fmax(hypot(lr, li), hypot(ur, ui)));
  return (diagonal ? fabs(li) : hypot(lr - ur, li + ui)) <= bound;
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
      if (!hermitian_pair(lower, upper, row == col)) {
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
        else if (!hermitian_pair(lower, upper, row == col)) flags |= 2;
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

/* ---- Matrix Market sector Hamiltonian file (<head>_ham_sector.mtx) ---- */

#define MATRIX_HEADER "%%MatrixMarket matrix coordinate complex hermitian"
#define MATRIX_FORMAT "sector-hamiltonian-v1"
#define MATRIX_SCALAR "complex-binary64"
#define METADATA_PREFIX "% HPHI "
enum { RECORD_CHUNK = 4096, LINE_LENGTH = 512, FIRST_IDENTITY_KEY = 2 };

/* format, scalar, then the identity words in pack_identity order. */
static const char *const required_keys[] = {
  "format", "scalar", "model", "nsite", "nup", "ndown", "ne",
  "local_site_mask", "fixed_flags", "phase", "full_dim", "sector_dim",
  "group_digest", "sector_count", "sector_xor", "sector_sum", "basis_order"
};
enum { REQUIRED_KEYS = sizeof(required_keys) / sizeof(required_keys[0]) };
typedef char matrix_keys_follow_identity
    [REQUIRED_KEYS == FIRST_IDENTITY_KEY + IDENTITY_WORDS ? 1 : -1];

static const struct { int model; const char *name; } matrix_models[] = {
  { Spin, "Spin" }, { SpinGC, "SpinGC" }, { SpinlessFermion, "SpinlessFermion" },
  { Hubbard, "Hubbard" }, { tJ, "tJ" }, { Kondo, "Kondo" },
  { KondoNConserved, "KondoNConserved" }, { KondoGC, "KondoGC" }
};

static const char *model_name(uint64_t model)
{
  for (size_t i = 0; i < sizeof(matrix_models) / sizeof(matrix_models[0]); ++i)
    if ((uint64_t)matrix_models[i].model == model) return matrix_models[i].name;
  return NULL;
}

/* Masks and digests are 16 hexadecimal digits, the other words decimal. */
static int hex_field(int field)
{
  return field == 5 || field == 10 || field == 12 || field == 13 || field == 14;
}

static void format_identity_value(int field, uint64_t value, char *text, size_t size)
{
  const char *name = field == 0 ? model_name(value) : NULL;
  if (name != NULL) snprintf(text, size, "%s", name);
  else if (hex_field(field)) snprintf(text, size, "%016" PRIx64, value);
  else snprintf(text, size, "%" PRIu64, value);
}

/* Every rank builds the same name, so a too long path fails everywhere. */
static int matrix_path(const struct DefineList *def, int temporary, char *path, size_t size)
{
  int length = snprintf(path, size, "%s%s_ham_sector.mtx%s", cParentOutputFolder,
                        def->CDataFileHead, temporary ? ".part" : "");
  return length < 0 || (size_t)length >= size ? -1 : 0;
}

/* Digits only: no sign, space, prefix or overflow. */
static int parse_decimal(const char *text, uint64_t *value)
{
  uint64_t result = 0;
  if (*text == '\0') return -1;
  for (; *text != '\0'; ++text) {
    unsigned int digit;
    if (*text < '0' || *text > '9') return -1;
    digit = (unsigned int)(*text - '0');
    if (result > (UINT64_MAX - digit) / 10) return -1;
    result = result * 10 + digit;
  }
  *value = result;
  return 0;
}

static int parse_hex16(const char *text, uint64_t *value)
{
  uint64_t result = 0;
  if (strlen(text) != 16) return -1;
  for (; *text != '\0'; ++text) {
    int digit;
    if (*text >= '0' && *text <= '9') digit = *text - '0';
    else if (*text >= 'a' && *text <= 'f') digit = *text - 'a' + 10;
    else if (*text >= 'A' && *text <= 'F') digit = *text - 'A' + 10;
    else return -1;
    result = (result << 4) | (uint64_t)digit;
  }
  *value = result;
  return 0;
}

static int parse_double(const char *text, double *value)
{
  char *end = NULL;
  *value = strtod(text, &end);
  return end == text || *end != '\0' ? -1 : 0;
}

/* Splits on spaces/tabs; returns max+1 when there are more tokens. */
static int split_tokens(char *text, char **tokens, int max)
{
  int count = 0;
  for (;;) {
    while (*text == ' ' || *text == '\t') *text++ = '\0';
    if (*text == '\0') return count;
    if (count == max) return max + 1;
    tokens[count++] = text;
    while (*text != '\0' && *text != ' ' && *text != '\t') ++text;
  }
}

static int blank(const char *text)
{
  for (; *text != '\0'; ++text)
    if (*text != ' ' && *text != '\t') return 0;
  return 1;
}

/* Root-only reader state. */
struct MatrixReader {
  FILE *fp;
  const char *path;
  unsigned long line;
  char text[LINE_LENGTH];
};

/* 1 for a line without its LF/CRLF, 0 at the end of the file, -1 for a
 * read error, embedded NUL or a line that does not fit the buffer. */
static int read_line(struct MatrixReader *reader)
{
  size_t length = 0;
  int byte = fgetc(reader->fp);
  if (byte == EOF)
    return ferror(reader->fp) ? -1 : 0;
  reader->line++;
  do {
    /* Count physical bytes, including the line end, before treating them as
     * a C string. A NUL in the file must not hide records or extra tokens. */
    if (byte == '\0' || length == sizeof(reader->text) - 1) return -1;
    reader->text[length++] = (char)byte;
    if (byte == '\n') break;
  } while ((byte = fgetc(reader->fp)) != EOF);
  if (ferror(reader->fp)) return -1;
  if (length > 0 && reader->text[length - 1] == '\n') --length;
  if (length > 0 && reader->text[length - 1] == '\r') --length;
  reader->text[length] = '\0';
  return 1;
}

static void reader_error(const struct MatrixReader *reader, const char *format, ...)
{
  va_list args;
  fprintf(stdoutMPI, "Error: sector Hamiltonian file %s line %lu: ", reader->path, reader->line);
  va_start(args, format);
  vfprintf(stdoutMPI, format, args);
  va_end(args);
  fprintf(stdoutMPI, ".\n");
}

static void unreadable_line(const struct MatrixReader *reader)
{
  reader_error(reader, "line contains NUL, is too long (at most %d bytes including its line end) or unreadable",
               LINE_LENGTH - 1);
}

/* "% HPHI key=value": each key once, values in their documented syntax. */
static int parse_metadata(const struct MatrixReader *reader, char *entry,
                          uint64_t *words, unsigned long *seen)
{
  char *value = strchr(entry, '=');
  int key, invalid;
  if (value == NULL) {
    reader_error(reader, "malformed metadata line (expected %skey=value)", METADATA_PREFIX);
    return -1;
  }
  *value++ = '\0';
  for (key = 0; key < REQUIRED_KEYS && strcmp(entry, required_keys[key]) != 0; ++key) ;
  if (key == REQUIRED_KEYS) {
    reader_error(reader, "unknown metadata key %s", entry);
    return -1;
  }
  if (*seen & (1UL << key)) {
    reader_error(reader, "duplicate metadata key %s", entry);
    return -1;
  }
  *seen |= 1UL << key;
  if (key == 0 || key == 1) {
    const char *expected = key == 0 ? MATRIX_FORMAT : MATRIX_SCALAR;
    if (strcmp(value, expected) == 0) return 0;
    reader_error(reader, "unsupported %s %s (expected %s)", entry, value, expected);
    return -1;
  }
  int field = key - FIRST_IDENTITY_KEY;
  if (field == 0) {
    invalid = 1;
    for (size_t i = 0; i < sizeof(matrix_models) / sizeof(matrix_models[0]); ++i)
      if (strcmp(value, matrix_models[i].name) == 0) {
        words[0] = (uint64_t)matrix_models[i].model;
        invalid = 0;
      }
  } else {
    invalid = hex_field(field) ? parse_hex16(value, &words[field])
                               : parse_decimal(value, &words[field]);
  }
  if (invalid) reader_error(reader, "invalid metadata value for %s: %s", entry, value);
  return invalid ? -1 : 0;
}

/* Root only: the first line, the HPHI keys and the dimension line. */
static int read_matrix_header(struct MatrixReader *reader, uint64_t *words, uint64_t *dims)
{
  char *tokens[3];
  unsigned long seen = 0;
  int status = read_line(reader);
  if (status != 1 || strcmp(reader->text, MATRIX_HEADER) != 0) {
    reader_error(reader, "first line must be %s", MATRIX_HEADER);
    return -1;
  }
  for (;;) {
    status = read_line(reader);
    if (status == 0) {
      reader_error(reader, "file ends before the dimension line");
      return -1;
    }
    if (status < 0) {
      unreadable_line(reader);
      return -1;
    }
    if (blank(reader->text)) continue;
    if (reader->text[0] != '%') break;
    if (strncmp(reader->text, METADATA_PREFIX, strlen(METADATA_PREFIX)) == 0 &&
        parse_metadata(reader, reader->text + strlen(METADATA_PREFIX), words, &seen) != 0)
      return -1;
  }
  for (int key = 0; key < REQUIRED_KEYS; ++key) {
    if (seen & (1UL << key)) continue;
    fprintf(stdoutMPI, "Error: sector Hamiltonian file %s: missing metadata key %s.\n",
            reader->path, required_keys[key]);
    return -1;
  }
  if (split_tokens(reader->text, tokens, 3) != 3 || parse_decimal(tokens[0], &dims[0]) != 0 ||
      parse_decimal(tokens[1], &dims[1]) != 0 || parse_decimal(tokens[2], &dims[2]) != 0) {
    reader_error(reader, "dimension line must hold three nonnegative integers (rows columns records)");
    return -1;
  }
  return 0;
}

/* Root only: the basis identity and the dimension line of the file. */
static int check_matrix_header(const struct MatrixReader *reader, const uint64_t *words,
                               const uint64_t *expected, const uint64_t *dims)
{
  uint64_t dim = expected[9];
  for (int field = 0; field < IDENTITY_WORDS; ++field) {
    char found[32], wanted[32];
    if (words[field] == expected[field]) continue;
    format_identity_value(field, words[field], found, sizeof(found));
    format_identity_value(field, expected[field], wanted, sizeof(wanted));
    fprintf(stdoutMPI, "Error: sector Hamiltonian file %s: metadata %s does not match this sector "
            "(file=%s, expected=%s).\n", reader->path,
            required_keys[FIRST_IDENTITY_KEY + field], found, wanted);
    return -1;
  }
  if (dims[0] != dim || dims[1] != dim) {
    reader_error(reader, "matrix %" PRIu64 " x %" PRIu64 " does not match sector dimension %" PRIu64,
                 dims[0], dims[1], dim);
    return -1;
  }
  if (dims[2] > dim * (dim + 1) / 2) {
    reader_error(reader, "%" PRIu64 " records declared, more records than the lower triangle "
                 "of dimension %" PRIu64 " holds", dims[2], dim);
    return -1;
  }
  return 0;
}

/* Root only: one "row column real imaginary" record of the lower triangle. */
static int read_record(struct MatrixReader *reader, uint64_t dim, uint64_t read, uint64_t nnz,
                       uint64_t *index, double *value)
{
  char *tokens[4];
  int status = read_line(reader);
  if (status == 0) {
    fprintf(stdoutMPI, "Error: sector Hamiltonian file %s: fewer records than declared "
            "(%" PRIu64 " of %" PRIu64 ").\n", reader->path, read, nnz);
    return -1;
  }
  if (status < 0) {
    unreadable_line(reader);
    return -1;
  }
  if (split_tokens(reader->text, tokens, 4) != 4 || parse_decimal(tokens[0], &index[0]) != 0 ||
      parse_decimal(tokens[1], &index[1]) != 0 || parse_double(tokens[2], &value[0]) != 0 ||
      parse_double(tokens[3], &value[1]) != 0) {
    reader_error(reader, "malformed record (expected: row column real imaginary)");
    return -1;
  }
  if (index[0] < 1 || index[0] > dim || index[1] < 1 || index[1] > dim) {
    reader_error(reader, "index out of range (%" PRIu64 ",%" PRIu64 ")", index[0], index[1]);
    return -1;
  }
  if (index[1] > index[0]) {
    reader_error(reader, "upper-triangle record (%" PRIu64 ",%" PRIu64 "); write row >= column only",
                 index[0], index[1]);
    return -1;
  }
  if (!isfinite(value[0]) || !isfinite(value[1])) {
    reader_error(reader, "non-finite coefficient");
    return -1;
  }
  if (index[0] == index[1] && value[1] != 0.0) {
    reader_error(reader, "nonzero diagonal imaginary part");
    return -1;
  }
  return 0;
}

/* Root only: blank lines and general comments may follow the records. */
static int read_trailer(struct MatrixReader *reader, uint64_t nnz)
{
  int status;
  while ((status = read_line(reader)) == 1) {
    if (blank(reader->text) ||
        (reader->text[0] == '%' &&
         strncmp(reader->text, METADATA_PREFIX, strlen(METADATA_PREFIX)) != 0)) continue;
    reader_error(reader, "unexpected content after the %" PRIu64 " declared records", nnz);
    return -1;
  }
  if (status < 0) {
    unreadable_line(reader);
    return -1;
  }
  return 0;
}

/* Zero the stored H: the file lists only its nonzero lower triangle. */
static void clear_storage(unsigned long dim)
{
#ifdef MPI
  if (iHamPanelActive) {
    if (HamColEnd >= HamColBegin)
      memset(Ham_local, 0, (size_t)dim * (size_t)(HamColEnd - HamColBegin + 1) * sizeof(*Ham_local));
    return;
  }
#endif
  for (unsigned long row = 0; row <= dim; ++row)
    memset(Ham[row], 0, (dim + 1) * sizeof(**Ham));
}

static void store_record(unsigned long row, unsigned long col, double complex value)
{
#ifdef MPI
  if (iHamPanelActive) {
    if (owns_column(col)) set_owned(row, col, value);
    if (row != col && owns_column(row)) set_owned(col, row, conj(value));
    return;
  }
#endif
  Ham[col][row] = conj(value);
  Ham[row][col] = value;
}

/* Collective: the root's status and count, then (on success) the records. */
static int share_chunk(int *share, uint64_t *index, double *value)
{
#ifdef MPI
  if (MPI_Bcast(share, 2, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS) return -1;
  if (share[0] != 0 || share[1] == 0) return 0;
  if (MPI_Bcast(index, 2 * share[1], MPI_UINT64_T, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(value, 2 * share[1], MPI_DOUBLE, 0, MPI_COMM_WORLD) != MPI_SUCCESS) return -1;
#else
  (void)share; (void)index; (void)value;
#endif
  return 0;
}

/* Every rank stores the records it owns. Duplicates are detected by the
 * unique column owner (col-1)/NC with a D x NC bitset, for every storage. */
static int apply_chunk(unsigned long dim, unsigned long chunk_columns, unsigned char *seen,
                       const uint64_t *index, const double *value, int count)
{
  unsigned long first = (unsigned long)myrank * chunk_columns + 1;
  int duplicate = 0;
  for (int i = 0; i < count; ++i) {
    unsigned long row = (unsigned long)index[2 * i], col = (unsigned long)index[2 * i + 1];
    if ((col - 1) / chunk_columns == (unsigned long)myrank) {
      uint64_t bit = (uint64_t)(col - first) * dim + (row - 1);
      unsigned char mask = (unsigned char)(1u << (bit & 7));
      if (seen[bit >> 3] & mask) {
        if (!duplicate)
          fprintf(stderr, "Error: rank %d sector Hamiltonian record (%lu,%lu) is duplicated.\n",
                  myrank, row, col);
        duplicate = 1;
        continue;
      }
      seen[bit >> 3] |= mask;
    }
    store_record(row, col, value[2 * i] + I * value[2 * i + 1]);
  }
  return duplicate;
}

int InputSymmetryHamiltonian(const struct BindStruct *X,
    struct SymmetryHamiltonianInfo *info)
{
  struct SymmetryHamiltonianIdentity identity;
  struct MatrixReader reader;
  char path[D_FileNameMax];
  uint64_t words[IDENTITY_WORDS] = {0}, expected[IDENTITY_WORDS], dims[3] = {0, 0, 0};
  uint64_t *index = NULL;
  double *value = NULL;
  unsigned char *seen = NULL;
  unsigned long dim, chunk_columns;
  int status = -1, failed;
  StartTimer(5800);
  memset(&reader, 0, sizeof(reader));
  failed = X == NULL || info == NULL || matrix_path(&X->Def, 0, path, sizeof(path)) != 0;
  if (SumMPI_i(failed) != 0) {
    fprintf(stdoutMPI, "Error: invalid sector Hamiltonian input arguments or path.\n");
    goto done;
  }
  if (BuildSymmetryHamiltonianIdentity(X, &identity) != 0) goto done;
  dim = X->Sym->dim;
  if (SumMPI_i(validate_storage(dim) != 0) != 0) {
    fprintf(stdoutMPI, "Error: sector Hamiltonian storage does not match the FullDiag sector dimension.\n");
    goto done;
  }
  chunk_columns = (dim + (unsigned long)nproc - 1) / (unsigned long)nproc;
  {
    uint64_t bits = (uint64_t)dim * chunk_columns;
    size_t bytes = (size_t)(bits / 8 + 1);
    seen = bits / 8 + 1 <= SIZE_MAX ? calloc(bytes, 1) : NULL;
  }
  index = malloc(2 * RECORD_CHUNK * sizeof(*index));
  value = malloc(2 * RECORD_CHUNK * sizeof(*value));
  if (SumMPI_i(seen == NULL || index == NULL || value == NULL) != 0) {
    fprintf(stdoutMPI, "Error: sector Hamiltonian input buffer allocation failed.\n");
    goto done;
  }
  clear_storage(dim);
  pack_identity(&identity, expected);
  failed = 0;
  if (myrank == 0) {
    reader.path = path;
    reader.fp = fopen(path, "r");
    if (reader.fp == NULL) {
      fprintf(stdoutMPI, "Error: could not open sector Hamiltonian input %s.\n", path);
      failed = 1;
    } else {
      failed = read_matrix_header(&reader, words, dims) != 0 ||
               check_matrix_header(&reader, words, expected, dims) != 0;
    }
  }
  if (SumMPI_i(failed) != 0) goto done;
#ifdef MPI
  if (MPI_Bcast(&dims[2], 1, MPI_UINT64_T, 0, MPI_COMM_WORLD) != MPI_SUCCESS) failed = 1;
  if (SumMPI_i(failed) != 0) goto done;
#endif
  for (uint64_t read = 0; read < dims[2];) {
    int share[2] = {0, 0}; /* root status, record count */
    if (myrank == 0) {
      uint64_t left = dims[2] - read;
      share[1] = left < RECORD_CHUNK ? (int)left : RECORD_CHUNK;
      for (int i = 0; i < share[1] && share[0] == 0; ++i)
        share[0] = read_record(&reader, dim, read + (uint64_t)i, dims[2],
                               &index[2 * i], &value[2 * i]) != 0;
    }
    failed = share_chunk(share, index, value) != 0;
    if (SumMPI_i(failed || share[0] != 0) != 0) goto done;
    if (SumMPI_i(apply_chunk(dim, chunk_columns, seen, index, value, share[1])) != 0) {
      fprintf(stdoutMPI, "Error: sector Hamiltonian file %s has a duplicated record.\n", path);
      goto done;
    }
    read += (uint64_t)share[1];
  }
  failed = 0;
  if (myrank == 0) {
    failed = read_trailer(&reader, dims[2]) != 0;
    if (fclose(reader.fp) != 0) failed = 1;
    reader.fp = NULL;
  }
  if (SumMPI_i(failed) != 0) goto done;
  /* The lower triangle and its conjugate are stored: canonicalize and hash. */
  if (PrepareSymmetryHamiltonian(X, 1, info) != 0) goto done;
  status = 0;
done:
  if (reader.fp != NULL) fclose(reader.fp);
  free(seen);
  free(index);
  free(value);
  StopTimer(5800);
  return status;
}

/* Collective. buffer[row-1] holds H(row,col) for row >= col on the root:
 * the owner sends a panel column, replicated storage reads the root's H. */
static int load_column(unsigned long dim, unsigned long col, double complex *buffer)
{
#ifdef MPI
  if (iHamPanelActive)
    return SumMPI_i(broadcast_column(dim, col, buffer) != 0) != 0 ? -1 : 0;
#endif
  if (myrank == 0)
    for (unsigned long row = col; row <= dim; ++row) buffer[row - 1] = Ham[row][col];
  return 0;
}

static void canonical_parts(double complex value, double *re, double *im)
{
  *re = creal(value);
  *im = cimag(value);
  if (*re == 0.0) *re = 0.0; /* signed zero is written as +0 */
  if (*im == 0.0) *im = 0.0;
}

static int write_matrix_header(FILE *fp, const struct SymmetryHamiltonianIdentity *identity,
                               unsigned long dim, uint64_t nnz)
{
  uint64_t words[IDENTITY_WORDS];
  char text[32];
  int failed = fprintf(fp, "%s\n%sformat=%s\n%sscalar=%s\n", MATRIX_HEADER,
                       METADATA_PREFIX, MATRIX_FORMAT, METADATA_PREFIX, MATRIX_SCALAR) < 0;
  pack_identity(identity, words);
  for (int field = 0; !failed && field < IDENTITY_WORDS; ++field) {
    format_identity_value(field, words[field], text, sizeof(text));
    failed = fprintf(fp, "%s%s=%s\n", METADATA_PREFIX,
                     required_keys[FIRST_IDENTITY_KEY + field], text) < 0;
  }
  if (!failed) failed = fprintf(fp, "%lu %lu %" PRIu64 "\n", dim, dim, nnz) < 0;
  return failed ? -1 : 0;
}

static int write_column(FILE *fp, unsigned long dim, unsigned long col,
                        const double complex *column)
{
  for (unsigned long row = col; row <= dim; ++row) {
    double re, im;
    canonical_parts(column[row - 1], &re, &im);
    if ((re != 0.0 || im != 0.0) &&
        fprintf(fp, "%lu %lu %.17g %.17g\n", row, col, re, im) < 0) return -1;
  }
  return 0;
}

int OutputSymmetryHamiltonian(const struct BindStruct *X,
    const struct SymmetryHamiltonianInfo *info)
{
  char path[D_FileNameMax], temporary[D_FileNameMax];
  FILE *fp = NULL;
  double complex *column = NULL;
  uint64_t nnz = 0;
  unsigned long dim = 0, col;
  int status = -1, failed, created = 0;
  StartTimer(5900);
  failed = X == NULL || info == NULL || X->Sym == NULL ||
           matrix_path(&X->Def, 0, path, sizeof(path)) != 0 ||
           matrix_path(&X->Def, 1, temporary, sizeof(temporary)) != 0;
  if (!failed) {
    dim = X->Sym->dim;
    failed = info->identity.sector_dim != dim || model_name(info->identity.model) == NULL ||
             validate_storage(dim) != 0;
  }
  if (SumMPI_i(failed) != 0) {
    fprintf(stdoutMPI, "Error: sector Hamiltonian output needs the prepared sector matrix and a valid path.\n");
    goto done;
  }
  column = malloc(dim * sizeof(*column));
  if (SumMPI_i(column == NULL) != 0) {
    fprintf(stdoutMPI, "Error: sector Hamiltonian output column buffer allocation failed.\n");
    goto done;
  }
  /* The open is agreed before any column is sent. */
  failed = 0;
  if (myrank == 0) {
    fp = fopen(temporary, "w");
    created = fp != NULL;
    failed = fp == NULL;
    if (failed) fprintf(stdoutMPI, "Error: could not open sector Hamiltonian output %s.\n", temporary);
  }
  if (SumMPI_i(failed) != 0) goto done;
  /* Pass 1: the root counts the nonzero lower-triangle elements. */
  for (col = 1; col <= dim && !failed; ++col) {
    if (load_column(dim, col, column) != 0) {
      fprintf(stdoutMPI, "Error: sector Hamiltonian column broadcast failed.\n");
      failed = 1;
      break;
    }
    if (myrank != 0) continue;
    for (unsigned long row = col; row <= dim; ++row) {
      double re, im;
      canonical_parts(column[row - 1], &re, &im);
      nnz += re != 0.0 || im != 0.0;
    }
  }
  if (myrank == 0 && !failed) failed = write_matrix_header(fp, &info->identity, dim, nnz) != 0;
  /* Pass 2: column-major lower triangle. Panel ranks agree before it and
   * after each column, so a failed root write stops the broadcasts; only
   * the root writes replicated storage. */
  if (iHamPanelActive) failed = SumMPI_i(failed) != 0;
  for (col = 1; col <= dim && !failed; ++col) {
    if (load_column(dim, col, column) != 0) {
      failed = 1;
      break;
    }
    if (myrank == 0) {
#ifdef HPHI_TEST_SYMMETRY_PROBE
      if (SymmetryFullDiagTestFail("matrix_write")) {
        fputs("1 ", fp); /* a real short write, reported like a failed fprintf */
        failed = 1;
      } else
#endif
      failed = write_column(fp, dim, col, column) != 0;
    }
    if (iHamPanelActive ? SumMPI_i(failed) != 0 : failed) break;
  }
  if (myrank == 0 && fp != NULL) {
    if (fflush(fp) != 0 || ferror(fp)) failed = 1;
    if (fclose(fp) != 0) failed = 1;
    fp = NULL;
#ifdef HPHI_TEST_SYMMETRY_PROBE
    if (SymmetryFullDiagTestFail("matrix_close")) failed = 1;
#endif
  }
  if (SumMPI_i(failed) != 0) {
    if (created) remove(temporary);
    fprintf(stdoutMPI, "Error: sector Hamiltonian output write failed for %s.\n", temporary);
    goto done;
  }
  if (myrank == 0) {
#ifdef HPHI_TEST_SYMMETRY_PROBE
    if (SymmetryFullDiagTestFail("matrix_rename")) failed = 1;
    else
#endif
    failed = rename(temporary, path) != 0;
    if (failed) remove(temporary);
  }
  if (SumMPI_i(failed) != 0) {
    fprintf(stdoutMPI, "Error: sector Hamiltonian output publish failed for %s.\n", path);
    goto done;
  }
  status = 0;
done:
  if (fp != NULL) {
    fclose(fp);
    remove(temporary);
  }
  free(column);
  StopTimer(5900);
  return status;
}
