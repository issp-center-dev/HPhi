#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include "symmetry_basis.h"
#include "symmetry_basis_io.h"
#include "symmetry_terms.h"
#include "readdef.h"
#include "struct.h"
#include "wrapperMPI.h"
#include "global.h"

static int read_header_count(FILE *fp, const char *defname, unsigned int *ntrans)
{
  char line[D_CharTmpReadDef + D_CharKWDMAX];
  char key[D_CharKWDMAX];
  int value = 0;
  int i;

  if (fgetsMPI(line, sizeof(line), fp) == NULL) return ReadDefFileError(defname);
  if (fgetsMPI(line, sizeof(line), fp) == NULL) return ReadDefFileError(defname);
  if (sscanf(line, "%199s %d", key, &value) != 2) return ReadDefFileError(defname);
  if (CheckWords(key, "NQPTrans") != 0 || value <= 0) return ReadDefFileError(defname);
  for (i = 0; i < 3; i++) {
    if (fgetsMPI(line, sizeof(line), fp) == NULL) return ReadDefFileError(defname);
  }
  *ntrans = (unsigned int)value;
  return 0;
}

int ReadTransSymNInt(const char *defname, struct DefineList *def)
{
  FILE *fp = fopenMPI(defname, "r");
  unsigned int ntrans = 0;
  if (fp == NULL) return ReadDefFileError(defname);
  if (read_header_count(fp, defname, &ntrans) != 0) {
    fclose(fp);
    return -1;
  }
  fclose(fp);
  def->iFlgSymmetryBasis = TRUE;
  def->NSymTrans = ntrans;
  return 0;
}

static int parse_character_line(const char *line,
                                unsigned int ntrans,
                                unsigned int *op,
                                double complex *ch)
{
  int iop = -1;
  double re = 0.0;
  double im = 0.0;
  int nread = sscanf(line, "%d %lf %lf", &iop, &re, &im);
  if (nread != 2 && nread != 3) return -1;
  if (iop < 0 || (unsigned int)iop >= ntrans) return -1;
  if (!isfinite(re) || !isfinite(im)) return -2;
  *op = (unsigned int)iop;
  *ch = re + im * I;
  return 0;
}

static int parse_perm_line(const char *line,
                           unsigned int ntrans,
                           unsigned int nsite,
                           unsigned int *op,
                           unsigned int *site,
                           int *target,
                           int *anti)
{
  int iop = -1;
  int isite = -1;
  int itarget = -1;
  int ianti = 0;
  if (sscanf(line, "%d %d %d %d", &iop, &isite, &itarget, &ianti) != 4) return -1;
  if (iop < 0 || (unsigned int)iop >= ntrans) return -1;
  if (isite < 0 || (unsigned int)isite >= nsite) return -1;
  if (itarget < 0 || (unsigned int)itarget >= nsite) return -1;
  *op = (unsigned int)iop;
  *site = (unsigned int)isite;
  *target = itarget;
  *anti = ianti;
  return 0;
}

/*
  Case-insensitive test for the metadata keyword. Returns 1 for an exact
  match, -1 when the word merely starts with the keyword (a typo such as
  "MomentumIndex=3" must not pass silently as an ordinary comment), 0 otherwise.
*/
static int match_momentum_index_key(const char *word)
{
  static const char keyword[] = "momentumindex";
  size_t i;
  for (i = 0; keyword[i] != '\0'; i++) {
    if (word[i] == '\0' || tolower((unsigned char)word[i]) != keyword[i]) return 0;
  }
  return word[i] == '\0' ? 1 : -1;
}

/*
  Optional metadata in comment lines of the TransSym file, for example
  "# MomentumIndex 3" written by Standard mode. Comment lines never reach
  the fixed-layout reader below because fgetsMPI() skips them, so older
  versions of HPhi ignore the metadata. It is scanned separately on rank 0
  and broadcast. *momentum_index is -1 when the file carries no metadata.
*/
static int read_transsym_metadata(const char *defname, int *momentum_index)
{
  int value = -1;
  int status = 0;
  if (myrank == 0) {
    FILE *fp = fopen(defname, "r");
    char line[D_CharTmpReadDef + D_CharKWDMAX];
    if (fp == NULL) {
      status = -1;
    } else {
      while (fgets(line, sizeof(line), fp) != NULL) {
        char key[D_CharKWDMAX];
        char *end;
        long parsed;
        int key_end = 0;
        int has_digits;
        int match;
        const char *p = line;
        if (*p != '#') continue; /* same rule as fgetsMPI(): '#' in column 0 */
        p++;
        if (sscanf(p, "%199s%n", key, &key_end) != 1) continue;
        match = match_momentum_index_key(key);
        if (match == 0) continue;
        p += key_end;
        errno = 0;
        parsed = strtol(p, &end, 10);
        has_digits = end != p;
        while (isspace((unsigned char)*end)) end++;
        /* Check the range before narrowing to int; scanf's %d can wrap. */
        if (match < 0 || !has_digits || errno == ERANGE ||
            parsed < 0 || parsed > INT_MAX || *end != '\0') {
          fprintf(stdoutMPI,
                  "Error: TransSym metadata must be \"# MomentumIndex <non-negative integer>\" "
                  "with value <= %d: %s", INT_MAX, line);
          status = -1;
          break;
        }
        if (value >= 0 && value != parsed) {
          fprintf(stdoutMPI,
                  "Error: TransSym metadata MomentumIndex is given twice with different values (%d and %d).\n",
                  value, (int)parsed);
          status = -1;
          break;
        }
        value = (int)parsed;
      }
      fclose(fp);
    }
  }
  status = BcastMPI_i(0, status);
  value = BcastMPI_i(0, value);
  if (status != 0) return -1;
  *momentum_index = value;
  return 0;
}

int ReadTransSymFile(const char *defname, struct DefineList *def)
{
  FILE *fp;
  char line[D_CharTmpReadDef + D_CharKWDMAX];
  unsigned int ntrans = 0;
  unsigned int i;
  int *seen_char;
  int *seen_perm;
  if (read_transsym_metadata(defname, &def->iSymMomentumIndex) != 0) {
    return ReadDefFileError(defname);
  }
  fp = fopenMPI(defname, "r");
  if (fp == NULL) return ReadDefFileError(defname);
  if (read_header_count(fp, defname, &ntrans) != 0) {
    fclose(fp);
    return -1;
  }
  if (ntrans != def->NSymTrans) {
    fclose(fp);
    return ReadDefFileError(defname);
  }

  seen_char = (int *)calloc(def->NSymTrans, sizeof(int));
  if (seen_char == NULL) {
    fclose(fp);
    return -1;
  }
  for (i = 0; i < def->NSymTrans; i++) {
    unsigned int op;
    double complex ch;
    if (fgetsMPI(line, sizeof(line), fp) == NULL) {
      free(seen_char);
      fclose(fp);
      return ReadDefFileError(defname);
    }
    {
      int parse_status = parse_character_line(line, def->NSymTrans, &op, &ch);
      if (parse_status != 0) {
        if (parse_status == -2) {
          fprintf(stdoutMPI, "Error: TransSym character must be finite.\n");
        }
        free(seen_char);
        fclose(fp);
        return ReadDefFileError(defname);
      }
    }
    if (seen_char[op] != 0) {
      fprintf(stdoutMPI, "Error: duplicate TransSym character entry op=%u.\n", op);
      free(seen_char);
      fclose(fp);
      return ReadDefFileError(defname);
    }
    seen_char[op] = 1;
    def->SymTransChar[op] = ch;
  }
  for (i = 0; i < def->NSymTrans; i++) {
    if (seen_char[i] == 0) {
      fprintf(stdoutMPI, "Error: missing TransSym character entry op=%u.\n", i);
      free(seen_char);
      fclose(fp);
      return ReadDefFileError(defname);
    }
  }
  free(seen_char);

  seen_perm = (int *)calloc(def->NSymTrans * def->Nsite, sizeof(int));
  if (seen_perm == NULL) {
    fclose(fp);
    return -1;
  }
  for (i = 0; i < def->NSymTrans * def->Nsite; i++) {
    unsigned int op;
    unsigned int site;
    int target;
    int anti;
    if (fgetsMPI(line, sizeof(line), fp) == NULL) {
      free(seen_perm);
      fclose(fp);
      return ReadDefFileError(defname);
    }
    if (parse_perm_line(line, def->NSymTrans, def->Nsite, &op, &site, &target, &anti) != 0) {
      free(seen_perm);
      fclose(fp);
      return ReadDefFileError(defname);
    }
    if (seen_perm[op * def->Nsite + site] != 0) {
      fprintf(stdoutMPI, "Error: duplicate TransSym permutation entry op=%u site=%u.\n", op, site);
      free(seen_perm);
      fclose(fp);
      return ReadDefFileError(defname);
    }
    seen_perm[op * def->Nsite + site] = 1;
    def->SymTrans[op][site] = target;
    def->SymTransAnti[op][site] = anti;
  }
  for (i = 0; i < def->NSymTrans * def->Nsite; i++) {
    if (seen_perm[i] == 0) {
      fprintf(stdoutMPI, "Error: missing TransSym permutation entry.\n");
      free(seen_perm);
      fclose(fp);
      return ReadDefFileError(defname);
    }
  }
  free(seen_perm);
  fclose(fp);
  if (ValidateSymmetryGroupInput(def) != 0) return -1;
  if (def->iSymMomentumIndex >= 0) {
    fprintf(stdoutMPI, "TransSym metadata: MomentumIndex=%d\n", def->iSymMomentumIndex);
  }
  return 0;
}

static int has_fixed_spin_sector(const struct DefineList *def)
{
  if (def->iFlgSzConserved == TRUE) return TRUE;
  if (def->iCalcModel == Spin &&
      def->Nsite > 0 &&
      def->Nup + def->Ndown == def->Nsite) {
    return TRUE;
  }
  return FALSE;
}

static int has_fixed_spinless_sector(const struct DefineList *def)
{
  return def->iCalcModel == SpinlessFermion && def->Ne <= def->Nsite;
}

static int has_fixed_hubbard_sector(const struct DefineList *def)
{
  return def->iCalcModel == Hubbard &&
         def->Nup <= def->Nsite &&
         def->Ndown <= def->Nsite &&
         def->Ne == def->Nup + def->Ndown;
}

/*
 * Runtime capability matrix of the TransSym symmetry basis.
 *
 * One row per CalcType.  `enabled` says whether the method may run with
 * TransSym at all; the other columns say which layout, options and output
 * families that method supports.  Enabling one more method or feature is a
 * change to this table (and its tests), not to the validation code below.
 *
 *   distributed_layout  the distributed symmetry basis layout may be used.
 *   correlation         OneBodyG/TwoBodyG/ThreeBodyG/FourBodyG/SixBodyG/NBodyG.
 *   spectrum            CalcSpec != CALCSPEC_NOT.
 *   restart             ReStart != RESTART_NOT.
 *   eigenvec_output     OutputEigenVec.
 *   eigenvec_input      InputEigenVec.
 *   ham_output          OutputHam.
 *   ham_input           InputHam.
 *
 * Every feature column is 0 for every method in this version, so TransSym
 * runs output energy/norm/convergence only.
 */
struct SymmetryMethodCapability {
  int calc_type;
  const char *name;
  int enabled;
  int distributed_layout;
  int correlation;
  int spectrum;
  int restart;
  int eigenvec_output;
  int eigenvec_input;
  int ham_output;
  int ham_input;
};

static const struct SymmetryMethodCapability symmetry_method_capabilities[] = {
  /* calc_type     name             enabled dist corr spec rest evout evin hout hin */
  { Lanczos,       "Lanczos",       1,      0,   0,   0,   0,   0,    0,   0,   0   },
  { TPQCalc,       "TPQ",           0,      0,   0,   0,   0,   0,    0,   0,   0   },
  { FullDiag,      "FullDiag",      0,      0,   0,   0,   0,   0,    0,   0,   0   },
  { CG,            "CG",            1,      1,   0,   0,   0,   0,    0,   0,   0   },
  { TimeEvolution, "TimeEvolution", 0,      0,   0,   0,   0,   0,    0,   0,   0   },
  { cTPQ,          "cTPQ",          0,      0,   0,   0,   0,   0,    0,   0,   0   },
};

/* Fails to compile when a CalcType is added without a row in the table. */
typedef char symmetry_method_capabilities_cover_every_calc_type
    [(sizeof(symmetry_method_capabilities) /
      sizeof(symmetry_method_capabilities[0])) == NUM_CALCTYPE ? 1 : -1];

/* One requested option (or output family) and whether the method supports it. */
struct SymmetryOptionGate {
  const char *option;
  const char *reason;
  int requested;
  int supported;
};

static const struct SymmetryMethodCapability *find_symmetry_method_capability(int calc_type)
{
  /* readdef already bounds CalcType, so this row is only a safe fallback. */
  static const struct SymmetryMethodCapability unknown_method = {
    -1, "an unknown CalcType", 0, 0, 0, 0, 0, 0, 0, 0, 0
  };
  size_t i;
  for (i = 0; i < sizeof(symmetry_method_capabilities) /
                      sizeof(symmetry_method_capabilities[0]); i++) {
    if (symmetry_method_capabilities[i].calc_type == calc_type) {
      return &symmetry_method_capabilities[i];
    }
  }
  return &unknown_method;
}

/* A feature column only counts when the method itself is enabled. */
static int symmetry_method_supports(const struct SymmetryMethodCapability *cap,
                                    int feature)
{
  return cap->enabled != 0 && feature != 0;
}

static int symmetry_method_is_listed(const struct SymmetryMethodCapability *cap,
                                     int distributed_only)
{
  return cap->enabled != 0 && (distributed_only == 0 || cap->distributed_layout != 0);
}

/* Writes "A", "A and B" or "A, B and C": the enabled methods, or with
 * distributed_only the enabled methods that may use the distributed layout. */
static void list_symmetry_methods(char *buf, size_t size, int distributed_only)
{
  const size_t nrow = sizeof(symmetry_method_capabilities) /
                      sizeof(symmetry_method_capabilities[0]);
  size_t i;
  size_t total = 0;
  size_t seen = 0;
  size_t len = 0;
  int written;
  if (size == 0) return;
  buf[0] = '\0';
  for (i = 0; i < nrow; i++) {
    if (symmetry_method_is_listed(&symmetry_method_capabilities[i], distributed_only)) total++;
  }
  for (i = 0; i < nrow; i++) {
    if (!symmetry_method_is_listed(&symmetry_method_capabilities[i], distributed_only)) continue;
    written = snprintf(buf + len, size - len, "%s%s",
                       seen == 0 ? "" : (seen + 1 == total ? " and " : ", "),
                       symmetry_method_capabilities[i].name);
    if (written < 0 || (size_t)written >= size - len) break;
    len += (size_t)written;
    seen++;
  }
}

static int reject_unsupported_method(const struct SymmetryMethodCapability *cap)
{
  char supported[128];
  list_symmetry_methods(supported, sizeof(supported), 0);
  fprintf(stdoutMPI,
          "Error: TransSym symmetry basis v1 does not support %s; "
          "it supports only %s in this version.\n",
          cap->name, supported);
  return -1;
}

/* Rejects the first requested option that the method does not support. */
static int reject_first_unsupported_option(const struct SymmetryMethodCapability *cap,
                                           const struct SymmetryOptionGate *gates,
                                           size_t ngate)
{
  size_t i;
  for (i = 0; i < ngate; i++) {
    if (gates[i].requested != 0 && gates[i].supported == 0) {
      fprintf(stdoutMPI,
              "Error: TransSym symmetry basis v1 does not support %s with %s; %s.\n",
              gates[i].option, cap->name, gates[i].reason);
      return -1;
    }
  }
  return 0;
}

/* 1. Model and conserved-quantity sector. */
static int validate_symmetry_model_sector(const struct DefineList *def)
{
  if (def->iCalcModel != Spin && def->iCalcModel != SpinlessFermion &&
      def->iCalcModel != Hubbard) {
    fprintf(stdoutMPI, "Error: TransSym symmetry basis supports only Spin, SpinlessFermion, and Hubbard canonical models.\n");
    return -1;
  }
  if (def->iCalcModel == Spin && def->iFlgGeneralSpin != FALSE) {
    fprintf(stdoutMPI, "Error: TransSym symmetry basis v1 supports only Spin-1/2.\n");
    return -1;
  }
  if (def->iCalcModel == Spin && has_fixed_spin_sector(def) != TRUE) {
    fprintf(stdoutMPI, "Error: TransSym symmetry basis requires fixed 2Sz.\n");
    return -1;
  }
  if (def->iCalcModel == SpinlessFermion && has_fixed_spinless_sector(def) != TRUE) {
    fprintf(stdoutMPI, "Error: TransSym SpinlessFermion symmetry basis requires fixed Ncond/Ne.\n");
    return -1;
  }
  if (def->iCalcModel == Hubbard && has_fixed_hubbard_sector(def) != TRUE) {
    fprintf(stdoutMPI, "Error: TransSym Hubbard symmetry basis requires fixed Nup/Ndown.\n");
    return -1;
  }
  return 0;
}

/* 2. Method, Hamiltonian/eigenvector I/O and restart.
 *
 * The evaluation order is part of the diagnostics: when an input violates
 * several rules, the first message printed is the first of FullDiag,
 * OutputHam, InputHam, OutputEigenVec, InputEigenVec, ReStart and then any
 * other unsupported method.  (These used to be one combined check followed by
 * the Lanczos/CG-only check.) */
static int validate_symmetry_method_capability(const struct DefineList *def)
{
  const struct SymmetryMethodCapability *cap =
      find_symmetry_method_capability(def->iCalcType);
  const struct SymmetryOptionGate io_gates[] = {
    { "OutputHam",
      "InputHam/OutputHam Hamiltonian I/O is not available in this version",
      def->iOutputHam != FALSE,
      symmetry_method_supports(cap, cap->ham_output) },
    { "InputHam",
      "InputHam/OutputHam Hamiltonian I/O is not available in this version",
      def->iInputHam != FALSE,
      symmetry_method_supports(cap, cap->ham_input) },
    { "OutputEigenVec",
      "EigenVec/ReStart vector I/O is not available in this version",
      def->iOutputEigenVec != FALSE,
      symmetry_method_supports(cap, cap->eigenvec_output) },
    { "InputEigenVec",
      "EigenVec/ReStart vector I/O is not available in this version",
      def->iInputEigenVec != FALSE,
      symmetry_method_supports(cap, cap->eigenvec_input) },
    { "ReStart",
      "EigenVec/ReStart vector I/O is not available in this version",
      def->iReStart != RESTART_NOT,
      symmetry_method_supports(cap, cap->restart) },
  };
  /* FullDiag led the old combined check, so it is still reported first. */
  if (def->iCalcType == FullDiag && cap->enabled == 0) {
    return reject_unsupported_method(cap);
  }
  if (reject_first_unsupported_option(
          cap, io_gates, sizeof(io_gates) / sizeof(io_gates[0])) != 0) {
    return -1;
  }
  if (cap->enabled == 0) return reject_unsupported_method(cap);
  return 0;
}

/* 3. Spectrum calculation and correlation functions. */
static int validate_symmetry_output_capability(const struct DefineList *def)
{
  const struct SymmetryMethodCapability *cap =
      find_symmetry_method_capability(def->iCalcType);
  static const char correlation_reason[] = "it outputs energy/norm/convergence only";
  const struct SymmetryOptionGate gates[] = {
    { "spectrum calculations",
      "CalcSpec is not available in this version",
      def->iFlgCalcSpec != CALCSPEC_NOT,
      symmetry_method_supports(cap, cap->spectrum) },
    { "correlation functions (OneBodyG)", correlation_reason,
      def->NCisAjt > 0,
      symmetry_method_supports(cap, cap->correlation) },
    { "correlation functions (TwoBodyG)", correlation_reason,
      def->NCisAjtCkuAlvDC > 0,
      symmetry_method_supports(cap, cap->correlation) },
    { "correlation functions (ThreeBodyG)", correlation_reason,
      def->NTBody > 0,
      symmetry_method_supports(cap, cap->correlation) },
    { "correlation functions (FourBodyG)", correlation_reason,
      def->NFBody > 0,
      symmetry_method_supports(cap, cap->correlation) },
    { "correlation functions (SixBodyG)", correlation_reason,
      def->NSBody > 0,
      symmetry_method_supports(cap, cap->correlation) },
    { "correlation functions (NBodyG)", correlation_reason,
      def->NNBodyG > 0,
      symmetry_method_supports(cap, cap->correlation) },
  };
  return reject_first_unsupported_option(
      cap, gates, sizeof(gates) / sizeof(gates[0]));
}

/* 4. Model-specific Hamiltonian term families. */
static int validate_symmetry_term_families(const struct DefineList *def)
{
  if (def->NNBodyInterAll || def->NAnomalousTerm || def->NPairLiftCoupling ||
      (def->iCalcModel != Hubbard && (def->NCoulombIntra || def->NPairHopping)) ||
      (def->iCalcModel == SpinlessFermion &&
       (def->NHundCoupling || def->NIsingCoupling || def->NExchangeCoupling))) {
    fprintf(stdoutMPI, "Error: TransSym rejects unsupported term families for this canonical model.\n");
    return -1;
  }
  return 0;
}

int ValidateSymmetryRuntimeOptions(const struct BindStruct *X)
{
  const struct DefineList *def = &X->Def;
  if (def->iFlgSymmetryBasis == FALSE) return 0;
  /* The order fixes which message is printed first for an input that
   * violates several rules; keep it when adding checks. */
  if (validate_symmetry_model_sector(def) != 0) return -1;
  if (validate_symmetry_method_capability(def) != 0) return -1;
  if (validate_symmetry_output_capability(def) != 0) return -1;
  if (validate_symmetry_term_families(def) != 0) return -1;
  return 0;
}

int ValidateSymmetryBasisLayoutOptions(
    const struct BindStruct *X,
    enum SymmetryBasisLayout layout)
{
  const struct DefineList *def;
  const struct SymmetryMethodCapability *cap;
  char methods[128];
  if (X == NULL) return -1;
  def = &X->Def;
  if (layout == SYMMETRY_BASIS_REPLICATED) return 0;
  if (layout != SYMMETRY_BASIS_DISTRIBUTED) {
    fprintf(stdoutMPI,
            "Error: invalid TransSym symmetry basis layout selection.\n");
    return -1;
  }
  cap = find_symmetry_method_capability(def->iCalcType);
  if (def->iFlgSymmetryBasis != TRUE ||
      cap->enabled == 0 || cap->distributed_layout == 0) {
    list_symmetry_methods(methods, sizeof(methods), 1);
    fprintf(stdoutMPI,
            "Error: distributed symmetry basis is supported for "
            "TransSym %s runs only.\n",
            methods);
    return -1;
  }
  return 0;
}

static int same_unordered_pair(int a0, int a1, int b0, int b1)
{
  return (a0 == b0 && a1 == b1) || (a0 == b1 && a1 == b0);
}

static double sum_exchange_pair(const struct DefineList *def, int site0, int site1)
{
  unsigned int i;
  double sum = 0.0;
  if (site0 == site1) return 0.0;
  for (i = 0; i < def->NExchangeCoupling; i++) {
    if (def->ExchangeCoupling[i][0] == def->ExchangeCoupling[i][1]) continue;
    if (same_unordered_pair(def->ExchangeCoupling[i][0], def->ExchangeCoupling[i][1], site0, site1)) {
      sum += def->ParaExchangeCoupling[i];
    }
  }
  return sum;
}

static unsigned int count_ising_hund_pair(const struct DefineList *def, int site0, int site1)
{
  unsigned int i;
  unsigned int count = 0;
  if (site0 == site1) return 0;
  for (i = 0; i < def->NHundCoupling; i++) {
    if (def->HundCoupling[i][0] == def->HundCoupling[i][1]) continue;
    if (same_unordered_pair(def->HundCoupling[i][0], def->HundCoupling[i][1], site0, site1)) {
      count++;
    }
  }
  return count;
}

static unsigned int count_ising_coulomb_pair(const struct DefineList *def, int site0, int site1)
{
  unsigned int i;
  unsigned int count = 0;
  if (site0 == site1) return 0;
  for (i = 0; i < def->NCoulombInter; i++) {
    if (def->CoulombInter[i][0] == def->CoulombInter[i][1]) continue;
    if (same_unordered_pair(def->CoulombInter[i][0], def->CoulombInter[i][1], site0, site1)) {
      count++;
    }
  }
  return count;
}

static double sum_ising_hund_pair(const struct DefineList *def, int site0, int site1)
{
  unsigned int i;
  double sum = 0.0;
  if (site0 == site1) return 0.0;
  for (i = 0; i < def->NHundCoupling; i++) {
    if (def->HundCoupling[i][0] == def->HundCoupling[i][1]) continue;
    if (same_unordered_pair(def->HundCoupling[i][0], def->HundCoupling[i][1], site0, site1)) {
      sum += def->ParaHundCoupling[i];
    }
  }
  return sum;
}

static double sum_ising_coulomb_pair(const struct DefineList *def, int site0, int site1)
{
  unsigned int i;
  double sum = 0.0;
  if (site0 == site1) return 0.0;
  for (i = 0; i < def->NCoulombInter; i++) {
    if (def->CoulombInter[i][0] == def->CoulombInter[i][1]) continue;
    if (same_unordered_pair(def->CoulombInter[i][0], def->CoulombInter[i][1], site0, site1)) {
      sum += def->ParaCoulombInter[i];
    }
  }
  return sum;
}

static int validate_exchange_invariance(const struct DefineList *def)
{
  unsigned int g, i;
  for (g = 0; g < def->NSymTrans; g++) {
    for (i = 0; i < def->NExchangeCoupling; i++) {
      int src0 = def->ExchangeCoupling[i][0];
      int src1 = def->ExchangeCoupling[i][1];
      if (src0 < 0 || src1 < 0 ||
          (unsigned int)src0 >= def->Nsite ||
          (unsigned int)src1 >= def->Nsite) {
        fprintf(stdoutMPI,
                "Error: TransSym Exchange term %u has a site outside [0, Nsite).\n",
                i);
        return -1;
      }
      if (src0 == src1) continue;
      int a = def->SymTrans[g][src0];
      int b = def->SymTrans[g][src1];
      double src_sum = sum_exchange_pair(def, src0, src1);
      double mapped_sum = sum_exchange_pair(def, a, b);
      if (fabs(mapped_sum - src_sum) > 1.0e-10) {
        fprintf(stdoutMPI,
                "Error: TransSym Hamiltonian invariance failed for Exchange term %u under op %u.\n",
                i, g);
        return -1;
      }
    }
  }
  return 0;
}

static int validate_ising_hund_invariance(const struct DefineList *def)
{
  unsigned int g, i;
  for (g = 0; g < def->NSymTrans; g++) {
    for (i = 0; i < def->NHundCoupling; i++) {
      int src0 = def->HundCoupling[i][0];
      int src1 = def->HundCoupling[i][1];
      if (src0 < 0 || src1 < 0 ||
          (unsigned int)src0 >= def->Nsite ||
          (unsigned int)src1 >= def->Nsite) {
        fprintf(stdoutMPI,
                "Error: TransSym Ising Hund term %u has a site outside [0, Nsite).\n",
                i);
        return -1;
      }
      if (src0 == src1) continue;
      {
        int a = def->SymTrans[g][src0];
        int b = def->SymTrans[g][src1];
        double src_sum = sum_ising_hund_pair(def, src0, src1);
        double mapped_sum = sum_ising_hund_pair(def, a, b);
        if (count_ising_hund_pair(def, a, b) == 0) {
          fprintf(stdoutMPI,
                  "Error: TransSym Ising Hund term %u maps to a missing pair under op %u.\n",
                  i, g);
          return -1;
        }
        if (fabs(mapped_sum - src_sum) > 1.0e-10) {
          fprintf(stdoutMPI,
                  "Error: TransSym Ising Hund invariance failed for term %u under op %u.\n",
                  i, g);
          return -1;
        }
      }
    }
  }
  return 0;
}

static int validate_ising_coulomb_invariance(const struct DefineList *def)
{
  unsigned int g, i;
  for (g = 0; g < def->NSymTrans; g++) {
    for (i = 0; i < def->NCoulombInter; i++) {
      int src0 = def->CoulombInter[i][0];
      int src1 = def->CoulombInter[i][1];
      if (src0 < 0 || src1 < 0 ||
          (unsigned int)src0 >= def->Nsite ||
          (unsigned int)src1 >= def->Nsite) {
        fprintf(stdoutMPI,
                "Error: TransSym Ising CoulombInter term %u has a site outside [0, Nsite).\n",
                i);
        return -1;
      }
      if (src0 == src1) continue;
      {
        int a = def->SymTrans[g][src0];
        int b = def->SymTrans[g][src1];
        double src_sum = sum_ising_coulomb_pair(def, src0, src1);
        double mapped_sum = sum_ising_coulomb_pair(def, a, b);
        if (count_ising_coulomb_pair(def, a, b) == 0) {
          fprintf(stdoutMPI,
                  "Error: TransSym Ising CoulombInter term %u maps to a missing pair under op %u.\n",
                  i, g);
          return -1;
        }
        if (fabs(mapped_sum - src_sum) > 1.0e-10) {
          fprintf(stdoutMPI,
                  "Error: TransSym Ising CoulombInter invariance failed for term %u under op %u.\n",
                  i, g);
          return -1;
        }
      }
    }
  }
  return 0;
}

static int validate_spinless_coulomb_invariance(const struct DefineList *def)
{
  unsigned int g, i;
  for (g = 0; g < def->NSymTrans; g++) {
    for (i = 0; i < def->NCoulombInter; i++) {
      int src0 = def->CoulombInter[i][0];
      int src1 = def->CoulombInter[i][1];
      if (src0 < 0 || src1 < 0 ||
          (unsigned int)src0 >= def->Nsite ||
          (unsigned int)src1 >= def->Nsite) {
        fprintf(stdoutMPI,
                "Error: TransSym SpinlessFermion CoulombInter term %u has a site outside [0, Nsite).\n",
                i);
        return -1;
      }
      if (src0 == src1) {
        fprintf(stdoutMPI,
                "Error: TransSym SpinlessFermion CoulombInter term %u is on-site (i==j), which is not supported by the symmetry basis.\n",
                i);
        return -1;
      }
      {
        int a = def->SymTrans[g][src0];
        int b = def->SymTrans[g][src1];
        double src_sum = sum_ising_coulomb_pair(def, src0, src1);
        double mapped_sum = sum_ising_coulomb_pair(def, a, b);
        if (count_ising_coulomb_pair(def, a, b) == 0) {
          fprintf(stdoutMPI,
                  "Error: TransSym SpinlessFermion CoulombInter term %u maps to a missing pair under op %u.\n",
                  i, g);
          return -1;
        }
        if (fabs(mapped_sum - src_sum) > 1.0e-10) {
          fprintf(stdoutMPI,
                  "Error: TransSym SpinlessFermion CoulombInter invariance failed for term %u under op %u.\n",
                  i, g);
          return -1;
        }
      }
    }
  }
  return 0;
}

static int validate_ising_diagonal_invariance(const struct DefineList *def)
{
  if (def->NIsingCoupling == 0) return 0;
  if (validate_ising_hund_invariance(def) != 0) return -1;
  if (validate_ising_coulomb_invariance(def) != 0) return -1;
  return 0;
}

static double complex sum_spinless_transfer(const struct DefineList *def,
                                            int src,
                                            int dst)
{
  unsigned int i;
  double complex sum = 0.0;
  for (i = 0; i < def->NTransfer; i++) {
    if (def->GeneralTransfer[i][0] == src &&
        def->GeneralTransfer[i][1] == 0 &&
        def->GeneralTransfer[i][2] == dst &&
        def->GeneralTransfer[i][3] == 0) {
      sum += def->ParaGeneralTransfer[i];
    }
  }
  return sum;
}

static int validate_spinless_transfer_invariance(const struct DefineList *def)
{
  unsigned int g, i;
  for (g = 0; g < def->NSymTrans; g++) {
    for (i = 0; i < def->NTransfer; i++) {
      int src = def->GeneralTransfer[i][0];
      int src_spin = def->GeneralTransfer[i][1];
      int dst = def->GeneralTransfer[i][2];
      int dst_spin = def->GeneralTransfer[i][3];
      double complex src_sum;
      double complex mapped_sum;
      if (src < 0 || dst < 0 ||
          (unsigned int)src >= def->Nsite ||
          (unsigned int)dst >= def->Nsite ||
          src_spin != 0 || dst_spin != 0) {
        fprintf(stdoutMPI,
                "Error: TransSym SpinlessFermion Transfer term %u is outside the supported spinless site range.\n",
                i);
        return -1;
      }
      src_sum = sum_spinless_transfer(def, src, dst);
      mapped_sum = sum_spinless_transfer(def, def->SymTrans[g][src], def->SymTrans[g][dst]);
      if (cabs(mapped_sum - src_sum) > 1.0e-10) {
        fprintf(stdoutMPI,
                "Error: TransSym SpinlessFermion Transfer invariance failed for term %u under op %u.\n",
                i, g);
        return -1;
      }
    }
  }
  return 0;
}

static double complex sum_hubbard_transfer(const struct DefineList *def,
                                           int src,
                                           int dst,
                                           int spin)
{
  unsigned int i;
  double complex sum = 0.0;
  for (i = 0; i < def->NTransfer; i++) {
    if (def->GeneralTransfer[i][0] == src &&
        def->GeneralTransfer[i][1] == spin &&
        def->GeneralTransfer[i][2] == dst &&
        def->GeneralTransfer[i][3] == spin) {
      sum += def->ParaGeneralTransfer[i];
    }
  }
  return sum;
}

static int validate_hubbard_transfer_invariance(const struct DefineList *def)
{
  unsigned int g, i;
  for (g = 0; g < def->NSymTrans; g++) {
    for (i = 0; i < def->NTransfer; i++) {
      int src = def->GeneralTransfer[i][0];
      int src_spin = def->GeneralTransfer[i][1];
      int dst = def->GeneralTransfer[i][2];
      int dst_spin = def->GeneralTransfer[i][3];
      double complex src_sum;
      double complex mapped_sum;
      if (src < 0 || dst < 0 ||
          (unsigned int)src >= def->Nsite ||
          (unsigned int)dst >= def->Nsite ||
          src_spin < 0 || dst_spin < 0 ||
          src_spin > 1 || dst_spin > 1 ||
          src_spin != dst_spin) {
        fprintf(stdoutMPI,
                "Error: TransSym Hubbard Transfer term %u is outside the supported fixed-spin hopping range.\n",
                i);
        return -1;
      }
      src_sum = sum_hubbard_transfer(def, src, dst, src_spin);
      mapped_sum = sum_hubbard_transfer(def, def->SymTrans[g][src], def->SymTrans[g][dst], src_spin);
      if (cabs(mapped_sum - src_sum) > 1.0e-10) {
        fprintf(stdoutMPI,
                "Error: TransSym Hubbard Transfer invariance failed for term %u under op %u.\n",
                i, g);
        return -1;
      }
    }
  }
  return 0;
}

static unsigned int count_hubbard_coulomb_intra_site(const struct DefineList *def, int site)
{
  unsigned int i;
  unsigned int count = 0;
  for (i = 0; i < def->NCoulombIntra; i++) {
    if (def->CoulombIntra[i][0] == site) count++;
  }
  return count;
}

static double sum_hubbard_coulomb_intra_site(const struct DefineList *def, int site)
{
  unsigned int i;
  double sum = 0.0;
  for (i = 0; i < def->NCoulombIntra; i++) {
    if (def->CoulombIntra[i][0] == site) sum += def->ParaCoulombIntra[i];
  }
  return sum;
}

static int validate_hubbard_coulomb_intra_invariance(const struct DefineList *def)
{
  unsigned int g, i;
  for (g = 0; g < def->NSymTrans; g++) {
    for (i = 0; i < def->NCoulombIntra; i++) {
      int src = def->CoulombIntra[i][0];
      int mapped;
      double src_sum;
      double mapped_sum;
      if (src < 0 || (unsigned int)src >= def->Nsite) {
        fprintf(stdoutMPI,
                "Error: TransSym Hubbard CoulombIntra term %u has a site outside [0, Nsite).\n",
                i);
        return -1;
      }
      mapped = def->SymTrans[g][src];
      src_sum = sum_hubbard_coulomb_intra_site(def, src);
      mapped_sum = sum_hubbard_coulomb_intra_site(def, mapped);
      if (count_hubbard_coulomb_intra_site(def, mapped) == 0) {
        fprintf(stdoutMPI,
                "Error: TransSym Hubbard CoulombIntra term %u maps to a missing site under op %u.\n",
                i, g);
        return -1;
      }
      if (fabs(mapped_sum - src_sum) > 1.0e-10) {
        fprintf(stdoutMPI,
                "Error: TransSym Hubbard CoulombIntra invariance failed for term %u under op %u.\n",
                i, g);
        return -1;
      }
    }
  }
  return 0;
}

int ValidateSymmetryHamiltonian(const struct BindStruct *X)
{
  const struct DefineList *def = &X->Def;
  if (def->iFlgSymmetryBasis == FALSE) return 0;
  /* Preserve established diagnostics for the original subset. Extended
   * inputs use a combined polynomial, including cross-family cancellations. */
  if (SymmetryUsesExtendedTerms(def))
    return ValidateSymmetryTerms(def);
  if (def->iCalcModel == Spin) {
    if (validate_exchange_invariance(def) != 0) return -1;
    if (validate_ising_diagonal_invariance(def) != 0) return -1;
  } else if (def->iCalcModel == SpinlessFermion) {
    if (validate_spinless_transfer_invariance(def) != 0) return -1;
    if (validate_spinless_coulomb_invariance(def) != 0) return -1;
  } else if (def->iCalcModel == Hubbard) {
    if (validate_hubbard_transfer_invariance(def) != 0) return -1;
    if (validate_hubbard_coulomb_intra_invariance(def) != 0) return -1;
  }
  return ValidateSymmetryTerms(def);
}
