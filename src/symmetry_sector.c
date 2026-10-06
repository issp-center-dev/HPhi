#include <float.h>
#include <inttypes.h>
#include <math.h>
#include "Common.h"
#include "FileIO.h"
#include "struct.h"
#include "symmetry_basis.h"
#include "symmetry_sector.h"
#include "wrapperMPI.h"
#ifdef MPI
#include <mpi.h>
#endif

#define FNV_OFFSET UINT64_C(14695981039346656037)
#define FNV_PRIME UINT64_C(1099511628211)

/* Explicit little-endian encoding; never hash struct padding/native endian. */
static uint64_t hash_integer(uint64_t hash, uint64_t value, unsigned int bytes)
{
  unsigned int b;
  for (b = 0; b < bytes; ++b) {
    hash = (hash ^ (value & UINT64_C(255))) * FNV_PRIME;
    value >>= 8;
  }
  return hash;
}

static uint64_t hash_tag(uint64_t hash, const char *tag)
{
  do {
    hash = (hash ^ (unsigned char)*tag) * FNV_PRIME;
  } while (*tag++ != '\0');
  return hash;
}

static int hash_double(uint64_t *hash, double value)
{
  uint64_t bits;
  if (sizeof(value) != sizeof(bits) || DBL_MANT_DIG != 53 ||
      DBL_MAX_EXP != 1024 || !isfinite(value)) return -1;
  memcpy(&bits, &value, sizeof(bits));
  *hash = hash_integer(*hash, bits, 8);
  return 0;
}

struct GroupDigestEntry {
  const int *perm;
  const int *anti;
  double complex character;
  unsigned int nsite;
};

static int compare_group_entries(const void *a, const void *b)
{
  const struct GroupDigestEntry *left = a, *right = b;
  unsigned int site;
  for (site = 0; site < left->nsite; ++site) {
    if (left->perm[site] < right->perm[site]) return -1;
    if (left->perm[site] > right->perm[site]) return 1;
  }
  return 0;
}

int ComputeSymmetryGroupDigest(const struct DefineList *def, uint64_t *digest)
{
  struct GroupDigestEntry *entries;
  unsigned int g, site;
  uint64_t hash = hash_tag(FNV_OFFSET, "hphi-group-fnv1a64-v1");
  if (digest == NULL) return -1;
  *digest = 0;
  if (def == NULL || def->Nsite == 0 || def->NSymTrans == 0 ||
      def->SymTrans == NULL || def->SymTransAnti == NULL ||
      def->SymTransChar == NULL) return -1;
  entries = calloc(def->NSymTrans, sizeof(*entries));
  if (entries == NULL) return -1;
  for (g = 0; g < def->NSymTrans; ++g) {
    double complex ch = def->SymTransChar[g];
    if (def->SymTrans[g] == NULL || def->SymTransAnti[g] == NULL ||
        !isfinite(creal(ch)) || !isfinite(cimag(ch)) ||
        fabs(creal(ch)) > 2.0 || fabs(cimag(ch)) > 2.0) {
      free(entries);
      return -1;
    }
    entries[g].perm = def->SymTrans[g];
    entries[g].anti = def->SymTransAnti[g];
    entries[g].character = ch;
    entries[g].nsite = def->Nsite;
  }
  qsort(entries, def->NSymTrans, sizeof(*entries), compare_group_entries);
  hash = hash_integer(hash, def->Nsite, 4);
  hash = hash_integer(hash, def->NSymTrans, 4);
  for (g = 0; g < def->NSymTrans; ++g) {
    for (site = 0; site < def->Nsite; ++site) {
      hash = hash_integer(hash, (uint32_t)entries[g].perm[site], 4);
      hash = hash_integer(hash, (uint32_t)entries[g].anti[site], 4);
    }
    hash = hash_integer(hash,
        (uint64_t)(int64_t)llround(1e10 * creal(entries[g].character)), 8);
    hash = hash_integer(hash,
        (uint64_t)(int64_t)llround(1e10 * cimag(entries[g].character)), 8);
  }
  free(entries);
  *digest = hash;
  return 0;
}

static int hash_terms(uint64_t *hash, const char *tag, unsigned int count,
                      unsigned int width, int *const *indices,
                      const double *real_values,
                      const double complex *complex_values)
{
  unsigned int term, col;
  *hash = hash_tag(*hash, tag);
  *hash = hash_integer(*hash, count, 4);
  if (count > 0 && (indices == NULL ||
                   (real_values == NULL && complex_values == NULL))) return -1;
  for (term = 0; term < count; ++term) {
    if (indices[term] == NULL) return -1;
    for (col = 0; col < width; ++col)
      *hash = hash_integer(*hash, (uint32_t)indices[term][col], 4);
    if (hash_double(hash, real_values != NULL ? real_values[term] :
                    creal(complex_values[term])) != 0) return -1;
    if (complex_values != NULL &&
        hash_double(hash, cimag(complex_values[term])) != 0) return -1;
  }
  return 0;
}

int ComputeSymmetryHamiltonianDigest(const struct DefineList *def, uint64_t *digest)
{
  uint64_t hash;
  if (digest == NULL) return -1;
  *digest = 0;
  if (def == NULL) return -1;
  hash = hash_tag(FNV_OFFSET, def->iCalcModel == SpinGC
      ? "hphi-parsed-hamiltonian-fnv1a64-v3"
      : "hphi-parsed-hamiltonian-fnv1a64-v2");
  /* Extend this scope when enabling additional term families. Never silently
   * omit an unsupported family from an otherwise successful fingerprint. */
  if (def->NNBodyInterAll || def->NAnomalousTerm ||
      (def->iCalcModel != SpinGC && def->NPairLiftCoupling))
    return -1;
  hash = hash_integer(hash, (uint32_t)def->iCalcModel, 4);
  hash = hash_integer(hash, def->Nsite, 4);
  hash = hash_integer(hash, def->NIsingCoupling, 4);
  if (hash_terms(&hash, "transfer", def->EDNTransfer, 4,
                 def->EDGeneralTransfer, NULL, def->EDParaGeneralTransfer) ||
      hash_terms(&hash, "coulomb_intra", def->NCoulombIntra, 1,
                 def->CoulombIntra, def->ParaCoulombIntra, NULL) ||
      hash_terms(&hash, "coulomb_inter", def->NCoulombInter, 2,
                 def->CoulombInter, def->ParaCoulombInter, NULL) ||
      hash_terms(&hash, "hund", def->NHundCoupling, 2,
                 def->HundCoupling, def->ParaHundCoupling, NULL) ||
      hash_terms(&hash, "exchange", def->NExchangeCoupling, 2,
                 def->ExchangeCoupling, def->ParaExchangeCoupling, NULL))
    return -1;
  if (hash_terms(&hash, "pair_hopping", def->NPairHopping, 2,
                 def->PairHopping, def->ParaPairHopping, NULL) ||
      hash_terms(&hash, "interall_diagonal", def->NInterAll_Diagonal, 4,
                 def->InterAll_Diagonal, def->ParaInterAll_Diagonal, NULL) ||
      hash_terms(&hash, "interall_offdiagonal", def->NInterAll_OffDiagonal, 8,
                 def->InterAll_OffDiagonal, NULL, def->ParaInterAll_OffDiagonal)) return -1;
  hash = hash_tag(hash, "chemical_potential");
  hash = hash_integer(hash, def->EDNChemi, 4);
  if (def->EDNChemi && (!def->EDChemi || !def->EDSpinChemi || !def->EDParaChemi)) return -1;
  for (unsigned int i = 0; i < def->EDNChemi; ++i) {
    hash = hash_integer(hash, (uint32_t)def->EDChemi[i], 4);
    hash = hash_integer(hash, (uint32_t)def->EDSpinChemi[i], 4);
    if (hash_double(&hash, def->EDParaChemi[i])) return -1;
  }
  if (def->iCalcModel == SpinGC &&
      hash_terms(&hash, "pair_lift", def->NPairLiftCoupling, 2,
                 def->PairLiftCoupling, def->ParaPairLiftCoupling, NULL)) return -1;
  *digest = hash;
  return 0;
}

int ComputeSymmetrySectorDigest(const struct SymmetryBasisRuntime *sym,
                                struct SymmetrySectorDigest *digest)
{
  unsigned long int i, count = 0;
  struct SymmetrySectorDigest local = {0, 0, 0};
  int error = digest == NULL || sym == NULL || sym->enabled != TRUE;
  if (digest != NULL) memset(digest, 0, sizeof(*digest));
  if (!error) {
    if (sym->basis_layout == SYMMETRY_BASIS_REPLICATED) count = sym->dim;
    else if (sym->basis_layout == SYMMETRY_BASIS_DISTRIBUTED) count = sym->local_dim;
    else error = 1;
  }
  for (i = 1; !error && i <= count; ++i) {
    const struct SymmetryBasisVector *entry =
        sym->basis_layout == SYMMETRY_BASIS_REPLICATED ?
        SymmetryBasisReplicatedGlobalEntry(sym, i) : SymmetryBasisLocalEntry(sym, i);
    uint64_t hash = hash_tag(FNV_OFFSET, "hphi-sector-multiset-v1");
    if (entry == NULL) { error = 1; break; }
    hash = hash_integer(hash, (uint64_t)entry->rep_state, 8);
    hash = hash_integer(hash, entry->orbit_size, 4);
    hash = hash_integer(hash, entry->stabilizer_size, 4);
    local.count++;
    local.xor_hash ^= hash;
    local.sum_hash += hash;
  }
  if (SumMPI_i(error) != 0) return -1;
#ifdef MPI
  if (sym->basis_layout == SYMMETRY_BASIS_DISTRIBUTED) {
    uint64_t sums[2] = {local.count, local.sum_hash}, global_sums[2], global_xor;
    if (MPI_Allreduce(sums, global_sums, 2, MPI_UINT64_T, MPI_SUM,
                      MPI_COMM_WORLD) != MPI_SUCCESS ||
        MPI_Allreduce(&local.xor_hash, &global_xor, 1, MPI_UINT64_T, MPI_BXOR,
                      MPI_COMM_WORLD) != MPI_SUCCESS) return -1;
    local.count = global_sums[0];
    local.sum_hash = global_sums[1];
    local.xor_hash = global_xor;
  }
#endif
  if (SumMPI_i(local.count != (uint64_t)sym->dim) != 0) return -1;
  *digest = local;
  return 0;
}

int WriteSymmetryCanonicalTPQSchedule(const struct BindStruct *X, int rows,
                                     const double *beta, const int *orders)
{
  char name[D_FileNameMax];
  FILE *fp = NULL;
  int error = 0;
  if (!X->Def.iFlgSymmetryBasis) return 0;
  if (myrank == 0) {
    int length = snprintf(name, sizeof(name), "%s%ssymmetry_sector.dat",
                          X->Def.iOutputDataHead ? X->Def.CDataFileHead : "",
                          X->Def.iOutputDataHead ? "_" : "");
    error = length < 0 || (size_t)length + strlen(cParentOutputFolder) >= sizeof(name) ||
            (rows > 0 && (beta == NULL || orders == NULL));
    if (!error && childfopenMPI(name, "a", &fp) != 0) error = 1;
    if (!error) {
      fprintf(fp, "canonical_tpq_steps=%u\nbeta_schedule=%s\n",
              rows > 0 ? (unsigned int)rows : X->Def.Lanczos_max,
              rows > 0 ? "explicit" : "uniform");
      if (rows > 0) {
        /* order[i] advances beta[i] to beta[i+1]; the final order is unused. */
        for (int i = 0; i < rows; ++i)
          fprintf(fp, "invtemp_row_%d=%.17g %d\n", i, beta[i], orders[i]);
      } else {
        fprintf(fp, "beta_initial=0\nbeta_step=%.17g\nexpand_coef=%d\n",
                1.0 / LargeValue, X->Def.Param.ExpandCoef);
      }
      if (ferror(fp)) error = 1;
      if (fclose(fp) != 0) error = 1;
    }
  }
  if (SumMPI_i(error) != 0) {
    fprintf(stdoutMPI, "Error: failed to record the symmetry cTPQ schedule.\n");
    return -1;
  }
  return 0;
}

int WriteSymmetrySectorManifest(const struct BindStruct *X)
{
  const struct DefineList *def = &X->Def;
  struct SymmetrySectorDigest sector;
  uint64_t group = 0, hamiltonian = 0;
  char name[D_FileNameMax];
  int error = 0, length, threads = 1;
  FILE *fp = NULL;
  const char *model, *method;
  if (!def->iFlgSymmetryBasis) return 0;
  if (ComputeSymmetrySectorDigest(X->Sym, &sector) != 0) return -1;
  if (myrank == 0) {
    error = ComputeSymmetryGroupDigest(def, &group) != 0 ||
            ComputeSymmetryHamiltonianDigest(def, &hamiltonian) != 0;
    length = snprintf(name, sizeof(name), "%s%ssymmetry_sector.dat",
                      def->iOutputDataHead == 1 ? def->CDataFileHead : "",
                      def->iOutputDataHead == 1 ? "_" : "");
    if (length < 0 || (size_t)length + strlen(cParentOutputFolder) >= sizeof(name))
      error = 1;
    if (!error && childfopenMPI(name, "w", &fp) != 0) error = 1;
    if (!error) {
#ifdef _OPENMP
      threads = omp_get_max_threads();
#endif
      model = def->iCalcModel == SpinGC ? "SpinGC" :
              def->iCalcModel == Spin ? "Spin" :
              def->iCalcModel == SpinlessFermion ? "SpinlessFermion" :
              def->iCalcModel == tJ ? "tJ" : "Hubbard";
      switch (def->iCalcType) {
      case Lanczos: method = "Lanczos"; break;
      case CG: method = "CG"; break;
      case TPQCalc: method = "TPQ"; break;
      case cTPQ: method = "cTPQ"; break;
      case FullDiag: method = "FullDiag"; break;
      case TimeEvolution: method = "TimeEvolution"; break;
      default: method = "unknown"; break;
      }
      fprintf(fp, "format=HPhiSymmetrySector version=1\ncalc_type=%s\nmodel=%s\n"
              "nsite=%u\nfull_dim=%lu\nsector_dim=%lu\ngroup_order=%u\n",
              method, model, def->Nsite, X->Sym->full_dim, X->Sym->dim, def->NSymTrans);
      if (def->iCalcModel == SpinGC)
        fprintf(fp, "fixed_quantities=none\npair_lift=%u\n", def->NPairLiftCoupling);
      else if (def->iCalcModel == Spin)
        fprintf(fp, "fixed_2sz=%d\n", (int)def->Nup - (int)def->Ndown);
      else {
        fprintf(fp, "fixed_ne=%u\n", def->Ne);
        if (def->iCalcModel == Hubbard || def->iCalcModel == tJ)
          fprintf(fp, "fixed_nup=%u\nfixed_ndown=%u\n", def->Nup, def->Ndown);
      }
      if (def->iSymMomentumIndex >= 0)
        fprintf(fp, "momentum_index=%d\n", def->iSymMomentumIndex);
      fprintf(fp, "group_digest=hphi-group-fnv1a64-v1:%016" PRIx64 "\n"
              "sector_digest=hphi-sector-multiset-v1:%" PRIu64 ":%016" PRIx64 ":%016" PRIx64 "\n"
              "hamiltonian_digest=hphi-parsed-hamiltonian-fnv1a64-v%d:%016" PRIx64 "\n",
              group, sector.count, sector.xor_hash, sector.sum_hash,
              def->iCalcModel == SpinGC ? 3 : 2, hamiltonian);
      fprintf(fp, "basis_layout=%s\nmpi_ranks=%d\nomp_threads=%d\n"
              "term_scope=transfer:%u coulomb_intra:%u coulomb_inter:%u hund:%u exchange:%u ising:%u pair_hopping:%u chemi:%u interall_diagonal:%u interall_offdiagonal:%u\n"
              "exct=%u\nlanczos_max=%u\nlanczos_eps=%d\ninitial_iv=%ld\n",
              X->Sym->basis_layout == SYMMETRY_BASIS_DISTRIBUTED ? "distributed" : "replicated",
              nproc, threads, def->EDNTransfer, def->NCoulombIntra, def->NCoulombInter,
              def->NHundCoupling, def->NExchangeCoupling, def->NIsingCoupling,
              def->NPairHopping, def->EDNChemi, def->NInterAll_Diagonal, def->NInterAll_OffDiagonal,
              def->k_exct, def->Lanczos_max, def->LanczosEps, def->initial_iv);
      if (def->iCalcType == TPQCalc || def->iCalcType == cTPQ)
        fprintf(fp, "ensemble=single_symmetry_sector\nnum_ave=%d\nlarge_value=%.17g\ninitial_vec_type=%d\n",
                NumAve, LargeValue, def->iInitialVecType);
      if (def->iCalcType == FullDiag)
        fprintf(fp, "output_scope=eigenvalues\nsolver_id=%d\neigenvalue_file=%s_energy_sector.dat\nmatrix_storage=%s\n",
                def->iSolver, def->CDataFileHead,
                def->iSolver == SOLVER_ELPA && nproc > 1 ? "column_panel" : "replicated");
      if (ferror(fp)) error = 1;
      if (fclose(fp) != 0) error = 1;
    }
  }
  if (SumMPI_i(error) != 0) {
    fprintf(stdoutMPI, "Error: could not write TransSym sector manifest.\n");
    return -1;
  }
  return 0;
}
