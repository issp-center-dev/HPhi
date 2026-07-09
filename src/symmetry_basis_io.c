#include <math.h>
#include <stdlib.h>
#include "symmetry_basis.h"
#include "symmetry_basis_io.h"
#include "readdef.h"
#include "struct.h"
#include "wrapperMPI.h"

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

int ReadTransSymFile(const char *defname, struct DefineList *def)
{
  FILE *fp = fopenMPI(defname, "r");
  char line[D_CharTmpReadDef + D_CharKWDMAX];
  unsigned int ntrans = 0;
  unsigned int i;
  int *seen_char;
  int *seen_perm;
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
  return ValidateSymmetryGroupInput(def);
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

int ValidateSymmetryRuntimeOptions(const struct BindStruct *X)
{
  const struct DefineList *def = &X->Def;
  if (def->iFlgSymmetryBasis == FALSE) return 0;
  if (def->iCalcModel != Spin && def->iCalcModel != SpinlessFermion) {
    fprintf(stdoutMPI, "Error: TransSym symmetry basis supports only Spin and SpinlessFermion canonical models.\n");
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
  if (def->iCalcType == FullDiag || def->iOutputHam != FALSE || def->iInputHam != FALSE ||
      def->iOutputEigenVec != FALSE || def->iInputEigenVec != FALSE ||
      def->iReStart != RESTART_NOT) {
    fprintf(stdoutMPI,
            "Error: TransSym symmetry basis v1 does not support FullDiag/InputHam/OutputHam/EigenVec/ReStart.\n");
    return -1;
  }
  if (def->iCalcType != Lanczos && def->iCalcType != CG) {
    fprintf(stdoutMPI, "Error: TransSym symmetry basis v1 supports only Lanczos and CG.\n");
    return -1;
  }
  if (def->iFlgCalcSpec != CALCSPEC_NOT) {
    fprintf(stdoutMPI, "Error: TransSym symmetry basis v1 does not support spectrum calculations.\n");
    return -1;
  }
  if (def->NCisAjt > 0 || def->NCisAjtCkuAlvDC > 0 ||
      def->NTBody > 0 || def->NFBody > 0 || def->NSBody > 0 ||
      def->NNBodyG > 0) {
    fprintf(stdoutMPI, "Error: TransSym symmetry basis v1 outputs energy/norm/convergence only.\n");
    return -1;
  }
  if (def->iCalcModel == Spin) {
    if (def->NTransfer > 0 || def->NPairHopping > 0 || def->NPairLiftCoupling > 0 ||
        def->NNBodyInterAll > 0 || def->NAnomalousTerm > 0) {
      fprintf(stdoutMPI, "Error: TransSym symmetry basis v1 rejects unsupported Spin term families.\n");
      return -1;
    }
    if (def->EDNChemi > 0 || def->NCoulombIntra > 0 || def->NInterAll > 0) {
      fprintf(stdoutMPI, "Error: TransSym symmetry basis v1 supports Exchange and Ising terms only.\n");
      return -1;
    }
    if (def->NCoulombInter != def->NIsingCoupling ||
        def->NHundCoupling != def->NIsingCoupling) {
      fprintf(stdoutMPI,
              "Error: TransSym symmetry basis v1 supports Exchange and Ising terms only; "
              "direct CoulombInter/Hund terms are not supported.\n");
      return -1;
    }
  } else if (def->iCalcModel == SpinlessFermion) {
    if (def->EDNChemi > 0 || def->NCoulombIntra > 0 ||
        def->NHundCoupling > 0 || def->NIsingCoupling > 0 || def->NExchangeCoupling > 0 ||
        def->NPairHopping > 0 || def->NPairLiftCoupling > 0 || def->NInterAll > 0 ||
        def->NNBodyInterAll > 0 || def->NAnomalousTerm > 0) {
      fprintf(stdoutMPI, "Error: TransSym SpinlessFermion symmetry basis supports Transfer and CoulombInter terms only.\n");
      return -1;
    }
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

int ValidateSymmetryHamiltonian(const struct BindStruct *X)
{
  const struct DefineList *def = &X->Def;
  if (def->iFlgSymmetryBasis == FALSE) return 0;
  if (def->iCalcModel == Spin) {
    if (validate_exchange_invariance(def) != 0) return -1;
    if (validate_ising_diagonal_invariance(def) != 0) return -1;
  } else if (def->iCalcModel == SpinlessFermion) {
    if (validate_spinless_transfer_invariance(def) != 0) return -1;
    if (validate_spinless_coulomb_invariance(def) != 0) return -1;
  }
  return 0;
}
