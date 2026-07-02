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
    if (parse_character_line(line, def->NSymTrans, &op, &ch) != 0) {
      free(seen_char);
      fclose(fp);
      return ReadDefFileError(defname);
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

int ValidateSymmetryRuntimeOptions(const struct BindStruct *X)
{
  const struct DefineList *def = &X->Def;
  if (def->iFlgSymmetryBasis == FALSE) return 0;
  if (nproc != 1) {
    fprintf(stdoutMPI, "Error: TransSym symmetry basis is serial-only in v1.\n");
    return -1;
  }
  if (def->iCalcModel != Spin) {
    fprintf(stdoutMPI, "Error: TransSym symmetry basis v1 supports only Spin canonical model.\n");
    return -1;
  }
  if (def->iFlgGeneralSpin != FALSE) {
    fprintf(stdoutMPI, "Error: TransSym symmetry basis v1 supports only Spin-1/2.\n");
    return -1;
  }
  if (def->iFlgSzConserved != TRUE) {
    fprintf(stdoutMPI, "Error: TransSym symmetry basis requires fixed 2Sz.\n");
    return -1;
  }
  if (def->iCalcType == FullDiag || def->iOutputHam != FALSE || def->iInputHam != FALSE) {
    fprintf(stdoutMPI, "Error: TransSym symmetry basis v1 does not support FullDiag/InputHam/OutputHam.\n");
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
  if (def->NTransfer > 0 || def->NPairHopping > 0 || def->NPairLiftCoupling > 0 ||
      def->NNBodyInterAll > 0 || def->NAnomalousTerm > 0) {
    fprintf(stdoutMPI, "Error: TransSym symmetry basis v1 rejects unsupported Spin term families.\n");
    return -1;
  }
  if (def->EDNChemi > 0 || def->NCoulombIntra > 0 || def->NCoulombInter > 0 ||
      def->NHundCoupling > 0 || def->NIsingCoupling > 0 ||
      def->NInterAll > 0) {
    fprintf(stdoutMPI, "Error: TransSym symmetry basis v1 initially supports Exchange terms only.\n");
    return -1;
  }
  return 0;
}

static int same_unordered_pair(int a0, int a1, int b0, int b1)
{
  return (a0 == b0 && a1 == b1) || (a0 == b1 && a1 == b0);
}

static int find_exchange_pair(const struct DefineList *def, int site0, int site1, double value)
{
  unsigned int i;
  for (i = 0; i < def->NExchangeCoupling; i++) {
    if (same_unordered_pair(def->ExchangeCoupling[i][0], def->ExchangeCoupling[i][1], site0, site1) &&
        fabs(def->ParaExchangeCoupling[i] - value) < 1.0e-10) {
      return TRUE;
    }
  }
  return FALSE;
}

static int validate_exchange_invariance(const struct DefineList *def)
{
  unsigned int g, i;
  for (g = 0; g < def->NSymTrans; g++) {
    for (i = 0; i < def->NExchangeCoupling; i++) {
      int a = def->SymTrans[g][def->ExchangeCoupling[i][0]];
      int b = def->SymTrans[g][def->ExchangeCoupling[i][1]];
      if (find_exchange_pair(def, a, b, def->ParaExchangeCoupling[i]) != TRUE) {
        fprintf(stdoutMPI,
                "Error: TransSym Hamiltonian invariance failed for Exchange term %u under op %u.\n",
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
  if (validate_exchange_invariance(def) != 0) return -1;
  return 0;
}
