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
  (void)X;
  return 0;
}

int ValidateSymmetryHamiltonian(const struct BindStruct *X)
{
  (void)X;
  return 0;
}
