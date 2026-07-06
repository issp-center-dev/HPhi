/* HPhi  -  Quantum Lattice Model Simulator */
#include "green_output.h"
#include "DefCommon.h"
#include "FileIO.h"
#include "global.h"
#include "wrapperMPI.h"

static int GreenOutputAggregateFamily(const struct BindStruct *X)
{
  if (X->Def.iOutputGreenFormat != OUTPUTGREENFORMAT_AGGREGATE) return 0;
  switch (X->Def.iCalcType) {
  case TPQCalc:
  case cTPQ:
    return 1;
  case TimeEvolution:
    return 2;
  case FullDiag:
  case CG:
    return 3;
  default:
    return 0;
  }
}

int GreenOutputUsesAggregate(const struct BindStruct *X)
{
  return GreenOutputAggregateFamily(X) != 0;
}

int GreenOutputKindUsesAggregate(const struct BindStruct *X, GreenOutputKind kind)
{
  if (!GreenOutputUsesAggregate(X)) return FALSE;
  return TRUE;
}

const char *GreenOutputOpenMode(const struct BindStruct *X)
{
  return GreenOutputUsesAggregate(X) ? "a" : "w";
}

static const char *GreenOutputFormatForFamily(int family, GreenOutputKind kind)
{
  if (family == 1) {
    switch (kind) {
    case GreenOutputOneBody: return cFileName1BGreen_TPQ_Aggregate;
    case GreenOutputTwoBody: return cFileName2BGreen_TPQ_Aggregate;
    case GreenOutputThreeBody: return cFileName3BGreen_TPQ_Aggregate;
    case GreenOutputFourBody: return cFileName4BGreen_TPQ_Aggregate;
    case GreenOutputSixBody: return cFileName6BGreen_TPQ_Aggregate;
    case GreenOutputNBody: return cFileNameNBodyG_TPQ_Aggregate;
    case GreenOutputAnomalous: return cFileNameAnomalousG_TPQ_Aggregate;
    }
  }
  if (family == 2) {
    switch (kind) {
    case GreenOutputOneBody: return cFileName1BGreen_TE_Aggregate;
    case GreenOutputTwoBody: return cFileName2BGreen_TE_Aggregate;
    case GreenOutputThreeBody: return cFileName3BGreen_TE_Aggregate;
    case GreenOutputFourBody: return cFileName4BGreen_TE_Aggregate;
    case GreenOutputSixBody: return cFileName6BGreen_TE_Aggregate;
    case GreenOutputNBody: return cFileNameNBodyG_TE_Aggregate;
    case GreenOutputAnomalous: return cFileNameAnomalousG_TE_Aggregate;
    }
  }
  if (family == 3) {
    switch (kind) {
    case GreenOutputOneBody: return cFileName1BGreen_Eigen_Aggregate;
    case GreenOutputTwoBody: return cFileName2BGreen_Eigen_Aggregate;
    case GreenOutputThreeBody: return cFileName3BGreen_Eigen_Aggregate;
    case GreenOutputFourBody: return cFileName4BGreen_Eigen_Aggregate;
    case GreenOutputSixBody: return cFileName6BGreen_Eigen_Aggregate;
    case GreenOutputNBody: return cFileNameNBodyG_Eigen_Aggregate;
    case GreenOutputAnomalous: return cFileNameAnomalousG_Eigen_Aggregate;
    }
  }
  return NULL;
}

int GreenOutputFileName(const struct BindStruct *X, GreenOutputKind kind, char *sdt)
{
  int family;
  const char *fmt = NULL;
  if (!GreenOutputKindUsesAggregate(X, kind)) return -1;
  family = GreenOutputAggregateFamily(X);
  fmt = GreenOutputFormatForFamily(family, kind);
  if (fmt == NULL) return -1;
  sprintf(sdt, fmt, X->Def.CDataFileHead);
  return 0;
}

int GreenOutputWriteIndexPrefix(FILE *fp, const struct BindStruct *X)
{
  switch (GreenOutputAggregateFamily(X)) {
  case 1:
    fprintf(fp, " %4d %4d ", X->Def.irand, X->Def.istep);
    return 0;
  case 2:
    fprintf(fp, " %4d ", X->Def.istep);
    return 0;
  case 3:
    fprintf(fp, " %4d ", X->Phys.eigen_num);
    return 0;
  default:
    return 0;
  }
}

static int GreenOutputTruncateFile(struct BindStruct *X, GreenOutputKind kind)
{
  FILE *fp = NULL;
  char sdt[D_FileNameMax];
  if (!GreenOutputKindUsesAggregate(X, kind)) return 0;
  if (GreenOutputFileName(X, kind, sdt) != 0) return -1;
  if (childfopenMPI(sdt, "w", &fp) != 0) return -1;
  fclose(fp);
  return 0;
}

int GreenOutputInitializeAggregateFiles(struct BindStruct *X)
{
  if (!GreenOutputUsesAggregate(X)) return 0;
  if (X->Def.NCisAjt > 0 &&
      GreenOutputTruncateFile(X, GreenOutputOneBody) != 0) return -1;
  if (X->Def.NCisAjtCkuAlvDC > 0 &&
      GreenOutputTruncateFile(X, GreenOutputTwoBody) != 0) return -1;
  if (X->Def.NTBody > 0 &&
      GreenOutputTruncateFile(X, GreenOutputThreeBody) != 0) return -1;
  if (X->Def.NFBody > 0 &&
      GreenOutputTruncateFile(X, GreenOutputFourBody) != 0) return -1;
  if (X->Def.NSBody > 0 &&
      GreenOutputTruncateFile(X, GreenOutputSixBody) != 0) return -1;
  if (X->Def.NNBodyG > 0 &&
      GreenOutputTruncateFile(X, GreenOutputNBody) != 0) return -1;
  if (X->Def.NAnomalousG > 0 &&
      GreenOutputTruncateFile(X, GreenOutputAnomalous) != 0) return -1;
  return 0;
}

int GreenOutputUsesTPQDataAggregate(const struct BindStruct *X)
{
  if (X->Def.iOutputGreenFormat != OUTPUTGREENFORMAT_AGGREGATE) return 0;
  return (X->Def.iCalcType == TPQCalc || X->Def.iCalcType == cTPQ);
}

static const char *GreenOutputTPQDataBaseName(GreenOutputTPQDataKind kind)
{
  switch (kind) {
  case GreenOutputTPQDataSS: return cFileNameSS_TPQ_Aggregate;
  case GreenOutputTPQDataNorm: return cFileNameNorm_TPQ_Aggregate;
  case GreenOutputTPQDataFlct: return cFileNameFlct_TPQ_Aggregate;
  }
  return NULL;
}

int GreenOutputTPQDataFileName(const struct BindStruct *X, GreenOutputTPQDataKind kind, char *sdt)
{
  const char *base = NULL;
  if (!GreenOutputUsesTPQDataAggregate(X)) return -1;
  base = GreenOutputTPQDataBaseName(kind);
  if (base == NULL) return -1;
  if (X->Def.iOutputDataHead == 1) {
    sprintf(sdt, "%s_%s", X->Def.CDataFileHead, base);
  } else {
    sprintf(sdt, "%s", base);
  }
  return 0;
}

static const char *GreenOutputTPQDataHeader(GreenOutputTPQDataKind kind)
{
  switch (kind) {
  case GreenOutputTPQDataSS:
    return " # set, step, inv_tmp, energy, phys_var, phys_doublon, phys_num\n";
  case GreenOutputTPQDataNorm:
    return " # set, step, inv_temp, global_norm, global_1st_norm\n";
  case GreenOutputTPQDataFlct:
    return " # set, step, inv_temp, N, N^2, D, D^2, Sz, Sz^2\n";
  }
  return NULL;
}

static int GreenOutputInitializeTPQDataFile(struct BindStruct *X, GreenOutputTPQDataKind kind)
{
  FILE *fp = NULL;
  char sdt[D_FileNameMax];
  const char *header = NULL;
  if (!GreenOutputUsesTPQDataAggregate(X)) return 0;
  if (GreenOutputTPQDataFileName(X, kind, sdt) != 0) return -1;
  header = GreenOutputTPQDataHeader(kind);
  if (header == NULL) return -1;
  if (childfopenMPI(sdt, "w", &fp) != 0) return -1;
  fprintf(fp, "%s", header);
  fclose(fp);
  return 0;
}

int GreenOutputInitializeTPQDataAggregateFiles(struct BindStruct *X)
{
  if (!GreenOutputUsesTPQDataAggregate(X)) return 0;
  if (GreenOutputInitializeTPQDataFile(X, GreenOutputTPQDataSS) != 0) return -1;
  if (GreenOutputInitializeTPQDataFile(X, GreenOutputTPQDataNorm) != 0) return -1;
  if (GreenOutputInitializeTPQDataFile(X, GreenOutputTPQDataFlct) != 0) return -1;
  return 0;
}

void GreenOutputWriteTPQSSRow(FILE *fp, const struct BindStruct *X, int step, double inv_temp)
{
  if (GreenOutputUsesTPQDataAggregate(X)) {
    fprintf(fp, " %4d %4d %.16lf  %.16lf %.16lf %.16lf %.16lf\n",
            X->Def.irand, step, inv_temp, X->Phys.energy, X->Phys.var,
            X->Phys.doublon, X->Phys.num);
  } else {
    fprintf(fp, "%.16lf  %.16lf %.16lf %.16lf %.16lf %d\n",
            inv_temp, X->Phys.energy, X->Phys.var, X->Phys.doublon,
            X->Phys.num, step);
  }
}

void GreenOutputWriteTPQNormRow(FILE *fp, const struct BindStruct *X, int step, double inv_temp,
                                double norm, double first_norm)
{
  if (GreenOutputUsesTPQDataAggregate(X)) {
    fprintf(fp, " %4d %4d %.16lf %.16lf %.16lf\n",
            X->Def.irand, step, inv_temp, norm, first_norm);
  } else {
    fprintf(fp, "%.16lf %.16lf %.16lf %d\n", inv_temp, norm, first_norm, step);
  }
}

void GreenOutputWriteTPQFlctRow(FILE *fp, const struct BindStruct *X, int step, double inv_temp)
{
  if (GreenOutputUsesTPQDataAggregate(X)) {
    fprintf(fp, " %4d %4d %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf\n",
            X->Def.irand, step, inv_temp, X->Phys.num, X->Phys.num2,
            X->Phys.doublon, X->Phys.doublon2, X->Phys.Sz, X->Phys.Sz2);
  } else {
    fprintf(fp, "%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %d\n",
            inv_temp, X->Phys.num, X->Phys.num2, X->Phys.doublon,
            X->Phys.doublon2, X->Phys.Sz, X->Phys.Sz2, step);
  }
}
