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
  if (kind == GreenOutputAnomalous && X->Def.iCalcType == CG) return FALSE;
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
