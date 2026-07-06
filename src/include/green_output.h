/* HPhi  -  Quantum Lattice Model Simulator */
#ifndef HPHI_GREEN_OUTPUT_H
#define HPHI_GREEN_OUTPUT_H

#include <stdio.h>
#include "struct.h"

typedef enum {
  GreenOutputOneBody = 0,
  GreenOutputTwoBody,
  GreenOutputThreeBody,
  GreenOutputFourBody,
  GreenOutputSixBody,
  GreenOutputNBody,
  GreenOutputAnomalous
} GreenOutputKind;

typedef enum {
  GreenOutputTPQDataSS = 0,
  GreenOutputTPQDataNorm,
  GreenOutputTPQDataFlct
} GreenOutputTPQDataKind;

int GreenOutputUsesAggregate(const struct BindStruct *X);
int GreenOutputKindUsesAggregate(const struct BindStruct *X, GreenOutputKind kind);
const char *GreenOutputOpenMode(const struct BindStruct *X);
int GreenOutputFileName(const struct BindStruct *X, GreenOutputKind kind, char *sdt);
int GreenOutputWriteIndexPrefix(FILE *fp, const struct BindStruct *X);
int GreenOutputInitializeAggregateFiles(struct BindStruct *X);
int GreenOutputUsesTPQDataAggregate(const struct BindStruct *X);
int GreenOutputTPQDataFileName(const struct BindStruct *X, GreenOutputTPQDataKind kind, char *sdt);
int GreenOutputInitializeTPQDataAggregateFiles(struct BindStruct *X);
void GreenOutputWriteTPQSSRow(FILE *fp, const struct BindStruct *X, int step, double inv_temp);
void GreenOutputWriteTPQNormRow(FILE *fp, const struct BindStruct *X, int step, double inv_temp,
                                double norm, double first_norm);
void GreenOutputWriteTPQFlctRow(FILE *fp, const struct BindStruct *X, int step, double inv_temp);

#endif
