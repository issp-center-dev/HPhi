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

int GreenOutputUsesAggregate(const struct BindStruct *X);
int GreenOutputKindUsesAggregate(const struct BindStruct *X, GreenOutputKind kind);
const char *GreenOutputOpenMode(const struct BindStruct *X);
int GreenOutputFileName(const struct BindStruct *X, GreenOutputKind kind, char *sdt);
int GreenOutputWriteIndexPrefix(FILE *fp, const struct BindStruct *X);
int GreenOutputInitializeAggregateFiles(struct BindStruct *X);

#endif
