#ifndef HPHI_SYMMETRY_BASIS_IO_H
#define HPHI_SYMMETRY_BASIS_IO_H

#include "Common.h"

struct DefineList;
struct BindStruct;

int ReadTransSymNInt(const char *defname, struct DefineList *def);
int ReadTransSymFile(const char *defname, struct DefineList *def);
int ValidateSymmetryRuntimeOptions(const struct BindStruct *X);
int ValidateSymmetryHamiltonian(const struct BindStruct *X);

#endif /* HPHI_SYMMETRY_BASIS_IO_H */
