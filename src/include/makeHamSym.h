#ifndef HPHI_MAKEHAMSYM_H
#define HPHI_MAKEHAMSYM_H

struct BindStruct;

/* Fill the existing one-based Ham, or its owned column-major Ham_local panel,
 * from replicated symmetry metadata. Clears the destination on each call.
 * No raw basis, raw diagonal list or vector communication is required.
 * Returns 0 on success, -1 on invalid storage or column-enumeration failure.
 * The caller propagates failures across ranks before diagonalization. */
int makeHamSym(const struct BindStruct *X);

#endif
