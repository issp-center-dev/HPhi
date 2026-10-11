#ifndef HPHI_SYMMETRY_FULLDIAG_H
#define HPHI_SYMMETRY_FULLDIAG_H

struct BindStruct;
struct SymmetryHamiltonianInfo;

/* Collective. After a successful sector FullDiag, takes eigenstates
 * 0..D-1 in energy order one at a time, gives every rank its owned block
 * of the state, and writes the requested sector checkpoints. Solver
 * dimensions and the distributed eigenvectors are left to the caller.
 * Returns 0 on every rank or -1 on every rank. */
int ProcessSymmetryFullDiagStates(const struct BindStruct *X,
    const struct SymmetryHamiltonianInfo *info);

#endif
