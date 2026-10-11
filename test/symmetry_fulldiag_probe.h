#ifndef HPHI_TEST_SYMMETRY_FULLDIAG_PROBE_H
#define HPHI_TEST_SYMMETRY_FULLDIAG_PROBE_H
struct BindStruct;
/* Test-only sector FullDiag matrix/state contracts. Returns 1 when the
 * action was handled successfully and -1 on failure (on every rank). */
int SymmetryFullDiagProbe(struct BindStruct *X, const char *action);
#endif
