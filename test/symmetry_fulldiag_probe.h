#ifndef HPHI_TEST_SYMMETRY_FULLDIAG_PROBE_H
#define HPHI_TEST_SYMMETRY_FULLDIAG_PROBE_H
struct BindStruct;
/* Test-only sector FullDiag matrix/state contracts. Returns 1 when the
 * action was handled successfully and -1 on failure (on every rank). */
int SymmetryFullDiagProbe(struct BindStruct *X, const char *action);
/* Fault injection: 1 only for the first call whose stage matches
 * HPHI_TEST_FD_FAIL on the rank named by HPHI_TEST_FD_FAIL_RANK. Production
 * sources reference it only under HPHI_TEST_SYMMETRY_PROBE. */
int SymmetryFullDiagTestFail(const char *stage);
#endif
