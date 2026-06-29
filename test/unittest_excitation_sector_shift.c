/* Unit test for the off-diagonal excitation sector-shift truth table
 * (GetExcitationSectorShift / GetExcitationOperatorSetShift).
 * Returns non-zero if any case disagrees with the expected behavior. */
#include "Common.h"
#include "ExcitationSectorShift.h"
#include <stdio.h>

static int g_fail = 0;

static void check_single(const char *name, int model, int gen, int op0, int op1, int op2,
                         int wantValid, int dNe, int dNup, int dNdown, int dT2Sz) {
  int op[3] = {op0, op1, op2};
  SectorShift s = GetExcitationSectorShift(model, gen, FALSE, op);
  int ok = (s.valid == wantValid);
  if (wantValid == TRUE)
    ok = ok && s.dNe == dNe && s.dNup == dNup && s.dNdown == dNdown && s.dTotal2Sz == dT2Sz;
  printf("%-40s : %s\n", name, ok ? "OK" : "FAIL");
  if (!ok) {
    printf("   got valid=%d (dNe=%d dNup=%d dNdown=%d dT2Sz=%d reason=%d)\n",
           s.valid, s.dNe, s.dNup, s.dNdown, s.dTotal2Sz, s.reason);
    g_fail = 1;
  }
}

static void check_pair(const char *name, int model, int gen,
                       int s1, int sp1, int s2, int sp2, int type,
                       int wantValid, int dNe, int dNup, int dNdown, int dT2Sz) {
  int op[5] = {s1, sp1, s2, sp2, type};
  SectorShift s = GetExcitationSectorShift(model, gen, TRUE, op);
  int ok = (s.valid == wantValid);
  if (wantValid == TRUE)
    ok = ok && s.dNe == dNe && s.dNup == dNup && s.dNdown == dNdown && s.dTotal2Sz == dT2Sz;
  printf("%-40s : %s\n", name, ok ? "OK" : "FAIL");
  if (!ok) {
    printf("   got valid=%d (dNe=%d dNup=%d dNdown=%d dT2Sz=%d reason=%d)\n",
           s.valid, s.dNe, s.dNup, s.dNdown, s.dTotal2Sz, s.reason);
    g_fail = 1;
  }
}

static void check_reason(const char *name, int model, int gen, int isPair, int wantReason,
                         int **op, int nOp) {
  SectorShift s = GetExcitationOperatorSetShift(model, gen, isPair, op, nOp);
  int ok = (s.valid == FALSE && s.reason == wantReason);
  printf("%-40s : %s\n", name, ok ? "OK" : "FAIL");
  if (!ok) {
    printf("   got valid=%d reason=%d (want reason=%d)\n", s.valid, s.reason, wantReason);
    g_fail = 1;
  }
}

int main(void) {
  /* --- Hubbard single (spin 0 = up). type 1 = creation, 0 = annihilation --- */
  check_single("Hubbard single c_up^dag (cis)", Hubbard, 0, 0, 0, 1, TRUE, +1, +1, 0, 0);
  check_single("Hubbard single c_up (ajt)",     Hubbard, 0, 0, 0, 0, TRUE, -1, -1, 0, 0);
  check_single("Hubbard single c_dn^dag (cis)", Hubbard, 0, 0, 1, 1, TRUE, +1, 0, +1, 0);
  check_single("Hubbard single c_dn (ajt)",     Hubbard, 0, 0, 1, 0, TRUE, -1, 0, -1, 0);

  /* --- HubbardNConserved single (Ne fixed, 2Sz free): shift tracks dNe ONLY; the Sz components
     are zero for BOTH spins, so a c_up and a c_dn of the same creation/annihilation sense share
     the Ne+-1 sector (this is what lets the cross-spin bra/ket projection work). --- */
  check_single("NConserved single c_up^dag", HubbardNConserved, 0, 0, 0, 1, TRUE, +1, 0, 0, 0);
  check_single("NConserved single c_up",     HubbardNConserved, 0, 0, 0, 0, TRUE, -1, 0, 0, 0);
  check_single("NConserved single c_dn^dag", HubbardNConserved, 0, 0, 1, 1, TRUE, +1, 0, 0, 0);
  check_single("NConserved single c_dn",     HubbardNConserved, 0, 0, 1, 0, TRUE, -1, 0, 0, 0);

  /* --- HubbardNConserved operator-set consistency: cross-spin same-sense operators agree on the
     Ne-only sector (valid), while creation+annihilation mix is rejected. --- */
  { /* up-annih + down-annih (cross spin, same sense) -> consistent, dNe=-1 */
    int r0[3] = {0, 0, 0}, r1[3] = {0, 1, 0}; int *op[2] = {(int*)r0, (int*)r1};
    SectorShift s = GetExcitationOperatorSetShift(HubbardNConserved, 0, FALSE, op, 2);
    int ok = (s.valid == TRUE && s.dNe == -1 && s.dNup == 0 && s.dNdown == 0 && s.dTotal2Sz == 0);
    printf("%-40s : %s\n", "NConserved set cross-spin (consistent)", ok ? "OK" : "FAIL");
    if (!ok) { printf("   got valid=%d dNe=%d dNup=%d dNdown=%d\n", s.valid, s.dNe, s.dNup, s.dNdown); g_fail = 1; }
  }
  { /* up-annih + down-creation (mixed sense) -> SET_INCONSISTENT */
    int r0[3] = {0, 0, 0}, r1[3] = {0, 1, 1}; int *op[2] = {(int*)r0, (int*)r1};
    check_reason("NConserved set mixed sense (inconsistent)", HubbardNConserved, 0, FALSE,
                 OFFDIAG_SHIFT_SET_INCONSISTENT, op, 2);
  }

  /* --- Hubbard pair: off-diagonal spin flips Nup/Ndown; diagonal -> no shift --- */
  check_pair("Hubbard pair n_up (diag)",        Hubbard, 0, 0, 0, 0, 0, 1, TRUE, 0, 0, 0, 0);
  check_pair("Hubbard pair c_up^dag c_dn",      Hubbard, 0, 0, 0, 0, 1, 1, TRUE, 0, +1, -1, 0);
  check_pair("Hubbard pair c_dn^dag c_up",      Hubbard, 0, 0, 1, 0, 0, 1, TRUE, 0, -1, +1, 0);

  /* --- SpinGC pair: full GC space, never shifts --- */
  check_pair("SpinGC pair (diag)",  SpinGC, 0, 0, 0, 0, 0, 1, TRUE, 0, 0, 0, 0);
  check_pair("SpinGC pair (flip)",  SpinGC, 0, 0, 0, 0, 1, 1, TRUE, 0, 0, 0, 0);

  /* --- canonical Spin (half): spin convention opposite to Hubbard --- */
  check_pair("Spin half pair flip(0->1)", Spin, 0, 0, 0, 0, 1, 1, TRUE, 0, -1, +1, 0);
  check_pair("Spin half pair flip(1->0)", Spin, 0, 0, 1, 0, 0, 1, TRUE, 0, +1, -1, 0);
  check_pair("Spin half pair (diag)",     Spin, 0, 0, 0, 0, 0, 1, TRUE, 0, 0, 0, 0);

  /* --- general Spin: Total2Sz = 2*(spin1-spin2) --- */
  check_pair("Spin general flip 0->2", Spin, 1, 0, 0, 0, 2, 1, TRUE, 0, 0, 0, -4);
  check_pair("Spin general diag",      Spin, 1, 0, 1, 0, 1, 1, TRUE, 0, 0, 0, 0);

  /* --- single excitation not allowed for spin --- */
  check_single("Spin single -> invalid",   Spin,   0, 0, 0, 1, FALSE, 0, 0, 0, 0);
  check_single("SpinGC single -> invalid", SpinGC, 0, 0, 0, 1, FALSE, 0, 0, 0, 0);

  /* --- deferred (non-allow-list) models -> invalid (MODEL_NOT_ALLOWED) --- */
  { int row[5] = {0, 0, 0, 0, 1}; int *op[1] = {row};
    check_reason("Kondo pair    -> MODEL_NOT_ALLOWED",   Kondo,   0, TRUE, OFFDIAG_SHIFT_MODEL_NOT_ALLOWED, op, 1);
    check_reason("KondoGC pair  -> MODEL_NOT_ALLOWED",   KondoGC, 0, TRUE, OFFDIAG_SHIFT_MODEL_NOT_ALLOWED, op, 1);
    check_reason("tJ pair       -> MODEL_NOT_ALLOWED",   tJ,      0, TRUE, OFFDIAG_SHIFT_MODEL_NOT_ALLOWED, op, 1);
    check_reason("tJGC pair     -> MODEL_NOT_ALLOWED",   tJGC,    0, TRUE, OFFDIAG_SHIFT_MODEL_NOT_ALLOWED, op, 1);
    check_reason("HubbardGC pair-> MODEL_NOT_ALLOWED",   HubbardGC,0, TRUE, OFFDIAG_SHIFT_MODEL_NOT_ALLOWED, op, 1);
  }

  /* --- operator-set internal consistency --- */
  { /* Hubbard single set: two up annihilations (same shift) -> valid */
    int r0[3] = {0, 0, 0}, r1[3] = {1, 0, 0}; int *op[2] = {(int*)r0, (int*)r1};
    SectorShift s = GetExcitationOperatorSetShift(Hubbard, 0, FALSE, op, 2);
    int ok = (s.valid == TRUE && s.dNe == -1);
    printf("%-40s : %s\n", "Hubbard single set (consistent)", ok ? "OK" : "FAIL");
    if (!ok) { printf("   got valid=%d dNe=%d\n", s.valid, s.dNe); g_fail = 1; }
  }
  { /* Hubbard single set: annihilation + creation (mixed shift) -> SET_INCONSISTENT */
    int r0[3] = {0, 0, 0}, r1[3] = {0, 0, 1}; int *op[2] = {(int*)r0, (int*)r1};
    check_reason("Hubbard single set (inconsistent)", Hubbard, 0, FALSE,
                 OFFDIAG_SHIFT_SET_INCONSISTENT, op, 2);
  }

  printf("\n%s\n", g_fail ? "UNIT TEST FAILED" : "UNIT TEST PASSED");
  return g_fail;
}
