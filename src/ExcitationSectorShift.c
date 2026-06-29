/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 Takahiro Misawa, Kazuyoshi Yoshimi, Mitsuaki Kawamura, Youhei Yamaji, Synge Todo, Naoki Kawashima */

/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>. */

/**@file
 * @brief Sector-shift truth table for off-diagonal dynamical Green functions.
 *
 * The shift returned here mirrors the *net* sector change applied by
 * MakeExcitedList() (src/CalcSpectrum.c). Note that MakeExcitedList's second
 * switch runs only when iFlgListModifed==TRUE, so several GC/NConserved
 * branches inside it are dead code. The net behavior is:
 *   - single: HubbardGC -> none; Spin/SpinGC -> N/A;
 *             Hubbard/Kondo/KondoGC/tJ/tJGC/(Hubbard|tJ)NConserved -> shift.
 *   - pair:   all GC (incl. KondoGC/tJGC) + NConserved -> none;
 *             Hubbard/Kondo/tJ -> off-diagonal-spin shift; Spin -> shift.
 *
 * v1 off-diagonal allow-list = {Hubbard, Spin, SpinGC}. Models outside the
 * allow-list return valid=FALSE here (the allow-list guard rejects them); their
 * exact shift is encoded together with a dedicated test when a model is
 * promoted into the allow-list.
 */
#include "Common.h"
#include "ExcitationSectorShift.h"

SectorShift GetExcitationSectorShift(int iCalcModel, int isGeneralSpin, int isPair, const int *op) {
  SectorShift s = {0, 0, 0, 0, TRUE, OFFDIAG_SHIFT_OK};

  if (isPair == FALSE) {
    const int spin = op[1];                /* 0 = up, 1 = down */
    const int isCreation = (op[2] == 1);   /* type==1 -> cis (creation), else ajt (annihilation) */
    switch (iCalcModel) {
    case Hubbard: /* v1 allow-list (CalcSpectrum.c:537-562; spin 0 == up) */
      if (isCreation) {
        s.dNe = +1;
        if (spin == 0) s.dNup = +1; else s.dNdown = +1;
      } else {
        s.dNe = -1;
        if (spin == 0) s.dNup = -1; else s.dNdown = -1;
      }
      break;
    case Spin:
    case SpinGC: /* single excitation N/A for spin (CalcSpectrum.c:563-565, SingleEx.c:54) */
      s.valid = FALSE;
      break;
    case HubbardGC:
      /* Grand canonical: the full Fock space is one Hilbert space, so a single c/c^dag
         excitation stays in it -- no sector shift (the doc's "single: HubbardGC -> none").
         s keeps the default {0,0,0,0, valid=TRUE}. This lets the bra/ket off-diagonal path
         (SpectrumNumBra) project ket and bra -- of either spin -- in the shared GC space. */
      break;
    case HubbardNConserved:
      /* Ne conserved, 2Sz free: c/c^dag shifts the sector by +-1 in Ne ONLY (the excited list
         is Ne+-1 with all 2Sz). Track dNe alone and leave the Sz components zero, so the
         bra/ket and cross-operator sector checks match same- AND cross-spin operators (both
         spins land in the same Ne+-1 space) -- the spin-orbit-capable off-diagonal route. */
      if (isCreation) s.dNe = +1;
      else            s.dNe = -1;
      break;
    case Kondo:
    case KondoGC:
    case tJ:
    case tJGC:
    case tJNConserved:
    case KondoNConserved:
      /* deferred: not in v1 allow-list (rejected upstream by the allow-list guard) */
      s.valid = FALSE;
      break;
    default:
      s.valid = FALSE;
      break;
    }
  } else {
    const int spin1 = op[1];
    const int spin2 = op[3];
    switch (iCalcModel) {
    case Hubbard: /* v1 allow-list (CalcSpectrum.c:576-590; spin 0 == up) */
      if (spin1 != spin2) {
        if (spin1 == 0) { s.dNup = +1; s.dNdown = -1; }
        else            { s.dNup = -1; s.dNdown = +1; }
      }
      break;
    case SpinGC: /* v1 allow-list: full grand-canonical space, no shift */
      break;
    case Spin: /* v1 allow-list; NOTE: spin convention is opposite to Hubbard (CalcSpectrum.c:591-606) */
      if (spin1 != spin2) {
        if (isGeneralSpin == FALSE) {
          if (spin1 == 0) { s.dNup = -1; s.dNdown = +1; }  /* CalcSpectrum.c:594 //down */
          else            { s.dNup = +1; s.dNdown = -1; }  /* CalcSpectrum.c:597 //up   */
        } else {
          s.dTotal2Sz = 2 * (spin1 - spin2);               /* CalcSpectrum.c:603 */
        }
      }
      break;
    case Kondo:
    case KondoGC:
    case tJ:
    case tJGC:
    case HubbardGC:
    case HubbardNConserved:
    case tJNConserved:
    case KondoNConserved:
      /* deferred: not in v1 allow-list (rejected upstream by the allow-list guard) */
      s.valid = FALSE;
      break;
    default:
      s.valid = FALSE;
      break;
    }
  }
  /* A single operator row can only be invalid because the model/operator is
     not in the allow-list (set-internal inconsistency needs >= 2 rows). */
  s.reason = (s.valid == TRUE) ? OFFDIAG_SHIFT_OK : OFFDIAG_SHIFT_MODEL_NOT_ALLOWED;
  return s;
}

SectorShift GetExcitationOperatorSetShift(int iCalcModel, int isGeneralSpin, int isPair,
                                          int **op, int nOp) {
  SectorShift acc = {0, 0, 0, 0, FALSE, OFFDIAG_SHIFT_MODEL_NOT_ALLOWED};
  int i;

  if (nOp <= 0) return acc; /* empty set: valid=FALSE */

  for (i = 0; i < nOp; i++) {
    SectorShift cur = GetExcitationSectorShift(iCalcModel, isGeneralSpin, isPair, op[i]);
    if (cur.valid == FALSE) {
      acc.valid = FALSE;
      acc.reason = cur.reason; /* model/operator not allowed */
      return acc;
    }
    if (i == 0) {
      acc = cur;
    } else if (acc.dNe != cur.dNe || acc.dNup != cur.dNup ||
               acc.dNdown != cur.dNdown || acc.dTotal2Sz != cur.dTotal2Sz) {
      acc.valid = FALSE; /* rows disagree on sector shift */
      acc.reason = OFFDIAG_SHIFT_SET_INCONSISTENT;
      return acc;
    }
  }
  return acc; /* valid=TRUE; shift is common to all rows */
}
