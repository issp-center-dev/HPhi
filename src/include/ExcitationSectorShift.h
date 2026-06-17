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
#pragma once

/**
 * @brief Sector shift (dNe, dNup, dNdown, dTotal2Sz) induced by an excitation
 *        operator, used to guard off-diagonal dynamical Green-function input so
 *        that the ket A|phi> and bra B|phi> land in the same excited Hilbert
 *        sector. @c valid is FALSE when the (model, operator) excitation is not
 *        supported in the current (v1) off-diagonal allow-list.
 */
/** @brief Why a SectorShift is invalid (for off-diagonal guard diagnostics). */
typedef enum {
  OFFDIAG_SHIFT_OK = 0,            /**< valid shift */
  OFFDIAG_SHIFT_MODEL_NOT_ALLOWED, /**< model/operator not in the off-diagonal allow-list */
  OFFDIAG_SHIFT_SET_INCONSISTENT   /**< rows in the operator set induce different shifts */
} SectorShiftReason;

typedef struct {
  int dNe;
  int dNup;
  int dNdown;
  int dTotal2Sz;
  int valid;
  SectorShiftReason reason;
} SectorShift;

/**
 * @brief Net sector shift of a single excitation operator row.
 * @param op single = {site, spin(0=up/1=down), type(1=creation 'cis', else annihilation 'ajt')}.
 *           pair   = {site1, spin1, site2, spin2, type}.
 * The table mirrors the *net* behavior of MakeExcitedList() in CalcSpectrum.c.
 */
SectorShift GetExcitationSectorShift(int iCalcModel, int isGeneralSpin, int isPair, const int *op);

/**
 * @brief Validate that every operator row in a set induces the same valid shift.
 *        Returns valid=FALSE if any row is unsupported or rows disagree.
 */
SectorShift GetExcitationOperatorSetShift(int iCalcModel, int isGeneralSpin, int isPair,
                                          int **op, int nOp);
