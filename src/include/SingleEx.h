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
#include "Common.h"

  double complex GC_Ajt(
          long unsigned int j,
          double complex *tmp_v0,
          double complex *tmp_v1,
          long unsigned int is1_spin,
          double complex tmp_V,
          long unsigned int *tmp_off
  );

  double complex GC_Cis(
          long unsigned int j,
          double complex *tmp_v0,
          double complex *tmp_v1,
          long unsigned int is1_spin,
          double complex tmp_V,
          long unsigned int *tmp_off
  );

  double complex X_GC_Cis_MPI(
				       int org_isite,
				       int org_ispin,
				       double complex tmp_trans,
  double complex *tmp_v0,
  double complex *tmp_v1,
  unsigned long int idim_max,
  double complex *tmp_v1buf,
  long int *Tpow 
  );

  double complex X_GC_Ajt_MPI(
				       int org_isite,
				       int org_ispin,
				       double complex tmp_trans,
  double complex *tmp_v0,
  double complex *tmp_v1,
  unsigned long int idim_max,
  double complex *tmp_v1buf,
  long int *Tpow 
  ); 


int GetSingleExcitedState
(
 struct BindStruct *X,
 double complex *tmp_v0, /**< [out] Result v0 = H v1*/
  double complex *tmp_v1 /**< [in] v0 = H v1*/
 );
