/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 The University of Tokyo */

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
/**
 * @file rearray_interactions.c
 *
 * @brief Two-body Green's-function operator reordering. The definition was
 * moved here VERBATIM from src/expec_cisajscktaltdc.c (phase 3b Task 2) so
 * that the trace-map unit test (test/unit/expec_trace_map_check.c) links the
 * real function; behavior is unchanged. This file performs no MPI and reads
 * only X->Def.{CisAjtCkuAlvDC,TBody,FBody,SBody}.
 */
#include "rearray_interactions.h"

///
/// \brief Rearray interactions
/// \param i
/// \param org_isite1 a site number on the site 1.
/// \param org_isite2 a site number on the site 2.
/// \param org_isite3 a site number on the site 3.
/// \param org_isite4 a site number on the site 4.
/// \param org_sigma1 a spin index on the site 1.
/// \param org_sigma2 a spin index on the site 2.
/// \param org_sigma3 a spin index on the site 3.
/// \param org_sigma4 a spin index on the site 4.
/// \param tmp_V a value of interaction
/// \param X  data list for calculation
/// \return 0 normally finished
/// \return -1 unnormally finished
int Rearray_Interactions(
                         int i,
                         long unsigned int *org_isite1,
                         long unsigned int *org_isite2,
                         long unsigned int *org_isite3,
                         long unsigned int *org_isite4,
                         long unsigned int *org_sigma1,
                         long unsigned int *org_sigma2,
                         long unsigned int *org_sigma3,
                         long unsigned int *org_sigma4,
                         double complex *tmp_V,
                         struct BindStruct *X,
                         int type
                         )
{
  long unsigned int tmp_org_isite1,tmp_org_isite2,tmp_org_isite3,tmp_org_isite4;
  long unsigned int tmp_org_sigma1,tmp_org_sigma2,tmp_org_sigma3,tmp_org_sigma4;

  if(type==6){
    tmp_org_isite1   = X->Def.SBody[i][0]+1;
    tmp_org_sigma1   = X->Def.SBody[i][1];
    tmp_org_isite2   = X->Def.SBody[i][2]+1;
    tmp_org_sigma2   = X->Def.SBody[i][3];
    tmp_org_isite3   = X->Def.SBody[i][4]+1;
    tmp_org_sigma3   = X->Def.SBody[i][5];
    tmp_org_isite4   = X->Def.SBody[i][6]+1;
    tmp_org_sigma4   = X->Def.SBody[i][7];
  }else if(type==4){
    tmp_org_isite1   = X->Def.FBody[i][0]+1;
    tmp_org_sigma1   = X->Def.FBody[i][1];
    tmp_org_isite2   = X->Def.FBody[i][2]+1;
    tmp_org_sigma2   = X->Def.FBody[i][3];
    tmp_org_isite3   = X->Def.FBody[i][4]+1;
    tmp_org_sigma3   = X->Def.FBody[i][5];
    tmp_org_isite4   = X->Def.FBody[i][6]+1;
    tmp_org_sigma4   = X->Def.FBody[i][7];
  }else if(type==3){
    tmp_org_isite1   = X->Def.TBody[i][0]+1;
    tmp_org_sigma1   = X->Def.TBody[i][1];
    tmp_org_isite2   = X->Def.TBody[i][2]+1;
    tmp_org_sigma2   = X->Def.TBody[i][3];
    tmp_org_isite3   = X->Def.TBody[i][4]+1;
    tmp_org_sigma3   = X->Def.TBody[i][5];
    tmp_org_isite4   = X->Def.TBody[i][6]+1;
    tmp_org_sigma4   = X->Def.TBody[i][7];
  }else{
    tmp_org_isite1   = X->Def.CisAjtCkuAlvDC[i][0]+1;
    tmp_org_sigma1   = X->Def.CisAjtCkuAlvDC[i][1];
    tmp_org_isite2   = X->Def.CisAjtCkuAlvDC[i][2]+1;
    tmp_org_sigma2   = X->Def.CisAjtCkuAlvDC[i][3];
    tmp_org_isite3   = X->Def.CisAjtCkuAlvDC[i][4]+1;
    tmp_org_sigma3   = X->Def.CisAjtCkuAlvDC[i][5];
    tmp_org_isite4   = X->Def.CisAjtCkuAlvDC[i][6]+1;
    tmp_org_sigma4   = X->Def.CisAjtCkuAlvDC[i][7];
  }

  if(tmp_org_isite1==tmp_org_isite2 && tmp_org_isite3==tmp_org_isite4){
    if(tmp_org_isite1 > tmp_org_isite3){
      *org_isite1   = tmp_org_isite3;
      *org_sigma1   = tmp_org_sigma3;
      *org_isite2   = tmp_org_isite4;
      *org_sigma2   = tmp_org_sigma4;
      *org_isite3   = tmp_org_isite1;
      *org_sigma3   = tmp_org_sigma1;
      *org_isite4   = tmp_org_isite2;
      *org_sigma4   = tmp_org_sigma2;
    }
    else{
      *org_isite1   = tmp_org_isite1;
      *org_sigma1   = tmp_org_sigma1;
      *org_isite2   = tmp_org_isite2;
      *org_sigma2   = tmp_org_sigma2;
      *org_isite3   = tmp_org_isite3;
      *org_sigma3   = tmp_org_sigma3;
      *org_isite4   = tmp_org_isite4;
      *org_sigma4   = tmp_org_sigma4;
    }
    *tmp_V = 1.0;

  }
  else if(tmp_org_isite1==tmp_org_isite4 && tmp_org_isite3==tmp_org_isite2){
    if(tmp_org_isite1 > tmp_org_isite3){
      *org_isite1   = tmp_org_isite3;
      *org_sigma1   = tmp_org_sigma3;
      *org_isite2   = tmp_org_isite2;
      *org_sigma2   = tmp_org_sigma2;
      *org_isite3   = tmp_org_isite1;
      *org_sigma3   = tmp_org_sigma1;
      *org_isite4   = tmp_org_isite4;
      *org_sigma4   = tmp_org_sigma4;
    }
    else{
      *org_isite1   = tmp_org_isite1;
      *org_sigma1   = tmp_org_sigma1;
      *org_isite2   = tmp_org_isite4;
      *org_sigma2   = tmp_org_sigma4;
      *org_isite3   = tmp_org_isite3;
      *org_sigma3   = tmp_org_sigma3;
      *org_isite4   = tmp_org_isite2;
      *org_sigma4   = tmp_org_sigma2;
    }
    *tmp_V =-1.0;
  }
  else{
    return -1;
  }
  return 0;
}
