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

#ifndef HPHI_MLTPLY_H
#define HPHI_MLTPLY_H

#include "Common.h"

#define M_MLTPLY 0
#define M_ENERGY 1
#define M_Ham 2
#define M_CORR 3
#define M_TOTALS 4
#define M_CALCSPEC 4

int mltply(struct BindStruct *X, double complex *tmp_v0,double complex *tmp_v1);

double complex child_general_hopp_element
(
 const long unsigned int j,
 double complex       *tmp_v0,
 double complex      *tmp_v1,
 struct BindStruct *X
 );

double complex GC_child_general_hopp
(
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 double complex trans
 );


double complex GC_child_general_int(
			 double complex *tmp_v0,
			 double complex *tmp_v1,
			 struct BindStruct *X
			 );


double complex child_general_int
(
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X
 );

double complex GC_child_general_int_spin
(
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X
 );

double complex child_general_int_spin
(
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X
 );


double complex child_general_hopp
(
 double complex       *tmp_v0,
 double complex       *tmp_v1,
 struct BindStruct *X,
 double complex trans
 );

double complex child_exchange
(
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X
 );

double complex child_pairhopp
(
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X
 );

double complex GC_child_exchange
(
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X
 );

double complex GC_child_pairlift
(
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X
 );

double complex GC_child_pairhopp
(
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X
 );


double complex GC_child_exchange_spin
(
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X
 );

double complex child_exchange_spin
(
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X
 );

double complex GC_child_pairlift_spin
(
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X
 );

double complex child_pairlift_spin
(
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X
 );

double complex child_exchange_element
(
 const long unsigned int j,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );

double complex child_CisAisCisAis_element
(
 long unsigned int j,
 long unsigned int isite1,
 long unsigned int isite3,
 double complex tmp_V,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );

double complex child_CisAisCjtAku_element
(
 long unsigned int j,
 long unsigned int isite1,
 long unsigned int isite3,
 long unsigned int isite4,
 long unsigned int Bsum,
 long unsigned int Bdiff,
 double complex tmp_V,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );

double complex child_CisAjtCkuAku_element
(
 long unsigned int j,
 long unsigned int isite1,
 long unsigned int isite2,
 long unsigned int isite3,
 long unsigned int Asum,
 long unsigned int Adiff,
 double complex tmp_V,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );

double complex child_CisAjtCkuAlv_element
(
 long unsigned int j,
 long unsigned int isite1,
 long unsigned int isite2,
 long unsigned int isite3,
 long unsigned int isite4,
 long unsigned int Asum,
 long unsigned int Adiff,
 long unsigned int Bsum,
 long unsigned int Bdiff,
 double complex tmp_V,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off_2
 );
//[s]Grand canonical 
double complex GC_child_CisAisCisAis_element
(
 long unsigned int j,
 long unsigned int isite1,
 long unsigned int isite3,
 double complex tmp_V,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );

double complex GC_child_CisAisCjtAku_element
(
 long unsigned int j,
 long unsigned int isite1,
 long unsigned int isite3,
 long unsigned int isite4,
 long unsigned int Bsum,
 long unsigned int Bdiff,
 double complex tmp_V,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );

double complex GC_child_CisAjtCkuAku_element
(
 long unsigned int j,
 long unsigned int isite1,
 long unsigned int isite2,
 long unsigned int isite3,
 long unsigned int Asum,
 long unsigned int Adiff,
 double complex tmp_V,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );

double complex GC_child_CisAjtCkuAlv_element
(
 long unsigned int j,
 long unsigned int isite1,
 long unsigned int isite2,
 long unsigned int isite3,
 long unsigned int isite4,
 long unsigned int Asum,
 long unsigned int Adiff,
 long unsigned int Bsum,
 long unsigned int Bdiff,
 double complex tmp_V,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off_2
 );

double complex GC_child_general_hopp_element_imp
(
 const long unsigned int j,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );
//[e]Grand canonical 


//[s]Spin
double complex child_CisAisCisAis_spin_element
(
 long unsigned int j,
 long unsigned int isA_up,
 long unsigned int isB_up,
 long unsigned int org_sigma2,
 long unsigned int org_sigma4,
 double complex tmp_V,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X
 );

double complex child_CisAjsCjtAit_spin_element
(
 long unsigned int j,
 long unsigned int isA_up,
 long unsigned int isB_up,
 long unsigned int org_sigma2,
 long unsigned int org_sigma4,
 double complex tmp_V,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );


double complex child_CisAisCitAiu_spin_element
(
 long unsigned int j,
 long unsigned int org_sigma2,
 long unsigned int org_sigma4,
 long unsigned int isA_up,
 long unsigned int isB_up,
 double complex tmp_V,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );

double complex child_CisAitCiuAiu_spin_element
(
 long unsigned int j,
 long unsigned int org_sigma2,
 long unsigned int org_sigma4,
 long unsigned int isA_up,
 long unsigned int isB_up,
 double complex tmp_V,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );

double complex child_CisAitCiuAiv_spin_element
(
 long unsigned int j,
 long unsigned int org_sigma2,
 long unsigned int org_sigma4,
 long unsigned int isA_up,
 long unsigned int isB_up,
 double complex tmp_V,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off_2
 );
//[e]Spin

//[s]GC Spin
double complex GC_child_CisAisCisAis_spin_element
(
 long unsigned int j,
 long unsigned int isA_up,
 long unsigned int isB_up,
 long unsigned int org_sigma2,
 long unsigned int org_sigma4,
 double complex tmp_V,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X
 );

double complex GC_child_CisAisCitAiu_spin_element
(
 long unsigned int j,
 long unsigned int org_sigma2,
 long unsigned int org_sigma4,
 long unsigned int isA_up,
 long unsigned int isB_up,
 double complex tmp_V,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );

double complex GC_child_CisAitCiuAiu_spin_element
(
 long unsigned int j,
 long unsigned int org_sigma2,
 long unsigned int org_sigma4,
 long unsigned int isA_up,
 long unsigned int isB_up,
 double complex tmp_V,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );

double complex GC_child_CisAitCiuAiv_spin_element
(
 long unsigned int j,
 long unsigned int org_sigma2,
 long unsigned int org_sigma4,
 long unsigned int isA_up,
 long unsigned int isB_up,
 double complex tmp_V,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off_2
 );
//[e]GC Spin




//[e]GC Spin
double complex child_pairhopp_element
(
 const long unsigned int j,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );

double complex GC_child_exchange_element
(
 const long unsigned int j,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );

double complex GC_child_pairhopp_element
(
 const long unsigned int j,
 double complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );

double complex child_exchange_spin_element
(
 const long unsigned int j,
 double complex       *tmp_v0,
 double complex       *tmp_v1,
 struct BindStruct *X,
  long unsigned int *tmp_off
 );

double complex GC_child_pairlift_spin_element
(
 const long unsigned int j,
 double complex       *tmp_v0,
 double complex       *tmp_v1,
 struct BindStruct *X, 
 long unsigned int *tmp_off
 );

double complex child_pairlift_spin_element
(
 const long unsigned int j,
 double complex       *tmp_v0,
 double complex       *tmp_v1,
 struct BindStruct *X,
 long unsigned int *tmp_off
 );

double complex GC_child_exchange_spin_element
(
 const long unsigned int j,
 double complex       *tmp_v0,
 double complex       *tmp_v1,
 struct BindStruct *X,  
  long unsigned int *tmp_off 
 );

int X_child_exchange_spin_element
(
 const long unsigned int j,
 struct BindStruct *X,
 const long unsigned int isA_up,
 const long unsigned int isB_up,
 const long unsigned int sigmaA,
 const long unsigned int sigmaB,
 long unsigned int *tmp_off
);



double complex CisAis
(		      
 long unsigned int j,
 double  complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int is1_spin,
 double complex tmp_trans
);


double complex GC_CisAis
(
 long unsigned int j,
 double  complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int is1_spin,
 double complex tmp_trans
);

double complex GC_AisCis(
        long unsigned int j,
        double complex *tmp_v0,
        double complex *tmp_v1,
        struct BindStruct *X,
        long unsigned int is1_spin,
        double complex tmp_trans
);

int X_CisAis
(
 long unsigned int list_1_j,
 struct BindStruct *X,
 long unsigned int is1_spin
 );

int X_CisAjt
(
 long unsigned int list_1_j,
 struct BindStruct *X,
 long unsigned int is1_spin,
 long unsigned int is2_spin,
 long unsigned int sum_spin,
 long unsigned int diff_spin,
 long unsigned int *tmp_off
 );


int X_GC_CisAjt
(
 long unsigned int list_1_j,
 struct BindStruct *X,
 long unsigned int is1_spin,
 long unsigned int is2_spin,
 long unsigned int sum_spin,
 long unsigned int diff_spin,
 long unsigned int *tmp_off
 );


double complex CisAjt
(
 long unsigned int j,
 double  complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int is1_spin,
 long unsigned int is2_spin,
 long unsigned int sum_spin,
 long unsigned int diff_spin,
 double complex tmp_V
 );


double complex GC_CisAjt
(
 long unsigned int j,
 double  complex *tmp_v0,
 double complex *tmp_v1,
 struct BindStruct *X,
 long unsigned int is1_spin,
 long unsigned int is2_spin,
 long unsigned int sum_spin,
 long unsigned int diff_spin,
 double complex tmp_V,
 long unsigned int *tmp_off
);

int child_general_hopp_GetInfo
(
 struct BindStruct *X,
 unsigned long int isite1,
 unsigned long int isite2,
 unsigned long int sigma1,
 unsigned long int sigma2
 );

int child_general_int_GetInfo
(
 const int iInterAll,
 struct BindStruct *X,
 long unsigned int isite1,
 long unsigned int isite2,
 long unsigned int isite3,
 long unsigned int isite4,
 long unsigned int sigma1,
 long unsigned int sigma2,
 long unsigned int sigma3,
 long unsigned int sigma4,
 double complex tmp_V
 );

int child_general_int_spin_GetInfo
(
 struct BindStruct *X,
 long unsigned int isite1,
 long unsigned int isite2,
 long unsigned int sigma1,
 long unsigned int sigma2,
 long unsigned int sigma3,
 long unsigned int sigma4,
 double complex tmp_V
 );


int child_pairhopp_GetInfo
(
 const int iPairHopp,
 struct BindStruct *X 
  );

int child_exchange_GetInfo
(
 const int iExchange,
 struct BindStruct *X 
 );

int child_exchange_spin_GetInfo
(
 const int iExchange,
 struct BindStruct *X 
 );

int child_pairlift_spin_GetInfo
(
 const int iPairLift,
 struct BindStruct *X 
 );

int X_SpinGC_CisAit(
long unsigned int j,
struct BindStruct *X,
long unsigned int is1_spin,
long unsigned int sigma2,
long unsigned int *tmp_off
);

int X_Spin_CisAit(
long unsigned int j,
struct BindStruct *X,
long unsigned int is1_spin,
long unsigned int sigma2,
long unsigned int *tmp_off
);

int X_Spin_CisAis(
long unsigned int j,
struct BindStruct *X,
long unsigned int is1_spin,
long unsigned int sigma1
);

int X_Spin_CisAjs(
long unsigned int j,
struct BindStruct *X,
long unsigned int is1_spin,
long unsigned int is2_spin,
long unsigned int sigma1,
long unsigned int *tmp_off
);


int X_SpinGC_CisAis(
long unsigned int j,
struct BindStruct *X,
long unsigned int is1_spin,
long unsigned int sigma1
);

#endif /* HPHI_MLTPLY_H */
