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
#include "CalcSpectrumByLanczos.h"
#include "Lanczos_EigenValue.h"
#include "FileIO.h"
#include "mfmemory.h"
#include "wrapperMPI.h"

/**
 * @file   CalcSpectrumByLanczos.c
 * @version 1.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 * @brief  File for givinvg functions of calculating spectrum by Lanczos
 * 
 * 
 */

/** 
 * @brief A main function to calculate spectrum by Lanczos
 * 
 * @param[in,out] X CalcStruct list for getting and pushing calculation information 
 * @retval 0 normally finished
 * @retval -1 unnormally finished
 *
 * @version 1.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 */
int CalcSpectrumByLanczos(		 
			  struct EDMainCalStruct *X,
			  double complex *tmp_v1,
			  double dnorm
				 )
{
  char sdt[D_FileNameMax];
  double diff_ene,var;
  unsigned long int i;
  unsigned long int i_max=0;
  FILE *fp;
  int iret=TRUE;
    unsigned long int liLanczosStp=X->Bind.Def.Lanczos_max;

  //ToDo: Nomega should be given as a parameter
  int Nomega=1000;
  double complex OmegaMax=20+0.01*I;
  double complex OmegaMin=-20+0.01*I;

  double complex *dcSpectrum;
  c_malloc1(dcSpectrum, Nomega);
  double complex *dcomega;
  c_malloc1(dcomega, Nomega);
  
  // calculate ai, bi
  // Functions in Lanczos_EigenValue
  iret= Lanczos_GetTridiagonalMatrixComponents( &(X->Bind), alpha, beta, v1,  &(liLanczosStp));
  if(iret != TRUE){
    //Error Message will be added.

    return FALSE;
  }

  // ToDo: Give dcomega
  for(i=0; i< Nomega; i++){
    dcomega[i]=(OmegaMax-OmegaMin)/Nomega+OmegaMin;
  }
  
  for( i = 0 ; i < Nomega; i++){
    //ToDo: calculate spectrum for a fixed omega
    iret = GetSpectrumByTridiagonalMatrixComponents(alpha, beta, dnorm, dcomega[i], &dcSpectrum[i], liLanczosStp);
    if(iret != TRUE){
      //Error Message will be added.
      return FALSE;
    }
  }
  
  //output spectrum
  sprintf(sdt, cFileNameLanczosStep, X->Bind.Def.CDataFileHead);   
  childfopenMPI(sdt,"w", &fp);
  for( i = 0 ; i < Nomega; i++){
    fprintf(fp, "%.10lf %.10lf %.10lf %.10lf \n",
	    creal(dcomega[i]), cimag(dcomega[i]),
	    creal(dcSpectrum[i]), cimag(dcSpectrum[i]));
  }
  fclose(fp);
  
  return TRUE;
}


int GetSpectrumByTridiagonalMatrixComponents(
		double *tmp_alpha, // alpha: 1,..., ilLanczosStp
		double *tmp_beta, // beta: 1, ..., ilLanczosStp
		double dnorm,
		double complex dcomega,
		double complex *dcSpectrum,
        unsigned long int ilLanczosStp
		)
{
  return TRUE;
}
