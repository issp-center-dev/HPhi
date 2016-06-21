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
			  double dnorm,
  int Nomega,
  double complex *dcSpectrum,
  double complex *dcomega
				 )
{
  char sdt[D_FileNameMax];
  double diff_ene,var;
  unsigned long int i;
  unsigned long int i_max=0;
  FILE *fp;
  int iret=TRUE;
  unsigned long int liLanczosStp=X->Bind.Def.Lanczos_max;
  
  // calculate ai, bi
  // Functions in Lanczos_EigenValue
  fprintf(stdoutMPI, "    Start: Calculate tridiagonal matrix components.\n");
  TimeKeeper(X, cFileNameTimeKeep, c_GetTridiagonalStart, "a");

  iret= Lanczos_GetTridiagonalMatrixComponents( &(X->Bind), alpha, beta, tmp_v1, &(liLanczosStp));
    if(iret != TRUE){
        //Error Message will be added.
        return FALSE;
    }

  fprintf(stdoutMPI, "    End:   Calculate tridiagonal matrix components.\n\n");
  TimeKeeper(X, cFileNameTimeKeep, c_GetTridiagonalEnd, "a");

    fprintf(stdoutMPI, "    Start: Caclulate spectrum from tridiagonal matrix components.\n");
    TimeKeeper(X, cFileNameTimeKeep, c_CalcSpectrumFromTridiagonalStart, "a");

    for( i = 0 ; i < Nomega; i++){
    iret = GetSpectrumByTridiagonalMatrixComponents(alpha, beta, dnorm, dcomega[i], &dcSpectrum[i], liLanczosStp);
    if(iret != TRUE){
      //ToDo: Error Message will be added.
      return FALSE;
    }
    dcSpectrum[i] = dnorm*dcSpectrum[i];
  }
    fprintf(stdoutMPI, "    End:   Caclulate spectrum from tridiagonal matrix components.\n\n");
    TimeKeeper(X, cFileNameTimeKeep, c_CalcSpectrumFromTridiagonalEnd, "a");


  /*
   Print alpha & beta to a file
  */
    sprintf(sdt, cFileNameTridiagonalMatrixComponents, X->Bind.Def.CDataFileHead);
    childfopenMPI(sdt,"w", &fp);
    fprintf(fp, "%ld \n",liLanczosStp);
    fprintf(fp, "%.10lf \n",dnorm);
    for( i = 1 ; i <= liLanczosStp; i++){
        fprintf(fp,"%.10lf %.10lf \n", creal(alpha[i]), beta[i]);
    }
    fclose(fp);

  return TRUE;
}


int GetSpectrumByTridiagonalMatrixComponents(
		double *tmp_alpha,
		double *tmp_beta,
        double dnorm,
		double complex dcomega,
		double complex *dcSpectrum,
        unsigned long int ilLanczosStp
		)
{
    unsigned long int istp=2;
    double complex dcDn;
    double complex dcb0;
    double complex dcbn, dcan;
    double complex dcDeltahn;
    double complex dch;

    if(ilLanczosStp < 1){
        //TODO: Add error message
        return FALSE;
    }

    dcb0 = dcomega-tmp_alpha[1];
    if(ilLanczosStp ==1) {
        if(cabs(dcb0)<eps_Energy){
            dcb0=eps_Energy;
        }
        *dcSpectrum = dnorm / (dcb0);
        return TRUE;
    }

    dcbn = dcomega-tmp_alpha[2];
    dcan = -pow(tmp_beta[1],2);
    dcDn=1.0/dcbn;
    dcDeltahn = dcan*dcDn;
    dch=dcb0+dcDeltahn;

    for(istp=2; istp<=ilLanczosStp; istp++){
        dcbn = dcomega-tmp_alpha[istp+1];
        dcan =-pow(tmp_beta[istp],2);
        dcDn = (dcbn+dcan*dcDn);
        if(cabs(dcDn)<eps_Energy){
            dcDn=eps_Energy;
        }
        dcDn =1.0/dcDn;
        dcDeltahn = (dcbn*dcDn-1.0)*dcDeltahn;
        dch += dcDeltahn;
        if(cabs(dcDeltahn/dch)<cabs(dcb0)*eps) break;
    }
    *dcSpectrum = dnorm/(dch);
  return TRUE;
}
