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
#include "mltply.h"
#include "CalcSpectrum.h"
#include "CalcSpectrumByLanczos.h"
#include "wrapperMPI.h"

/**
 * @file   CalcSpectrum.c
 * @version 1.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 * @brief  File for givinvg functions of calculating spectrum 
 * 
 * 
 */

/** 
 * @brief A main function to calculate spectrum 
 * 
 * @param[in,out] X CalcStruct list for getting and pushing calculation information 
 * @retval 0 normally finished
 * @retval -1 unnormally finished
 *
 * @version 1.1
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 */
int CalcSpectrum(		 
		 struct EDMainCalStruct *X
				 )
{
  char sdt[D_FileNameMax];
  double diff_ene,var;
  double complex cdnorm;
  unsigned long int i;
  unsigned long int i_max=0;
  FILE *fp;
  double dnorm;

  //input eigen vector
  fprintf(stdoutMPI, "An Eigenvector is inputted.\n");
  sprintf(sdt, cFileNameInputEigen, X->Bind.Def.CDataFileHead, X->Bind.Def.k_exct-1, myrank);
  fp = fopen(sdt, "rb");
  if(fp==NULL){
    fprintf(stderr, "Error: A file of Inputvector does not exist.\n");
    fclose(fp);
    exitMPI(-1);
  }
  fread(&i_max, sizeof(long int), 1, fp);
  if(i_max != X->Bind.Check.idim_max){
    fprintf(stderr, "Error: A file of Inputvector is incorrect.\n");
    fclose(fp);
    exitMPI(-1);
  }
  fread(v1, sizeof(complex double),X->Bind.Check.idim_max+1, fp);
  fclose(fp);
  
  //mltply Operator by using mltply.c (not yet modified) 
  mltply(&(X->Bind), v0, v1);
  
  //calculate norm
  dnorm = NormMPI_dc(i_max, v0);
  
  //normalize vector
#pragma omp parallel for default(none) private(i) shared(v1, v0) firstprivate(i_max, dnorm)
  for(i=1;i<=i_max;i++){
    v1[i] = v0[i]/dnorm;
  }

  int CalcSpecByLanczos=0;
  int iCalcSpecType=CalcSpecByLanczos;
  int iret=TRUE;
  switch (iCalcSpecType){
  case 0:
    iret = CalcSpectrumByLanczos(X, v1, dnorm);
    if(iret != TRUE){
      //Error Message will be added.
      return -1;
    }    
    break;    
    // case CalcSpecByShiftedKlyrov will be added
  default:
    break;
  }

  return 0;
}
