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
#include <expec_cisajs.h>
#include <expec_cisajscktaltdc.h>
#include <expec_totalspin.h>
#include "CalcByLanczos.h"
#include "wrapperMPI.h"

/**
 * @file   CalcByLanczos.c
 * @version 0.1, 0.2
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 * @brief  File for givinvg functions of calculating eigenvalues and eigenvectors by Lanczos method 
 * 
 * 
 */


/** 
 * @brief A main function to calculate eigenvalues and eigenvectors by Lanczos method 
 * 
 * @param[in,out] X CalcStruct list for getting and pushing calculation information 
 * @retval 0 normally finished
 * @retval -1 unnormally finished
 *
 * @version 0.2
 * @date 2015/10/20 add function of using a flag of iCalcEigenVec
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 */
int CalcByLanczos(		 
		  struct EDMainCalStruct *X
				 )
{
  char sdt[D_FileNameMax];
  double diff_ene,var;
  int iconv=0;
  int flag=0;
  // this part will be modified
  switch(X->Bind.Def.iCalcModel){
  case HubbardGC:
  case SpinGC:
  case KondoGC:
    initial_mode = 1; // 1 -> random initial vector
    break;
  case Hubbard:
  case Kondo:
  case Spin:
    if(X->Bind.Def.iFlgGeneralSpin ==TRUE){
      initial_mode=1;
    }
    else{
      if(X->Bind.Def.initial_iv>0){
	initial_mode = 0; // 0 -> only v[iv] = 1
      }else{
	initial_mode = 1; // 1 -> random initial vector
      }
    }
    break;
  default:
    return -1;
  }
 
  if(Lanczos_EigenValue(&(X->Bind))!=0){
    fprintf(stdoutMPI, "Lanczos Eigenvalue is not converged in this process.\n");
    return -1;
  }  
  Lanczos_EigenVector(&(X->Bind));
  expec_energy(&(X->Bind));
//check for the accuracy of the eigenvector
  var      = fabs(X->Bind.Phys.var-X->Bind.Phys.energy*X->Bind.Phys.energy)/fabs(X->Bind.Phys.var);
  diff_ene = fabs(X->Bind.Phys.Target_energy-X->Bind.Phys.energy)/fabs(X->Bind.Phys.Target_energy);

  fprintf(stdoutMPI, "\n");
  fprintf(stdoutMPI, "Accuracy check !!!\n");
  fprintf(stdoutMPI, "%.14e %.14e: diff_ere=%.14e var=%.14e \n ",X->Bind.Phys.Target_energy,X->Bind.Phys.energy,diff_ene,var);
  if(diff_ene < eps_Energy && var< eps_Energy){
    fprintf(stdoutMPI, "Accuracy of Lanczos vectors is enough\n");
    fprintf(stdoutMPI, "\n");
  }else{
    fprintf(stdoutMPI, "Accuracy of Lanczos vectors is NOT enough\n");
    iconv=1;
    fprintf(stdoutMPI, "Eigenvector is improved by power Lanczos method \n");
    fprintf(stdoutMPI, "Power Lanczos starts\n");
    flag=PowerLanczos(&(X->Bind));
    fprintf(stdoutMPI, "Power Lanczos ends\n");
    if(flag==1){
      var      = fabs(X->Bind.Phys.var-X->Bind.Phys.energy*X->Bind.Phys.energy)/fabs(X->Bind.Phys.var);
      diff_ene = fabs(X->Bind.Phys.Target_energy-X->Bind.Phys.energy)/fabs(X->Bind.Phys.Target_energy);
      fprintf(stdoutMPI,"\n");
      fprintf(stdoutMPI,"Power Lanczos Accuracy check !!!\n");
      fprintf(stdoutMPI,"%.14e %.14e: diff_ere=%.14e var=%.14e \n ",X->Bind.Phys.Target_energy,X->Bind.Phys.energy,diff_ene,var);
      fprintf(stdoutMPI,"\n");
    }
    else if(X->Bind.Def.iCalcEigenVec==CALCVEC_LANCZOSCG && iconv==1){
      fprintf(stdoutMPI, "Accuracy of Lanczos vectors is NOT enough\n");
      fprintf(stdoutMPI, "Eigenvector is improved by CG method \n");
      X->Bind.Def.St=1;
      CG_EigenVector(&(X->Bind));
      expec_energy(&(X->Bind));
      var      = fabs(X->Bind.Phys.var-X->Bind.Phys.energy*X->Bind.Phys.energy)/fabs(X->Bind.Phys.var);
      diff_ene = fabs(X->Bind.Phys.Target_energy-X->Bind.Phys.energy)/fabs(X->Bind.Phys.Target_energy);
      fprintf(stdoutMPI, "\n");
      fprintf(stdoutMPI, "CG Accuracy check !!!\n");
      fprintf(stdoutMPI, "%.14e %.14e: diff_ere=%.14e var=%.14e \n ",X->Bind.Phys.Target_energy,X->Bind.Phys.energy,diff_ene,var);
      fprintf(stdoutMPI, "\n");
    }
  }
  /*
  expec_cisajs(&(X->Bind),v1);
  // v1 is eigen vector
  expec_cisajscktaltdc(&(X->Bind), v1);
  */
  if(!expec_cisajs(&(X->Bind), v1)==0){
    fprintf(stderr, "Error: calc OneBodyG.\n");
    exitMPI(-1);
  }
  if(!expec_cisajscktaltdc(&(X->Bind), v1)==0){
    fprintf(stderr, "Error: calc TwoBodyG.\n");
    exitMPI(-1);
  }
  if(!expec_totalspin(&(X->Bind), v1)==0){
    fprintf(stderr, "Error: calc TotalSpin.\n");
    exitMPI(-1);
  }
  
    if(X->Bind.Def.St==0){
      sprintf(sdt, cFileNameEnergy_Lanczos, X->Bind.Def.CDataFileHead);
    }else if(X->Bind.Def.St==1){
      sprintf(sdt, cFileNameEnergy_CG, X->Bind.Def.CDataFileHead);
    }
	
    FILE *fp;    
    if(childfopenMPI(sdt, "w", &fp)!=0){
      return -1;
    }  

    fprintf(fp,"Energy  %.10lf \n",X->Bind.Phys.energy);
    fprintf(fp,"Doublon  %.10lf \n",X->Bind.Phys.doublon);
    fprintf(fp,"Sz  %.10lf \n",X->Bind.Phys.sz);
    fprintf(fp,"total S^2  %.10lf \n",X->Bind.Phys.s2);
    
    fclose(fp);
  
  return 0;
}
