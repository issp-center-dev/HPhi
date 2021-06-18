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
#include "expec_cisajs.h"
#include "expec_cisajscktaltdc.h"
#include "expec_totalspin.h"
#include "CG_EigenVector.h"
#include "expec_energy_flct.h"
#include "Lanczos_EigenValue.h"
#include "Lanczos_EigenVector.h"
#include "CalcByLanczos.h"
#include "FileIO.h"
#include "wrapperMPI.h"
#include "CalcTime.h"
#include "mltply.h"

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
  long int i_max=0;
  long int i;
  double dnorm;
  FILE *fp;
  size_t byte_size;

  if(X->Bind.Def.iInputEigenVec==FALSE){
    // this part will be modified
    switch(X->Bind.Def.iCalcModel){
    case HubbardGC:
    case SpinGC:
    case KondoGC:
    case SpinlessFermionGC:
      initial_mode = 1; // 1 -> random initial vector
      break;
    case Hubbard:
    case Kondo:
    case Spin:
    case SpinlessFermion:

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
      exitMPI(-1);
    }

    StartTimer(4100);
    int iret=0;
    iret=Lanczos_EigenValue(&(X->Bind));
    StopTimer(4100);
    if(iret == 1) return(TRUE);//restart mode.
    else if(iret != 0) return(FALSE);

    if(X->Bind.Def.iCalcEigenVec==CALCVEC_NOT){
       fprintf(stdoutMPI, "  Lanczos EigenValue = %.10lf \n ",X->Bind.Phys.Target_energy);
       return(TRUE);
    }

    fprintf(stdoutMPI, "%s", cLogLanczos_EigenVecStart);

    if(X->Bind.Check.idim_maxMPI != 1){
      StartTimer(4200);
      Lanczos_EigenVector(&(X->Bind));
      StopTimer(4200);

      StartTimer(4300);
      iret=expec_energy_flct(&(X->Bind));
      StopTimer(4300);
      if(iret != 0) return(FALSE);

      //check for the accuracy of the eigenvector
      var      = fabs(X->Bind.Phys.var-X->Bind.Phys.energy*X->Bind.Phys.energy)/fabs(X->Bind.Phys.var);
      diff_ene = fabs(X->Bind.Phys.Target_CG_energy-X->Bind.Phys.energy)/fabs(X->Bind.Phys.Target_CG_energy);

      fprintf(stdoutMPI, "\n");
      fprintf(stdoutMPI, "  Accuracy check !!!\n");
      fprintf(stdoutMPI, "  LanczosEnergy = %.14e \n  EnergyByVec   = %.14e \n  diff_ene      = %.14e \n  var           = %.14e \n",X->Bind.Phys.Target_CG_energy,X->Bind.Phys.energy,diff_ene,var);
      if(diff_ene < eps_Energy && var< eps_Energy){
        fprintf(stdoutMPI, "  Accuracy of Lanczos vectors is enough.\n");
        fprintf(stdoutMPI, "\n");
      }
      else if(X->Bind.Def.iCalcEigenVec==CALCVEC_LANCZOSCG){
        fprintf(stdoutMPI, "  Accuracy of Lanczos vectors is NOT enough\n\n");
        X->Bind.Def.St=1;
        StartTimer(4400);
        iret=CG_EigenVector(&(X->Bind));
        StopTimer(4400);
        if(iret != 0) return(FALSE);

        StartTimer(4300);
        iret=expec_energy_flct(&(X->Bind));
        StopTimer(4300);
        if(iret != 0) return(FALSE);

        var      = fabs(X->Bind.Phys.var-X->Bind.Phys.energy*X->Bind.Phys.energy)/fabs(X->Bind.Phys.var);
        diff_ene = fabs(X->Bind.Phys.Target_CG_energy-X->Bind.Phys.energy)/fabs(X->Bind.Phys.Target_CG_energy);
        fprintf(stdoutMPI, "\n");
        fprintf(stdoutMPI, "  CG Accuracy check !!!\n");
        fprintf(stdoutMPI, "  LanczosEnergy = %.14e\n  EnergyByVec   = %.14e\n  diff_ene      = %.14e\n  var           = %.14e \n ",X->Bind.Phys.Target_CG_energy,X->Bind.Phys.energy,diff_ene,var);
        fprintf(stdoutMPI, "\n");
        //}
      }
    }
    else{//idim_max=1
      v0[1]=1;
      StartTimer(4300);
      iret=expec_energy_flct(&(X->Bind));
      StopTimer(4300);
      if(iret != 0) return(FALSE);
    }
  }
  else if(X->Bind.Def.iInputEigenVec==1){// X->Bind.Def.iInputEigenVec=true :input v1:
    fprintf(stdoutMPI, "An Eigenvector is inputted.\n");
    StartTimer(4800);
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, cReadEigenVecStart, "a");
    StartTimer(4801);
    sprintf(sdt, cFileNameInputEigen, X->Bind.Def.CDataFileHead, X->Bind.Def.k_exct-1, myrank);
    childfopenALL(sdt, "rb", &fp);
    if(fp==NULL){
      fprintf(stderr, "Error: A file of Inputvector does not exist.\n");
      exitMPI(-1);
    }
    byte_size = fread(&step_i, sizeof(int), 1, fp);
    byte_size = fread(&i_max, sizeof(long int), 1, fp);
    if(i_max != X->Bind.Check.idim_max){
      fprintf(stderr, "Error: A file of Inputvector is incorrect.\n");
      exitMPI(-1);
    }
    byte_size = fread(v1, sizeof(complex double),X->Bind.Check.idim_max+1, fp);
    fclose(fp);
    StopTimer(4801);
    StopTimer(4800);
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, cReadEigenVecFinish, "a");
    if (byte_size == 0) printf("byte_size: %d \n", (int)byte_size);
  }
  else if(X->Bind.Def.iInputEigenVec==2){// X->Bind.Def.iInputEigenVec=true :input v1:
    fprintf(stdoutMPI, "An Eigenvector is inputted.\n");
    StartTimer(4800);
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, cReadEigenVecStart, "a");
    StartTimer(4801);
    sprintf(sdt, cFileNameInputEigen, X->Bind.Def.CDataFileHead, X->Bind.Def.k_exct-1, myrank);
    childfopenALL(sdt, "rb", &fp);
    if(fp==NULL){
      fprintf(stderr, "Error: A file of Inputvector does not exist.\n");
      exitMPI(-1);
    }
    byte_size = fread(&step_i, sizeof(int), 1, fp);
    byte_size = fread(&i_max, sizeof(long int), 1, fp);
    if(i_max != X->Bind.Check.idim_max){
      fprintf(stderr, "Error: A file of Inputvector is incorrect.\n");
      exitMPI(-1);
    }
    byte_size = fread(v1, sizeof(complex double),X->Bind.Check.idim_max+1, fp);
    fclose(fp);
    StopTimer(4801);
    StopTimer(4800);
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, cReadEigenVecFinish, "a");
    if (byte_size == 0) printf("byte_size: %d \n", (int)byte_size);
    /*[s]mltply*/
    mltply(&(X->Bind), v0, v1); // v0+=H*v1
    /*[e]mltply*/
    dnorm = 0.0;
#pragma omp parallel for reduction(+:dnorm) default(none) private(i) shared(v0) firstprivate(i_max)
    for (i = 1; i <= i_max; i++) {
       dnorm += conj(v0[i]) * (v0[i]);
    }
    dnorm=SumMPI_dc(dnorm);
    dnorm = sqrt(dnorm);
#pragma omp parallel for default(none) private(i) shared(v0,v1,dnorm) firstprivate(i_max)
    for (i = 1; i <= i_max; i++) {
        v1[i]=v0[i]/dnorm;
    }
  }


  fprintf(stdoutMPI, "%s", cLogLanczos_EigenVecEnd);
  // v1 is eigen vector

  StartTimer(4500);
  if(expec_cisajs(&(X->Bind), v1)!=0){
    fprintf(stderr, "Error: calc OneBodyG.\n");
    exitMPI(-1);
  }
  StopTimer(4500);
  StartTimer(4600);  
  if(expec_cisajscktaltdc(&(X->Bind), v1)!=0){
    fprintf(stderr, "Error: calc TwoBodyG.\n");
    exitMPI(-1);
  }
  StopTimer(4600);

  if(expec_totalSz(&(X->Bind), v1)!=0){
    fprintf(stderr, "Error: calc TotalSz.\n");
    exitMPI(-1);
  }

  if(X->Bind.Def.St==0){
    sprintf(sdt, cFileNameEnergy_Lanczos, X->Bind.Def.CDataFileHead);
  }else if(X->Bind.Def.St==1){
    sprintf(sdt, cFileNameEnergy_CG, X->Bind.Def.CDataFileHead);
  }
  
  if(childfopenMPI(sdt, "w", &fp)!=0){
    exitMPI(-1);
  }

  fprintf(fp,"Energy  %.16lf \n",X->Bind.Phys.energy);
  fprintf(fp,"Doublon  %.16lf \n",X->Bind.Phys.doublon);
  fprintf(fp,"Sz  %.16lf \n",X->Bind.Phys.Sz);
  //    fprintf(fp,"total S^2  %.10lf \n",X->Bind.Phys.s2);    
  fclose(fp);

  if(X->Bind.Def.iOutputEigenVec==TRUE){
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, cOutputEigenVecStart, "a");
    sprintf(sdt, cFileNameOutputEigen, X->Bind.Def.CDataFileHead, X->Bind.Def.k_exct-1, myrank);
    if(childfopenALL(sdt, "wb", &fp)!=0){
      exitMPI(-1);
    }
    fwrite(&X->Bind.Large.itr, sizeof(X->Bind.Large.itr),1,fp);
    fwrite(&X->Bind.Check.idim_max, sizeof(X->Bind.Check.idim_max),1,fp);
    fwrite(v1, sizeof(complex double),X->Bind.Check.idim_max+1, fp);
    fclose(fp);
    TimeKeeper(&(X->Bind), cFileNameTimeKeep, cOutputEigenVecFinish, "a");
  }

  return TRUE;
}
