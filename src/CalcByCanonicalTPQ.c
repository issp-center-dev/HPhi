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
#include "FirstMultiply.h"
#include "Multiply.h"
#include "expec_energy_flct.h"
#include "expec_cisajs.h"
#include "expec_cisajscktaltdc.h"
#include "MakeIniVec.h"
#include "CalcByCanonicalTPQ.h"
#include "FileIO.h"
#include "wrapperMPI.h"
#include "CalcTime.h"


/**
 * @file   CalcByCanonicalTPQ.c 
 * @version 3.4
 * @author Takahiro Misawa (BAQIS)
 * @brief  File for giving functions of the canonical TPQ (cTPQ) method
 *
 */

/** 
 * 
 * @brief A main function to calculate physical quqntities by the cTPQ method
 *
 * @param [in] NumAve  Number of samples
 * @param [in] ExpecInterval interval steps between the steps to calculate physical quantities
 * @param [in,out] X CalcStruct list for getting and giving calculation information
 * 
 * @author Takahiro Misawa (BAQIS)
 *
 * @retval 0 normally finished
 * @retval -1 unnormally finished
 */
int CalcByCanonicalTPQ(
      const int NumAve,
      const int ExpecInterval,
      struct EDMainCalStruct *X
)
{
  char sdt[D_FileNameMax];
  char sdt_phys[D_FileNameMax];
  char sdt_norm[D_FileNameMax];
  char sdt_flct[D_FileNameMax];
  int rand_i, rand_max, iret;
  unsigned long int i_max;
  int step_iO=0;
  FILE *fp;
  double inv_temp, Ns,delta_tau;
  struct TimeKeepStruct tstruct;
  size_t byte_size;

  tstruct.tstart=time(NULL);
  
  rand_max = NumAve;
  step_spin = ExpecInterval;
  if (X->Bind.Def.Param.ExpandCoef==0){
     X->Bind.Def.Param.ExpandCoef=10;   
     fprintf(stdout, "In cTPQ calc., the default value of ExpandCoef (=10) is used. \n");
  }else{
     fprintf(stdout, "In cTPQ calc., ExpandCoef is specified as %d. \n",X->Bind.Def.Param.ExpandCoef);
  }
  X->Bind.Def.St=0;
  fprintf(stdoutMPI, "%s", cLogTPQ_Start);
  for (rand_i = 0; rand_i<rand_max; rand_i++){
    sprintf(sdt_phys, cFileNameSSRand, rand_i);      
    sprintf(sdt_norm, cFileNameNormRand, rand_i);
    sprintf(sdt_flct, cFileNameFlctRand, rand_i);
    Ns = 1.0 * X->Bind.Def.NsiteMPI;
    fprintf(stdoutMPI, cLogTPQRand, rand_i+1, rand_max);
    iret=0;
    X->Bind.Def.irand=rand_i;

  //Make or Read initial vector
    if(X->Bind.Def.iReStart==RESTART_INOUT || X->Bind.Def.iReStart==RESTART_IN) {
      StartTimer(3600);
      TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecStart, "a", rand_i, step_i);
      fprintf(stdoutMPI, "%s", cLogInputVecStart);
      sprintf(sdt, cFileNameInputVector, rand_i, myrank);
      childfopenALL(sdt, "rb", &fp);
      if(fp==NULL){
        fprintf(stdout, "A file of Inputvector does not exist.\n");
        fprintf(stdout, "Start to calculate in normal procedure.\n");
        iret=1;
      }
      byte_size = fread(&step_i, sizeof(step_i), 1, fp);
      byte_size = fread(&i_max, sizeof(long int), 1, fp);
      if(i_max != X->Bind.Check.idim_max){
        fprintf(stderr, "Error: A file of Inputvector is incorrect.\n");
        exitMPI(-1);
      }
      byte_size = fread(v0, sizeof(complex double), X->Bind.Check.idim_max+1, fp);
      TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecFinish, "a", rand_i, step_i);
      fprintf(stdoutMPI, "%s", cLogInputVecFinish);
      fclose(fp);
      StopTimer(3600);
      X->Bind.Def.istep=step_i;
      StartTimer(3200);
      iret=expec_energy_flct(&(X->Bind)); //v1 <- v0 and v0 = H*v1
      StopTimer(3200);
      if(iret != 0) return -1;

      step_iO=step_i-1;
      if (byte_size == 0) printf("byte_size: %d \n", (int)byte_size);
    }
    
    if(X->Bind.Def.iReStart==RESTART_NOT || X->Bind.Def.iReStart==RESTART_OUT || iret ==1) {
      StartTimer(3600);
      if (childfopenMPI(sdt_phys, "w", &fp) != 0) {
        return -1;
      }
      fprintf(fp, "%s", cLogSSRand);
      fclose(fp);
// for norm
      if (childfopenMPI(sdt_norm, "w", &fp) != 0) {
        return -1;
      }
      fprintf(fp, "%s", cLogNormRand);
      fclose(fp);
// for fluctuations
      if (childfopenMPI(sdt_flct, "w", &fp) != 0) {
        return -1;
      }
      fprintf(fp, "%s", cLogFlctRand);
      fclose(fp);

      StopTimer(3600);

      step_i = 0;

      StartTimer(3100);
      if(rand_i==0){
        TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cTPQStep, "w", rand_i, step_i);
      }
      else{
        TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cTPQStep, "a", rand_i, step_i);
      }
      /**@brief
      Initialize v1 and v0 = v1
      */
      MakeIniVec(rand_i, &(X->Bind)); 
      /*[s] tau*/
      inv_temp  = 0.0;
      delta_tau = 1.0/LargeValue;
      /*[e] tau*/
      StopTimer(3100);
      // for norm
      if (childfopenMPI(sdt_norm, "a", &fp) != 0) {
        return -1;
      }
      fprintf(fp, "%.16lf %.16lf %.16lf %d\n", inv_temp, global_1st_norm, global_1st_norm, step_i);
      fclose(fp);
      /**@brief
      Compute expectation value at infinite temperature
      */
      X->Bind.Def.istep = 0;
      StartTimer(3300);
      iret=expec_cisajs(&(X->Bind), v1);
      StopTimer(3300);
      if(iret !=0) return -1;

      StartTimer(3400);
      iret=expec_cisajscktaltdc(&(X->Bind), v1);
      StopTimer(3400);
      if(iret !=0) return -1;

      
      StartTimer(3200);
      iret=expec_energy_flct(&(X->Bind)); //v1 <- v0 and v0 = H*v1
      StopTimer(3200);
      if(iret !=0) return -1;
      //inv_temp = 0; /* (2.0 / Ns) / (LargeValue - X->Bind.Phys.energy / Ns);*/
      if (childfopenMPI(sdt_phys, "a", &fp) != 0) {
        return -1;
      }
      fprintf(fp, "%.16lf  %.16lf %.16lf %.16lf %.16lf %d\n", inv_temp, X->Bind.Phys.energy, X->Bind.Phys.var,
        X->Bind.Phys.doublon, X->Bind.Phys.num, step_i);
      fclose(fp);

      StartTimer(3600);
// for fluctuations
      if (childfopenMPI(sdt_flct, "a", &fp) != 0) {
        return -1;
      }
      fprintf(fp, "%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %d\n", inv_temp,X->Bind.Phys.num,X->Bind.Phys.num2, X->Bind.Phys.doublon,X->Bind.Phys.doublon2, X->Bind.Phys.Sz,X->Bind.Phys.Sz2,step_i);
      fclose(fp);
//
      StopTimer(3600);
      step_i += 1;
      X->Bind.Def.istep = step_i;
      step_iO=0;
    }

    for (step_i = X->Bind.Def.istep; step_i<X->Bind.Def.Lanczos_max; step_i++){
      X->Bind.Def.istep=step_i;
      if(step_i %((X->Bind.Def.Lanczos_max-step_iO)/10)==0){
        fprintf(stdoutMPI, cLogTPQStep, step_i, X->Bind.Def.Lanczos_max);
      }
      X->Bind.Def.istep=step_i;
      StartTimer(3600);
      TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cTPQStep, "a", rand_i, step_i);
      StopTimer(3600);
      StartTimer(3500);
      MultiplyForCanonicalTPQ(&(X->Bind),delta_tau); // v0=exp[-delta_tau*H/2]*v1 in 4th order
      StopTimer(3500);

      StartTimer(3200);
      iret=expec_energy_flct(&(X->Bind)); //v1 <- v0 and v0 = H*v1
      StopTimer(3200);
      if(iret !=0) return -1;
//
      //inv_temp = (2.0*step_i / Ns) / (LargeValue - X->Bind.Phys.energy / Ns);
      inv_temp  += delta_tau;
      //temp      = 1.0/inv_temp;

      StartTimer(3600);
      if(childfopenMPI(sdt_phys, "a", &fp)!=0){
        return FALSE;
      }
      fprintf(fp, "%.16lf  %.16lf %.16lf %.16lf %.16lf %d\n", inv_temp, X->Bind.Phys.energy, X->Bind.Phys.var, X->Bind.Phys.doublon, X->Bind.Phys.num ,step_i);
// for
      fclose(fp);

      if(childfopenMPI(sdt_norm, "a", &fp)!=0){
        return FALSE;
      }
      fprintf(fp, "%.16lf %.16lf %.16lf %d\n", inv_temp, global_norm, global_1st_norm, step_i);
      fclose(fp);

// for fluctuations
      if (childfopenMPI(sdt_flct, "a", &fp) != 0) {
        return -1;
      }
      fprintf(fp, "%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %d\n", inv_temp,X->Bind.Phys.num,X->Bind.Phys.num2, X->Bind.Phys.doublon,X->Bind.Phys.doublon2, X->Bind.Phys.Sz,X->Bind.Phys.Sz2,step_i);
      fclose(fp);
//
      StopTimer(3600);


      if (step_i%step_spin == 0){
        StartTimer(3300);
        iret=expec_cisajs(&(X->Bind),v1);
        StopTimer(3300);
        if(iret !=0) return -1;

        StartTimer(3400);
        iret=expec_cisajscktaltdc(&(X->Bind), v1);
        StopTimer(3400);
        if(iret !=0) return -1;
      }
    }

    if(X->Bind.Def.iReStart== RESTART_OUT || X->Bind.Def.iReStart==RESTART_INOUT){
      TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecStart, "a", rand_i, step_i);
      fprintf(stdoutMPI, "%s", cLogOutputVecStart);
      sprintf(sdt, cFileNameOutputVector, rand_i, myrank);
      if(childfopenALL(sdt, "wb", &fp)!=0){
        exitMPI(-1);
      }
      fwrite(&step_i, sizeof(step_i), 1, fp);
      fwrite(&X->Bind.Check.idim_max, sizeof(X->Bind.Check.idim_max),1,fp);
      fwrite(v1, sizeof(complex double),X->Bind.Check.idim_max+1, fp);
      fclose(fp);
      TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecFinish, "a", rand_i, step_i);
      fprintf(stdoutMPI, "%s", cLogOutputVecFinish);
    }
  }
  fprintf(stdoutMPI, "%s", cLogTPQ_End);

  tstruct.tend=time(NULL);
  fprintf(stdoutMPI, cLogTPQEnd, (int)(tstruct.tend-tstruct.tstart));
  return TRUE;
}
