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
#include "CalcByTPQ.h"
#include "FileIO.h"
#include "wrapperMPI.h"
#include "CalcTime.h"
#include "common/setmemory.h"
/**
 * @file   CalcByTPQ.c
 * @version 0.1, 0.2
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 *
 * @brief  File for givinvg functions of TPQ method
 */
/** 
 * @brief A main function to calculate physical quqntities by TPQ method
 *
 * @param [in] NumAve  Number of samples
 * @param [in] ExpecInterval interval steps between the steps to calculate physical quantities
 * @param [in,out] X CalcStruct list for getting and giving calculation information
 * 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 *
 * @retval 0 normally finished
 * @retval -1 unnormally finished
 */
int CalcByTPQ(
  const int NumAve,
  const int ExpecInterval,
  struct EDMainCalStruct *X
) {
  char sdt[D_FileNameMax];
  char **sdt_phys, **sdt_norm, **sdt_flct;
  int rand_i, iret;
  unsigned long int i_max;
  int step_iO = 0;
  FILE *fp;
  double *inv_temp, Ns;
  struct TimeKeepStruct tstruct;
  size_t byte_size;

  tstruct.tstart = time(NULL);
  inv_temp = d_1d_allocate(NumAve);

  step_spin = ExpecInterval;
  X->Bind.Def.St = 0;
  fprintf(stdoutMPI, "%s", cLogTPQ_Start);

  //for rand_i =0, rand_i<NumAve
  sdt_phys = (char**)malloc(sizeof(char*)*NumAve);
  sdt_norm = (char**)malloc(sizeof(char*)*NumAve);
  sdt_flct = (char**)malloc(sizeof(char*)*NumAve);
  for (rand_i = 0; rand_i < NumAve; rand_i++) {
    sdt_phys[rand_i] = (char*)malloc(sizeof(char)*D_FileNameMax);
    sdt_norm[rand_i] = (char*)malloc(sizeof(char)*D_FileNameMax);
    sdt_flct[rand_i] = (char*)malloc(sizeof(char)*D_FileNameMax);
    sprintf(sdt_phys[rand_i], cFileNameSSRand, rand_i);
    sprintf(sdt_norm[rand_i], cFileNameNormRand, rand_i);
    sprintf(sdt_flct[rand_i], cFileNameFlctRand, rand_i);
  }
  Ns = 1.0 * X->Bind.Def.NsiteMPI;
  fprintf(stdoutMPI, cLogTPQRand, 1, NumAve);
  iret = 0;

  //Make or Read initial vector
  if (X->Bind.Def.iReStart == RESTART_INOUT || X->Bind.Def.iReStart == RESTART_IN) {
    StartTimer(3600);
    TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecStart, "a", 0, step_i);
    fprintf(stdoutMPI, "%s", cLogInputVecStart);
    sprintf(sdt, cFileNameInputVector, rand_i, myrank);
    childfopenALL(sdt, "rb", &fp);
    if (fp == NULL) {
      fprintf(stdout, "A file of Inputvector does not exist.\n");
      fprintf(stdout, "Start to calculate in normal procedure.\n");
      iret = 1;
    }
    byte_size = fread(&step_i, sizeof(step_i), 1, fp);
    byte_size = fread(&i_max, sizeof(long int), 1, fp);
    if (i_max != X->Bind.Check.idim_max) {
      fprintf(stderr, "Error: A file of Inputvector is incorrect.\n");
      exitMPI(-1);
    }
    byte_size = fread(v0, sizeof(complex double), (X->Bind.Check.idim_max + 1)*NumAve, fp);
    TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecFinish, "a", 0, step_i);
    fprintf(stdoutMPI, "%s", cLogInputVecFinish);
    fclose(fp);
    StopTimer(3600);
    X->Bind.Def.istep = step_i;
    StartTimer(3200);
    iret = expec_energy_flct(&(X->Bind), NumAve, v0, v1);
    StopTimer(3200);
    if (iret != 0) return -1;

    step_iO = step_i - 1;
    if (byte_size == 0) printf("byte_size: %d \n", (int)byte_size);
  }/*if (X->Bind.Def.iReStart == RESTART_INOUT || X->Bind.Def.iReStart == RESTART_IN)*/

  if (X->Bind.Def.iReStart == RESTART_NOT || X->Bind.Def.iReStart == RESTART_OUT || iret == 1) {
    StartTimer(3600);
    for (rand_i = 0; rand_i < NumAve; rand_i++) {
      if (childfopenMPI(sdt_phys[rand_i], "w", &fp) == 0) {
        fprintf(fp, "%s", cLogSSRand);
        fclose(fp);
      }
      else return -1;
      // for norm
      if (childfopenMPI(sdt_norm[rand_i], "w", &fp) == 0) {
        fprintf(fp, "%s", cLogNormRand);
        fclose(fp); 
      }
      else return -1;
      // for fluctuations
      if (childfopenMPI(sdt_flct[rand_i], "w", &fp) == 0) {
        fprintf(fp, "%s", cLogFlctRand);
        fclose(fp);
      }
      else return -1;
    }
    StopTimer(3600);

    step_i = 0;

    StartTimer(3100);
    if (rand_i == 0) {
      TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cTPQStep, "w", 0, step_i);
    }
    else {
      TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cTPQStep, "a", 0, step_i);
    }
    /**@brief
    Initialize v1 and compute v0 = H*v1
    */
    FirstMultiply(&(X->Bind));
    StopTimer(3100);
    for (rand_i = 0; rand_i < NumAve; rand_i++) {
      inv_temp[rand_i] = 0.0;
      if (childfopenMPI(sdt_phys[rand_i], "a", &fp) == 0) {
        fprintf(fp, "%.16lf  %.16lf %.16lf %.16lf %.16lf %d\n", 
          inv_temp[rand_i], X->Bind.Phys.energy[rand_i], X->Bind.Phys.var[rand_i],
          X->Bind.Phys.doublon[rand_i], X->Bind.Phys.num[rand_i], step_i);
        fclose(fp);
      }
      else return -1;
      // for norm
      if (childfopenMPI(sdt_norm[rand_i], "a", &fp) == 0) {
        fprintf(fp, "%.16lf %.16lf %.16lf %d\n", 
          inv_temp[rand_i], global_1st_norm[rand_i], global_1st_norm[rand_i], step_i);
        fclose(fp);
      }
      else return -1;
    }
    
     /**@brief
    Compute expectation value at infinite temperature
    */
    X->Bind.Def.istep = 0;
    StartTimer(3300);
    iret = expec_cisajs(&(X->Bind), NumAve, v0, v1);
    StopTimer(3300);
    if (iret != 0) return -1;

    StartTimer(3400);
    iret = expec_cisajscktaltdc(&(X->Bind), NumAve, v0, v1);
    StopTimer(3400);
    if (iret != 0) return -1;

    /**@brief
    Compute v1=0, and compute v0 = H*v1
    */
    StartTimer(3200);
    iret = expec_energy_flct(&(X->Bind), NumAve, v0, v1); //v0 = H*v1
    StopTimer(3200);
    if (iret != 0) return -1;
    step_i += 1;
    StartTimer(3600);
    for (rand_i = 0; rand_i < NumAve; rand_i++) {
      inv_temp[rand_i] = (2.0 / Ns) / (LargeValue - X->Bind.Phys.energy[rand_i] / Ns);
      if (childfopenMPI(sdt_phys[rand_i], "a", &fp) == 0) {
        fprintf(fp, "%.16lf  %.16lf %.16lf %.16lf %.16lf %d\n", 
          inv_temp[rand_i], X->Bind.Phys.energy[rand_i], X->Bind.Phys.var[rand_i],
          X->Bind.Phys.doublon[rand_i], X->Bind.Phys.num[rand_i], step_i);
        fclose(fp);
      }
      else return -1;     
      // for norm
      if (childfopenMPI(sdt_norm[rand_i], "a", &fp) == 0) {
        fprintf(fp, "%.16lf %.16lf %.16lf %d\n", 
          inv_temp[rand_i], global_norm[rand_i], global_1st_norm[rand_i], step_i);
        fclose(fp);
      }
      else return -1;
      // for fluctuations
      if (childfopenMPI(sdt_flct[rand_i], "a", &fp) == 0) {
        fprintf(fp, "%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %d\n", 
          inv_temp[rand_i], X->Bind.Phys.num[rand_i], X->Bind.Phys.num2[rand_i], 
          X->Bind.Phys.doublon[rand_i], X->Bind.Phys.doublon2[rand_i], 
          X->Bind.Phys.Sz[rand_i], X->Bind.Phys.Sz2[rand_i], step_i);
        fclose(fp);
      }
      else return -1;      
    }/*for (rand_i = 0; rand_i < NumAve; rand_i++)*/
    StopTimer(3600);
    step_i += 1;
    X->Bind.Def.istep = step_i;
    step_iO = 0;
  }/*if (X->Bind.Def.iReStart == RESTART_NOT || X->Bind.Def.iReStart == RESTART_OUT || iret == 1)*/

  for (step_i = X->Bind.Def.istep; step_i < X->Bind.Def.Lanczos_max; step_i++) {
    X->Bind.Def.istep = step_i;
    if (step_i % ((X->Bind.Def.Lanczos_max - step_iO) / 10) == 0) {
      fprintf(stdoutMPI, cLogTPQStep, step_i, X->Bind.Def.Lanczos_max);
    }
    X->Bind.Def.istep = step_i;
    StartTimer(3600);
    TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cTPQStep, "a", 0, step_i);
    StopTimer(3600);
    StartTimer(3500);
    Multiply(&(X->Bind));
    StopTimer(3500);

    StartTimer(3200);
    iret = expec_energy_flct(&(X->Bind), NumAve, v0, v1);
    StopTimer(3200);
    if (iret != 0) return -1;

    StartTimer(3600);
    for (rand_i = 0; rand_i < NumAve; rand_i++) {
      inv_temp[rand_i] = (2.0*step_i / Ns) / (LargeValue - X->Bind.Phys.energy[rand_i] / Ns);
      if (childfopenMPI(sdt_phys[rand_i], "a", &fp) == 0) {
        fprintf(fp, "%.16lf  %.16lf %.16lf %.16lf %.16lf %d\n", 
          inv_temp[rand_i], X->Bind.Phys.energy[rand_i], X->Bind.Phys.var[rand_i],
          X->Bind.Phys.doublon[rand_i], X->Bind.Phys.num[rand_i], step_i);
        // for
        fclose(fp);
      }
      else return FALSE;

      if (childfopenMPI(sdt_norm[rand_i], "a", &fp) == 0) {
        fprintf(fp, "%.16lf %.16lf %.16lf %d\n", 
          inv_temp[rand_i], global_norm[rand_i], global_1st_norm[rand_i], step_i);
        fclose(fp);
      }
      else return FALSE;

      // for fluctuations
      if (childfopenMPI(sdt_flct[rand_i], "a", &fp) == 0) {
        fprintf(fp, "%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %d\n", 
          inv_temp[rand_i], X->Bind.Phys.num[rand_i], X->Bind.Phys.num2[rand_i],
          X->Bind.Phys.doublon[rand_i], X->Bind.Phys.doublon2[rand_i],
          X->Bind.Phys.Sz[rand_i], X->Bind.Phys.Sz2[rand_i], step_i);
        fclose(fp);
      }
      else return -1;
    }/*for (rand_i = 0; rand_i < NumAve; rand_i++)*/
    StopTimer(3600);

    if (step_i%step_spin == 0) {
      StartTimer(3300);
      iret = expec_cisajs(&(X->Bind), NumAve, v0, v1);
      StopTimer(3300);
      if (iret != 0) return -1;

      StartTimer(3400);
      iret = expec_cisajscktaltdc(&(X->Bind), NumAve, v0, v1);
      StopTimer(3400);
      if (iret != 0) return -1;
    }
  }/*for (step_i = X->Bind.Def.istep; step_i < X->Bind.Def.Lanczos_max; step_i++)*/

  if (X->Bind.Def.iReStart == RESTART_OUT || X->Bind.Def.iReStart == RESTART_INOUT) {
    TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecStart, "a", 0, step_i);
    fprintf(stdoutMPI, "%s", cLogOutputVecStart);
    sprintf(sdt, cFileNameOutputVector, 0, myrank);
    if (childfopenALL(sdt, "wb", &fp) != 0) {
      exitMPI(-1);
    }
    fwrite(&step_i, sizeof(step_i), 1, fp);
    fwrite(&X->Bind.Check.idim_max, sizeof(X->Bind.Check.idim_max), 1, fp);
    fwrite(v1, sizeof(complex double), (X->Bind.Check.idim_max + 1)*NumAve, fp);
    fclose(fp);
    TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecFinish, "a", 0, step_i);
    fprintf(stdoutMPI, "%s", cLogOutputVecFinish);
  }

  fprintf(stdoutMPI, "%s", cLogTPQ_End);

  tstruct.tend = time(NULL);
  fprintf(stdoutMPI, cLogTPQEnd, (int)(tstruct.tend - tstruct.tstart));
  free_d_1d_allocate(inv_temp);

  for (rand_i = 0; rand_i < NumAve; rand_i++) {
    free(sdt_phys[rand_i]);
    free(sdt_norm[rand_i]);
    free(sdt_flct[rand_i]);
  }
  free(sdt_phys);
  free(sdt_norm);
  free(sdt_flct);

  return TRUE;
}
