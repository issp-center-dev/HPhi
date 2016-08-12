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
#include "expec_energy.h"
#include "expec_cisajs.h"
#include "expec_cisajscktaltdc.h"
#include "CalcByTPQ.h"
#include "FileIO.h"
#include "wrapperMPI.h"

/** 
 * 
 * 
 * @param NumAve 
 * @param ExpecInterval 
 * @param X 
 * 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 *
 * @return 
 */
int CalcByTPQ(
	    const int NumAve,
	    const int ExpecInterval,
	    struct EDMainCalStruct *X
)
{
  char sdt[D_FileNameMax];
  char sdt_phys[D_FileNameMax];
  char  sdt_norm[D_FileNameMax];
  int rand_i, rand_max, iret;
  long int i_max;
  FILE *fp;
  double inv_temp, Ns;
  struct TimeKeepStruct tstruct;
  tstruct.tstart=time(NULL);

  rand_max = NumAve;
  step_spin = ExpecInterval;
  X->Bind.Def.St=0;
  fprintf(stdoutMPI, cLogTPQ_Start);
  for (rand_i = 0; rand_i<rand_max; rand_i++){

    sprintf(sdt_phys, cFileNameSSRand, rand_i);
    sprintf(sdt_norm, cFileNameNormRand, rand_i);      
    Ns = 1.0 * X->Bind.Def.NsiteMPI;
    fprintf(stdoutMPI, cLogTPQRand, rand_i+1, rand_max);
    iret=0;
    X->Bind.Def.irand=rand_i;

    if(rand_i==0){
      TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cTPQStep, "w", rand_i, step_i);
    }
    else{
      TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cTPQStep, "a", rand_i, step_i);
    }

  //Make or Read initial vector
    if(X->Bind.Def.iReStart==RESTART_INOUT || X->Bind.Def.iReStart==RESTART_IN) {
      TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecStart, "a", rand_i, step_i);
      fprintf(stdoutMPI, cLogInputVecStart);
      sprintf(sdt, cFileNameInputVector, rand_i, myrank);
      childfopenALL(sdt, "rb", &fp);
      if(fp==NULL){
        fprintf(stdout, "A file of Inputvector does not exist.\n");
        fprintf(stdout, "Start to calculate in normal procedure.\n");
        iret=1;
      }
      fread(&step_i, sizeof(step_i), 1, fp);
      fread(&i_max, sizeof(long int), 1, fp);
      //fprintf(stdoutMPI, "Debug: i_max=%ld, step_i=%d\n", i_max, step_i);
      if(i_max != X->Bind.Check.idim_max){
        fprintf(stderr, "Error: A file of Inputvector is incorrect.\n");
        exitMPI(-1);
      }
      fread(v0, sizeof(complex double), X->Bind.Check.idim_max+1, fp);
      TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecFinish, "a", rand_i, step_i);
      fprintf(stdoutMPI, cLogInputVecFinish);
      fclose(fp);
      X->Bind.Def.istep=step_i;
      expec_energy(&(X->Bind));
    }
    
    if(X->Bind.Def.iReStart==RESTART_NOT || X->Bind.Def.iReStart==RESTART_OUT || iret ==1) {
      
      if (!childfopenMPI(sdt_phys, "w", &fp) == 0) {
        return -1;
      }
      fprintf(fp, cLogSSRand);
      fclose(fp);
      if (!childfopenMPI(sdt_norm, "w", &fp) == 0) {
        return -1;
      }
      fprintf(fp, cLogNormRand);
      fclose(fp);

      step_i = 1;
      FirstMultiply(rand_i, &(X->Bind));
      expec_energy(&(X->Bind)); //v0 = H*v1
      inv_temp = (2.0 / Ns) / (LargeValue - X->Bind.Phys.energy / Ns);

      X->Bind.Def.istep = step_i;
      expec_cisajs(&(X->Bind), v1);
      expec_cisajscktaltdc(&(X->Bind), v1);

      if (!childfopenMPI(sdt_phys, "a", &fp) == 0) {
        return -1;
      }
      fprintf(fp, "%.16lf  %.16lf %.16lf %.16lf %.16lf %d\n", inv_temp, X->Bind.Phys.energy, X->Bind.Phys.var,
              X->Bind.Phys.doublon, X->Bind.Phys.num, step_i);
      fclose(fp);

      if (!childfopenMPI(sdt_norm, "a", &fp) == 0) {
        return -1;
      }
      fprintf(fp, "%.16lf %.16lf %.16lf %d\n", inv_temp, global_norm, global_1st_norm, step_i);
      fclose(fp);
    }

    fprintf(stdoutMPI, "step_i=%d\n", X->Bind.Def.istep);

    for (step_i = X->Bind.Def.istep; step_i<X->Bind.Def.Lanczos_max; step_i++){
      X->Bind.Def.istep=step_i;
      if(step_i %(X->Bind.Def.Lanczos_max/10)==0){
        fprintf(stdoutMPI, cLogTPQStep, step_i, X->Bind.Def.Lanczos_max);
      }
      X->Bind.Def.istep=step_i;
      TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cTPQStep, "a", rand_i, step_i);
      Multiply(&(X->Bind));
      expec_energy(&(X->Bind));
      //expec(&(X->Bind));
      inv_temp = (2.0*step_i / Ns) / (LargeValue - X->Bind.Phys.energy / Ns);
      if(!childfopenMPI(sdt_phys, "a", &fp)==0){
        return FALSE;
      }
      fprintf(fp, "%.16lf  %.16lf %.16lf %.16lf %.16lf %d\n", inv_temp, X->Bind.Phys.energy, X->Bind.Phys.var, X->Bind.Phys.doublon, X->Bind.Phys.num ,step_i);
      fclose(fp);

      if(!childfopenMPI(sdt_norm, "a", &fp)==0){
        return FALSE;
      }
      fprintf(fp, "%.16lf %.16lf %.16lf %d\n", inv_temp, global_norm, global_1st_norm, step_i);
      fclose(fp);

      if (step_i%step_spin == 0){
        expec_cisajs(&(X->Bind),v1);
        expec_cisajscktaltdc(&(X->Bind), v1);
      }
    }

    if(X->Bind.Def.iReStart== RESTART_OUT || X->Bind.Def.iReStart==RESTART_INOUT){
      TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecStart, "a", rand_i, step_i);
      fprintf(stdoutMPI, cLogOutputVecStart);
      sprintf(sdt, cFileNameOutputVector, rand_i, myrank);
      if(childfopenALL(sdt, "wb", &fp)!=0){
        exitMPI(-1);
      }
      fwrite(&step_i, sizeof(step_i), 1, fp);
      fwrite(&X->Bind.Check.idim_max, sizeof(X->Bind.Check.idim_max),1,fp);
      fwrite(v1, sizeof(complex double),X->Bind.Check.idim_max+1, fp);
      fclose(fp);
      TimeKeeperWithRandAndStep(&(X->Bind), cFileNameTPQStep, cOutputVecFinish, "a", rand_i, step_i);
      fprintf(stdoutMPI, cLogOutputVecFinish);
    }
  }

  fprintf(stdoutMPI, cLogTPQ_End);
  tstruct.tend=time(NULL);
  fprintf(stdoutMPI, cLogTPQEnd, (int)(tstruct.tend-tstruct.tstart));
  return TRUE;
}
