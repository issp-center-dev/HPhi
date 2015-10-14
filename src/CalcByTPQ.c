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
#include "CalcByTPQ.h"
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
int CalcBySSM(
	    const int NumAve,
	    const int ExpecInterval,
	    struct EDMainCalStruct *X
)
{
  long unsigned int u_long_i;
  dsfmt_t dsfmt;
  char sdt_phys[D_FileNameMax];
  char  sdt_norm[D_FileNameMax];
  int rand_i, rand_max;
  FILE *fp;
  double inv_temp, Ns;
  struct TimeKeepStruct tstruct;
  tstruct.tstart=time(NULL);
  
  
  rand_max = NumAve;
  step_spin = ExpecInterval;
  X->Bind.Def.St=0;
  for (rand_i = 0; rand_i<rand_max; rand_i++){
    fprintf(stdoutMPI, "rand_i=%d \n", rand_i);
    u_long_i = 123432 + rand_i*X->Bind.Def.initial_iv;
    dsfmt_init_gen_rand(&dsfmt, u_long_i);    
    sprintf(sdt_phys, cFileNameSSRand, rand_i);
    if(!childfopenMPI(sdt_phys, "w", &fp)==0){
      return -1;
    }
    fprintf(fp, " # inv_tmp, energy, phys_var, phys_doublon, phys_num, step_i\n");
    fclose(fp);
    
    sprintf(sdt_norm, cFileNameNormRand, rand_i);
    if(!childfopenMPI(sdt_norm, "w", &fp)==0){
      return -1;
    }
    fprintf(fp, " # inv_temp, global_norm, global_1st_norm, step_i \n");
    fclose(fp);
    
    FirstMultiply(&dsfmt, &(X->Bind));
    
    //    expec(&(X->Bind));
    expec_energy(&(X->Bind));
    Ns = 1.0*X->Bind.Def.Nsite;
    inv_temp = (2.0 / Ns) / (LargeValue - X->Bind.Phys.energy / Ns);
    step_i = 1;
    X->Bind.Def.istep=step_i;
    X->Bind.Def.irand=rand_i;
    
    if(X->Bind.Def.iCalcModel!=Spin){
      expec_cisajs(&(X->Bind),v1);
    }

    expec_cisajscktaltdc(&(X->Bind), v1);
    if(!childfopenMPI(sdt_phys, "a", &fp)==0){
      return -1;
    }
    fprintf(fp, "%lf  %lf %lf %lf %lf %d\n", inv_temp, X->Bind.Phys.energy, X->Bind.Phys.var, X->Bind.Phys.doublon, X->Bind.Phys.num ,step_i);
    fclose(fp);

    if(!childfopenMPI(sdt_norm, "a", &fp)==0){
      return -1;
    }
    fprintf(fp, "%lf %lf %lf %d\n", inv_temp, global_norm, global_1st_norm, step_i);
    fclose(fp);
    for (step_i = 2; step_i<X->Bind.Def.Lanczos_max; step_i++){

      if(step_i %(X->Bind.Def.Lanczos_max/10)==0){
	fprintf(stdoutMPI, "step_i/total_step=%d/%d \n", step_i, X->Bind.Def.Lanczos_max);
      }
      X->Bind.Def.istep=step_i;
      TimeKeeperWithStep(&(X->Bind), cFileNameTPQStep, cTPQStep, "a", step_i);
      Multiply(&(X->Bind));
      TimeKeeperWithStep(&(X->Bind), cFileNameTimeKeep, cTPQStepEnd, "a", step_i);
      expec_energy(&(X->Bind));
      //expec(&(X->Bind));
      inv_temp = (2.0*step_i / Ns) / (LargeValue - X->Bind.Phys.energy / Ns);
      if(!childfopenMPI(sdt_phys, "a", &fp)==0){
	return -1;
      }
      fprintf(fp, "%lf  %lf %lf %lf %lf %d\n", inv_temp, X->Bind.Phys.energy, X->Bind.Phys.var, X->Bind.Phys.doublon, X->Bind.Phys.num ,step_i);
      fclose(fp);

      if(!childfopenMPI(sdt_norm, "a", &fp)==0){
	return -1;
      }
      fprintf(fp, "%lf %lf %lf %d\n", inv_temp, global_norm, global_1st_norm, step_i);
      fclose(fp);

      if (step_i%step_spin == 0){
	if(X->Bind.Def.iCalcModel!=Spin){
	  expec_cisajs(&(X->Bind),v1);
	}
	expec_cisajscktaltdc(&(X->Bind), v1);
      }
    }
  }

  tstruct.tend=time(NULL);
  fprintf(stdoutMPI, "Finish: Elapsed time is %d [s].\n", (int)(tstruct.tend-tstruct.tstart));
  return 0;
}
