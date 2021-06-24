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
/*-------------------------------------------------------------*/
/*-------------------------------------------------------------
 * HPhi
 * timer program
 * "-lrt" option is needed for clock_gettime().
 *-------------------------------------------------------------
 * original code written by Satoshi Morita
 *-------------------------------------------------------------*/

#include "Common.h"
#include "FileIO.h"
#include "CalcTime.h"

#ifdef MPI
#include <mpi.h>
#endif
/** 
 * 
 * @brief function for displaying elapse time
 * 
 * @version 2.0
 */

void StampTime(FILE *fp, char *str, int num){
#ifdef MPI
  char str1[256];
  sprintf(str1, "%-50s [%04d] %12.5lf\n", str, num, Timer[num]);
  fprintf(fp, "%s", str1);
#endif
}

/** 
 * 
 * @brief function for initializing Timer[]
 * 
 * @version 2.0
 */
void InitTimer() {
#ifdef MPI
  int i;
  int NTimer=10000;
  Timer       = (double*)malloc((NTimer)*sizeof(double));
  TimerStart  = (double*)malloc((NTimer)*sizeof(double));
  for(i=0;i<NTimer;i++) Timer[i]=0.0;
  for(i=0;i<NTimer;i++) TimerStart[i]=0.0;
#endif
  return;
}
/** 
 * 
 * @brief function for initializing elapse time [start]
 * 
 * @version 2.0
 */

void StartTimer(int n) {
#ifdef MPI
  TimerStart[n]=MPI_Wtime();
#endif
  return;
}
/** 
 * 
 * @brief function for calculating elapse time [elapse time=StartTimer-StopTimer]
 * 
 * @version 2.0
 */
void StopTimer(int n) {
#ifdef MPI
  Timer[n] += MPI_Wtime() - TimerStart[n];
#endif
  return;
}
/** 
 * 
 * @brief function for outputting elapse time for each function
 * 
 * @version 2.0
 */
void OutputTimer(struct BindStruct *X) {

#ifdef MPI
  char fileName[D_FileNameMax];
  FILE *fp;
  sprintf(fileName, "CalcTimer.dat"); //TBC
  childfopenMPI(fileName,"w", &fp);
  //fp = fopen(fileName, "w");
  //fp = childfopenMPI(fileName, "w");
  StampTime(fp, "All", 0);
  StampTime(fp, "  sz", 1000);
  StampTime(fp, "  diagonalcalc", 2000);
  if(X->Def.iFlgCalcSpec == CALCSPEC_NOT){
    if(X->Def.iCalcType==TPQCalc || X->Def.iCalcType==cTPQ) {
      StampTime(fp, "  CalcByTPQ", 3000);
      StampTime(fp, "    FirstMultiply", 3100);
      StampTime(fp, "      rand   in FirstMultiply", 3101);
      StampTime(fp, "      mltply in FirstMultiply", 3102);
      StampTime(fp, "    expec_energy_flct        ", 3200);
      StampTime(fp, "      calc flctuation in expec_energy_flct ", 3201);
      StampTime(fp, "      mltply in expec_energy_flct ", 3202);
      StampTime(fp, "    expec_onebody            ", 3300);
      StampTime(fp, "    expec_twobody            ", 3400);
      StampTime(fp, "    Multiply                 ", 3500);
      StampTime(fp, "    FileIO                   ", 3600);
    }
    else if(X->Def.iCalcType==Lanczos){
      StampTime(fp, "  CalcByLanczos", 4000);
      StampTime(fp, "    LanczosEigenValue", 4100);
      StampTime(fp, "      mltply      in LanczosEigenValue", 4101);
      StampTime(fp, "      vec12       in LanczosEigenValue", 4102);
      StampTime(fp, "      DSEVvalue   in LanczosEigenValue", 4103);
      StampTime(fp, "    LanczosEigenVector", 4200);
      StampTime(fp, "      mltply      in LanczosEigenVector", 4201);
      StampTime(fp, "    expec_energy_flct", 4300);
      StampTime(fp, "      calc flctuation in expec_energy_flct ", 4301);
      StampTime(fp, "      mltply in expec_energy_flct ", 4302);
      StampTime(fp, "    CGEigenVector", 4400);
      StampTime(fp, "      mltply in CGEigenVector ", 4401);
      StampTime(fp, "    expec_onebody            ", 4500);
      StampTime(fp, "    expec_twobody            ", 4600);
      StampTime(fp, "    expec_TotalSz            ", 4700);
      StampTime(fp, "    FileIO                   ", 4800);
      StampTime(fp, "      Read Input Eigenvec ", 4801);    
    }
    else if(X->Def.iCalcType==FullDiag){
      StampTime(fp, "  CalcByFullDiag", 5000);
      StampTime(fp, "    MakeHam", 5100);
      StampTime(fp, "    LapackDiag", 5200);
      StampTime(fp, "    CalcPhys", 5300);
    StampTime(fp, "      calc flctuation in expec_energy_flct ", 5301);
    StampTime(fp, "      mltply in expec_energy_flct ", 5302);
        StampTime(fp, "    Output", 5400);
      StampTime(fp, "    OutputHam", 5500);
    }
  }
  else{ 
    StampTime(fp, "  CalcSpectrum by Lanczos method", 6000);
    StampTime(fp, "    Make excited state", 6100);
    StampTime(fp, "      Read origin state", 6101);
    StampTime(fp, "      Multiply excited operator", 6102);
    StampTime(fp, "    Calculate spectrum", 6200);
    if(X->Def.iCalcType==Lanczos){
      StampTime(fp, "      Read vector for recalculation", 6201);
      StampTime(fp, "      Read tridiagonal components for recalculation", 6202);
      StampTime(fp, "      Calculate tridiagonal components", 6203);
      StampTime(fp, "      Output tridiagonal components", 6204);
      StampTime(fp, "      Calculate spectrum by Lanczos method", 6205);
      StampTime(fp, "      Output vectors for recalculation", 6206);
    }
    else if(X->Def.iCalcType==FullDiag){
      StampTime(fp, "      MakeHam", 6301);
      StampTime(fp, "      lapackdiag", 6302);
      StampTime(fp, "      Calculate v1", 6303);
      StampTime(fp, "      Calculate spectrum", 6304);
    }
  }
  
  fprintf(fp,"================================================\n");
  
  StampTime(fp,"All mltply",1);
  StampTime(fp,"  diagonal", 100);

  switch(X->Def.iCalcModel){
  case HubbardGC:
    StampTime(fp,"  HubbardGC", 200);
    StampTime(fp,"    trans    in HubbardGC", 210);
    StampTime(fp,"      double", 211);
    StampTime(fp,"      single", 212);
    StampTime(fp,"      inner", 213);
    StampTime(fp,"    interall in HubbardGC", 220);
    StampTime(fp,"      interPE", 221);
    StampTime(fp,"      inner", 222);
    StampTime(fp,"    pairhopp in HubbardGC", 230);
    StampTime(fp,"      interPE", 231);
    StampTime(fp,"      inner", 232);
    StampTime(fp,"    exchange in HubbardGC", 240);
    StampTime(fp,"      interPE", 241);
    StampTime(fp,"      inner", 242);
    break;
    
  case Hubbard:
    StampTime(fp,"  Hubbard", 300);
    StampTime(fp,"    trans    in Hubbard", 310);
    StampTime(fp,"      double", 311);
    StampTime(fp,"      single", 312);
    StampTime(fp,"      inner", 313);
    StampTime(fp,"    interall in Hubbard", 320);
    StampTime(fp,"      interPE", 321);
    StampTime(fp,"      inner", 322);
    StampTime(fp,"    pairhopp in Hubbard", 330);
    StampTime(fp,"      interPE", 331);
    StampTime(fp,"      inner", 332);
    StampTime(fp,"    exchange in Hubbard", 340);
    StampTime(fp,"      interPE", 341);
    StampTime(fp,"      inner", 342);
    break;
    
  case Spin:
    fprintf(fp,"\n");
    StampTime(fp,"  Spin", 400);
    StampTime(fp,"    interall in Spin", 410);
    StampTime(fp,"      double", 411);
    StampTime(fp,"      single1", 412);
    StampTime(fp,"      single2", 413);
    StampTime(fp,"      inner", 414);
    StampTime(fp,"    exchange in Spin", 420);
    StampTime(fp,"      double", 421);
    StampTime(fp,"      single1", 422);
    StampTime(fp,"      single2", 423);
    StampTime(fp,"      inner", 424);
    break;
    
  case SpinGC:
    StampTime(fp,"  SpinGC", 500);
    StampTime(fp,"    trans    in SpinGC", 510);
    StampTime(fp,"      double", 511);
    StampTime(fp,"      inner", 512);
    StampTime(fp,"    interall in SpinGC", 520);
    StampTime(fp,"      double", 521);
    StampTime(fp,"      single", 522);
    StampTime(fp,"      inner", 523);
    StampTime(fp,"    exchange in SpinGC", 530);
    StampTime(fp,"      double", 531);
    StampTime(fp,"      single", 532);
    StampTime(fp,"      inner", 533);
    StampTime(fp,"    pairlift in SpinGC", 540);
    StampTime(fp,"      double", 541);
    StampTime(fp,"      single", 542);
    StampTime(fp,"      inner", 543);
    break;

  default:
    break;
  }
  fprintf(fp,"================================================\n");

  fclose(fp);
  free(Timer);
  free(TimerStart);
#endif
  return;
}

/**
@page page_time Compute elapsed time for new functions

 Using StartTimer and StopTimer functions defined in time.c, we can measure the elapsed time for computation.

 1. Define an index and an output message by using StampTime function in time.c.
 For example, the index and the output message for the elapsed time of TPQ calculation is defined as follows.
 ```
       StampTime(fp, "  CalcByTPQ", 3000);
 ```

 2. Include CalcTime.h in the target source file.

 3. Set StartTimer and StopTimer functions in the region where you want to measure the time.
    It is noted that both functions must have the same index defined in time.c.
 For example, the elapsed time of TPQ calculation can be measured as follows.
 ```
       case TPQCalc:
        StartTimer(3000);
        if (CalcByTPQ(NumAve, X.Bind.Def.Param.ExpecInterval, &X) != TRUE) {
          FinalizeMPI();
          StopTimer(3000);
          return 0;
        }
        StopTimer(3000);
      break;

 ```

When above procedures were done, after calculation, you can see the elapsed time in CalcTimer.dat file as follows. The time unit is second.

```
 All                                                [0000]     37.94046
  sz                                               [1000]      0.00058
  diagonalcalc                                     [2000]      0.00046
  CalcByTPQ                                        [3000]     37.93129
```


*/
