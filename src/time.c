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

//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#include <time.h>
//#include <unistd.h>
//#include <errno.h>
//#include <limits.h>
//#include <string.h>
//#include <complex.h>
//#include "struct.h"
#ifdef MPI
  #include <mpi.h>
#endif

/*
void InitTimer();
inline void StartTimer(int n);
inline void StopTimer(int n);
void OutputTimer();
*/

void InitTimer() {
  int i;
  int NTimer=1000;
  Timer       = (double*)malloc((NTimer)*sizeof(double));
  TimerStart  = (double*)malloc((NTimer)*sizeof(double));
  for(i=0;i<NTimer;i++) Timer[i]=0.0;
  for(i=0;i<NTimer;i++) TimerStart[i]=0.0;
  return;
}

void StartTimer(int n) {
#ifdef MPI
  TimerStart[n]=MPI_Wtime();
#else
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME,&ts);
  TimerStart[n]=ts.tv_sec + ts.tv_nsec*1.0e-9;
#endif
  return;
}

void StopTimer(int n) {
#ifdef MPI
  Timer[n] += MPI_Wtime() - TimerStart[n];
#else
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME,&ts);
  Timer[n] += ts.tv_sec + ts.tv_nsec*1.0e-9 - TimerStart[n];
#endif
  return;
}

void OutputTimer(struct BindStruct *X) {
  char fileName[D_FileNameMax];
  FILE *fp;
  sprintf(fileName, "zvo_CalcTimer.dat"); //TBC
  childfopenMPI(fileName,"w", &fp);
  //fp = fopen(fileName, "w");
  //fp = childfopenMPI(fileName, "w");
  StampTime(fp, "All", 0);
  StampTime(fp, "  sz", 1);
  StampTime(fp, "  diagonalcalc", 2);
  if(X->Def.iCalcType==TPQCalc) {
    StampTime(fp, "  CalcByTPQ", 3);
    StampTime(fp, "    FirstMultiply", 101);
    StampTime(fp, "      rand   in FirstMultiply", 106);
    StampTime(fp, "      mltply in FirstMultiply", 107);
    StampTime(fp, "    expec_energy             ", 102);
    StampTime(fp, "      mltply in expec_energy ", 108);
    StampTime(fp, "    expec_onebody            ", 103);
    StampTime(fp, "    expec_twobody            ", 104);
    StampTime(fp, "    Multiply                 ", 105);
    StampTime(fp, "    FileIO                   ", 109);
  }
  else if(X->Def.iCalcType==Lanczos){
    StampTime(fp, "  CalcByLanczos", 3);
  }
  else if(X->Def.iCalcType==FullDiag){
    StampTime(fp, "  CalcByFullDiag", 3);
  }
  else if(X->Def.iFlgCalcSpec != CALCSPEC_NOT){
    StampTime(fp, "  CalcSpectrum", 3);
  }
  
  fprintf(fp,"================================================\n");
  
  StampTime(fp,"All mltply",200);
  StampTime(fp,"  diagonal", 201);
  StampTime(fp,"  trans    in HubbardGC", 202);
  StampTime(fp,"  interall in HubbardGC", 203);
  StampTime(fp,"  pairhopp in HubbardGC", 204);
  StampTime(fp,"  exchange in HubbardGC", 205);
  fprintf(fp,"\n");
  StampTime(fp,"  trans    in Hubbard", 206);
  StampTime(fp,"  interall in Hubbard", 207);
  StampTime(fp,"  pairhopp in Hubbard", 208);
  StampTime(fp,"  exchange in Hubbard", 209);
  fprintf(fp,"\n");
  StampTime(fp,"  interall in Spin", 210);
  StampTime(fp,"    double", 220);
  StampTime(fp,"    single1", 221);
  StampTime(fp,"    single2", 222);
  StampTime(fp,"    else", 223);
  StampTime(fp,"  exchange in Spin", 211);
  StampTime(fp,"    double", 224);
  StampTime(fp,"    single1", 225);
  StampTime(fp,"    single2", 226);
  StampTime(fp,"    else", 227);
  fprintf(fp,"\n");
  StampTime(fp,"  trans    in SpinGC", 212);
  StampTime(fp,"  interall in SpinGC", 213);
  StampTime(fp,"  exchange in SpinGC", 214);
  StampTime(fp,"  pairlift in SpinGC", 215);
  fprintf(fp,"================================================\n");

  fclose(fp);
  free(Timer);
  free(TimerStart);
}

void StampTime(FILE *fp, char *str, int num){
  char str1[256];
  sprintf(str1, "%-40s [%03d] %12.5lf\n", str, num, Timer[num]);
  fprintf(fp, str1);
}
