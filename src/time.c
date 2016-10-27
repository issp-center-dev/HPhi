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

void OutputTimer() {
  char fileName[D_FileNameMax];
  FILE *fp;
  sprintf(fileName, "zvo_CalcTimer.dat"); //TBC
  childfopenMPI(fileName,"w", &fp);
  //fp = fopen(fileName, "w");
  //fp = childfopenMPI(fileName, "w");

  fprintf(fp,"All                           [0]   %12.5lf\n",Timer[0]);
  fprintf(fp,"  sz                          [1]   %12.5lf\n",Timer[1]);
  fprintf(fp,"  diagonalcalc                [2]   %12.5lf\n",Timer[2]);
  fprintf(fp,"  CalcByTPQ                   [100] %12.5lf\n",Timer[100]);
  fprintf(fp,"    FirstMultiply             [101] %12.5lf\n",Timer[101]);
  fprintf(fp,"      rand   in FirstMultiply [106] %12.5lf\n",Timer[106]);
  fprintf(fp,"      mltply in FirstMultiply [107] %12.5lf\n",Timer[107]);
  fprintf(fp,"    expec_energy              [102] %12.5lf\n",Timer[102]);
  fprintf(fp,"      mltply in expec_energy  [108] %12.5lf\n",Timer[108]);
  fprintf(fp,"    expec_onebody             [103] %12.5lf\n",Timer[103]);
  fprintf(fp,"    expec_twobody             [104] %12.5lf\n",Timer[104]);
  fprintf(fp,"    Multiply                  [105] %12.5lf\n",Timer[105]);
  fprintf(fp,"    FileIO                    [109] %12.5lf\n",Timer[109]);
  fprintf(fp,"================================================\n");
  fprintf(fp,"All mltply                    [200] %12.5lf\n",Timer[200]);
  fprintf(fp,"  diagonal                    [201] %12.5lf\n",Timer[201]);
  fprintf(fp,"  trans    in HubbardGC       [202] %12.5lf\n",Timer[202]);
  fprintf(fp,"  interall in HubbardGC       [203] %12.5lf\n",Timer[203]);
  fprintf(fp,"  pairhopp in HubbardGC       [204] %12.5lf\n",Timer[204]);
  fprintf(fp,"  exchange in HubbardGC       [205] %12.5lf\n",Timer[205]);
  fprintf(fp,"\n");
  fprintf(fp,"  trans    in Hubbard         [206] %12.5lf\n",Timer[206]);
  fprintf(fp,"  interall in Hubbard         [207] %12.5lf\n",Timer[207]);
  fprintf(fp,"  pairhopp in Hubbard         [208] %12.5lf\n",Timer[208]);
  fprintf(fp,"  exchange in Hubbard         [209] %12.5lf\n",Timer[209]);
  fprintf(fp,"\n");
  fprintf(fp,"  interall in Spin            [210] %12.5lf\n",Timer[210]);
  fprintf(fp,"    double                    [220] %12.5lf\n",Timer[220]);
  fprintf(fp,"    single1                   [221] %12.5lf\n",Timer[221]);
  fprintf(fp,"    single2                   [222] %12.5lf\n",Timer[222]);
  fprintf(fp,"    else                      [223] %12.5lf\n",Timer[223]);
  fprintf(fp,"  exchange in Spin            [211] %12.5lf\n",Timer[211]);
  fprintf(fp,"    double                    [224] %12.5lf\n",Timer[224]);
  fprintf(fp,"    single1                   [225] %12.5lf\n",Timer[225]);
  fprintf(fp,"    single2                   [226] %12.5lf\n",Timer[226]);
  fprintf(fp,"    else                      [227] %12.5lf\n",Timer[226]);
  fprintf(fp,"\n");
  fprintf(fp,"  trans    in SpinGC          [212] %12.5lf\n",Timer[212]);
  fprintf(fp,"  interall in SpinGC          [213] %12.5lf\n",Timer[213]);
  fprintf(fp,"  exchange in SpinGC          [214] %12.5lf\n",Timer[214]);
  fprintf(fp,"  pairlift in SpinGC          [215] %12.5lf\n",Timer[215]);
  fprintf(fp,"================================================\n");

  fclose(fp);
  free(Timer);
  free(TimerStart);
}
