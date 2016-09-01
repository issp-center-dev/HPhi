/*-------------------------------------------------------------
 * HPhi
 * timer program
 * "-lrt" option is needed for clock_gettime().
 *-------------------------------------------------------------
 * original code written by Satoshi Morita
 *-------------------------------------------------------------*/

#include "Common.h"
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
//#ifdef _mpi_use
  #include <mpi.h>
//#endif


void InitTimer();
inline void StartTimer(int n);
inline void StopTimer(int n);
void OutputTimer();

void InitTimer() {
  int i;
  int NTimer=100;
  Timer       = (double*)malloc((NTimer)*sizeof(double));
  TimerStart  = (double*)malloc((NTimer)*sizeof(double));
  for(i=0;i<NTimer;i++) Timer[i]=0.0;
  for(i=0;i<NTimer;i++) TimerStart[i]=0.0;
  return;
}

void StartTimer(int n) {
#ifdef _mpi_use
  TimerStart[n]=MPI_Wtime();
#else
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME,&ts);
  TimerStart[n]=ts.tv_sec + ts.tv_nsec*1.0e-9;
#endif
  return;
}

void StopTimer(int n) {
#ifdef _mpi_use
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
  fp = fopen(fileName, "w");

  fprintf(fp,"All                         [0] %12.5lf\n",Timer[0]);

  fclose(fp);
  free(Timer);
  free(TimerStart);
}
