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

void current_utc_time(struct timespec *ts) {
#ifdef _OSX // OS X does not have clock_gettime, use clock_get_time
  clock_serv_t clock_mac;
  mach_timespec_t mts_mac;
  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &clock_mac);
  clock_get_time(clock_mac, &mts_mac);
  mach_port_deallocate(mach_task_self(), clock_mac);
  ts->tv_sec = mts_mac.tv_sec;
  ts->tv_nsec = mts_mac.tv_nsec;
#else
  clock_gettime(CLOCK_REALTIME, ts);
#endif

}


void StampTime(FILE *fp, char *str, int num){
  char str1[256];
  sprintf(str1, "%-40s [%04d] %12.5lf\n", str, num, Timer[num]);
  fprintf(fp, str1);
}

void InitTimer() {
  int i;
  int NTimer=10000;
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
  current_utc_time(&ts);
  //clock_gettime(CLOCK_REALTIME,&ts);
  TimerStart[n]=ts.tv_sec + ts.tv_nsec*1.0e-9;
#endif
  return;
}

void StopTimer(int n) {
#ifdef MPI
  Timer[n] += MPI_Wtime() - TimerStart[n];
#else
  struct timespec ts;
  current_utc_time(&ts);
  //clock_gettime(CLOCK_REALTIME,&ts);
  Timer[n] += ts.tv_sec + ts.tv_nsec*1.0e-9 - TimerStart[n];
#endif
  return;
}

void OutputTimer(struct BindStruct *X) {
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
    if(X->Def.iCalcType==TPQCalc) {
      StampTime(fp, "  CalcByTPQ", 3000);
      StampTime(fp, "    FirstMultiply", 3100);
      StampTime(fp, "      rand   in FirstMultiply", 3101);
      StampTime(fp, "      mltply in FirstMultiply", 3102);
      StampTime(fp, "    expec_energy             ", 3200);
      StampTime(fp, "      mltply in expec_energy ", 3201);
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
      StampTime(fp, "    expec_energy", 4300);
      StampTime(fp, "      mltply in expec_energy ", 4301);
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
}

