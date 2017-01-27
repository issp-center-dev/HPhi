/*
HPhi-mVMC-StdFace - Common input generator
Copyright (C) 2015 The University of Tokyo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <complex.h>

struct StdIntList {
  /*
  Parameters for LATTICE
  */
  char lattice[256];
  double a; /**< The lattice constant */
  double a0;
  double a1;
  int L;
  int W;
  double Lx;
  double Ly;
  double Wx;
  double Wy;
  int a0L;
  int a0W;
  int a1L;
  int a1W;
  int bW0;
  int bW1;
  int bL0;
  int bL1;
  int NCell;
  int **Cell;
  int NsiteUC;
  double **tau;
  /*
  Parameters for MODEL
  */
  char model[256];
  double mu;
  double complex t;
  double complex tp;
  double complex t0;
  double complex t1;
  double complex t1p;
  double complex t2;
  double complex t2p;
  double U;
  double V;
  double Vp;
  double V0;
  double V1;
  double V1p;
  double V2;
  double V2p;
  /**/
  double JAll;
  double JpAll;
  double J0All;
  double J1All;
  double J1pAll;
  double J2All;
  double J2pAll;
  double J[3][3];
  double Jp[3][3];
  double J0[3][3];
  double J1[3][3];
  double J1p[3][3];
  double J2[3][3];
  double J2p[3][3];
  double D[3][3];
  double h;
  double Gamma;
  double K;
  /*
   Phase for the boundary
  */
  double pi180;
  double phase0;
  double phase1;
  double complex ExpPhase0;
  double complex ExpPhase1;
  int AntiPeriod0;
  int AntiPeriod1;
  /*
   Transfer, Interaction, Locspin
  */
  int nsite;
  int *locspinflag;
  int ntrans;
  int Ltrans;
  int **transindx;
  double complex *trans;
  int nintr;
  int Lintr;
  int **intrindx;
  double complex *intr;
  int NCintra;
  int LCintra;
  int **CintraIndx;
  double *Cintra;
  int NCinter;
  int LCinter;
  int **CinterIndx;
  double *Cinter;
  int NHund;
  int LHund;
  int **HundIndx;
  double *Hund;
  int NEx;
  int LEx;
  int **ExIndx;
  double *Ex;
  int NPairLift;
  int LPairLift;
  int **PLIndx;
  double *PairLift;
  int lBoost;
  /*
   Calculation conditions
  */
  int lGC;
  int nelec;
  int S2;
  char outputmode[256];
  char CDataFileHead[256];
  int Sz2;
  int ioutputmode;
  /*HPhi modpara*/
  char method[256];
  char Restart[256];
  char InitialVecType[256];
  char EigenVecIO[256];
  int FlgTemp;
  int Lanczos_max;
  int initial_iv; 
  int nvec;
  int exct;
  int LanczosEps;
  int LanczosTarget;
  int NumAve;
  int ExpecInterval;
  double LargeValue;
  /*Boost*/
  int ***list_6spin_pair;
  int **list_6spin_star;
  int num_pivot;
  int ishift_nspin;
  /*Spectrum*/
  char CalcSpec[256];
  char SpectrumType[256];
  int Nomega;
  double OmegaMax;
  double OmegaMin;
  double OmegaIm;
  double SpectrumQL;
  double SpectrumQW;
  int SpectrumBody;

};
