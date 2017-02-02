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
  double length[3];
  int L;
  int W;
  int Height;
  double direct[3][3];
  int box[3][3];
  int rbox[3][3];
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
  double complex t0p;
  double complex t1;
  double complex t1p;
  double complex t2;
  double complex t2p;
  double complex tpp;
  double U;
  double V;
  double Vp;
  double V0;
  double V0p;
  double V1;
  double V1p;
  double V2;
  double V2p;
  double Vpp;
  /**/
  double JAll;
  double JpAll;
  double J0All;
  double J0pAll;
  double J1All;
  double J1pAll;
  double J2All;
  double J2pAll;
  double JppAll;
  double J[3][3];
  double Jp[3][3];
  double J0[3][3];
  double J0p[3][3];
  double J1[3][3];
  double J1p[3][3];
  double J2[3][3];
  double J2p[3][3];
  double Jpp[3][3];
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
  double phase2;
  double complex ExpPhase0;
  double complex ExpPhase1;
  double complex ExpPhase2;
  int AntiPeriod0;
  int AntiPeriod1;
  int AntiPeriod2;
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

#if defined(_HPhi)
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
#elif defined(_mVMC)
  /*mVMC modpara*/
  char CParaFileHead[256];
  int NVMCCalMode;
  int NLanczosMode;
  int NDataIdxStart;
  int NDataQtySmp;
  int NSPGaussLeg;
  int NMPTrans;
  int NSROptItrStep;
  int NSROptItrSmp;
  int NSROptFixSmp;
  double DSROptRedCut;
  double DSROptStaDel;
  double DSROptStepDt;
  int NVMCWarmUp;
  int NVMCInterval;
  int NVMCSample;
  int NExUpdatePath;
  int RndSeed;
  int NSplitSize;
  int NStore;
  int ComplexType;
  /*
   Sub-lattice
  */
  int Lsub;
  int Wsub;
  int a0Lsub;
  int a0Wsub;
  int a1Lsub;
  int a1Wsub;
  int bW0sub;
  int bW1sub;
  int bL0sub;
  int bL1sub;
  int NCellsub;
  /*
   2-body part of the trial wavefunction
  */
  int **Orb;
  int **AntiOrb;
  int NOrb;
  int NSym;
#endif
};
