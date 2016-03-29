/*
HPhi  -  Quantum Lattice Model Simulator
Copyright (C) 2015 Takahiro Misawa, Kazuyoshi Yoshimi, Mitsuaki Kawamura, Youhei Yamaji, Synge Todo, Naoki Kawashima

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
  /*
  Parameters for MODEL
  */
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
   Calculation conditions
  */
  int FlgTemp;
  int Lanczos_max;
  int initial_iv; 
  int nvec;
  int exct;
  int LanczosEps;
  int LanczosTarget;
  int NumAve;
  int ExpecInterval;
  int Sz2;
  int nelec;
  int ioutputmode;
  double LargeValue;
  int S2;
  /*
   Input strings
  */
  char model[256];
  char lattice[256];
  char method[256];
  char outputmode[256];
  char filehead[256];
  /*
   Parameter for lattice
  */
  double bW0;
  double bW1;
  double bL0;
  double bL1;
  int NCell;
  int **Cell;
  int NsiteUC;
  /*
   Transfer, Interaction, Locspin
  */
  int nsite;
  int *locspinflag;
  int ntrans;
  int **transindx;
  double complex *trans;
  int nintr;
  int **intrindx;
  double complex *intr;
  /*
   Boost
  */
  int lBoost;
  int ***list_6spin_pair;
  int **list_6spin_star;
  int num_pivot;
  int ishift_nspin;
  /**/
  int lGC;
};
